#!/usr/bin/env python3
"""
Boundary convergence module.

  For every individual sample it repeats, up to the coverage walk, what the main script does for
  the merged (cumulative) coverage: load the sample's coverage, get the intergenic background,
  check the expression level of each gene, and run the coverage walk. This gives one walk TSS per
  (gene, sample). The walk TSS positions of all samples of a gene are then summarised
  (centre, spread, clusters, closeness to the TSS points found by the main script).

Two stop rules (they apply to the whole run, not to single genes)
  * fewer than `min_samples` individual samples            -> nothing is run
  * walk TSS found in < `min_fraction` of all the samples  -> the gene is flagged low_confidence,
                                                              no statistics, no plot should be drawn
"""
import bisect
import contextlib
import io
import math
import multiprocessing
import os
import statistics

import numpy as np
import pandas as pd

_CTX = {}  # read-only context shared with forked worker processes (avoids pickling the genome etc.)

CUMULATIVE_FILES = ('Cumulative_aligned_reads.cov', 'Cumulative_spanning_reads.cov')  # fixed names given by the main script

def list_sample_cov_files(folder):
	"""List the per-sample .cov files in a folder, leaving out the merged cumulative file the main
	script writes into the same folder (its name is fixed, see CUMULATIVE_FILES)."""
	return sorted(
		os.path.join(folder, f)
		for f in os.listdir(folder)
		if f.endswith('.cov') and f not in CUMULATIVE_FILES
	)

def pair_sample_files(cov_files, scov_files):
	"""Pair the aligned and the spanning cov file of each sample. Only two cases are accepted:
	   NAME.cov          (aligned folder)  <->  NAME.cov           (spanning folder)
	   NAME_aligned.cov  (aligned folder)  <->  NAME_spanning.cov  (spanning folder)
	A mix such as NAME_aligned.cov with NAME.cov is not accepted. If any file finds no partner, it
	stops and lists the files without a partner."""
	scov_by_file_name = {os.path.basename(f): f for f in scov_files}
	pairs = []
	partners_used = set()
	aligned_without_partner = []
	for cov_file in sorted(cov_files):
		file_name = os.path.basename(cov_file)
		if file_name.endswith('_aligned.cov'):
			name = file_name[:-len('_aligned.cov')]
			partner = name + '_spanning.cov'
		else:
			name = file_name[:-len('.cov')]
			partner = file_name
		if partner in scov_by_file_name:
			pairs.append((name, cov_file, scov_by_file_name[partner]))
			partners_used.add(partner)
		else:
			aligned_without_partner.append(file_name)
	spanning_without_partner = sorted(set(scov_by_file_name) - partners_used)
	if aligned_without_partner or spanning_without_partner:
		raise ValueError('Aligned and spanning cov files do not match (NAME.cov with NAME.cov, or NAME_aligned.cov with '
						 'NAME_spanning.cov). Aligned files without a partner: ' + str(aligned_without_partner[:5])
						 + ' Spanning files without a partner: ' + str(spanning_without_partner[:5]))
	return pairs

# intergenic windows: positions only, so they are worked out once for all samples
def get_intergenic_windows(genes_per_chromosome, gene_infos, intergenic_buffer, intergenic_region_size):
	"""Start position of every intergenic window used for the background.
	Same rules as get_intergenic_background_seqs of the main script (region between two neighbouring
	genes, buffer left out on both sides, rest cut into windows of intergenic_region_size), but only the
	positions are worked out here. They do not depend on the sample, so this runs once."""
	windows = {}
	for contig in genes_per_chromosome:
		gene_list = genes_per_chromosome[contig]
		window_starts = []
		for i, gene in enumerate(gene_list):
			if i + 1 >= len(gene_list):  #skip the last gene in the contig
				break
			next_gene = gene_list[i + 1]
			intergenic_start = gene_infos[gene]['end']
			intergenic_end = gene_infos[next_gene]['start']
			intergenic_length = intergenic_end - intergenic_start + 1
			if intergenic_length <= 0:
				continue#skip if the intergenic region belongs to the overlapping gene region
			if intergenic_length > intergenic_buffer and intergenic_length > intergenic_region_size:
				usable_length = max(0, intergenic_length - 2 * intergenic_buffer)#leave a buffer near the gene regions for both the genes
				chops = int(usable_length / intergenic_region_size)#only whole windows are used
				first_window_start = intergenic_start + intergenic_buffer
				for k in range(chops):
					window_starts.append(first_window_start + k * intergenic_region_size)
		windows[contig] = np.array(window_starts, dtype=np.int64)
	return windows

# sample-dependent helper functions
def load_coverage_arrays(cov_file):
	"""Read a cov file (bedtools -d output, plain or .gz) into one array of coverage values per contig.
	Same layout as the lists returned by load_coverage of the main script: the value of position p
	is at index p - 1. Nothing else is built (no tuples, no dictionaries)."""
	table = pd.read_csv(cov_file, sep='\t', header=None, usecols=[0, 2], dtype={0: 'category', 2: 'float64'})
	contig_codes = table[0].cat.codes.to_numpy()
	contig_names = table[0].cat.categories
	values = table[2].to_numpy()
	first_rows = np.concatenate(([0], np.flatnonzero(contig_codes[1:] != contig_codes[:-1]) + 1))#row where each contig starts
	last_rows = np.concatenate((first_rows[1:], [len(values)]))
	return {contig_names[contig_codes[a]]: values[a:b] for a, b in zip(first_rows, last_rows)}

def get_intergenic_background(coverage, intergenic_windows, intergenic_region_size):
	"""Intergenic background of ONE sample = average coverage of every intergenic window, as returned
	(first return value) by get_intergenic_background_seqs of the main script: windows with zero
	average are left out and the outliers above q3 + 1.5 * IQR are removed. Returned sorted, ready for
	the bisect inside find_gene_exp_level. No plots and no extra statistics are made here."""
	window_averages = []
	for contig, window_starts in intergenic_windows.items():
		if len(window_starts) == 0:
			continue
		running_sum = np.concatenate(([0.0], np.cumsum(coverage[contig])))#running_sum[n] = sum of the first n positions
		window_sums = running_sum[window_starts - 1 + intergenic_region_size] - running_sum[window_starts - 1]
		window_averages.append(window_sums / intergenic_region_size)
	cov_arr = np.concatenate(window_averages) if window_averages else np.array([])
	cov_arr = cov_arr[cov_arr != 0.0]#adding only non-zero average coverage values to the background noise
	if cov_arr.size == 0:
		return [0.0]
	q1, q3 = np.percentile(cov_arr, [25, 75])
	iqr = q3 - q1
	intergenic_window_coverages_wo_outliers = cov_arr[cov_arr <= q3 + 1.5*iqr]
	return sorted(intergenic_window_coverages_wo_outliers.tolist())

class PositionLookup:
	"""allows find_gene_exp_level (written for a dictionary position -> coverage) read from an array:
	coverage_lookup[position] gives the value stored at index position - 1."""
	def __init__(self, coverage_array):
		self.coverage_array = coverage_array
	def __getitem__(self, position):
		return self.coverage_array[position - 1]


def find_gene_exp_level(intergenic_window_coverages, genes_per_chromosome, coverage_lookup, start, end, gene,intergenic_region_size,pvalue):
	passed_windows = 0
	tot_windows = 0
	if abs(end - start) > intergenic_region_size:
		B = len(intergenic_window_coverages)
		cumulative_coverage = sum(coverage_lookup[pos] for pos in range(start, start + intergenic_region_size))  #build once
		while start <= end - intergenic_region_size:
			# calculating average coverage of the window
			avg_coverage = float(cumulative_coverage / intergenic_region_size)
			idx = bisect.bisect_left(intergenic_window_coverages, avg_coverage)
			hits = B - idx
			psig = float(hits / B)
			if psig < pvalue:
				passed_windows+=1
			tot_windows+=1
			cumulative_coverage += coverage_lookup[start + intergenic_region_size] - coverage_lookup[start]  #sliding the sum instead of rebuilding
			start +=1
		percentage_passed_windows = (float(passed_windows) / float(tot_windows))*100
		if percentage_passed_windows <= 10:
			exp_status = 'low'
		elif 10 < percentage_passed_windows <= 60:
			exp_status = 'moderate'
		elif  percentage_passed_windows > 60:
			exp_status = 'high'
	else:
		exp_status = None
	return exp_status

# sample-independent part, done once
def build_gene_table(goi, shared, params, helpers):
	"""Walk origin (start/end), orientation, hard cutoff for every gene of interest.
	Genes skipped by the main script's own rules (not protein encoding, blocking overlap) are
	returned in `skipped`."""
	gene_infos = shared['gene_infos']
	genes_per_chromosome = shared['genes_per_chromosome']
	protein_coding_genes_per_chromosome = shared['protein_coding_genes_per_chromosome']
	transcripts_per_gene = shared['transcripts_per_gene']
	five_utr_infos = shared['five_utr_infos']
	cds_infos = shared['cds_infos']
	genome_seq = shared['genome_seq']
	protein_encoding = params['protein_encoding']
	coverage_walk_origin = params['coverage_walk_origin']
	window = params['window']
	find_flanking_genes = helpers['find_flanking_genes']
	get_overlap_type = helpers['get_overlap_type']

	gene_table = {}
	skipped = {}
	five_utr_dic = {}
	for gene in goi:
		if protein_encoding == "yes":
			if gene not in transcripts_per_gene:
				print(f'You chose only protein encoding gene analysis mode. But your GOI {gene} not found in protein-encoding gene group. Skipping TSS analysis for {gene}\n')
				skipped[gene] = 'not in protein-encoding gene group'
				continue
		try:
			seq_per_contig = genome_seq[gene_infos[gene]['chromosome']]  # get the sequence of the contig/pseudochromosome that harbours the gene of interest
			# code block to check if the goi has annotated 5'UTR and if yes take the most upstream/ downstream 5'UTR start/ end as start or end according to + or - strand orientation
			# Initialize with gene coordinates as default with 5'UTR checks downstream
			gene_start, gene_end, orientation, contig = gene_infos[gene]['start'], gene_infos[gene]['end'], gene_infos[gene]['orientation'], gene_infos[gene]['chromosome']  # get information about gene of interest if it does not have 5'UTRs annotated
			if gene in transcripts_per_gene:
				transcript_list = transcripts_per_gene[gene]
				for each in transcript_list:
					if coverage_walk_origin == 'utr':
						if each in five_utr_infos and gene_infos[gene]['orientation'] == '+' and each == transcript_list[0]:
							start = min(utr_start for utr_start, utr_end in five_utr_infos[each])   # in case a gene has transcripts with 5'UTR annotated, the most upstream 5'UTR start will be taken as the walk start for the + strand gene
							end, orientation = gene_infos[gene]['end'], gene_infos[gene]['orientation']  # get information about gene of interest
							five_utr_dic[gene] = f"5'UTR start of {each} used for TSS prediction."
							break
						elif each in five_utr_infos and gene_infos[gene]['orientation'] == '-' and each == transcript_list[-1]:
							end = max(utr_end for utr_start, utr_end in five_utr_infos[each])  # in case a gene has transcripts with 5'UTR annotated, the most downstream 5'UTR end will be taken as the start for the - strand gene coverage walk
							start, orientation = gene_infos[gene]['start'], gene_infos[gene]['orientation']  # get information about gene of interest
							five_utr_dic[gene] = f"5'UTR end of {each} used for TSS prediction."
							break
					elif coverage_walk_origin == 'cds':
						if each in cds_infos and gene_infos[gene]['orientation'] == '+' and each == transcript_list[0]:
							cds_infos[each].sort(key=lambda cds: cds[0])
							start = cds_infos[each][0][0]  # the most upstream CDS start will be taken as the coverage walking start for the + strand gene
							end, orientation = gene_infos[gene]['end'], gene_infos[gene]['orientation']  # get information about gene of interest
							five_utr_dic[gene] = f"CDS start of {each} used for TSS prediction."
							break
						elif each in cds_infos and gene_infos[gene]['orientation'] == '-' and each == transcript_list[-1]:
							cds_infos[each].sort(key=lambda cds: cds[1])
							end = cds_infos[each][-1][1]  # the most upstream CDS end will be taken as the coverage walking start for the - strand gene
							start, orientation = gene_infos[gene]['start'], gene_infos[gene]['orientation']  # get information about gene of interest
							five_utr_dic[gene] = f"CDS end of {each} used for TSS prediction."
							break
			if coverage_walk_origin == 'utr':
				if gene not in five_utr_dic and orientation == '+':
					five_utr_dic[gene] = f"No 5'UTR annotated. {gene} start used for TSS prediction."
					start = gene_start
					end = gene_end
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No 5'UTR annotated. {gene} end used for TSS prediction."
					start = gene_start
					end = gene_end
			elif coverage_walk_origin == 'cds':#gene not in five_utr_dic is not the right check here, nevertheless it is harmless/ benign and may need a rewrite later
				if gene not in five_utr_dic and orientation == '+':
					five_utr_dic[gene] = f"No CDS annotated. {gene} start used for TSS prediction."
					start = gene_start
					end = gene_end
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No CDS annotated. {gene} end used for TSS prediction."
					start = gene_start
					end = gene_end
			if protein_encoding == "no":
				upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes( gene, gene_infos, genes_per_chromosome, window )
			else:
				upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes(gene, gene_infos,protein_coding_genes_per_chromosome,window)
			#check if goi is an overlapping gene
			#TSS-blocking overlap types
			SKIP_TSS_TYPES = {'head_head', 'head_into_neighbor', 'same_strand', 'nested'}#'tail_tail' and 'tail_into_neighbor' overlap types are safe for TSS analysis
			goi_strand = gene_infos[gene]['orientation']
			blocking_overlaps = 0
			for ugene in upstream_gene_list:
				nbr = gene_infos[ugene]

				if gene_start <= nbr['end']:  # positional overlap exists
					ov_type = get_overlap_type(goi_strand, gene_start, gene_end,nbr['orientation'], nbr['start'], nbr['end'])
					print(f"  {gene} - {ugene}: {ov_type} overlap")
					if ov_type in SKIP_TSS_TYPES:
						blocking_overlaps += 1
				#'tail_tail' → do NOT increment; TSS analysis can still run
			#downstream check is re-enabled now because strand matters:
			#a downstream neighbor can cause a head_head overlap for − strand GOIs
			for dgene in downstream_gene_list:
				nbr = gene_infos[dgene]
				if gene_end >= nbr['start']:  # positional overlap exists
					ov_type = get_overlap_type( goi_strand, gene_start, gene_end,nbr['orientation'], nbr['start'], nbr['end'])
					print(f"  {gene} - {dgene}: {ov_type} overlap")
					if ov_type in SKIP_TSS_TYPES:
						blocking_overlaps += 1
			if blocking_overlaps > 0:
				print(f"{gene} has {blocking_overlaps} overlap(s). TSS analysis skipped.")
				skipped[gene] = f'{blocking_overlaps} blocking overlap(s)'
				continue
			if orientation == "+":	#only works on forward strand
				if upstream_gene:
					hard_cutoff = gene_infos[ upstream_gene ]['end']

				else:
					hard_cutoff = 1
			else:	#solution for reverse strand genes
				if downstream_gene:
					hard_cutoff = gene_infos[ downstream_gene ]['start']
				else:
					hard_cutoff = len( seq_per_contig )
			print(f"hardcutoff of {gene} is set to {hard_cutoff}")
			gene_table[gene] = {
				'contig': contig, 'orientation': orientation,
				'start': start, 'end': end,
				'gene_start': gene_start, 'gene_end': gene_end,
				'hard_cutoff': hard_cutoff,
			}
		except KeyError as e:
			print(f"Missing gene error: {gene}, missing key: {e}")
			skipped[gene] = f'missing key {e}'
	return gene_table, skipped

# the coverage walks, copied from run_fwd_analysis / run_rev_analysis (up to walk TSS)
def fwd_walk(gene, cov_per_contig, scov_per_contig, seq_per_contig, start, mincov, min_exon_size,
			 hard_cutoff, tolerated_gap, splicesites, percentdiff_threshold):
	"""! @brief run analysis on forward strand """
	hard_cutoff_reached = False
	intron_boundary_marker = {}#dictionary to collect intron boundary positions where key is donor splice site and value is acceptor splice site
	most_upstream_pos = start
	final_pos_status = False
	pending_acceptor = None
	while not final_pos_status:
		# --- walk coverage upstream of transcription start while there is coverage --- #
		while cov_per_contig[most_upstream_pos - 2] >= mincov:  # index = genomic position -1 but coverage of gene that is next to the current position needs to be assessed before moving in there
			if scov_per_contig[most_upstream_pos - 2] !=0:
				cov_percent_diff = ((scov_per_contig[most_upstream_pos - 2] - cov_per_contig[most_upstream_pos - 2])/scov_per_contig[most_upstream_pos - 2])*100
			elif scov_per_contig[most_upstream_pos - 2] ==0:
				cov_percent_diff = percentdiff_threshold #if spanning cov is zero then it is obviously not an intron. So shortcircuiting the check by assigning cov percent diff a value equal to the percent diff
			boundary_marker = 0#this is a trick flag to mark the intron acceptor splice site position when the while loop to check percent diff begins
			intron_acceptor_splice_site_position = None
			hit_hard_cutoff = False#this is also a trick flag to turn on hard cutoff hit if hardcutoff is reached while the walk is happening in the bounds of a likely intron territory
			low_coverage_regime = False #trick flag to track aligned read coverage value relative to mincov in the percent diff while loop
			while cov_percent_diff > percentdiff_threshold:#in-line intron boundary marking during the coverage walk itself
				if boundary_marker == 0:
					intron_acceptor_splice_site_position = most_upstream_pos - 1 # most_upstream_pos - 1 because always the coverage walk starts from CDS start for + and CDS end for - gene, so the first current position is definitiely not an intron-exon boundary. So check starts from next position onwards.
					boundary_marker = 1
				most_upstream_pos -= 1 # move one step upstream
				if most_upstream_pos == hard_cutoff:
					hit_hard_cutoff = True
					break
				if cov_per_contig[most_upstream_pos - 2] < mincov:#guard to check if aligned read coverage falls less than mincov within this percent diff checking while loop since mincov is checked only in the outer while loop
					low_coverage_regime = True
					break
				if scov_per_contig[most_upstream_pos - 2] != 0:#compute the percentage cov diff in the next upstream position for the next iteration of the while loop
					cov_percent_diff = ((scov_per_contig[most_upstream_pos - 2] - cov_per_contig[most_upstream_pos - 2]) / scov_per_contig[most_upstream_pos - 2]) * 100
				elif scov_per_contig[most_upstream_pos - 2] == 0:
					cov_percent_diff = percentdiff_threshold  # if spanning cov is zero then it is obviously not an intron. So shortcircuiting the check by assigning cov percent diff a value equal to the percent diff
			if hit_hard_cutoff == True:#skip crossing the intron block when hard cutoff is already reached while marking the intron boundaries
				hard_cutoff_reached = 'True;'+str(hard_cutoff)
				final_pos_status = True
			elif intron_acceptor_splice_site_position is not None and low_coverage_regime == False:
				intron_donor_splice_site_position = most_upstream_pos# the intron marking loop already checks percent difference in the next upstream pos. so the while loop must have exited when the next upstream pos did not satisfy the percent criteria. therefore the intron boundary should be most_upstream_pos and not most_upstream_pos -1 which is the position that did not satisfy the percent diff condition to stay in the while loop
				if intron_donor_splice_site_position != intron_acceptor_splice_site_position:
					avg_gap_coverage = sum(scov_per_contig[intron_donor_splice_site_position - 1:intron_acceptor_splice_site_position - 1]) / (intron_acceptor_splice_site_position - intron_donor_splice_site_position)
					if most_upstream_pos > min_exon_size and avg_gap_coverage > mincov:#avg spanning gap coverage should be greater than mincov and most upstream pos should be greater than min exon size. matters for first contig end
						# --- check coverage gaps for (canonical) splice sites to continue across introns --- #
						donor_splice_site = seq_per_contig[intron_donor_splice_site_position - 1:intron_donor_splice_site_position + 1].upper()  # this should be GT for it to be canonical splice site
						acceptor_splice_site = seq_per_contig[intron_acceptor_splice_site_position - 2:intron_acceptor_splice_site_position].upper()  # this should be AG for it to be canonical splice site; changed the slicing from -3 and -1 respectively to -2 and nothing as only then the intron acceptor boundaries are correctly considered inclusively
						if splicesites == "off":  # ignore check for canonical splice sites #tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
							if intron_acceptor_splice_site_position - intron_donor_splice_site_position < tolerated_gap:
								intron_boundary_marker[intron_donor_splice_site_position] = intron_acceptor_splice_site_position
								most_upstream_pos = most_upstream_pos - 1
							else:#if splice sites off and tolerated gap check fails, don't mark intron boundaries. just walk
								most_upstream_pos = most_upstream_pos - 1
						elif donor_splice_site == "GT" and acceptor_splice_site == "AG":
							print(f"Canonical donor splice site starting at {intron_donor_splice_site_position}: " + donor_splice_site + '\n')
							print(f"Canonical acceptor splice site ending at {intron_acceptor_splice_site_position}: " + acceptor_splice_site + '\n')
							intron_boundary_marker[intron_donor_splice_site_position] = intron_acceptor_splice_site_position
							most_upstream_pos = most_upstream_pos - 1
						elif donor_splice_site != "GT" or acceptor_splice_site != "AG":#tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
							print(f"Warning! Either one or both the donor and acceptor splice sites are non-canonical."+ '\n')
							print(f"Donor splice site starting at {intron_donor_splice_site_position}: " + donor_splice_site+ '\n')
							print(f"Acceptor splice site ending at {intron_acceptor_splice_site_position}: " + acceptor_splice_site+ '\n')
							if intron_acceptor_splice_site_position - intron_donor_splice_site_position < tolerated_gap:
								intron_boundary_marker[intron_donor_splice_site_position] = intron_acceptor_splice_site_position
								most_upstream_pos = most_upstream_pos - 1
							else:#if non canonical splice sites show up and tolerated gap check fails, don't mark intron boundaries. just walk
								most_upstream_pos = most_upstream_pos - 1
					else:
						final_pos_status = True
				elif intron_donor_splice_site_position == intron_acceptor_splice_site_position:#to account for one base noise that looks like intron but is just noise and walk pas it
					most_upstream_pos = most_upstream_pos - 1

			else:
				"""
				#the following additional check is placed here since the percent diff while loop exits with a break when the next upstream position is less than mincov; if i walk upstream here again then i would be walking one more base than 
				I should have and the outer while loop will not look at the same next upstream position that was looked at by the percent diff and would be looking at the next upstream position instead
				"""
				if low_coverage_regime == False:
					most_upstream_pos -= 1  # move one step upstream
				else:
					pending_acceptor = intron_acceptor_splice_site_position # to facilitate correct carry over of marked acceptor from high to low coverage regime
			if most_upstream_pos == hard_cutoff:  # stop if end of upstream contig/pseudochromosome is reached
				hard_cutoff_reached = 'True;'+str(hard_cutoff)
				break

		# --- try to cross intron --- #
		"""
		#this block to fires when block one while loop exits due to cov being less than mincov
		#the intron crossing logic is repeated here again similar to the above block with the moving variable current_position and fixed variable most_upstream_position
		to specifically catch cases where the cov is exactly 0 in the intron region while spanning cov is finite for a default mincov of 1
		"""
		current_position = most_upstream_pos - 1
		if current_position > hard_cutoff:
			# code block to handle noisy cases where spurious high cov regime like reads show up in low cov regime and vice versa
			# this is being added to make the transition marking pending_acceptor from the high cov regime to the low cov one mark actual intron boundaries and not intermediate points
			while True:
				aligned_cov = cov_per_contig[current_position - 2] #-2 here becuase for this block, current_position is the running variable. I look for the coverage in the next position upstream (current_position - 1) beforing stepping on to it
				spanning_cov = scov_per_contig[current_position -2]
				if spanning_cov != 0:
					divergence_percent_diff = (((spanning_cov - aligned_cov)) / spanning_cov) * 100
				else:
					divergence_percent_diff = percentdiff_threshold  # trick to equate the divergence to the threshold so that it fails the next condition check as spanning cov is 0
				if aligned_cov >= mincov and divergence_percent_diff <= percentdiff_threshold:
					break
				current_position -= 1
				if current_position == hard_cutoff:
					print(f'current position is {current_position}')
					hard_cutoff_reached = 'True;' + str(hard_cutoff)
					break
			acceptor_pos = pending_acceptor if pending_acceptor is not None else most_upstream_pos
			"""
			while cov_per_contig[ current_position - 2 ] < mincov:#check if upstream position has low coverage
				current_position -= 1
				if current_position == hard_cutoff:
					hard_cutoff_reached = 'True;'+str(hard_cutoff)
					break
			acceptor_pos = pending_acceptor if pending_acceptor is not None else most_upstream_pos
			"""
			avg_gap_coverage = sum( scov_per_contig[ current_position-1:acceptor_pos-1 ] )/(acceptor_pos-current_position)
			#average coverage in intron should be very low
			if current_position > min_exon_size and avg_gap_coverage > mincov:
			# --- check coverage gaps for (canonical) splice sites to continue across introns --- #
				donor_splice_site = seq_per_contig[current_position - 1:current_position + 1].upper()  # this should be GT for it to be canonical splice site
				acceptor_splice_site = seq_per_contig[acceptor_pos - 2:acceptor_pos].upper()  # this should be AG for it to be canonical splice site; changed the slicing from -3 and -1 respectively to -2 and nothing as only then the intron acceptor boundaries are correctly considered inclusively
				if splicesites == "off":  # ignore check for canonical splice sites #tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
					if acceptor_pos - current_position < tolerated_gap:
						intron_boundary_marker[current_position] = acceptor_pos
						most_upstream_pos = current_position - 1
					else:#if splice sites off and tolerated gap check fails, make final pos status true; making the check stringent for the loop with cov less than mincov
						final_pos_status = True
						most_upstream_pos = acceptor_pos
				elif donor_splice_site == "GT" and acceptor_splice_site == "AG":
					print(f"Low coverage regime: Canonical donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Low coverage regime: Canonical acceptor splice site ending at {acceptor_pos}: " + acceptor_splice_site + '\n')
					intron_boundary_marker[current_position] = acceptor_pos
					most_upstream_pos = current_position - 1
				elif donor_splice_site != "GT" or acceptor_splice_site != "AG":#tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
					print(f"Low coverage regime: Warning! Either one or both the donor and acceptor splice sites are non-canonical."+ '\n')
					print(f"Low coverage regime: Donor splice site starting at {current_position}: " + donor_splice_site+ '\n')
					print(f"Low coverage regime: Acceptor splice site ending at {acceptor_pos}: " + acceptor_splice_site+ '\n')
					if acceptor_pos - current_position < tolerated_gap:
						intron_boundary_marker[current_position] = acceptor_pos
						most_upstream_pos = current_position - 1
					else:#if non canonical splice sites show up and tolerated gap check fails, make final pos status true; making the check stringent for the loop with cov less than mincov
						final_pos_status = True
						most_upstream_pos = acceptor_pos
			else:
				final_pos_status = True
				most_upstream_pos = acceptor_pos
			pending_acceptor = None
		else:
			final_pos_status = True

	print("Walk TSS position of " + gene + ": " + str(most_upstream_pos))
	walk_tss = most_upstream_pos
	return walk_tss, hard_cutoff_reached

def rev_walk(gene, cov_per_contig, scov_per_contig, seq_per_contig, end, mincov, min_exon_size,
			 hard_cutoff, tolerated_gap, splicesites, percentdiff_threshold):
	"""! @brief run analysis on reverse strand """
	hard_cutoff_reached = False
	intron_boundary_marker = {}
	most_downstream_pos = end
	final_pos_status = False
	pending_acceptor = None
	while not final_pos_status:
		# --- walk coverage upstream of transcription start while there is coverage --- #
		while cov_per_contig[most_downstream_pos] >= mincov:  # index = next genomic position but the coverage of the next successive position needs to be looked at and not the current position
			if scov_per_contig[most_downstream_pos] != 0:
				cov_percent_diff = ((scov_per_contig[most_downstream_pos] - cov_per_contig[most_downstream_pos])/scov_per_contig[most_downstream_pos])*100
			else:
				cov_percent_diff = percentdiff_threshold
			boundary_marker = 0
			intron_acceptor_splice_site_position = None
			hit_hard_cutoff = False
			low_coverage_regime = False
			while cov_percent_diff > percentdiff_threshold:
				if boundary_marker == 0:
					intron_acceptor_splice_site_position = most_downstream_pos + 1
					boundary_marker = 1
				most_downstream_pos += 1
				if most_downstream_pos == hard_cutoff:
					hit_hard_cutoff = True
					break
				if cov_per_contig[most_downstream_pos] < mincov:
					low_coverage_regime = True
					break
				if scov_per_contig[most_downstream_pos] != 0:
					cov_percent_diff = ((scov_per_contig[most_downstream_pos] - cov_per_contig[most_downstream_pos]) / scov_per_contig[most_downstream_pos]) * 100
				elif scov_per_contig[most_downstream_pos] == 0:
					cov_percent_diff = percentdiff_threshold
			if hit_hard_cutoff == True:
				hard_cutoff_reached = 'True;'+str(hard_cutoff)
				final_pos_status = True
			elif intron_acceptor_splice_site_position is not None and low_coverage_regime == False:
				intron_donor_splice_site_position = most_downstream_pos
				if intron_donor_splice_site_position != intron_acceptor_splice_site_position:
					avg_gap_coverage = sum(scov_per_contig[intron_acceptor_splice_site_position - 1:intron_donor_splice_site_position - 1]) / (intron_donor_splice_site_position - intron_acceptor_splice_site_position)
					if most_downstream_pos < (len(seq_per_contig) - min_exon_size) and avg_gap_coverage > mincov:
						donor_splice_site = seq_per_contig[intron_donor_splice_site_position - 2:intron_donor_splice_site_position].upper()
						acceptor_splice_site = seq_per_contig[intron_acceptor_splice_site_position - 1:intron_acceptor_splice_site_position + 1].upper()
						if splicesites == "off":
							if intron_donor_splice_site_position - intron_acceptor_splice_site_position < tolerated_gap:
								intron_boundary_marker[intron_acceptor_splice_site_position] = intron_donor_splice_site_position #the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
								most_downstream_pos += 1
							else:
								most_downstream_pos += 1
						elif donor_splice_site == 'AC' and acceptor_splice_site == "CT":
							print(f" Canonical donor splice site starting at {intron_donor_splice_site_position}: " + donor_splice_site + '\n')
							print(f" Canonical acceptor splice site ending at {intron_acceptor_splice_site_position}: " + acceptor_splice_site + '\n')
							intron_boundary_marker[intron_acceptor_splice_site_position] = intron_donor_splice_site_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
							most_downstream_pos += 1
						elif donor_splice_site != 'AC' or acceptor_splice_site != "CT":
							print(f"Warning! Either one or both the donor and acceptor splice sites are non-canonical." + '\n')
							print(f"Donor splice site starting at {intron_donor_splice_site_position}: " + donor_splice_site + '\n')
							print(f"Acceptor splice site ending at {intron_acceptor_splice_site_position}: " + acceptor_splice_site + '\n')
							if intron_donor_splice_site_position - intron_acceptor_splice_site_position < tolerated_gap:
								intron_boundary_marker[intron_acceptor_splice_site_position] = intron_donor_splice_site_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
								most_downstream_pos += 1
							else:
								most_downstream_pos += 1
					else:
						final_pos_status = True
				elif intron_donor_splice_site_position == intron_acceptor_splice_site_position:  # to account for one base noise that looks like intron but is just noise and walk pas it
					most_downstream_pos = most_downstream_pos + 1
			else:
				if low_coverage_regime == False:
					most_downstream_pos += 1
				else:
					pending_acceptor = intron_acceptor_splice_site_position #to carry over correct acceptor spliced marked from high cov region over to low coverage regime
			if most_downstream_pos == hard_cutoff:  # stop if end of contig/pseudochromosome is reached
				hard_cutoff_reached = 'True;'+str(hard_cutoff)
				break
		print(f'most downstream pos is {most_downstream_pos}')
		# --- try to cross intron --- #
		current_position = most_downstream_pos + 1  # most_downstream_pos has coverage above cutoff (position, not index!)
		if current_position < hard_cutoff:
			#code block to handle noisy cases where spurious high cov regime like reads show up in low cov regime and vice versa
			#this is being added to make the transition marking pending_acceptor from the high cov regime to the low cov one mark actual intron boundaries and not intermediate points
			while True:
				aligned_cov = cov_per_contig[current_position]
				spanning_cov = scov_per_contig[current_position]
				if spanning_cov != 0:
					divergence_percent_diff = (((spanning_cov - aligned_cov))/spanning_cov) * 100
				else:
					divergence_percent_diff = percentdiff_threshold#trick to equate the divergence to the threshold so that it fails the next condition check as spanning cov is 0
				if aligned_cov >= mincov and divergence_percent_diff <= percentdiff_threshold:
					break
				current_position += 1
				if current_position == hard_cutoff:
					print(f'current position is {current_position}')
					hard_cutoff_reached = 'True;'+str(hard_cutoff)
					break
			acceptor_pos = pending_acceptor if pending_acceptor is not None else most_downstream_pos
			"""
			while cov_per_contig[current_position] < mincov:  # check if downstream position has low coverage
				current_position += 1  # move one step downstream
				if current_position == hard_cutoff:
					print(f'current position is {current_position}')
					hard_cutoff_reached = 'True;'+str(hard_cutoff)
					break
			acceptor_pos = pending_acceptor if pending_acceptor is not None else most_downstream_pos
			"""
			avg_gap_coverage = sum( scov_per_contig[ acceptor_pos-1 : current_position-1 ] )/(current_position - acceptor_pos)
			if current_position < (len(seq_per_contig) - min_exon_size) and avg_gap_coverage > mincov:
				donor_splice_site = seq_per_contig[current_position - 2:current_position].upper()
				acceptor_splice_site = seq_per_contig[acceptor_pos - 1:acceptor_pos + 1].upper()
				if splicesites == "off":
					if current_position - acceptor_pos < tolerated_gap:
						intron_boundary_marker[acceptor_pos] = current_position #the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
						most_downstream_pos = current_position + 1
					else:
						final_pos_status = True
						most_downstream_pos = acceptor_pos
				elif donor_splice_site == 'AC' and acceptor_splice_site == "CT":
					print(f"Low coverage regime: Canonical donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Low coverage regime: Canonical acceptor splice site ending at {acceptor_pos}: " + acceptor_splice_site + '\n')
					intron_boundary_marker[acceptor_pos] = current_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
					most_downstream_pos = current_position + 1
				elif donor_splice_site != 'AC' or acceptor_splice_site != "CT":
					print(f"Low coverage regime: Warning! Either one or both the donor and acceptor splice sites are non-canonical." + '\n')
					print(f"Low coverage regime: Donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Low coverage regime: Acceptor splice site ending at {acceptor_pos}: " + acceptor_splice_site + '\n')
					if current_position - acceptor_pos < tolerated_gap:
						intron_boundary_marker[acceptor_pos] = current_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
						most_downstream_pos = current_position + 1
					else:
						final_pos_status = True
						most_downstream_pos = acceptor_pos
			else:
				final_pos_status = True
				most_downstream_pos = acceptor_pos
			pending_acceptor = None
		else:
			final_pos_status = True

	print("Walk TSS position of " + gene + ": " + str(most_downstream_pos))
	walk_tss = most_downstream_pos
	return walk_tss, hard_cutoff_reached

# sample-dependent part, once per sample
def process_sample(sample, cov_file, scov_file, gene_table, shared, params, intergenic_windows, verbose):
	"""The steps of the main script for ONE sample: load the coverage, get the intergenic background,
	then for every gene check the expression level and run the coverage walk.
	Returns the walk TSS of every gene that was walked in this sample."""
	genes_per_chromosome = shared['genes_per_chromosome']
	genome_seq = shared['genome_seq']
	sink = io.StringIO()
	out = {'sample': sample, 'walks': {}, 'expression': {}, 'not_expressed': [], 'errors': {}}
	with (contextlib.nullcontext() if verbose else contextlib.redirect_stdout(sink)):
		coverage = load_coverage_arrays(cov_file)
		scoverage = load_coverage_arrays(scov_file)

		intergenic_window_coverages = get_intergenic_background(coverage, intergenic_windows, params['intergenic_region_size'])

		for gene, g in gene_table.items():
			try:
				cov_per_contig = coverage[g['contig']]  # get coverage of the sequence that harbours the gene of interest
				scov_per_contig = scoverage[g['contig']]  # get spanning read coverage of the sequence that harbours the gene of interest
				seq_per_contig = genome_seq[g['contig']]  # get the sequence of the contig/pseudochromosome that harbours the gene of interest
				coverage_lookup = PositionLookup(cov_per_contig)
				gene_exp_status = find_gene_exp_level(intergenic_window_coverages, genes_per_chromosome, coverage_lookup, g['gene_start'], g['gene_end'], gene, params['intergenic_region_size'], params['background_percentage'])
				if gene_exp_status == None or gene_exp_status == 'low':
					out['not_expressed'].append(gene)
					continue
				if g['orientation'] == "+":
					walk_tss, _ = fwd_walk(gene, cov_per_contig, scov_per_contig, seq_per_contig, g['start'], params['mincov'], params['min_exon_size'], g['hard_cutoff'], params['tolerated_gap'], params['splicesites'], params['percentdiff_threshold'])
				else:
					walk_tss, _ = rev_walk(gene, cov_per_contig, scov_per_contig, seq_per_contig, g['end'], params['mincov'], params['min_exon_size'], g['hard_cutoff'], params['tolerated_gap'], params['splicesites'], params['percentdiff_threshold'])
				out['walks'][gene] = walk_tss
				out['expression'][gene] = gene_exp_status  # 'high' or 'moderate' (low / None never reach the walk)
			except Exception as e:  # one failing (sample, gene) must not stop the whole run; it is reported, not hidden
				out['errors'][gene] = repr(e)
	return out


def _worker(item):
	sample, cov_file, scov_file = item
	c = _CTX
	return process_sample(sample, cov_file, scov_file, c['gene_table'], c['shared'], c['params'], c['intergenic_windows'], c['verbose'])

# statistics of one gene's boundary positions
def summarize_boundaries(positions, orientation, tiers=None, cluster_gap=10, tolerance=10,
						 cluster_fraction=0.1, sharpness=50):
	"""Simple summary of the walk positions of one gene.

	centre       median of the samples inside the narrowest window that holds half of them
	w50          width (bp) of that window = spread of the dominant group
	clusters     groups separated by gaps > cluster_gap; only groups of at least
	             cluster_fraction * n samples are kept  (--cluster_fraction of the main script)
	tiers        for every pass-1 tier point: upstream_offset (bp, strand-aware: positive = the
	             boundaries lie further UPSTREAM than the tier, negative = downstream of it) and the
	             fraction of samples within +-tolerance
	promoter     provisional label: 'multi-site' (>=2 clusters), 'focused' (w50 <= sharpness, --sharpness of the main script), else 'broad'
	"""
	if not positions:
		return None
	s = sorted(positions)
	n = len(s)
	k = math.ceil(n / 2)
	best_i = min(range(n - k + 1), key=lambda i: s[i + k - 1] - s[i])
	window = s[best_i:best_i + k]
	w50 = window[-1] - window[0]
	centre = statistics.median(window)

	groups = [[s[0]]]
	for v in s[1:]:
		if v - groups[-1][-1] > cluster_gap:
			groups.append([v])
		else:
			groups[-1].append(v)
	min_size = math.ceil(cluster_fraction * n)
	clusters = [{'centre': statistics.median(c), 'size': len(c), 'fraction': len(c) / n, 'low': c[0], 'high': c[-1]}
				for c in groups if len(c) >= min_size]
	clusters.sort(key=lambda c: -c['size'])

	if len(clusters) >= 2:
		promoter = 'multi-site'
	elif w50 <= sharpness:
		promoter = 'focused'
	else:
		promoter = 'broad'

	tier_stats = {}
	for name, value in (tiers or {}).items():
		if value is None:
			continue
		offset = (value - centre) if orientation == '+' else (centre - value)
		within = sum(1 for p in s if abs(p - value) <= tolerance) / n
		tier_stats[name] = {'value': value, 'upstream_offset': offset, 'fraction_within': within}
	nearest = min(tier_stats, key=lambda t: abs(tier_stats[t]['upstream_offset'])) if tier_stats else None

	return {
		'n': n, 'median': statistics.median(s), 'centre': centre, 'w50': w50, 'tolerance': tolerance,
		'clusters': clusters, 'promoter': promoter,
		'tiers': tier_stats, 'nearest_tier': nearest,
		'nearest_within_tolerance': (abs(tier_stats[nearest]['upstream_offset']) <= tolerance) if nearest else None,
	}

# which tier point do the samples support most
TIER_PREFERENCE = ('accelerated', 'elevated', 'basal')  # only used when two tiers are at exactly the same position

def pick_dominant_tier(stats, candidates=('basal', 'elevated', 'accelerated'), min_fraction=0.1):
	"""Ranks the tier points of one gene by how many samples agree with them.

	stats        result of summarize_boundaries(...) called WITH the tier points of the gene
	candidates   the tiers that can become dominant (the other tool is only shown, never chosen; the
	             cumulative walk TSS needs no entry of its own, it is the basal or the elevated tier)
	min_fraction a tier is dominant only if at least this fraction of the samples (--cluster_fraction)
	             have their walk TSS within +-tolerance of it

	Ranking: the larger fraction of samples within +-tolerance first; if equal, the tier closer to the
	sample centre; if still equal (same position), accelerated before elevated before basal.
	Tiers that lie within +-tolerance of the dominant tier share its samples, so the order among them is
	not evidence: they are listed in `tied_with`.

	Returns None when there are no statistics (low confidence gene), otherwise a dictionary with
	  dominant_tier / dominant_tss   (None, None when no tier reaches min_fraction)
	  fraction_within                of the dominant tier (or of the best tier when none is dominant)
	  supported                      True / False
	  ranking                        [{tier, value, fraction_within, upstream_offset}, ...] best first
	  tied_with                      other candidate tiers within +-tolerance of the dominant tier
	"""
	if not stats:
		return None
	ranking = []
	for tier in candidates:
		info = stats['tiers'].get(tier)
		if info is None:
			continue
		ranking.append({'tier': tier, 'value': info['value'], 'fraction_within': info['fraction_within'],
						'upstream_offset': info['upstream_offset']})
	ranking.sort(key=lambda r: (-r['fraction_within'], abs(r['upstream_offset']),
								TIER_PREFERENCE.index(r['tier']) if r['tier'] in TIER_PREFERENCE else len(TIER_PREFERENCE)))
	result = {'dominant_tier': None, 'dominant_tss': None, 'fraction_within': None, 'supported': False,
			  'ranking': ranking, 'tied_with': []}
	if not ranking:
		return result
	best = ranking[0]
	result['fraction_within'] = best['fraction_within']
	if best['fraction_within'] >= min_fraction:
		result['dominant_tier'] = best['tier']
		result['dominant_tss'] = best['value']
		result['supported'] = True
		result['tied_with'] = [r['tier'] for r in ranking[1:] if abs(r['value'] - best['value']) <= stats['tolerance']]
	return result

# boundary convergence plot of one gene
TIER_COLOURS = {'basal': 'brown', 'elevated': 'blue', 'accelerated': 'pink'}  # same colours as the coverage plot
OTHER_TOOL_COLOUR = 'red'
TIER_LABELS = {'basal': 'Edge point', 'elevated': 'Elevation point', 'accelerated': 'Acceleration point'}  # names used in the coverage plot
EXPRESSION_COLOURS = {'high': '#2ca02c', 'moderate': '#ff8c00'}  # green / dark orange, as in the earlier batch script
EXPRESSION_NAMES = {'high': 'High expression', 'moderate': 'Moderate expression'}

def plot_boundary_convergence(gene, orientation, bc_gene, stats, ranking, out_file, dpi=300):
	"""One dot per sample at its walk TSS (samples at the same position are stacked upwards), green for
	samples where the gene is highly expressed and dark orange where it is moderately expressed.
	The tier points are vertical lines; the dominant tier is the thick solid one. The legend gives, for
	every tier, the fraction of samples within +-tolerance of it and its distance from the sample centre."""
	import matplotlib.pyplot as plt
	from matplotlib.lines import Line2D
	from matplotlib.ticker import MaxNLocator

	positions = bc_gene['by_sample']
	expression = bc_gene['by_sample_expression']
	tolerance = stats['tolerance']
	centre = stats['centre']

	stack = {}
	dots = {'high': ([], []), 'moderate': ([], [])}
	for sample, value in sorted(positions.items(), key=lambda item: item[1]):
		key = round(value)
		stack[key] = stack.get(key, 0) + 1
		category = expression[sample]
		dots[category][0].append(value)
		dots[category][1].append(stack[key])
	max_height = max(stack.values())

	fig, ax = plt.subplots(figsize=(10, 6.5))
	for category in ('high', 'moderate'):
		xs, ys = dots[category]
		if xs:
			ax.scatter(xs, ys, s=28, color=EXPRESSION_COLOURS[category], alpha=0.85, edgecolors='none', zorder=3,
					   label=f'{EXPRESSION_NAMES[category]} (n={len(xs)})')

	dominant_tier = ranking['dominant_tier'] if ranking else None
	tier_handles = []
	ax.axvline(centre, color='grey', linestyle='-.', linewidth=1, zorder=1)
	tier_handles.append(Line2D([0], [0], color='grey', linestyle='-.', linewidth=1,
							   label=f'sample centre {centre:.12g}  (n={stats["n"]}, w50={stats["w50"]:.12g} bp)'))
	for tier, info in stats['tiers'].items():
		colour = TIER_COLOURS.get(tier, OTHER_TOOL_COLOUR)
		is_dominant = tier == dominant_tier
		ax.axvline(info['value'], color=colour, linestyle='-' if is_dominant else ':',
				   linewidth=2.6 if is_dominant else 1.4, zorder=2)
		offset = info['upstream_offset']  # positive: the samples lie further upstream than the tier
		if offset > 0:
			where = f'{offset:.12g} bp downstream of centre'
		elif offset < 0:
			where = f'{-offset:.12g} bp upstream of centre'
		else:
			where = 'at centre'
		label = f'{TIER_LABELS.get(tier, tier)} {info["value"]:.12g}: {info["fraction_within"]:.2f} of samples within \u00b1{tolerance:g} bp, {where}'
		if is_dominant:
			label += '   <- dominant'
		tier_handles.append(Line2D([0], [0], color=colour, linestyle='-' if is_dominant else ':',
								   linewidth=2.6 if is_dominant else 1.4, label=label))

	ax.set_ylim(0, (max_height + 1) / 0.5)  # empty space above the dots for the legend
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	ax.ticklabel_format(axis='x', useOffset=False, style='plain')
	every_x = list(positions.values()) + [info['value'] for info in stats['tiers'].values()]
	pad = max(5, 0.05 * (max(every_x) - min(every_x)))
	ax.set_xlim(min(every_x) - pad, max(every_x) + pad)
	ax.set_xlabel('Walk TSS' + ('   (minus strand gene: upstream is to the right)' if orientation == '-' else ''))
	ax.set_ylabel('Number of samples (stacked)')
	if ranking and ranking['supported']:
		verdict = f'dominant point: {TIER_LABELS.get(dominant_tier, dominant_tier)} ({ranking["fraction_within"]:.2f} of samples within \u00b1{tolerance:g} bp)'
		if ranking['tied_with']:
			verdict += f'; within \u00b1{tolerance:g} bp of it: {", ".join(TIER_LABELS.get(t, t) for t in ranking["tied_with"])}'
	elif ranking and ranking['ranking']:
		verdict = f'no point supported by the samples (best: {TIER_LABELS.get(ranking["ranking"][0]["tier"], ranking["ranking"][0]["tier"])}, {ranking["fraction_within"]:.2f})'
	else:
		verdict = 'no point to compare with'
	ax.set_title(f'{gene} ({orientation}): {stats["promoter"]}\n{verdict}', fontsize=11)

	first_legend = ax.legend(handles=tier_handles, loc='upper left', fontsize=8, frameon=True, framealpha=0.9)
	ax.add_artist(first_legend)
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
	fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
	plt.close(fig)


# entry point
def run_boundary_convergence(sample_cov_files, sample_scov_files, goi, shared, params, helpers, out_folder,
							 min_samples=50, min_fraction=1/3, n_workers=1, verbose=False,
							 tier_points=None, **summary_kwargs):
	"""
	sample_cov_files / sample_scov_files : lists of per-sample aligned / spanning cov files
	                                       (see list_sample_cov_files; paired as described in pair_sample_files)
	goi                                  : gene ids, as in the main script
	shared  : dict with gene_infos, genes_per_chromosome, protein_coding_genes_per_chromosome,
	          transcripts_per_gene, five_utr_infos, cds_infos, genome_seq  (all already loaded by main)
	params  : dict with protein_encoding, coverage_walk_origin, window, mincov, min_exon_size,
	          tolerated_gap, splicesites, percentdiff_threshold, intergenic_buffer,
	          intergenic_region_size, background_percentage
	helpers : dict of the main script's own functions: find_flanking_genes and get_overlap_type
	          (the coverage loading, intergenic background and expression level are inside this module)
	tier_points : optional {gene: {'basal': .., 'elevated': .., 'accelerated': .., 'EPD': ..}} from pass 1.
	          Can be left out and summarize_boundaries() called later, when pass 1 is finished.
	n_workers : samples processed in parallel. Every worker holds the genome-wide cov/scov of one
	          sample in memory (as the main script does for the merged files), so keep this low.

	Returns {'status', 'n_samples', 'genes': {gene: {...}}, 'skipped_genes', 'errors'}.
	"""
	pairs = pair_sample_files(sample_cov_files, sample_scov_files)
	n_samples = len(pairs)
	if n_samples < min_samples:
		msg = f'{n_samples} individual samples found; boundary convergence needs at least {min_samples}. Not run.'
		print(msg)
		return {'status': 'skipped_too_few_samples', 'message': msg, 'n_samples': n_samples,
				'genes': {}, 'skipped_genes': {}, 'errors': {}}

	os.makedirs(out_folder, exist_ok=True)
	gene_table, skipped_genes = build_gene_table(goi, shared, params, helpers)
	intergenic_windows = get_intergenic_windows(shared['genes_per_chromosome'], shared['gene_infos'], params['intergenic_buffer'], params['intergenic_region_size'])

	_CTX.clear()
	_CTX.update({'gene_table': gene_table, 'shared': shared, 'params': params,
				 'intergenic_windows': intergenic_windows, 'verbose': verbose})
	results = []
	if n_workers > 1:
		with multiprocessing.get_context('fork').Pool(n_workers) as pool:
			for i, r in enumerate(pool.imap_unordered(_worker, pairs), start=1):
				results.append(r)
				print(f'[boundary convergence] {r["sample"]} done ({i}/{n_samples})')
	else:
		for i, item in enumerate(pairs, start=1):
			r = _worker(item)
			results.append(r)
			print(f'[boundary convergence] {r["sample"]} done ({i}/{n_samples})')
	results.sort(key=lambda r: r['sample'])

	genes_out = {}
	errors = {}
	rows = []
	for gene, g in gene_table.items():
		by_sample = {}
		by_sample_expression = {}
		for r in results:
			if gene in r['walks']:
				walk_tss = r['walks'][gene]
				by_sample[r['sample']] = walk_tss
				by_sample_expression[r['sample']] = r['expression'][gene]
				rows.append((gene, r['sample'], walk_tss, r['expression'][gene]))
			if gene in r['errors']:
				errors[(r['sample'], gene)] = r['errors'][gene]
		n_valid = len(by_sample)
		fraction = n_valid / n_samples
		low_conf = fraction < min_fraction
		stats = None
		if not low_conf:
			stats = summarize_boundaries(list(by_sample.values()), g['orientation'], (tier_points or {}).get(gene),
										 **summary_kwargs)
		genes_out[gene] = {
			'orientation': g['orientation'], 'hard_cutoff': g['hard_cutoff'],
			'n_samples': n_samples, 'n_valid': n_valid, 'fraction_valid': fraction,
			'low_confidence': low_conf,
			'positions': sorted(by_sample.values()), 'by_sample': by_sample,
			'by_sample_expression': by_sample_expression,
			'statistics': stats,
		}

	with open(os.path.join(out_folder, 'boundary_walk_positions.tsv'), 'w') as f:
		f.write('gene\tsample\twalk_tss\texpression_level\n')
		for row in rows:
			f.write('\t'.join(str(x) for x in row) + '\n')

	if errors:
		print(f'[boundary convergence] {len(errors)} (sample, gene) analyses failed and were left out; see result["errors"]')
	return {'status': 'ok', 'n_samples': n_samples, 'genes': genes_out,
			'skipped_genes': skipped_genes, 'errors': errors}
