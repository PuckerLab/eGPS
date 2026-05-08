### Boas Pucker ###
### Shakunthala Natarajan ###
### pucker@uni-bonn.de ###
__version__ = "v0.1"

__reference__ = "https://github.com/bpucker/TSS_finder"

__usage__ = """
					TSS_finder """ + __version__ + """("""+ __reference__ +""")
					
					Usage:
					python3 TSS_finder.py
					--fasta <GENOMIC_FASTA_FILE>
					--gff <GFF_FILE>
					--goi <TXT_FILE_WITH_LIST_OF_GENES_OF_INTEREST_ONE_PER_LINE>
					--out <OUTPUT_FOLDER>
					[--sra_folder <READ_FILES> | --bam <BAM_FILE>|--cov <COV_FILE> --scov <SCOV_FILE>]
					
					optional:
					--coverage_walk_origin <cds or utr> DEFAULT is cds
					--mincov <MINIMAL_COVERAGE>[1]
					--bam_is_sorted <PREVENTS_BAM_FILE_SORTING>
					--samtools <FULL_PATH_TO_SAMTOOLS>[samtools]
					--bedtools <FULL_PATH_TO_genomeCoverageBed>[genomeCoverageBed]
					--m <MEM_FOR_SAMTOOLS_SORTING>[5000000000]
					--threads <NUMBER_THREADS_FOR_SAMTOOLS_SORTING>[4]
					--minexon <MINIMAL_EXON_SIZE>[10]
					--flanksize <FLANKING_REGION_SIZE>[50]
					--gapsize <COVERAGE_GAP_SIZE>[5]
					--splicesites <HANDLING_OF_SPLICE_SITES>[strict](strict|off)
					--gff_config <TXT_CONFIG_FILE_WITH_OPTIONS_TO_DEFINE_THE_GFF_ATTRIBUTE_FIELDS>
					--intron_percentile_cutoff <INTRON_SIZE_CUTOFF_FOR_STAR+HISAT2_ALIGNMENT>
					--neighbourhood <GENE_NEIGHBOURHOOD_WINDOW_FOR_OVERLAPPING_GENE_ANALYSIS>
					--min_promoter_size <MINIMUM_PROMOTER_SIZE>
					--max_promoter_size <MAXIMUM_PROMOTER_SIZE>
					--background < NUMBER_OF_RANDOM_BACKGROUND_SEQUENCES_TO_BE_CONSIDERED_FOR_MOTIF_SCORING>
					--downstream_size <LENGTH_OF_ANALYSIS_REGION_DOWNSTREAM_TO_TSS>
					--upstream_slice <LENGTH_OF_PROMOTER_REGION_TO_BE_SLICED_FOR_MOTIF_ANALYSIS>
					--downstream_slice <LENGTH_OF_DOWNSTREAM_REGION_TO_BE_SLICED_FOR_MOTIF_ANALYSIS>
					--aligner <STAR or HISAT2> STAR is default
					--HISAT2 <FULL_PATH_TO_HISAT2_FOR_RNA_Seq_MAPPING>
					--STAR <FULL_PATH_TO_STAR_FOR_RNA_Seq_MAPPING>
					--index_bases <PARAMETER_FOR_GENOME_INDEX_GENERATION_IN_STAR>
					--run_mode < find_tss or make_bam> find_tss is default
					--fastq_pattern <FASTQ_FILE_NAME_PATTERN_SEPARATED_BY_COMMA><EG: pass1,pass2 >
					--analyse_promoter <yes or no> default is no
					--moods <FULL_PATH_TO_MOODS_SCRIPT>
					--PFM <TXT_FILE_FOR_PROMOTER_MOTIF_ANALYSIS; TAB_SEPARATED; MOTIF_NAME	PATH_TO_PFM_FILE	UPSTREAM_BOUNDARY	DOWNSTREAM_BOUNDARY	DIRECTION_SENSITIVITY>
					--pval <MOODS_MOTIF_ANALYSIS_THRESHOLD>
					"""


import re, os, sys, subprocess, gzip
import random
import math
import tempfile
import traceback
import numpy as np
from scipy.stats import percentileofscore
from decimal import Decimal, ROUND_HALF_DOWN
from collections import defaultdict
try:
	import matplotlib.pyplot as plt
	from matplotlib.lines import Line2D
	import matplotlib.ticker as ticker
	from matplotlib.patches import FancyArrow
except ImportError:
	pass

# --- end of imports --- #

def construct_cov_file( bam_file, cov_file, bedtools ):
	""" @brief calculate read coverage depth per position """
	
	print ( "calculating coverage per position ...." )
	cmd = bedtools + " -d -split -ibam " + bam_file + " > " + cov_file	#-split ignored spanning reads when calculating depth
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()
	return cov_file


def construct_scov_file( bam_file, scov_file, bedtools ):
	""" @brief calculate read coverage depth per position """
	
	print ( "calculating coverage per position (spanning) ...." )
	cmd = bedtools + " -d -ibam " + bam_file + " > " + scov_file	#-include spanning reads when calculating depth
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()
	return scov_file


def load_coverage( cov_file, input_mode ):
	"""! @brief load coverage per genomic position """
	coverage_per_seq = {}
	if not cov_file.endswith("gz"):#uncompressed file
		with open( cov_file, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					coverage_per_seq[ parts[0] ].append( float( parts[2] ) )
				except KeyError:
					coverage_per_seq.update( { parts[0]: [ float( parts[2] ) ] } )
				line = f.readline()
	else:#compressed file
		with gzip.open( cov_file, "rt" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				try:
					coverage_per_seq[ parts[0] ].append( float( parts[2] ) )
				except KeyError:
					coverage_per_seq.update( { parts[0]: [ float( parts[2] ) ] } )
				line = f.readline()
	return coverage_per_seq


def load_gene_infos( gff_file, child_attribute, child_parent_linker, parent_attribute):
	"""! @brief load gene ID, position, and orientation from GFF3 file """
	
	gene_infos = {}
	mrna_infos = {}
	five_utr_infos={}
	cds_infos = {}
	genes_per_chromosome = {}
	transcripts_per_gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[2].upper() == "GENE" or parts[2].upper() == "TRANSPOSABLE_ELEMENT_GENE" or parts[2].upper() == "PSEUDOGENE" or parts[2].upper() == "PSEUDO_GENE":	#could be extended to other feature types
					ID = parts[-1].split(f'{parent_attribute}=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					gene_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						genes_per_chromosome[ parts[0] ].append( ID )
					except KeyError:
						genes_per_chromosome.update( { parts[0]: [ ID ] } )
				if parts[2].upper() == "MRNA":
					ID = parts[-1].split(f'{child_attribute}=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					Parent = parts[-1].split(f'{child_parent_linker}=')[-1]
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					mrna_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						transcripts_per_gene[ Parent ].append( ID )
					except KeyError:
						transcripts_per_gene.update( { Parent: [ ID ] } )
				if parts[2].upper() == 'FIVE_PRIME_UTR':
					Parent = parts[-1].split(f'{child_parent_linker}=')[-1]#Parent of 5'UTR is transcript
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					five_utr_infos.update({ Parent: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } })# key of this nested dictionary is the transcript name
				if parts[2].upper() == 'CDS':
					cds_parents = parts[-1].split(f'{child_parent_linker}=')[-1]
					if ";" in cds_parents:
						cds_parents = cds_parents.split(';')[0]
					for cds_parent in cds_parents.split(','):  # handle multiple parents
						cds_parent = cds_parent.strip()
						cds_tuple = (int(parts[3]), int(parts[4]))
						try:
							cds_infos[cds_parent].append(cds_tuple)
						except KeyError:
							cds_infos[cds_parent] = [cds_tuple]
			line = f.readline()
	for chromosome in genes_per_chromosome: #sort the genes in each contig/ chromosome in the ascending order of start positions
		genes_per_chromosome[chromosome].sort(key=lambda gene: (gene_infos[gene]['start'], gene))
	for gene in transcripts_per_gene:#sort the transcripts per gene in the ascending order of mRNA start positions for + strand and mRNA end positions for - strand
		if gene_infos[gene]['orientation']=='+':
			transcripts_per_gene[gene].sort(key=lambda transcript: (mrna_infos[transcript]['start'], transcript))
		elif gene_infos[gene]['orientation']=='-':
			transcripts_per_gene[gene].sort(key=lambda transcript: (mrna_infos[transcript]['end'], transcript))
	gene_atg_dic = {}
	for gene in transcripts_per_gene:
		if gene not in gene_infos:
			continue
		orientation = gene_infos[gene]['orientation']
		transcript_list = transcripts_per_gene[gene]

		# select most upstream transcript for +, most downstream for -
		if orientation == '+':
			selected_transcript = transcript_list[0]  # already sorted by start ascending
		else:
			selected_transcript = transcript_list[-1]  # already sorted by end ascending, last = highest end

		if selected_transcript not in cds_infos:
			continue  # no CDS annotated for this transcript

		cds_list = cds_infos[selected_transcript]

		if orientation == '+':
			cds_list.sort(key=lambda x: x[0])  # sort by start ascending
			gene_atg_dic[gene] = cds_list[0][0]  # start of most upstream CDS = ATG
		else:
			cds_list.sort(key=lambda x: x[1])  # sort by end ascending
			gene_atg_dic[gene] = cds_list[-1][1]  # end of most downstream CDS = ATG
	return gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	seq_counter = 1
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				seq_counter+=1
				sequences.update( { header: "".join( seq ) } )
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
	return sequences, seq_counter

def get_gene_orientation_for_transcript(transcripts_per_gene,gene_infos, genome_seq):
	transcript_orientation_dic = {}
	transcript_contig_dic = {}
	for gene, transcript_list in transcripts_per_gene.items():
		transcript_orientation = gene_infos[gene]['orientation']
		transcript_containing_contig = genome_seq[ gene_infos[ gene ]['chromosome'] ]
		for transcript_id in transcript_list:
			transcript_orientation_dic[transcript_id] = transcript_orientation
			transcript_contig_dic[transcript_id] = transcript_containing_contig
	return transcript_orientation_dic, transcript_contig_dic

def get_background_seqs (promoter_slice, downstream_slice, cds_infos, fraclen, max_promoter_size, gff_file, downstream_size,genome_seq, gene_infos, transcripts_per_gene, child_parent_linker, background_strength):
	seen_introns = set()
	seen_cds = set()
	background_positive_strand_seqs =[]
	background_negative_strand_seqs = []
	background_positive_strand_cds_seqs = []
	background_negative_strand_cds_seqs = []
	trans_exon_map = {}
	transcript_orientation_dic, transcript_contig_dic, transcript_index_dic = get_gene_orientation_for_transcript(transcripts_per_gene, gene_infos, genome_seq)
	with open(gff_file, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				parts = line.strip().split('\t')
				if parts[2].lower() == 'exon':
					fields = parts[8].strip().split(';')
					for each in fields:
						if f'{child_parent_linker}=' in each:
							parents = each.replace(f'{child_parent_linker}=', '').split(',')
							for transcript in parents:
								if parents.index(transcript) == 0:
									transcript = transcript.strip()
									start = int(parts[3])
									end = int(parts[4])
									if transcript not in trans_exon_map:
										trans_exon_map[transcript] = [(start, end)]
									else:
										trans_exon_map[transcript].append((start, end))

	max_bg_len = promoter_slice + downstream_slice
	cds_coords = []
	# calculate intron sizes
	for transcript in trans_exon_map.keys():
		if transcript in transcript_orientation_dic:
			orientation = transcript_orientation_dic[transcript]
			seq = transcript_contig_dic[transcript]
			exons = trans_exon_map[transcript]
			if len(exons) < 2:
				continue
			counter = 0
			if transcript in cds_infos:
				cds_coords= cds_infos[transcript]
			if orientation == '+':
				# sort by ascending order of start coordinates
				exons.sort(key=lambda x: x[0])
				while counter <= (len(exons) - 2):
					intron_start = exons[counter][1]  # end of current exon
					intron_end = exons[counter + 1][0]  # start of next exon
					intron_key = (intron_start, intron_end)
					if intron_key in seen_introns:
						counter += 1
						continue
					seen_introns.add(intron_key)
					intron_seq = seq[intron_start:intron_end]
					if len(intron_seq) <= 0:
						pass
					else:
						if len(intron_seq) >= max_bg_len:
							# intron alone is long enough — truncate and use directly
							background_seq = intron_seq[:max_bg_len]
							if len(background_positive_strand_seqs) < background_strength:
									background_positive_strand_seqs.append(background_seq)
					counter += 1
				#collecting cds background seqs
				if transcript_index_dic[transcript] != 0:  # ensuring that the cds of the most upstream transcript is not taken for the background since it might have the TSS
					if cds_coords:
						cds_coords.sort(key=lambda x: x[0]) # sort by ascending order of start coordinates
						for each in cds_coords:
							if cds_coords.index(each) == 0: #avoid taking the most upstream CDS
								pass
							else:
								cds_start = each[0]
								cds_end = each[1]
								cds_seq = seq[cds_start:cds_end]
								if len(cds_seq) >= max_bg_len:
									background_seq = cds_seq[:max_bg_len]
									if len(background_positive_strand_cds_seqs) < background_strength:
										background_positive_strand_cds_seqs.append(background_seq)

			elif orientation == '-':
				exons.sort(key=lambda x: x[1], reverse=True)  # sort by end descending
				while counter <= (len(exons) - 2):
					intron_start = exons[counter + 1][1]   # end of downstream exon
					intron_end = exons[counter][0]          # start of upstream exon
					intron_key = (intron_start, intron_end)
					if intron_key in seen_introns:
						counter += 1
						continue
					seen_introns.add(intron_key)
					intron_seq = seq[intron_start:intron_end]
					if len(intron_seq) <= 0:
						pass
					else:
						if len(intron_seq) >= max_bg_len:
							# intron alone is long enough — truncate and use directly
							background_seq = intron_seq[-max_bg_len:]
							if len(background_negative_strand_seqs) < background_strength:
								background_negative_strand_seqs.append(background_seq)
					counter += 1
				# collecting cds background seqs
				if transcript_index_dic[transcript] != -1:  # ensuring that the cds of the most upstream transcript is not taken for the background since it might have the TSS
					if cds_coords:
						cds_coords.sort(key=lambda x: x[0])  # sort by ascending order of start coordinates
						for each in cds_coords:
							if cds_coords.index(each) == -1:  # avoid taking the most upstream CDS
								pass
							else:
								cds_start = each[0]
								cds_end = each[1]
								cds_seq = seq[cds_start:cds_end]
								if len(cds_seq) >= max_bg_len:
									background_seq = cds_seq[:max_bg_len]
									if len(background_negative_strand_cds_seqs) < background_strength:
										background_negative_strand_cds_seqs.append(background_seq)

	return background_positive_strand_seqs, background_positive_strand_cds_seqs, background_negative_strand_seqs, background_negative_strand_cds_seqs

#function to randomly subsample a fragment of defined length from a given contig
def get_random_fragment(contig, seq_len):
	start = random.randint(0,(len(contig) - seq_len))
	random_seq = contig[start:start+seq_len]
	return random_seq

#function to randomly keep sampling one random seq from each contig in the genome iteratively until the background strength cap is reached
def get_random_background_seqs(genome_seq, seq_len, background_strength):
	background = []
	while len(background) < background_strength:
		for contig in genome_seq:
			if len(genome_seq[contig]) > seq_len:
				if len(background) >= background_strength:
					break
				random_seq = get_random_fragment(genome_seq[contig], seq_len)
				background.append(random_seq)
	return background

def generate_plot( values, svalues, fig_file, atg_pos, cov_walk_start, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot ):
	"""! @brief generate a coverage plot """
	fig, (ax1, ax_features) = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios':[4,1]}, sharex=True)
	ax1.plot(values, color="black", linestyle="solid")  # coverage of aligned bases
	ax2 = ax1.twinx()
	ax2.plot(svalues, color="red", linestyle="dotted")  # coverage of spanning reads
	ax2.plot([atg_pos, atg_pos], [0, max(svalues + values)], color="green", linestyle="dotted", label="ATG")  # ATG position
	ax2.plot([cov_walk_start, cov_walk_start], [0, max(svalues + values)], color="orange", linestyle="dotted", label="Coverage walk origin")  # 5'UTR start or rnd or gene start or end position depending on strandedness and 5'UTR annotation being present for the gene's most upstream or downstream transcripts
	ax2.plot([tss_pos, tss_pos], [0, max(svalues + values)], color="blue", linestyle="dotted", label="TSS")  # TSS position
	ax2.legend(loc='best')
	"""
	# feature track plotting with features represented as rectangles
	for mrna_start, mrna_end in mrnas_to_plot:
		ax_features.broken_barh([(mrna_start, mrna_end - mrna_start)], (0.1, 0.2), facecolors='steelblue')
	for cds_start, cds_end in cds_to_plot:
		ax_features.broken_barh([(cds_start, cds_end - cds_start)], (0.4, 0.2), facecolors='lightgreen')
	for five_utr_start, five_utr_end in mrnas_to_plot:
		ax_features.broken_barh([(five_utr_start, five_utr_end - five_utr_start)], (0.7, 0.2), facecolors='salmon')
	"""
	# feature track plotting with features represented as arrows
	for mrna_start, mrna_end in mrnas_to_plot:
		feature_length = mrna_end - mrna_start
		if orientation == '+':
			dx = feature_length
			x_origin = mrna_start
			xdot_origin = tss_pos
			predicted_tss_mrna_feature = mrna_end - tss_pos
		elif orientation == '-':
			dx = -feature_length
			x_origin = mrna_end
			xdot_origin = tss_pos
			predicted_tss_mrna_feature = -(tss_pos - mrna_start)
		head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
		arrow = FancyArrow(x=x_origin, y=0.2, dx=dx, dy=0, width=0.2, head_width=0.3, head_length=20, length_includes_head=True, facecolor='steelblue', alpha=0.5, edgecolor='black', linewidth=0.5)
		ax_features.add_patch(arrow)
		dotted_arrow = FancyArrow(x=xdot_origin, y=0.2, dx=predicted_tss_mrna_feature, dy=0, width=0.2, head_width=0.3, head_length=20, length_includes_head=True, facecolor='none', alpha=0.5, edgecolor='steelblue', linewidth=0.5, linestyle='dotted')
		ax_features.add_patch(dotted_arrow)
	for cds_start, cds_end in cds_to_plot:
		feature_length = cds_end - cds_start
		if orientation == '+':
			dx = feature_length
			x_origin = cds_start
		elif orientation == '-':
			dx = -feature_length
			x_origin = cds_end
		head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
		arrow = FancyArrow(x=x_origin, y=0.5, dx=dx, dy=0, width=0.2, head_width=0.3, head_length=20, length_includes_head=True, facecolor='lightgreen', alpha=0.5, edgecolor='black', linewidth=0.5)
		ax_features.add_patch(arrow)
	for five_utr_start, five_utr_end in five_utr_to_plot:
		feature_length = five_utr_end - five_utr_start
		if orientation == '+':
			dx = feature_length
			x_origin = five_utr_start
			xdot_origin = tss_pos
			predicted_tss_five_utr_feature = five_utr_end - tss_pos
		elif orientation == '-':
			dx = -feature_length
			x_origin = five_utr_end
			xdot_origin = tss_pos
			predicted_tss_five_utr_feature = -(tss_pos - five_utr_start)
		head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
		arrow = FancyArrow(x=x_origin, y=0.8, dx=dx, dy=0, width=0.2, head_width=0.3, head_length=20, length_includes_head=True, facecolor='salmon', alpha=0.5, edgecolor='black', linewidth=0.5)
		ax_features.add_patch(arrow)
		dotted_arrow = FancyArrow(x=xdot_origin, y=0.8, dx=predicted_tss_five_utr_feature, dy=0, width=0.2, head_width=0.3,head_length=20, length_includes_head=True, facecolor='none', alpha=0.5,edgecolor='salmon', linewidth=0.5, linestyle='dotted')
		ax_features.add_patch(dotted_arrow)

	# extend vertical lines into feature axis
	ax_features.axvline(atg_pos, color="green", linestyle="dotted")
	ax_features.axvline(cov_walk_start, color="orange", linestyle="dotted")
	ax_features.axvline(tss_pos, color="blue", linestyle="dotted")

	ax1.set_title(gene + "   (" + orientation + ")")
	ax_features.set_yticks([0.2, 0.5, 0.8])
	ax_features.set_yticklabels(['mRNA', 'CDS', "5'UTR"], fontsize=8)
	ax1.set_xlabel("Position in genomic region from " + str(genomic_start) + " to " + str(genomic_end))
	ax1.set_ylabel("Aligned RNA-seq coverage")
	ax1.yaxis.label.set_color('black')
	ax2.set_ylabel("Spanning RNA-seq coverage")
	ax2.yaxis.label.set_color('red')
	"""
	#adaptive figure sizing based on coverage intervals and genomic range
	genomic_range = genomic_end - genomic_start
	fig_width = max(10, genomic_range / 100)  # 1 inch per 100 bp, minimum 10 inches

	y_max = max(max(svalues + values), 1)
	fig_height = max(5, y_max / 500)  # 1 inch per 500 coverage units, minimum 5 inches
	
	fig, ax1 = plt.subplots(figsize=(fig_width, fig_height))
	ax1.plot( values, color="black", linestyle="solid" )	#coverage of aligned bases
	ax2 = ax1.twinx()
	ax2.plot( svalues, color="red", linestyle="dotted" )	#coverage of spanning reads
	ax2.plot( [ atg_pos, atg_pos ], [ 0, max( svalues+values ) ], color="green", linestyle="dotted", label="ATG")	#ATG position
	ax2.plot([tss_pos, tss_pos], [0, max(svalues + values)], color="blue", linestyle="dotted", label="TSS")  # TSS position
	ax2.legend()
	
	#replacing the above commented code block with the code block below to overlay with and display the nucleotide sequence on the plot
	genomic_range = genomic_end - genomic_start
	fig_width = max(10, genomic_range / 100)
	y_max = max(max(svalues + values), 1)
	fig_height = max(5, y_max / 500)

	# base colors for sequence track
	base_colors = {'A': '#2ecc71', 'T': '#e74c3c', 'G': '#3498db', 'C': '#f39c12',
				   'a': '#2ecc71', 't': '#e74c3c', 'g': '#3498db', 'c': '#f39c12',
				   'N': '#cccccc', 'n': '#cccccc'}

	# two-panel layout: coverage on top (85%), sequence track on bottom (15%)
	fig = plt.figure(figsize=(fig_width, fig_height + 1))
	gs = fig.add_gridspec(2, 1, height_ratios=[0.85, 0.15], hspace=0.05)

	ax1 = fig.add_subplot(gs[0])
	ax_seq = fig.add_subplot(gs[1], sharex=ax1)  # shares x-axis with coverage plot

	ax2 = ax1.twinx()
	ax1.plot(values, color="black", linestyle="solid")
	ax2.plot(svalues, color="red", linestyle="dotted")

	y_max = max(max(svalues + values), 1)
	ax2.plot([atg_pos, atg_pos], [0, y_max], color="green", linestyle="dotted", label="ATG")
	ax2.plot([tss_pos, tss_pos], [0, y_max], color="blue", linestyle="dotted", label="TSS")
	ax2.legend()

	# --- sequence track ---
	ax_seq.set_xlim(0, len(dna_sequence_for_plot))
	ax_seq.set_ylim(0, 1)
	ax_seq.axis('off')  # hide axes frame and ticks for sequence track

	if genomic_range <= 300:
		# show individual base letters with colored background rectangles
		for i, base in enumerate(dna_sequence_for_plot):
			color = base_colors.get(base, '#cccccc')
			ax_seq.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, alpha=0.6))
			ax_seq.text(i + 0.5, 0.5, base.upper(), ha='center', va='center',
						fontsize=max(4, min(8, fig_width / len(dna_sequence_for_plot) * 10)),
						fontweight='bold', color='black')
	else:
		# show color strip only — no text since bases would be illegible
		for i, base in enumerate(dna_sequence_for_plot):
			color = base_colors.get(base, '#cccccc')
			ax_seq.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, alpha=0.8))
		# add a compact legend for base colors
		for base, color in [('A', '#2ecc71'), ('T', '#e74c3c'), ('G', '#3498db'), ('C', '#f39c12')]:
			ax_seq.plot([], [], color=color, linewidth=6, label=base)
		ax_seq.legend(loc='upper right', ncol=4, fontsize=7,
					  framealpha=0.7, borderpad=0.3, handlelength=1)

	ax1.set_title( gene + "   (" + orientation + ")" )
	ax_seq.set_xlabel( "Position in genomic region from " + str( genomic_start ) + " to " + str( genomic_end ) )
	ax1.set_ylabel( "Aligned RNA-seq coverage", labelpad=15 )
	ax1.yaxis.label.set_color('black')
	ax2.set_ylabel( "Spanning RNA-seq coverage", labelpad=15 )
	ax2.yaxis.label.set_color('red')
	"""
	fig.savefig( fig_file, dpi=600, bbox_inches='tight' )
	plt.close(fig)

def get_overlap_type(goi_strand, goi_start, goi_end, nbr_strand, nbr_start, nbr_end):
	#define head/tail coordinates based on strand
	if goi_strand == '+':
		goi_head, goi_tail = goi_start, goi_end
	else:# '- strand'
		goi_head, goi_tail = goi_end,   goi_start

	if nbr_strand == '+':
		nbr_head, nbr_tail = nbr_start, nbr_end
	else:
		nbr_head, nbr_tail = nbr_end,   nbr_start

	#find the actual overlap window
	ov_start = max(goi_start, nbr_start)
	ov_end   = min(goi_end,   nbr_end)

	if ov_start > ov_end:
		return 'no_overlap'

	if goi_strand == nbr_strand:
		return 'same_strand'# tandem or nested, promoter may be affected

	#check which functional ends fall inside the overlap window
	goi_head_in_ov = ov_start <= goi_head <= ov_end
	goi_tail_in_ov = ov_start <= goi_tail <= ov_end
	nbr_head_in_ov = ov_start <= nbr_head <= ov_end
	nbr_tail_in_ov = ov_start <= nbr_tail <= ov_end

	if goi_head_in_ov and nbr_head_in_ov:
		return 'head_head'
	if goi_tail_in_ov and nbr_tail_in_ov:
		return 'tail_tail'
	if goi_head_in_ov:
		return 'head_into_neighbor'
	if goi_tail_in_ov:
		return 'tail_into_neighbor'
	return 'nested'



def run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos ):
	"""! @brief run analysis on forward strand """
	most_upstream_pos = start
	final_pos_status = False
	while not final_pos_status:
		
		# --- walk coverage upstream of transcription start while there is coverage --- #
		while cov_per_contig[ most_upstream_pos-2 ] >= mincov:	#index = genomic position -1 but coverage of gene that is next to the current position needs to be assessed before moving in there
			most_upstream_pos -= 1	#move one step upstream
			if most_upstream_pos == hard_cutoff:	#stop if end of upstream contig/pseudochromosome is reached
				break
		
		# --- try to cross intron --- #
		current_position = most_upstream_pos - 1	#most_upstream_pos has coverage above cutoff (position, not index!)
		if current_position > hard_cutoff:
			while cov_per_contig[ current_position - 2 ] < mincov:	#check if upstream position has low coverage
				current_position -= 1	#move one step upstream
				if current_position == hard_cutoff:
					break
		else:
			final_pos_status = True
		
		avg_gap_coverage = sum( scov_per_contig[ current_position-1:most_upstream_pos-1 ] )/(most_upstream_pos-current_position)
		#average coverage in intron should be very low
		if current_position > min_exon_size and avg_gap_coverage > mincov:
			# --- check coverage gaps for (canonical) splice sites to continue across introns --- #
			donor_splice_site = seq_per_contig[current_position-1:current_position+1].upper()	#this should be GT
			acceptor_splice_site = seq_per_contig[most_upstream_pos-3:most_upstream_pos-1].upper()	#this should be AG
			if donor_splice_site == "GT" and acceptor_splice_site == "AG":
				print( "donor splice site: " + donor_splice_site )
				print( "acceptor splice site: " + acceptor_splice_site )
				most_upstream_pos = current_position - 1
			elif donor_splice_site == "GC" and acceptor_splice_site == "AG":
				print( "Non-canonical donor splice site: " + donor_splice_site )
				print( "Non-canonical acceptor splice site: " + acceptor_splice_site )
				most_upstream_pos = current_position - 1
			elif splicesites == "off":	#ignore check for canonical splice sites
				most_upstream_pos = current_position - 1
			elif most_upstream_pos - current_position < tolerated_gap:
				most_upstream_pos = current_position - 1
			else:
				final_pos_status = True
		else:
			final_pos_status = True
	print( "TSS position of " + gene + ": " + str( most_upstream_pos ) )
	
	# --- generate figures to visualize coverage around the TSS for manual inspection --- #
	transcript_list = transcripts_per_gene[gene]
	mrna_dic = {}
	cds_dic={}
	five_utr_dic={}
	for each in transcript_list:
		if each in mrna_infos:
			mrna_dic[each]= (mrna_infos[each]['start'], mrna_infos[each]['end'])
		if each in cds_infos:
			cds_dic[each] = cds_infos[each]
		if each in five_utr_infos:
			five_utr_dic[each] = (five_utr_infos[each]['start'], five_utr_infos[each]['end'])
	if most_upstream_pos > flank_region_for_plot:
		plot_start_region = most_upstream_pos - flank_region_for_plot
	else:
		plot_start_region = 0
	plot_end_region = max(start + flank_region_for_plot, atg_genomic_pos + flank_region_for_plot)
	mrnas_to_plot = []
	cds_to_plot = []
	five_utr_to_plot = []
	for transcript, (mrna_start, mrna_end) in mrna_dic.items():
		if mrna_end >= plot_start_region and mrna_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
			if mrna_start >= plot_start_region and mrna_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (mrna_end - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end <= plot_end_region:#case2 where feature start is out of the plot bounds
				mrnas_to_plot.append((0, (mrna_end - plot_start_region)))
			elif mrna_start >= plot_start_region and mrna_end >= plot_end_region:#case3 where feature end is out of the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region),(plot_end_region - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				mrnas_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, cds_list in cds_dic.items():
		for (cds_start, cds_end) in cds_list:#two level loop for cds_dic alone since cds_dic structure is a list of tuples per transcript similar to the cds_infos structure from which it is derived
			if cds_end >= plot_start_region and cds_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
				if cds_start >= plot_start_region and cds_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (cds_end - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end <= plot_end_region:#case2 where feature start is out of the plot bounds
					cds_to_plot.append((0, (cds_end - plot_start_region)))
				elif cds_start >= plot_start_region and cds_end >= plot_end_region:#case3 where feature end is out of the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region),(plot_end_region - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					cds_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, (five_utr_start, five_utr_end) in five_utr_dic.items():
		if five_utr_end >= plot_start_region and five_utr_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
			if five_utr_start >= plot_start_region and five_utr_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
				five_utr_to_plot.append(((five_utr_start - plot_start_region), (five_utr_end - plot_start_region)))
			elif five_utr_start <= plot_start_region and five_utr_end <= plot_end_region:#case2 where feature start is out of the plot bounds
				five_utr_to_plot.append((0, (five_utr_end - plot_start_region)))
			elif five_utr_start >= plot_start_region and five_utr_end >= plot_end_region:#case3 where feature end is out of the plot bounds
				five_utr_to_plot.append(((five_utr_start - plot_start_region),(plot_end_region - plot_start_region)))
			elif five_utr_start <= plot_start_region and five_utr_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				five_utr_to_plot.append((0, (plot_end_region - plot_start_region)))

	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = atg_genomic_pos - plot_start_region
	tss_pos = most_upstream_pos - plot_start_region
	cov_walk_start = start - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "+"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end+1]
	try:
		generate_plot( values, svalues, fig_file, atg_pos, cov_walk_start, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot)
	except:
		print( "ERROR: plot failed" + gene )
		
	return { 'TSS': most_upstream_pos, 'start': start, 'end': end }


def run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos ):
	"""! @brief run analysis on reverse strand """
	
	most_downstream_pos = end
	final_pos_status = False
	while not final_pos_status:
		
		# --- walk coverage upstream of transcription start while there is coverage --- #
		while cov_per_contig[ most_downstream_pos ] >= mincov:	#index = next genomic position but the coverage of the next successive position needs to be looked at and not the current position
			most_downstream_pos += 1	#move one step downstream
			if most_downstream_pos == hard_cutoff:	#stop if end of contig/pseudochromosome is reached
				break
		
		# --- try to cross intron --- #
		current_position = most_downstream_pos + 1	#most_downstream_pos has coverage above cutoff (position, not index!)
		if current_position < hard_cutoff:
			while cov_per_contig[ current_position] < mincov:	#check if downstream position has low coverage
				current_position += 1	#move one step downstream
				if current_position == hard_cutoff:
					break
		else:
			final_pos_status = True
		
		avg_gap_coverage = sum( scov_per_contig[ most_downstream_pos-1:current_position-1] )/( current_position-most_downstream_pos )
		#average coverage in intron should be very low
		if current_position < ( len( seq_per_contig ) - min_exon_size ) and avg_gap_coverage > mincov:
			# --- check coverage gaps for (canonical) splice sites to continue across introns --- #
			acceptor_splice_site = seq_per_contig[most_downstream_pos-1:most_downstream_pos+1] # CT
			donor_splice_site = seq_per_contig[current_position-3:current_position-1] # AC
			if donor_splice_site == "AC" and acceptor_splice_site == "CT":	#reverse sequences of GT-AG
				print("donor splice site: " + donor_splice_site)
				print("acceptor splice site: " + acceptor_splice_site)
				most_downstream_pos = current_position + 1
			elif donor_splice_site == "GC" and acceptor_splice_site == "CT":
				print("Non-canonical donor splice site: " + donor_splice_site)
				print("Non-canonical acceptor splice site: " + acceptor_splice_site)
				most_downstream_pos = current_position + 1
			elif splicesites == "off":	#ignore check for canonical splice sites
				most_downstream_pos = current_position + 1
			elif current_position - most_downstream_pos < tolerated_gap:
				most_downstream_pos = current_position + 1
			else:
				final_pos_status = True
		else:
			final_pos_status = True
	print( "TSS position of " + gene + ": " + str( most_downstream_pos ) )

	# --- generate figures to visualize coverage around the TSS for manual inspection --- #
	transcript_list = transcripts_per_gene[gene]
	mrna_dic = {}
	cds_dic={}
	five_utr_dic={}
	for each in transcript_list:
		if each in mrna_infos:
			mrna_dic[each]= (mrna_infos[each]['start'], mrna_infos[each]['end'])
		if each in cds_infos:
			cds_dic[each] = cds_infos[each]
		if each in five_utr_infos:
			five_utr_dic[each] = (five_utr_infos[each]['start'], five_utr_infos[each]['end'])
	plot_start_region = min(end - flank_region_for_plot, atg_genomic_pos - flank_region_for_plot)
	if most_downstream_pos < ( len( seq_per_contig ) - flank_region_for_plot ):
		plot_end_region = most_downstream_pos + flank_region_for_plot
	else:
		plot_end_region = len( seq_per_contig )

	mrnas_to_plot = []
	cds_to_plot = []
	five_utr_to_plot = []
	for transcript, (mrna_start, mrna_end) in mrna_dic.items():
		if mrna_end >= plot_start_region and mrna_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
			if mrna_start >= plot_start_region and mrna_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (mrna_end - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end <= plot_end_region:#case2 where feature start is out of the plot bounds
				mrnas_to_plot.append((0, (mrna_end - plot_start_region)))
			elif mrna_start >= plot_start_region and mrna_end >= plot_end_region:#case3 where feature end is out of the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region),(plot_end_region - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				mrnas_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, cds_list in cds_dic.items():
		for (cds_start, cds_end) in cds_list:#two level loop for cds_dic alone since cds_dic structure is a list of tuples per transcript similar to the cds_infos structure from which it is derived
			if cds_end >= plot_start_region and cds_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
				if cds_start >= plot_start_region and cds_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (cds_end - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end <= plot_end_region:#case2 where feature start is out of the plot bounds
					cds_to_plot.append((0, (cds_end - plot_start_region)))
				elif cds_start >= plot_start_region and cds_end >= plot_end_region:#case3 where feature end is out of the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region),(plot_end_region - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					cds_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, (five_utr_start, five_utr_end) in five_utr_dic.items():
		if five_utr_end >= plot_start_region and five_utr_start <= plot_end_region:#primary check to see if feature is within the bounds of the plot
			if five_utr_start >= plot_start_region and five_utr_end <= plot_end_region:#case1 where feature is entirely within the plot bounds
				five_utr_to_plot.append(((five_utr_start - plot_start_region), (five_utr_end - plot_start_region)))
			elif five_utr_start <= plot_start_region and five_utr_end <= plot_end_region:#case2 where feature start is out of the plot bounds
				five_utr_to_plot.append((0, (five_utr_end - plot_start_region)))
			elif five_utr_start >= plot_start_region and five_utr_end >= plot_end_region:#case3 where feature end is out of the plot bounds
				five_utr_to_plot.append(((five_utr_start - plot_start_region),(plot_end_region - plot_start_region)))
			elif five_utr_start <= plot_start_region and five_utr_end >= plot_end_region:#case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				five_utr_to_plot.append((0, (plot_end_region - plot_start_region)))

	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = atg_genomic_pos - plot_start_region
	tss_pos = most_downstream_pos - plot_start_region
	cov_walk_start = end - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "-"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end + 1]
	try:
		generate_plot( values, svalues, fig_file, atg_pos, cov_walk_start, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot)
	except:
		print( "ERROR: plot failed" + gene )
		
	return { 'TSS': most_downstream_pos, 'start': start, 'end': end }


def find_flanking_genes( gene, gene_infos, genes_per_chromosome, window ):
	"""! @brief find upstream and downstream genes """
	
	just_up_gene = False
	just_down_gene = False
	chromosome = gene_infos[ gene ]['chromosome']
	gene_order = genes_per_chromosome[ chromosome ]
	index = gene_order.index( gene )
	counter = 1
	up_genes=[]
	if index > 0:
		if index < window:
			cutoff = index
		else:
			cutoff = window
		while counter <= cutoff:
			up_gene = gene_order[ index - counter ]
			if counter == 1:
				just_up_gene = up_gene
			up_genes.append(up_gene)
			counter+=1
	counter=1
	down_genes=[]
	if index < len( gene_order )-1:
		if ((len(gene_order)) - (index + 1)) < window:
			cutoff = ((len(gene_order)) - (index + 1))
		else:
			cutoff = window
		while counter <= cutoff:
			down_gene = gene_order[index + counter]
			if counter==1:
				just_down_gene = down_gene
			down_genes.append(down_gene)
			counter+=1
	return just_up_gene, just_down_gene, up_genes, down_genes
	

def extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size ):
	"""! @brief extract promoter region """
	
	tss = result['TSS']
	gene_start = result['start']
	gene_end = result['end']
	if orientation == "+":	#forward strand
		if (tss - hard_cutoff) > min_promoter_size:
			if (tss - hard_cutoff) > max_promoter_size:
				promoter = seq_per_contig[ tss-max_promoter_size:tss ]
				promoter_status = True
			else:
				promoter = seq_per_contig[ hard_cutoff:tss ]
				promoter_status = True
		else:
			#no promoter detected (returning everything upstream of start codon
			promoter = seq_per_contig[ hard_cutoff:gene_start ]
			promoter_status = False
		# retrieving downstream sequence
		downstream_to_tss = seq_per_contig[tss:tss + downstream_size]
		#retrieving a sequence of fixed length (250 bp around the TSS with 200 bp upstream and 50 bp downstream)
		full_seq = promoter[-(upstream_slice):] + downstream_to_tss[:downstream_slice]
	else:	#reverse strand
		if (hard_cutoff - tss) > min_promoter_size:
			if (hard_cutoff - tss) > max_promoter_size:
				promoter = seq_per_contig[ tss:tss+max_promoter_size ]
				promoter_status = True
			else:
				promoter = seq_per_contig[ tss:hard_cutoff ]
				promoter_status = True
		else:
			#no promoter detected (returning everything upstream of start codon
			promoter = seq_per_contig[ gene_end:hard_cutoff ]
			promoter_status = False
		# retrieving downstream sequence
		downstream_to_tss = seq_per_contig[tss - downstream_size:tss]
		# retrieving a sequence of fixed length (250 bp around the TSS with 200 bp upstream and 50 bp downstream)
		full_seq = downstream_to_tss[-(downstream_slice):] + promoter[:upstream_slice]
	return promoter_status, promoter, downstream_to_tss, full_seq

#function to sort sequences from MOODS results on the basis of strandedness and proximity to TSS
def sort_key(row, orientation):
	position = int(row[2])
	strand=row[3]
	if orientation == '+':
		priority = 0 if strand == '+' else 1
	elif orientation == '-':
		priority = 0 if strand =='-' else 1
	if strand == '+':
		return(priority,-position)
	elif strand == '-':
		return(priority, position)

def run_moods (sequence, pfm_folder, output_file, moods, pvalue):
	cmd = 'python3 ' + moods + ' -m ' + pfm_folder + '/*.pfm ' + '-s ' + sequence + ' -p ' + str(pvalue) + ' -o ' + output_file
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()
	return output_file

def compute_moods_score(moods_file,seq_len):
	score = 0
	rows = []
	with open(moods_file, 'r') as f:
		line = f.readline()
		while line:
			row = line.strip().rstrip(',').split(',')
			rows.append(row)
			line = f.readline()
	for hit in rows:
		score += float(hit[4])
	final_score = score/seq_len
	return final_score, rows

def compute_background_moods_scores (background_seqs, pfm_config_dic, tmp_folder, output_folder, moods, pvalue):
	print(f"total no. of background seqs is {len(background_seqs)}")
	background_scores = []
	counter = 1
	for seq in background_seqs:
		seq_len = len(seq)
		tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', dir=tmp_folder, delete=False)
		try:
			tmp.write(f">seq_{counter}\n{seq}\n")
			tmp.flush()
			tmp.close()
			cumulative_promoter_motif_score = 0
			for elements in pfm_config_dic:
				pfm_folder = pfm_config_dic[elements][0]
				output_file = os.path.join(tmp_folder, f"seq_{counter}_moods_{elements}.txt")
				moods_file = run_moods(tmp.name, pfm_folder, output_file, moods, pvalue)
				score, sorted_hits = compute_moods_score(moods_file, seq_len)
				cumulative_promoter_motif_score += score
			background_scores.append(cumulative_promoter_motif_score)
		finally:
			if os.path.exists(tmp.name):
				os.remove(tmp.name)
			counter += 1
	return background_scores

#function to scan extracted promoter sequences for user-specified promoter motif elements
def promoter_motif_analysis (numcols, upstream_slice, downstream_slice, bg_scores, result, gene, orientation, promoter_seq, downstream_to_tss, moods, pvalue, top_motifs, pfm_config_dic, tmp_folder, output_folder):
	tss_neighbourhood = os.path.join(output_folder,f'{gene}_tss_neighbourhood.png')
	tss = result['TSS']
	tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', dir=tmp_folder, delete=False)
	try:
		if orientation == '+':
			full_seq = promoter_seq[-(upstream_slice):] + downstream_to_tss[:downstream_slice]
			promoter_seq = promoter_seq[-(upstream_slice):]
			downstream_to_tss = downstream_to_tss[:downstream_slice]
		elif orientation == '-':
			full_seq = downstream_to_tss[-(downstream_slice):] + promoter_seq[:upstream_slice]
			promoter_seq = promoter_seq[:upstream_slice]
			downstream_to_tss = downstream_to_tss[-(downstream_slice):]
		tmp.write(f">{gene}\n{full_seq}\n")
		tmp.flush()
		tmp.close()
		seq_len = len(full_seq)
		downstream_len = len(downstream_to_tss)

		sorted_hits_dic = {}
		cumulative_promoter_motif_score = 0
		for elements in pfm_config_dic:
			pfm_folder = pfm_config_dic[elements][0]
			output_file = os.path.join(tmp_folder, f"{gene}_moods_{elements}.txt")
			moods_file = run_moods(tmp.name, pfm_folder, output_file, moods, pvalue)
			score, sorted_hits = compute_moods_score( moods_file, seq_len)
			cumulative_promoter_motif_score += score
			sorted_hits_dic[elements]=sorted_hits
		#computing the percentile of cumulative promoter score with respect to the background scores based on the MOODS scoring
		percentile_of_promoter_score = percentileofscore(bg_scores, cumulative_promoter_motif_score, kind='rank')
		#positional scanning of promoter seqs for canonical hits
		canonical_hits = 0
		position_dic = {}
		for elements in sorted_hits_dic:
			position_dic[elements] = defaultdict(list)  # defaultdict method is used here to account for cases where + and - hits occur at the same genomic position in which since genomic position is the key, it will be overwritten to retain just the last entry
			# checking for direction sensitive motif presence both upstream and downstream
			for hit in sorted_hits_dic[elements]:
				if pfm_config_dic[elements][3] == 'yes': # if direction sensitivity is yes in config
					if orientation == "+":# for the code block with direction sensitivity check, collecting hits in all positions irrespective of direction before the direction check is performed
						plot_pos = tss - (len(promoter_seq) - int(hit[2]))
						position_dic[elements][plot_pos].append(hit[3])
					elif orientation == '-':
						plot_pos = tss + (int(hit[2]))
						position_dic[elements][plot_pos].append(hit[3])
					if hit[3] == orientation:
						if orientation == '+':
							if int(hit[2]) <= len(promoter_seq):
								if ((len(promoter_seq) - int(hit[2])) <= pfm_config_dic[elements][1]):
									canonical_hits += 1
							else:
								if int(hit[2]) - len(promoter_seq) <= pfm_config_dic[elements][2]:
									canonical_hits += 1
						elif orientation == '-':
							if int(hit[2]) >= downstream_len:  # upstream
								distance = int(hit[2]) - downstream_len
								if distance <= pfm_config_dic[elements][1]:
									canonical_hits += 1
							elif int(hit[2]) < downstream_len:  # downstream
								distance = downstream_len - int(hit[2])
								if distance <= pfm_config_dic[elements][2]:
									canonical_hits += 1
				else:# if direction sensitivity is no in config
					if orientation == '+':
						plot_pos = tss - (len(promoter_seq) - int(hit[2]))
						position_dic[elements][plot_pos].append(hit[3])
						if int(hit[2]) <= len(promoter_seq):
							if ((len(promoter_seq) - int(hit[2])) <= pfm_config_dic[elements][1]):
								canonical_hits += 1
						else:
							if int(hit[2]) - len(promoter_seq) <= pfm_config_dic[elements][2]:
								canonical_hits += 1
					elif orientation == '-':
						plot_pos = tss + (int(hit[2]))
						position_dic[elements][plot_pos].append(hit[3])
						if int(hit[2]) >= downstream_len:  # upstream
							distance = int(hit[2]) - downstream_len
							if distance <= pfm_config_dic[elements][1]:
								canonical_hits += 1
						elif int(hit[2]) < downstream_len:  # downstream
							distance = downstream_len - int(hit[2])
							if distance <= pfm_config_dic[elements][2]:
								canonical_hits += 1
		motifs = []
		for elements in position_dic:
			motifs.append({
				"label": f"{elements} motif",
				"positions": position_dic[elements],
			})
		n_motifs = len(motifs)
		ncols = numcols
		nrows = math.ceil(n_motifs / ncols)
		fig, axes = plt.subplots(nrows, ncols,figsize=(6 * ncols, 3 * nrows),squeeze=False)
		#lollipop plot of up and downstream regions of predicted TSS
		fixed_heights = [0.4, 0.6, 0.8, 1.0]  # fixed height levels since no histogram to scale to
		legend_elements = [
			Line2D([0], [0], color='black', linewidth=2, linestyle='solid', label='TSS'),
			Line2D([0], [0], color='green', linewidth=1, linestyle='dotted', label='Same orientation motif'),
			Line2D([0], [0], color='red', linewidth=1, linestyle='dotted', label='Reverse orientation motif')
		]
		for idx, motif in enumerate(motifs):
			row, col = divmod(idx, ncols)
			ax = axes[row][col]

			for i, (pos, signs) in enumerate(motif["positions"].items()):
				for j, sign in enumerate(signs):
					colour = 'green' if sign == orientation else 'red'
					h = fixed_heights[(i + j) % len(fixed_heights)]
					ax.vlines(pos, ymin=0, ymax=h,
							  linestyle='dotted', color=colour, linewidth=1.5, zorder=2)
					ax.plot(pos, h, 'o', color=colour, markersize=4, zorder=2)

			ax.vlines(tss, ymin=0, ymax=1.2,
					  linestyle='solid', color='black', linewidth=2, zorder=3)
			ax.axhline(0, color='black', linewidth=0.5)
			#setting the plot boundary to the length of the sequence used for the promoter analysis
			if orientation == '+':
				ax.set_xlim(tss - upstream_slice, tss + downstream_slice)
			elif orientation == '-':
				ax.set_xlim(tss - downstream_slice, tss + upstream_slice)
			ax.set_yticks([])
			ax.set_xlabel('Genomic position')
			ax.set_title(motif["label"], fontsize=10, fontweight='bold')
			ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
			ax.ticklabel_format(style='plain', axis='x')
			plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
			ax.legend(handles=legend_elements, loc='best', fontsize=7)

		# Hide any unused subplot panels
		for idx in range(n_motifs, nrows * ncols):
			row, col = divmod(idx, ncols)
			axes[row][col].set_visible(False)

		plt.tight_layout()
		plt.savefig(tss_neighbourhood, dpi=600)
		plt.close()
	finally:
		if os.path.exists(tmp.name):
			os.remove(tmp.name)
	return cumulative_promoter_motif_score, percentile_of_promoter_score, canonical_hits


def main( arguments ):
	"""! @brief run everything """
	
	#input: GFF file with gene positions (translation start sites)
	gff_file = arguments[ arguments.index('--gff')+1 ]
	
	#input FASTA with genome sequence (check splite sites)
	fasta_file = arguments[ arguments.index('--fasta')+1 ]
	
	#gene of interest (string or file)
	goi_input = arguments[ arguments.index('--goi')+1 ]
	if not os.path.isfile( goi_input ):
		if "," in goi_input:
			goi = goi_input.strip().split(',')
		else:
			goi = [ goi_input.strip() ]
	else:
		with open( goi_input, "r" ) as f:
			lines = f.readlines()
			goi = []
			for line in lines:
				if len( line.strip() ) > 3:
					goi.append( line.strip() )
	print( "Number of detected genes of interest: " + str( len( goi ) ) )

	#gff file config params
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config')+1]
		with open (gff_config_file, 'r') as f:
			for line in f:
				parts = line.strip().split()
				child_attribute = parts[0]
				child_parent_linker = parts[1]
				parent_attribute = parts[2]
	else:
		child_attribute = 'ID'
		child_parent_linker = 'Parent'
		parent_attribute = 'ID'

	#output folder
	output_folder = arguments[ arguments.index('--out')+1 ]
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	#tmp folder
	tmp_folder = os.path.join(output_folder,'Tmp')
	if not os.path.exists(tmp_folder):
		os.makedirs(tmp_folder)

	if '--samtools' in arguments:
		samtools = arguments[arguments.index('--samtools') + 1]
	else:
		samtools = "samtools"

	if '--bedtools' in arguments:
		bedtools = arguments[arguments.index('--bedtools') + 1]
	else:
		bedtools = "genomeCoverageBed"

	if '--m' in arguments:
		m = arguments[arguments.index('--m') + 1]
	else:
		m = "5000000"

	if '--threads' in arguments:
		t = arguments[arguments.index('--threads') + 1]
	else:
		t = "4"


	if '--intron_percentile_cutoff' in arguments:
		percentile_cut = int(arguments[arguments.index('--intron_percentile_cutoff') + 1])
	else:
		percentile_cut = 99

	if '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]

		if '--bam_is_sorted' in arguments:
			bam_sorted_status = True
		else:
			bam_sorted_status = False
		
		if not bam_sorted_status:	#sorting the BAM file if it was not sorted already
			print ("sorting BAM file ...")
			sorted_bam_file = output_folder + "sorted.bam"
			cmd = samtools + " sort -m " + str(m) + " --threads " + t + " " + bam_file + " > " + sorted_bam_file
			p = subprocess.Popen( args= cmd, shell=True )
			p.communicate()
		
		cov_file = output_folder + "reads_aligned.cov"
		scov_file = output_folder + "reads_spanning.cov"
		
		if bam_sorted_status:
			if not os.path.isfile( cov_file ):
				construct_cov_file( bam_file, cov_file, bedtools )
			if not os.path.isfile( scov_file ):
				construct_scov_file( bam_file, scov_file, bedtools )
		else:
			if not os.path.isfile( cov_file ):
				construct_cov_file( sorted_bam_file, cov_file, bedtools )
			if not os.path.isfile( scov_file ):
				construct_scov_file( sorted_bam_file, scov_file, bedtools )
		input_mode = "cov"
		
	elif '--sra_folder' not in arguments:
		cov_file = arguments[ arguments.index('--cov')+1 ]
		scov_file = arguments[ arguments.index('--scov')+1 ]
		if cov_file.split('.')[-1].lower() == "gz":
			input_mode = "gz"
		else:
			input_mode = "cov"

	#coverage cutoff
	if '--mincov' in arguments:
		mincov = int( arguments[ arguments.index('--mincov')+1 ] )
	else:
		mincov = 1
	
	#minimal exon size
	if '--minexon' in arguments:
		min_exon_size = int( arguments[ arguments.index('--minexon')+1 ] )
	else:
		min_exon_size = 10
	
	#flanking region for plot
	if '--flanksize' in arguments:
		flank_region_for_plot = int( arguments[ arguments.index('--flanksize')+1 ] )
	else:
		flank_region_for_plot = 50

	#gene neighbourhood window for overlapping gene analysis
	if '--neighbourhood' in arguments:
		window = int( arguments[arguments.index('--neighbourhood')+1])
	else:
		window = 5
	
	#tolerated coverage gap size (due to sequence variant)
	if '--gapsize' in arguments:
		tolerated_gap = int( arguments[ arguments.index('--gapsize')+1 ] )
	else:
		tolerated_gap = 5
	
	if '--splicesites' in arguments:
		splicesites = arguments[ arguments.index('--splicesites')+1 ]
		if splicesites not in [ "strict", "off" ]:
			splicesites = "strict"
			print( "WARNING: splice site handling set to strict due to unrecognized user input." )
	else:
		splicesites = "strict"
	
	if '--min_promoter_size' in arguments:
		min_promoter_size = arguments[arguments.index('--min_promoter_size')+1]
	else:
		min_promoter_size = 50

	if '--max_promoter_size' in arguments:
		max_promoter_size = arguments[arguments.index('--max_promoter_size')+1]
	else:
		max_promoter_size = 1000

	if '--downstream_size' in arguments:
		downstream_size = arguments[arguments.index('--downstream_size')+1]
	else:
		downstream_size = 300

	if '--background' in arguments:
		background_strength = int(arguments[arguments.index('--background')+1])
	else:
		background_strength = 1000

	if '--aligner' in arguments:# <STAR/HISAT2> option for user to use star or HISAT2 aligners
		aligner = arguments[arguments.index('--aligner')+1]
	else:
		aligner = 'STAR'

	if '--HISAT2' in arguments:#full path to HISAT2 for RNAseq mapping
		hisat2 = arguments[arguments.index('--HISAT2')+1]
	else:
		hisat2 = 'hisat2'

	if '--STAR'in arguments:#full path to STAR for RNAseq mapping
		star = arguments[arguments.index('--STAR')+1]
	else:#recommended to use STARlong as it is applicable to both short and long reads.
		star = 'STARlong'

	if '--index_bases' in arguments:#parameter for the genomeSAindexNbases flag in STAR
		index_bases = arguments[arguments.index('--index_bases')+1]
	else:
		index_bases = 12

	if '--run_mode' in arguments: #to provide mode options for the user; make_bam starts from STAR indexing; find_tss starts from bam or coverage file.
		run_mode = arguments[arguments.index('--run_mode')+1]
	else:
		run_mode = 'find_tss'

	if '--sra_folder' in arguments:
		sra_folder = arguments[arguments.index('--sra_folder')+1]
	else:
		sra_folder = ''

	if '--fastq_pattern' in arguments:
		pattern_names=arguments[arguments.index('--fastq_patern')+1]#specify fastq read file name pattern for the paired end files separated by commas without spaces like - _pass_1,_pass_2_
		pattern_names_list = pattern_names.split(',')
	else:
		pattern_names_list = ["_pass_1", "_pass_2"]

	if '--analyse_promoter' in arguments: #yes or no for promoter analysis with MOODS
		promoter_analysis = arguments[arguments.index('--analyse_promoter')+1]
	else:
		promoter_analysis = 'no'

	if '--moods' in arguments: #full path to moods python script
		moods = arguments[arguments.index('--moods')+1]
	else:
		moods = 'moods-dna.py'

	pfm_config_dic = {}
	if promoter_analysis == 'yes':
		#full path to config file with full path to folder(s) with PFM matrices, upstream, downstream motif position boundaries, directional sensitivity of the promoter motif (yes or no) for different motif elements respectively
		pfm_config = arguments[arguments.index('--PFM')+1]
		with open (pfm_config, 'r') as f:
			for line in f:
				parts = line.strip().split()
				pfm_config_dic[parts[0]] = (parts[1], int(parts[2]), int(parts[3]), parts[4])

		#number of columns for plotting the motif hits
		if '--numcols' in arguments:
			numcols = int(arguments[arguments.index('--numcols')+1])
		else:
			numcols = 2

	#p-value threshold for MOODS analysis
	if '--pval' in arguments:
		pvalue = float(arguments[arguments.index('--pval')+1])
	else:
		pvalue = 0.001

	if '--top_motif_hit' in arguments:
		top_motifs = int(arguments[arguments.index('--top_motif_hit')+1])
	else:
		top_motifs = 1

	if '--coverage_walk_origin' in arguments: #cds or utr
		coverage_walk_origin = arguments[arguments.index('--coverage_walk_origin')+1]
	else:
		coverage_walk_origin = 'cds'

	if '--upstream_slice' in arguments:
		upstream_slice = int(arguments[arguments.index('--upstream_slice')+1])
	else:
		upstream_slice = 200

	if '--downstream_slice' in arguments:
		downstream_slice = int(arguments[arguments.index('--downstream_slice')+1])
	else:
		downstream_slice = 50

	#code block to do RNAseq mapping of SRA files and then produce the cov files for the tss analysis
	if run_mode == 'make_bam':
		#Calculating the optimal intron size for the species to be analysed
		trans_exon_map = {}
		intron_sizes = []

		# parsing the gff file and making a dictionary where keys are transcript names and values are lists of tuples where each tuple is the start, end coordinate of exons in that transcript
		with open(gff_file,'r') as f:
			for line in f:
				if not line.startswith('#'):
					parts = line.strip().split('\t')
					if parts[2].lower() == 'exon':
						fields = parts[8].strip().split(';')
						for each in fields:
							if f'{child_parent_linker}=' in each:
								transcript = each.replace(f'{child_parent_linker}=', '')
								start = int(parts[3])
								end = int(parts[4])
								if transcript not in trans_exon_map:
									trans_exon_map[transcript] = [(start, end)]
								else:
									trans_exon_map[transcript].append((start, end))

		# calculate intron sizes
		for transcript in trans_exon_map.keys():
			exons = trans_exon_map[transcript]
			if len(exons) < 2:
				continue
			counter = 0
			# sort by ascending order of start coordinates
			exons.sort(key=lambda x: x[0])
			while counter <= (len(exons) - 2):
				intron = exons[counter + 1][0] - exons[counter][1] - 1
				if intron <= 0:
					pass
				else:
					intron_sizes.append(intron)
				counter += 1

		# Calculate statistics
		median_size = np.median(intron_sizes)
		mean_size = np.mean(intron_sizes)
		intron_cutoff = np.percentile(intron_sizes, percentile_cut)
		intron_cutoff = int(Decimal(str(intron_cutoff)).quantize(0, rounding=ROUND_HALF_DOWN))
		intron_min = np.min(intron_sizes)
		intron_max = np.max(intron_sizes)
		print(f"Total number of introns: {len(intron_sizes)}")
		print(f"Minimum intron size is {intron_min} and maximum intron size is {intron_max}")
		print(f"Median intron size is {median_size} and mean intron size is {mean_size}")
		print(f"Intron size for the --alignIntronMax flag after rounding is {intron_cutoff}")
		# plotting the intron size distribution
		intron_plot=os.path.join(output_folder,'Intron_size_distribution.png')
		plt.figure(figsize=(8, 5))
		plt.hist(intron_sizes, bins=100, color='steelblue', edgecolor='black')
		plt.xlabel('Intron length (bp)')
		plt.ylabel('Frequency')
		plt.title('Distribution of intron sizes')
		plt.grid(True, linestyle='--', alpha=0.6)
		plt.tight_layout()
		plt.savefig(intron_plot,dpi=600)
		star_indexing_folder=os.path.join(output_folder,'Genome_index')
		os.mkdir(star_indexing_folder)
		if aligner == 'STAR':
			# Indexing with STAR
			cmd = star+' --runMode genomeGenerate --genomeDir '+star_indexing_folder+' --genomeFastaFiles '+fasta_file+' --sjdbGTFfile '+gff_file+' --sjdbGTFtagExonParentTranscript Parent --runThreadN '+t+' --genomeSAindexNbases '+str(index_bases)
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
		elif aligner == 'HISAT2':
			#Indexing with HISAT2
			index_file_name = os.path.join(star_indexing_folder,'Index')
			cmd = hisat2 + '-build -p ' + t + ' ' + fasta_file + ' ' + index_file_name
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
		star_mapping_folder = os.path.join(output_folder, 'RNA-seq_map')
		os.mkdir(star_mapping_folder)
		if sra_folder:
			pairs = defaultdict(dict)  # create a default dictionary for holding the paired end files

			# Recursively walk through all subdirectories
			for root, dirs, files in os.walk(sra_folder):
				for f in files:
					# Only consider FASTQ files
					if not f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
						continue

					# Store the full path instead of just filename
					full_path = os.path.join(root, f)

					if pattern_names_list[0] in f:
						sample = f.split(f"{pattern_names_list[0]}.")[0]
						pairs[sample]["R1"] = full_path

					elif pattern_names_list[1] in f:
						sample = f.split(f"{pattern_names_list[1]}.")[0]
						pairs[sample]["R2"] = full_path
			# Add check to ensure pairs were found
			if not pairs:
				raise ValueError(f"No paired-end files found in {sra_folder}. Check file naming pattern.")

			print(f"Found {len(pairs)} sample(s) to process:")
			for sample in pairs.keys():
				print(f"  - {sample}")

			for sample, reads in pairs.items():
				if "R1" in reads and "R2" in reads:
					r1 = os.path.join(sra_folder, reads["R1"])
					r2 = os.path.join(sra_folder, reads["R2"])
					if aligner == 'STAR':
						#RNAseq mapping with STAR
						prefix=os.path.join(star_mapping_folder,sample+'_')
						cmd = 'ulimit -n 4096 && '+star+' --runMode alignReads --genomeDir '+star_indexing_folder+' --outSAMtype BAM SortedByCoordinate --readFilesIn '+r1+' '+r2+' --runThreadN '+t+' --outFileNamePrefix '+prefix+' --readFilesCommand zcat --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --alignIntronMax '+ str(intron_cutoff)
						p = subprocess.Popen(args=cmd, shell=True)
						p.communicate()
					elif aligner == 'HISAT2':
						sorted_bam = os.path.join(star_mapping_folder, sample + '_sorted.bam')
						cmd = hisat2 + ' --max-intronlen ' + str(intron_cutoff) + ' -p ' + str(t) + ' -x ' + index_file_name + ' -1 ' + r1 + ' -2 ' + r2 + ' | ' + samtools + ' sort --threads ' + str(t) + ' -O BAM -o ' + sorted_bam
						p = subprocess.Popen(args=cmd, shell=True)
						p.communicate()
			#merging all the sorted BAM files obtained from STARlong mapping
			bam_files=os.path.join(output_folder,'bam_files.txt')
			bam_list = [
				os.path.join(star_mapping_folder, f)
				for f in os.listdir(star_mapping_folder)
				if f.endswith('.bam')
			]
			if len(bam_list) > 1:
				with open(bam_files, 'w') as out:
					for f in bam_list:
						if f.endswith('bam'):
							out.write(f + '\n')
				merged_bam=os.path.join(output_folder,'merged.bam')
				cmd = samtools + ' merge --threads '+ t + ' -o ' + merged_bam + ' -b ' + bam_files
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()
				# Check if merge was successful
				if not os.path.exists(merged_bam):
					raise FileNotFoundError(f"Merge failed: {merged_bam} was not created")

			elif len(os.listdir(star_mapping_folder)) == 1:
				for f in (os.listdir(star_mapping_folder)):
					sorted_merged_bam_file = f
			cov_file = output_folder + "reads_aligned.cov"
			scov_file = output_folder + "reads_spanning.cov"

			if not os.path.isfile(cov_file):
				construct_cov_file(merged_bam, cov_file, bedtools)
			if not os.path.isfile(scov_file):
				construct_scov_file(merged_bam, scov_file, bedtools)
			input_mode = "cov"

	# --- load data --- #
	coverage = load_coverage( cov_file, input_mode )
	scoverage = load_coverage( scov_file, input_mode )
	
	gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos = load_gene_infos( gff_file, child_attribute, child_parent_linker, parent_attribute )
	genome_seq, seq_counter = load_sequences( fasta_file )
	if promoter_analysis == 'yes':
		print("Retrieving background seqs.")
		background_seq_len = upstream_slice + downstream_slice
		background_seqs = get_random_background_seqs(genome_seq, background_seq_len, background_strength)
		print("Completed retrieving background seqs.")
		print("Starting MOODS scoring of background seqs.")
		bg_scores = compute_background_moods_scores(background_seqs, pfm_config_dic, tmp_folder, output_folder, moods,pvalue)
		print("Completed MOODS scoring of background seqs.")
	#run analysis per gene of interest
		# run analysis per gene of interest
		results = {}
		motifs = []
		motif_scores = {}
		pvalue_dic = {}
		full_seq_pos_strand_gene = {}
		full_seq_neg_strand_gene = {}
		# confidence_score_dic = {}
		tss_confidence_dic = {}
		isoforms_dic = {}
		distance_dic = {}
		coverage_dic = {}
		five_utr_dic = {}
		percentile_dic = {}
		canonical_hits_dic = {}
	for gene in goi:
		try:
			cov_per_contig = coverage[gene_infos[gene]['chromosome']]  # get coverage of the sequence that harbours the gene of interest
			scov_per_contig = scoverage[gene_infos[gene]['chromosome']]  # get spanning read coverage of the sequence that harbours the gene of interest
			seq_per_contig = genome_seq[gene_infos[gene]['chromosome']]  # get the sequence of the contig/pseudochromosome that harbours the gene of interest
			# code block to check if the goi has annotated 5'UTR and if yes take the most upstream/ downstream 5'UTR start/ end as start or end according to + or - strand orientation
			# Initialize with gene coordinates as default with 5'UTR checks downstream
			start, end, orientation, contig = gene_infos[gene]['start'], gene_infos[gene]['end'], gene_infos[gene]['orientation'], gene_infos[gene]['chromosome']  # get information about gene of interest if it does not have 5'UTRs annotated
			# store ATG position as fixed reference before start may be modified by 5'UTR code block
			if gene in gene_atg_dic:
				atg_pos = gene_atg_dic[gene]
			else:
				atg_pos = start if orientation == '+' else end  # fallback to gene boundary if no CDS found
			if gene in transcripts_per_gene:
				transcript_list = transcripts_per_gene[gene]
				for each in transcript_list:
					if coverage_walk_origin == 'utr':
						if each in five_utr_infos and gene_infos[gene]['orientation'] == '+' and each == transcript_list[0]:
							start = five_utr_infos[each]['start']  # in case a gene has transcripts with 5'UTR annotated, the most upstream 5'UTR will be taken as the start for the + strand gene
							end, orientation = gene_infos[gene]['end'], gene_infos[gene]['orientation']  # get information about gene of interest
							five_utr_dic[gene] = f"5'UTR start of {each} used for TSS prediction."
							break
						elif each in five_utr_infos and gene_infos[gene]['orientation'] == '-' and each == transcript_list[-1]:
							end = five_utr_infos[each]['end']  # in case a gene has transcripts with 5'UTR annotated, the most downstream 5'UTR will be taken as the start for the - strand gene
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
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No 5'UTR annotated. {gene} end used for TSS prediction."
			elif coverage_walk_origin == 'cds':
				if gene not in five_utr_dic and orientation == '+':
					five_utr_dic[gene] = f"No CDS annotated. {gene} start used for TSS prediction."
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No CDS annotated. {gene} end used for TSS prediction."
			upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes( gene, gene_infos, genes_per_chromosome, window )
			avg_cov_gene = sum(cov_per_contig[start - 1:end]) / (end - (start - 1))  # get average coverage of the gene of interest for confidence thresholding
			#check if goi is an overlapping gene
			#TSS-blocking overlap types
			SKIP_TSS_TYPES = {'head_head', 'head_into_neighbor', 'same_strand', 'nested'}#'tail_tail' and 'tail_into_neighbor' overlap types are safe for TSS analysis
			goi_strand = gene_infos[gene]['orientation']
			blocking_overlaps = 0
			for ugene in upstream_gene_list:
				nbr = gene_infos[ugene]

				if start <= nbr['end']:  # positional overlap exists
					ov_type = get_overlap_type(goi_strand, start, end,nbr['orientation'], nbr['start'], nbr['end'])
					print(f"  {gene} - {ugene}: {ov_type} overlap")
					if ov_type in SKIP_TSS_TYPES:
						blocking_overlaps += 1
				#'tail_tail' → do NOT increment; TSS analysis can still run
			#downstream check is re-enabled now because strand matters:
			#a downstream neighbor can cause a head_head overlap for − strand GOIs
			for dgene in downstream_gene_list:
				nbr = gene_infos[dgene]
				if end >= nbr['start']:  # positional overlap exists
					ov_type = get_overlap_type(goi_strand, start, end,nbr['orientation'], nbr['start'], nbr['end'])
					print(f"  {gene} - {dgene}: {ov_type} overlap")
					if ov_type in SKIP_TSS_TYPES:
						blocking_overlaps += 1
			if blocking_overlaps > 0:
				print(f"{gene} has {blocking_overlaps} overlap(s). TSS analysis skipped.")
				continue
			if orientation == "+":	#only works on forward strand
				fig_file = output_folder + gene + ".png"
				if upstream_gene:
					hard_cutoff = gene_infos[ upstream_gene ]['end']

				else:
					hard_cutoff = 1
				result = run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos )
				promoter_status, promoter, downstream_to_tss, full_seq = extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size )
				full_seq_pos_strand_gene[gene] = full_seq
				if promoter_analysis == 'yes':
					if os.path.exists(moods):
						cumulative_promoter_motif_score, percentile_of_promoter_score, canonical_hits = promoter_motif_analysis(numcols, upstream_slice, downstream_slice, bg_scores, result, gene, orientation,promoter, downstream_to_tss, moods, pvalue, top_motifs, pfm_config_dic, tmp_folder,output_folder)
						#motifs.append(best_motif_hits)
						motif_scores[gene] = cumulative_promoter_motif_score
						percentile_dic[gene] = percentile_of_promoter_score
						canonical_hits_dic[gene] = canonical_hits
					else:
						print('MOODS not found. Promoter analysis not possible.')
				result.update( { 'promoter_status': promoter_status } )
				result.update( { 'promoter': promoter } )
				results.update( { gene: result } )
			else:	#solution for reverse strand genes
				fig_file = output_folder + gene + ".png"
				if downstream_gene:
					hard_cutoff = gene_infos[ downstream_gene ]['start']
				else:
					hard_cutoff = len( seq_per_contig )
				result = run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos )
				promoter_status, promoter, downstream_to_tss, full_seq = extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size )
				full_seq_neg_strand_gene[gene] = full_seq
				if promoter_analysis == 'yes':
					if os.path.exists(moods):
						cumulative_promoter_motif_score, percentile_of_promoter_score, canonical_hits = promoter_motif_analysis(numcols, upstream_slice, downstream_slice, bg_scores, result, gene, orientation,promoter, downstream_to_tss, moods, pvalue, top_motifs, pfm_config_dic, tmp_folder,output_folder)
						motif_scores[gene] = cumulative_promoter_motif_score
						percentile_dic[gene] = percentile_of_promoter_score
						canonical_hits_dic[gene] = canonical_hits
					else:
						print('MOODS not found. Promoter analysis not possible.')
				result.update( { 'promoter_status': promoter_status } )
				result.update( { 'promoter': promoter } )
				results.update( { gene: result } )
			# calculating confidence score
			isoforms = len(transcripts_per_gene[gene])  # no. of isoforms per gene
			distance = abs((results[gene]['TSS']) - atg_pos)  # distance between predicted TSS and upstream CDS start (+ gene)/ downstream CDS end (- gene)
			"""
			total_contribution = avg_cov_gene + distance + isoforms
			coverage_score = (avg_cov_gene / total_contribution)  # contribution of coverage to total score
			distance_score = 1 - (distance / total_contribution)  # contribution of distance to total score; 1 - is used since distance is inversely related to total score
			isoform_score = 1 - (isoforms / total_contribution)  # contribution of no. of isoforms to total score; 1 - is used since no. of isoforms is inversely related to total score
			TSS_confidence_score = (coverage_score * distance_score * isoform_score) ** (1 / 3) # confidence score is a geometric mean of the individual factor scores; each factor is independent of one another and must contribute well for the overal confidence making gemoetric mean and the multiplicative approach preferred over arithmetic mean and the additive approach
			"""
			isoforms_dic[gene]=isoforms
			distance_dic[gene]=distance
			coverage_dic[gene]=avg_cov_gene
		except KeyError as e:
			print(f"Missing gene error: {gene}, missing key: {e}")
			traceback.print_exc()
			"""
			# If gene was added to results but confidence calculation failed, add fall back
			if gene in results and gene not in confidence_score_dic:
				confidence_score_dic[gene] = "NA"
			"""
	# --- report TSS in output file --- #
	final_output_file = output_folder + "Results.tsv"
	pos_strand_tss_neighbourhood_file = output_folder + "Positive_strand_gene_promoter_downstream_seqs.fasta"
	neg_strand_tss_neighbourhood_file = output_folder + "Negative_strand_gene_promoter_downstream_seqs.fasta"
	with open (pos_strand_tss_neighbourhood_file, 'w') as out:
		for each in full_seq_pos_strand_gene:
			out.write(f">{each}\n{full_seq_pos_strand_gene[each]}\n")
	with open (neg_strand_tss_neighbourhood_file, 'w') as out:
		for each in full_seq_neg_strand_gene:
			out.write(f">{each}\n{full_seq_neg_strand_gene[each]}\n")
	with open( final_output_file, "w" ) as out:
		if promoter_analysis != 'yes':
			out.write( "\t".join( [ "GeneID", "TSS", "Average gene coverage", "Number of isoforms", "Start", "End", "PromoterStatus", "Promoter", "Additional comments" ] ) + "\n" )
			for gene in list( results.keys() ):
				out.write( "\t".join( [ 	gene,
													str( results[ gene ]['TSS'] ),
													str(coverage_dic[gene]),
													str(isoforms_dic[gene]),
													str( results[ gene ]['start'] ),
													str( results[ gene ]['end'] ),
													str( results[ gene ]['promoter_status'] ),
													str( results[ gene ]['promoter'] ),
													str(five_utr_dic[gene])
											] ) + "\n" )
		elif promoter_analysis == 'yes':
			out.write( "\t".join( [ "GeneID", "TSS", "Promoter motif score", "Promoter motif percentile", "Canonical motif hits", "Average gene coverage", "Number of isoforms", "Start", "End", "PromoterStatus", "Promoter", "Additional comments" ] ) + "\n" )
			for gene in list( results.keys() ):
				if gene in results:
					out.write( "\t".join( [ 	gene,
														str( results[ gene ]['TSS'] ),
														str(motif_scores[gene]),
														str(percentile_dic[gene]),
														str(canonical_hits_dic[gene]),
														str(coverage_dic[gene]),
														str(isoforms_dic[gene]),
														str( results[ gene ]['start'] ),
														str( results[ gene ]['end'] ),
														str( results[ gene ]['promoter_status'] ),
														str( results[ gene ]['promoter'] ),
														str(five_utr_dic[gene])
												] ) + "\n" )


if '--bam' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--scov' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--run_mode' in sys.argv and '--sra_folder' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main(sys.argv)
else:
	sys.exit( __usage__ )
