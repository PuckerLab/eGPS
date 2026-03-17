### Boas Pucker ###
### Shakunthala Natarajan ###
### pucker@uni-bonn.de ###
__version__ = "v0.03"

__reference__ = "Pucker et al., 2025: https://github.com/bpucker/TSS_finder"

__usage__ = """
					TSS_finder """ + __version__ + """("""+ __reference__ +""")
					
					Usage:
					python3 TSS_finder.py
					--fasta <GENOMIC_FASTA_FILE>
					--gff <GFF_FILE>
					--out <OUTPUT_FOLDER>
					[--bam <BAM_FILE>|--cov <COV_FILE> --scov <SCOV_FILE>]
					
					optional:
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
					"""


import re, os, sys, subprocess, gzip
import tempfile
import numpy as np
from decimal import Decimal, ROUND_HALF_DOWN
from collections import defaultdict
try:
	import matplotlib.pyplot as plt
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


def load_gene_infos( gff_file ):
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
					ID = parts[-1].split('ID=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					gene_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						genes_per_chromosome[ parts[0] ].append( ID )
					except KeyError:
						genes_per_chromosome.update( { parts[0]: [ ID ] } )
				if parts[2].upper() == "MRNA":
					ID = parts[-1].split('ID=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					Parent = parts[-1].split('Parent=')[-1]
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					mrna_infos.update( { ID: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
					try:
						transcripts_per_gene[ Parent ].append( ID )
					except KeyError:
						transcripts_per_gene.update( { Parent: [ ID ] } )
				if parts[2].upper() == 'FIVE_PRIME_UTR':
					Parent = parts[-1].split('Parent=')[-1]#Parent of 5'UTR is transcript
					if ";" in Parent:
						Parent = Parent.split(';')[0]
					five_utr_infos.update({ Parent: { 'chromosome': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } })# key of this nested dictionary is the transcript name
				if parts[2].upper() == 'CDS':
					cds_parents = parts[-1].split('Parent=')[-1]
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
	return gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def generate_plot( values, svalues, fig_file, atg_pos, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot ):
	"""! @brief generate a coverage plot """
	fig, ax1 = plt.subplots()
	ax1.plot(values, color="black", linestyle="solid")  # coverage of aligned bases
	ax2 = ax1.twinx()
	ax2.plot(svalues, color="red", linestyle="dotted")  # coverage of spanning reads
	ax2.plot([atg_pos, atg_pos], [0, max(svalues + values)], color="green", linestyle="dotted", label="ATG")  # ATG position
	ax2.plot([tss_pos, tss_pos], [0, max(svalues + values)], color="blue", linestyle="dotted", label="TSS")  # TSS position
	ax2.legend()

	ax1.set_title(gene + "   (" + orientation + ")")
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



def run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig, genome_seq ):
	"""! @brief run analysis on forward strand """
	most_upstream_pos = start
	final_pos_status = False
	while not final_pos_status:
		
		# --- walk coverage upstream of transcription start while there is coverage --- #
		while cov_per_contig[ most_upstream_pos-2 ] >= mincov:	#index = genomic position -1 but coverage of gene that is next to the current position needs to be assessed before moving in there
			most_upstream_pos -= 1	#move one step upstream
			if most_upstream_pos == hard_cutoff:	#stop if start of contig/pseudochromosome is reached
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
			print("current position is " + str( current_position ) )
			print ("most upstream position is " + str( most_upstream_pos ) )
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
	if most_upstream_pos > flank_region_for_plot:
		plot_start_region = most_upstream_pos - flank_region_for_plot
	else:
		plot_start_region = 0
	plot_end_region = max(start + flank_region_for_plot, atg_genomic_pos + flank_region_for_plot)
	
	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = atg_pos = atg_genomic_pos - plot_start_region
	tss_pos = most_upstream_pos - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "+"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end+1]
	try:
		generate_plot( values, svalues, fig_file, atg_pos, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot )
	except:
		print( "ERROR: plot failed" + gene )
		
	return { 'TSS': most_upstream_pos, 'start': start, 'end': end }


def run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig, genome_seq ):
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
	plot_start_region = min(end - flank_region_for_plot, atg_genomic_pos - flank_region_for_plot)
	if most_downstream_pos < ( len( seq_per_contig ) - flank_region_for_plot ):
		plot_end_region = most_downstream_pos + flank_region_for_plot
	else:
		plot_end_region = len( seq_per_contig )
		
	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = atg_genomic_pos - plot_start_region
	tss_pos = most_downstream_pos - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "-"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end + 1]
	try:
		generate_plot( values, svalues, fig_file, atg_pos, tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot )
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
	

def extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size ):
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
		"""
		if 'TATAAAA' in promoter.upper() or 'TATAAAT' in promoter.upper() or 'TATATAA' in promoter.upper() or 'TATATAT' in promoter.upper():
			tata_status = "TATA box found!"
		else:
			tata_status = 'No TATA box detected ...'
		"""
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
		"""
		if 'TTTTATA' in promoter.upper() or 'ATTTATA' in promoter.upper() or 'TTATATA' in promoter.upper() or 'ATATATA' in promoter.upper():#searching for reverse complements of TATA consensus sequences in the reverse strand genes
			tata_status = "TATA box found!"
		else:
			tata_status = 'No TATA box detected ...'
		"""
	return promoter_status, promoter

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

#function to scan extracted promoter sequences for TATA binding motifs
def promoter_motif_analysis (gene, orientation, promoter_seq, moods, pvalue, top_motifs, tss_prox, pfm_folder, tmp_folder, output_folder):
	tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', dir=tmp_folder, delete=False)
	try:
		tmp.write(f">{gene}\n{promoter_seq}\n")
		tmp.flush()
		tmp.close()
		moods_output_file = os.path.join(tmp_folder,f"{gene}_moods.txt")
		cmd = 'python3 ' + moods + ' -m ' + pfm_folder+'/*.pfm ' + '-s ' + tmp.name + ' -p ' + str(pvalue) + ' -o ' + moods_output_file
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
		rows = []
		with open (moods_output_file, 'r') as f:
			line = f.readline()
			while line:
				row = line.strip().rstrip(',').split(',')
				rows.append(row)
				line=f.readline()
		sorted_rows = sorted(rows,key=lambda row:sort_key(row,orientation))
		best_motif_hits=sorted_rows[:(top_motifs)]
		motif_closeness = None
		for hit in best_motif_hits:
			if orientation=='+':
				pos = int(hit[2])
				if len(promoter_seq) - pos <= tss_prox:
					motif_closeness = 'proximal motif'
				else:
					motif_closeness = 'distal motif'
			elif orientation=='-':
				pos=int(hit[2])
				if pos <= tss_prox:
					motif_closeness = 'proximal motif'
				else:
					motif_closeness = 'distal motif'
			hit.append(motif_closeness)
	finally:
		if os.path.exists(tmp.name):
			os.remove(tmp.name)
	return best_motif_hits

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

	if '--analyse_promoter' in arguments: #yes or no for promoter analysis with MOODS
		promoter_analysis = arguments[arguments.index('--analyse_promoter')+1]
	else:
		promoter_analysis = 'no'

	if '--moods' in arguments: #full path to moods python script
		moods = arguments[arguments.index('--moods')+1]
	else:
		moods = 'moods-dna.py'

	if promoter_analysis == 'yes':
		#full path to folder with PFM matrices
		pfm_folder = arguments[arguments.index('--PFM')+1]

	#p-value threshold for MOODS analysis
	if '--pval' in arguments:
		pvalue = float(arguments[arguments.index('--pval')+1])
	else:
		pvalue = 0.0001

	if '--top_motif_hit' in arguments:
		top_motifs = int(arguments[arguments.index('--top_motif_hit')+1])
	else:
		top_motifs = 1

	if '--tss_proximity' in arguments:
		tss_prox = int(arguments[arguments.index('--tss_proximity')+1])
	else:
		tss_prox = 30


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
							if 'Parent=' in each:
								transcript = each.replace('Parent=', '')
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
			cmd = star+' --runMode genomeGenerate --genomeDir '+star_indexing_folder+' --genomeFastaFiles '+fasta_file+' --sjdbGTFfile '+gff_file+' --sjdbGTFtagExonParentTranscript Parent --runThreadN '+t+' --genomeSAindexNbases '+index_bases
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

					if "_pass_1" in f:
						sample = f.split("_pass_1.")[0]
						pairs[sample]["R1"] = full_path

					elif "_pass_2" in f:
						sample = f.split("_pass_2.")[0]
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
						cmd = 'ulimit -n 4096 && '+star+' --runMode alignReads --genomeDir '+star_indexing_folder+' --outSAMtype BAM SortedByCoordinate --readFilesIn '+r1+' '+r2+' --runThreadN '+t+' --outFileNamePrefix '+prefix+' --readFilesCommand zcat --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --alignIntronMax '+intron_cutoff
						p = subprocess.Popen(args=cmd, shell=True)
						p.communicate()
					elif aligner == 'HISAT2':
						sorted_bam = os.path.join(star_mapping_folder, sample + '_sorted.bam')
						cmd = hisat2 + ' --max-intronlen ' + str(intron_cutoff) + ' -p ' + str(t) + ' -x ' + index_file_name + ' -1 ' + r1 + ' -2 ' + r2 + ' | ' + samtools + ' sort --threads ' + str(t) + ' -O BAM -o ' + sorted_bam
						p = subprocess.Popen(args=cmd, shell=True)
						p.communicate()
			#merging all the sorted BAM files obtained from STARlong mapping
			bam_files=os.path.join(output_folder,'bam_files.txt')
			if len(os.listdir(star_mapping_folder)) > 1:
				with open(bam_files, 'w') as out:
					for f in os.listdir(star_mapping_folder):
						if f.endswith('bam'):
							out.write(f + '\n')
				merged_bam=os.path.join(output_folder,sample+'_merged.bam')
				cmd = samtools+' merge -o '+merged_bam+' -b '+bam_files
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
	
	gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic = load_gene_infos( gff_file )
	
	genome_seq = load_sequences( fasta_file )
	
	#run analysis per gene of interest
	results = {}
	motifs = []
	#confidence_score_dic = {}
	isoforms_dic = {}
	distance_dic = {}
	coverage_dic = {}
	five_utr_dic = {}
	for gene in goi:
		try:
			cov_per_contig = coverage[ gene_infos[ gene ]['chromosome'] ]	#get coverage of the sequence that harbours the gene of interest
			scov_per_contig = scoverage[ gene_infos[ gene ]['chromosome'] ]	#get spanning read coverage of the sequence that harbours the gene of interest
			seq_per_contig = genome_seq[ gene_infos[ gene ]['chromosome'] ]	#get the sequence of the contig/pseudochromosome that harbours the gene of interest
			#code block to check if the goi has annotated 5'UTR and if yes take the most upstream/ downstream 5'UTR start/ end as start or end according to + or - strand orientation
			# Initialize with gene coordinates as default with 5'UTR checks downstream
			start, end, orientation, contig = gene_infos[gene]['start'], gene_infos[gene]['end'], gene_infos[gene]['orientation'], gene_infos[ gene ]['chromosome'] # get information about gene of interest if it does not have 5'UTRs annotated
			# store ATG position as fixed reference before start may be modified by 5'UTR code block
			if gene in gene_atg_dic:
				atg_pos = gene_atg_dic[gene]
			else:
				atg_pos = start if orientation == '+' else end  # fallback to gene boundary if no CDS found
			if gene in transcripts_per_gene:
				transcript_list = transcripts_per_gene[ gene ]
				for each in transcript_list:
					if each in five_utr_infos and gene_infos[ gene ]['orientation']=='+' and each==transcript_list[0]:
						start=five_utr_infos[each]['start'] #in case a gene has transcripts with 5'UTR annotated, the most upstream 5'UTR will be taken as the start for the + strand gene
						end, orientation = gene_infos[gene]['end'], gene_infos[gene]['orientation']  # get information about gene of interest
						five_utr_dic[gene]=f"5'UTR start of {each} used for TSS prediction."
						break
					elif each in five_utr_infos and gene_infos[ gene ]['orientation']=='-' and each==transcript_list[-1]:
						end = five_utr_infos[each]['end']  # in case a gene has transcripts with 5'UTR annotated, the most downstream 5'UTR will be taken as the start for the - strand gene
						start, orientation = gene_infos[gene]['start'], gene_infos[gene]['orientation']  # get information about gene of interest
						five_utr_dic[gene] = f"5'UTR end of {each} used for TSS prediction."
						break
			if gene not in five_utr_dic and orientation == '+':
				five_utr_dic[gene] = f"No 5'UTR annotated. {gene} start used for TSS prediction."
			if gene not in five_utr_dic and orientation == '-':
				five_utr_dic[gene] = f"No 5'UTR annotated. {gene} end used for TSS prediction."
			upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes( gene, gene_infos, genes_per_chromosome, window )
			avg_cov_gene = sum(cov_per_contig[start-1:end])/(end-(start-1)) #get average coverage of the gene of interest for confidence thresholding
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
				result = run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq )
				promoter_status, promoter = extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size )
				if promoter_analysis == 'yes':
					best_motif_hits = promoter_motif_analysis(gene, orientation, promoter, moods, pvalue,top_motifs, tss_prox, pfm_folder, tmp_folder, output_folder)
					motifs.append(best_motif_hits)
				result.update( { 'promoter_status': promoter_status } )
				result.update( { 'promoter': promoter } )
				results.update( { gene: result } )
			else:	#solution for reverse strand genes
				fig_file = output_folder + gene + ".png"
				if downstream_gene:
					hard_cutoff = gene_infos[ downstream_gene ]['start']
				else:
					hard_cutoff = len( seq_per_contig )
				result = run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq )
				promoter_status, promoter = extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size )
				if promoter_analysis == 'yes':
					best_motif_hits = promoter_motif_analysis(gene, orientation, promoter, moods, pvalue, top_motifs,tss_prox, pfm_folder, tmp_folder, output_folder)
					motifs.append(best_motif_hits)
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
		except KeyError:
			print( "Missing gene error: " + gene )
			"""
			# If gene was added to results but confidence calculation failed, add fall back
			if gene in results and gene not in confidence_score_dic:
				confidence_score_dic[gene] = "NA"
			"""
	# --- report TSS in output file --- #
	final_output_file = output_folder + "Results.txt"
	promoter_motif_output_file = os.path.join(output_folder,'Top_promoter_motifs.txt')
	with open( final_output_file, "w" ) as out:
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
	# --- write top promoter motif hits in a separate output file --- #
	if promoter_analysis == 'yes':
		if motifs:
			with open (promoter_motif_output_file,'w')as out:
				for hit_list in motifs:
					for hits in hit_list:
						out.write('\t'.join(hits)+'\n')
		else:
			print('No valid top motifs list found. Promoter motif analysis file not written.')


if '--bam' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--scov' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--run_mode' in sys.argv and '--sra_folder' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main(sys.argv)
else:
	sys.exit( __usage__ )
