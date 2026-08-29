### Boas Pucker ###
### Shakunthala Natarajan ###
### pucker@uni-bonn.de ###
### s64snata@uni-bonn.de ###
__version__ = "v0.1.0"

__reference__ = "https://github.com/bpucker/eGPS"

__usage__ = """
					eGPS """ + __version__ + """("""+ __reference__ +""")
					
					Usage:
					python3 eGPS.py
					--fasta <GENOMIC_FASTA_FILE>
					--gff <GFF_FILE>
					--goi <TXT_FILE_WITH_LIST_OF_GENES_OF_INTEREST_ONE_PER_LINE>
					--out <OUTPUT_FOLDER>
					[--sra_folder <READ_FILES> | --bam <BAM_FILE|BAM_FOLDER>|--cov <COV_FILE|COV_FOLDER> --scov <SCOV_FILE|SCOV_FOLDER>]
					
					optional:
					--gff_config <TXT_CONFIG_FILE_WITH_OPTIONS_TO_DEFINE_THE_GFF_ATTRIBUTE_FIELDS>
					--protein_encoding <YES OR NO TO CONSIDER ONLY PROTEIN ENCODING GENES FOR TSS ANALYSIS> default is no
					--compare_tss <SPECIFY NAME OF THE ALTERNATE TSS SOURCE TO COMPARE AGAINST EGPS>
					--compare_tss_input <FULL PATH TO TWO COLUMN TAB SEPARATED CONFIG FILE WHERE COLUMN ONE IS GOI AND COLUMN TWO IS THE TSS FOR IT FROM THE SOURCE BEING COMPARED AGAINST>
					--sample_support <yes or no FOR SAMPLE SUPPORT-BASED TSS SCORING> default is no
					--coverage_difference <THRESHOLD TO BE USED AS PERCENTAGE DIFFERENCE THRESHOLD BETWEEN ALIGNED AND SPANNING READ COVERAGE> default is 10%
					--coverage_walk_origin <cds or utr> DEFAULT is cds
					--mincov <MINIMAL_COVERAGE>[1]
					--bam_is_sorted <PREVENTS_BAM_FILE_SORTING>
					--samtools <FULL_PATH_TO_SAMTOOLS>[samtools]
					--bedtools <FULL_PATH_TO_genomeCoverageBed>[genomeCoverageBed]
					--m <MEM_FOR_SAMTOOLS_SORTING>[5000000000]
					--threads <NUMBER_THREADS_FOR_SAMTOOLS_SORTING>[4]
					--minexon <MINIMAL_EXON_SIZE>[10]
					--flanksize <FLANKING_REGION_SIZE_FOR_COVERAGE_PLOT_VISUALIZATION>[50]
					--gapsize <COVERAGE_GAP_SIZE>[5]
					--splicesites <HANDLING_OF_SPLICE_SITES>[strict](strict|off)
					--intron_percentile_cutoff <INTRON_SIZE_CUTOFF_FOR_STAR+HISAT2_ALIGNMENT>
					--neighbourhood <GENE_NEIGHBOURHOOD_WINDOW_FOR_OVERLAPPING_GENE_ANALYSIS>[5]
					--min_promoter_size <MINIMUM_PROMOTER_SIZE>[50]
					--max_promoter_size <MAXIMUM_PROMOTER_SIZE>[1000]
					--background < NUMBER_OF_RANDOM_BACKGROUND_SEQUENCES_TO_BE_CONSIDERED_FOR_MOTIF_SCORING_IN_PROMOTER_ANALYSIS>[1000]
					--downstream_size <LENGTH_OF_ANALYSIS_REGION_DOWNSTREAM_TO_TSS>[300]
					--upstream_slice <LENGTH_OF_PROMOTER_REGION_TO_BE_SLICED_FOR_MOTIF_ANALYSIS>[200]
					--downstream_slice <LENGTH_OF_DOWNSTREAM_REGION_TO_BE_SLICED_FOR_MOTIF_ANALYSIS>[50]
					--aligner <STAR or HISAT2> STAR is default
					--HISAT2 <FULL_PATH_TO_HISAT2_FOR_RNA_Seq_MAPPING>
					--STAR <FULL_PATH_TO_STAR_FOR_RNA_Seq_MAPPING>
					--index_bases <PARAMETER_FOR_GENOME_INDEX_GENERATION_IN_STAR>
					--run_mode < find_tss or make_bam> find_tss is default
					--fastq_pattern <FASTQ_FILE_NAME_PATTERN_SEPARATED_BY_COMMA><EG: pass1,pass2 >
					--analyse_promoter <yes or no> [no]
					--moods <FULL_PATH_TO_MOODS_SCRIPT>
					--PFM <TXT_FILE_FOR_PROMOTER_MOTIF_ANALYSIS; TAB_SEPARATED; MOTIF_NAME	PATH_TO_PFM_FILE	UPSTREAM_BOUNDARY	DOWNSTREAM_BOUNDARY	DIRECTION_SENSITIVITY>
					--moods_pval <MOODS_MOTIF_ANALYSIS_THRESHOLD>
					--background_percentage <PERCENTAGE_OF_BACKGROUND_SEQUENCES_BELOW_THE_AVG_COVERAGE_PER_ANALYSIS_SLIDING_WINDOW_TO_DETERMINE_BASAL_ELEVATED_TRANSCRIPTION_REGIONS>[0.05]
					--ks_pval <P-VALUE_OF_KOLOMOGOROV-SMIRNOV_DISTRIBUTION_FIT_TEST_USED_FOR_aCCELERATED_TSS_FINDING>[0.01]
					--background_unit <SIZE_OF_INTERGENIC_INTRONIC_EXONIC_WINDOWS_USED_IN_SLIDING_WINDOW_ANALYSIS_TO_EVALUATE_AVERAGE_COVERAGE_PER_WINDOW_FOR_RNA-SEQ_BACKGROUND_BASELINE_ANALYSIS>[10]
					--intron_trim <NUMBER_OF_BASES_TO_TRIM_FROM_INTRON_ENDS_TO_AVOID_SPLICE_JUNCTION_INCLUSION_IN_INTRON_COVERAGE_CALCULATION> [5]
					--slide <NUMBER_OF_SLIDING_WINDOW_PROGRESSION_INTERVAL_FOR_BACKGROUND_BASELINE_ANALYSIS_AND_DISTINCTION>[1]
					--signal_strength <NUMBER_OF_CONSECUTIVE_WINDOWS_WITH_MEAN_COVERAGE_VALUES_ABOVE_FIXED_PERCENTAGE_OF_BACKGROUND_UNITS>[3]
					--lookahead <NUMBER_OF_BASES_OR_POSITIONS_TO_BE_LOOKED_AHEAD_TO_DETERMINE_THE_MINIMA_POINT_CORRESPONDING_TO_THE_STEEPEST_POSITIVE_TRANSITION_IN_COVERAGE_FOR_ACCELERATED_TSS_DETECTION>
					--buffer <SIZE_OF_BUFFER_REGIONS_TO_BE_IGNORED_FOR_SLICING_INTERGENIC_REGIONS_TO_AVOID_READTHROUGH_SIGNAL_INTERRIPTIONS_FROM_ACTUAL_EXONIC_REGIONS>
					"""


import re, os, sys, subprocess, gzip
from pathlib import Path
import bisect
import random
import math
import tempfile
import time
import traceback
import numpy as np
from scipy import stats
from scipy.stats import percentileofscore
import seaborn as sns
from decimal import Decimal, ROUND_HALF_DOWN
from collections import defaultdict
import shutil
from build_offsets import run_build_offsets_parallel
from sum_chunk_seek import run_sum_chunk_seek_parallel
from concatenate_chunks import concatenate_chunks
from parallel_coverage_generation import parallelize_coverage_file_generation
try:
	import matplotlib.pyplot as plt
	from matplotlib.lines import Line2D
	import matplotlib.ticker as ticker
	from matplotlib.patches import FancyArrow
	from statistics import mode
	from copy import deepcopy
except ImportError:
	pass

# --- end of imports --- #

def construct_coverage_files ( bam_list, bam_files, read_coverage_folder, bedtools, parallelize_cov_generation, parallel, cores, tss_scoring, coverage_type ):
	""" @brief calculate read coverage depth per position """

	if parallelize_cov_generation:
		parallelize_coverage_file_generation(bam_files, cores, read_coverage_folder, parallel, bedtools, coverage_type)
		offset_dir = os.path.join(read_coverage_folder, 'Offsets')
		if not os.path.exists(offset_dir):
			os.makedirs(offset_dir)
		chunk_dir = os.path.join(read_coverage_folder, 'Cov_chunks')
		if not os.path.exists(chunk_dir):
			os.makedirs(chunk_dir)
		cov_files_list = [
			  os.path.join(read_coverage_folder, f)
						for f in os.listdir(read_coverage_folder)
						if f.endswith('.cov')
					]
		single_cov_file = cov_files_list[0]
		cmd = f'wc -l {single_cov_file}'#compute number of lines in a cov file to determine size of a chunk
		p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True)
		output, _ = p.communicate()
		n = int(output.split()[0])
		chunk_size = ((n + cores )-1)//cores

		#build offsets to mark the indices of the chunk starts for parallelized chunk sum approach via seek
		run_build_offsets_parallel(cov_files_list, chunk_size, cores, offset_dir, parallel)
		"""
		take a chunk from all cov files parallely and run the summing as n chunk sum processes where n is the number of cores;
		eg for 100 cov files, and 28 cores; split each file into 28 chunks. Now take 100 chunks across each of the 100 cov files and sum them;
		do this summation as 28 parallel processes where each process has 100 chunks to be summed across
		this chunk per sample and sum across samples leveraging parallelism is implemented to speed up cumulative cov file generation
		"""
		run_sum_chunk_seek_parallel(read_coverage_folder, offset_dir, cores, chunk_size, chunk_dir, parallel)
		if coverage_type == 'aligned':
			out_path = os.path.join(read_coverage_folder,"Cumulative_aligned_reads.cov")
		elif coverage_type == 'spanning':
			out_path = os.path.join(read_coverage_folder,"Cumulative_spanning_reads.cov")
		coverage_output_file = concatenate_chunks(chunk_dir, out_path)
		if tss_scoring:
			pass
		else:
			shutil.rmtree(offset_dir)
			shutil.rmtree(chunk_dir)

	else:
		if coverage_type == 'aligned':
			coverage_output_file = os.path.join(read_coverage_folder, 'Reads_aligned.cov')
			print("calculating coverage per position (aligned) ....")
			cmd = bedtools + " -d -split -ibam " + bam_list[0] + " > " + coverage_output_file  # -include spanning reads when calculating depth
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
		elif coverage_type == 'spanning':
			coverage_output_file = os.path.join(read_coverage_folder, 'Reads_spanning.cov')
			print("calculating coverage per position (spanning) ....")
			cmd = bedtools + " -d -ibam " + bam_list[0] + " > " + coverage_output_file  # -include spanning reads when calculating depth
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()

	return coverage_output_file

def sum_coverage_files ( cov_files_list, input_read_coverage_folder, read_coverage_folder, parallel, cores, tss_scoring, coverage_type ):
	""" @brief calculate read coverage depth per position """

	offset_dir = os.path.join(read_coverage_folder, 'Offsets')
	if not os.path.exists(offset_dir):
		os.makedirs(offset_dir)
	chunk_dir = os.path.join(read_coverage_folder, 'Cov_chunks')
	if not os.path.exists(chunk_dir):
		os.makedirs(chunk_dir)
	single_cov_file = cov_files_list[0]
	cmd = f'wc -l {single_cov_file}'#compute number of lines in a cov file to determine size of a chunk
	p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True)
	output, _ = p.communicate()
	n = int(output.split()[0])
	chunk_size = ((n + cores )-1)//cores

	#build offsets to mark the indices of the chunk starts for parallelized chunk sum approach via seek
	run_build_offsets_parallel(cov_files_list, chunk_size, cores, offset_dir, parallel)
	"""
	take a chunk from all cov files parallely and run the summing as n chunk sum processes where n is the number of cores;
	eg for 100 cov files, and 28 cores; split each file into 28 chunks. Now take 100 chunks across each of the 100 cov files and sum them;
	do this summation as 28 parallel processes where each process has 100 chunks to be summed across
	this chunk per sample and sum across samples leveraging parallelism is implemented to speed up cumulative cov file generation
	"""
	run_sum_chunk_seek_parallel(input_read_coverage_folder, offset_dir, cores, chunk_size, chunk_dir, parallel)
	if coverage_type == 'aligned':
		out_path = os.path.join(read_coverage_folder,"Cumulative_aligned_reads.cov")
	elif coverage_type == 'spanning':
		out_path = os.path.join(read_coverage_folder,"Cumulative_spanning_reads.cov")
	coverage_output_file = concatenate_chunks(chunk_dir, out_path)
	if tss_scoring:
		pass
	else:
		shutil.rmtree(offset_dir)
		shutil.rmtree(chunk_dir)

	return coverage_output_file


def load_coverage( cov_file):
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


def load_gene_infos( protein_encoding, gff_file, child_attribute, child_parent_linker, parent_attribute):
	"""! @brief load gene ID, position, and orientation from GFF3 file """
	
	gene_infos = {}
	mrna_infos = {}
	five_utr_infos={}
	cds_infos = {}
	genes_per_chromosome = {}
	protein_coding_genes_per_chromosome = defaultdict(list)
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
					utr_parents = parts[-1].split(f'{child_parent_linker}=')[-1]#Parent of 5'UTR is transcript
					if ";" in utr_parents:
						utr_parents = utr_parents.split(';')[0]
					for utr_parent in utr_parents.split(','): #handle multiple parents
						utr_parent = utr_parent.strip()
						utr_tuple = (int(parts[3]), int(parts[4]))
						try:
							five_utr_infos[ utr_parent ].append(utr_tuple)
						except KeyError:
							five_utr_infos[ utr_parent ] = [utr_tuple]
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
	if protein_encoding == "yes":
		for chromosome, genes in genes_per_chromosome.items():
			for gene in genes:
				if gene in transcripts_per_gene:
					protein_coding_genes_per_chromosome[chromosome].append(gene)
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
	return gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos, protein_coding_genes_per_chromosome

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

#function to get intergenic background sequence coverages
def get_intergenic_background_seqs (outdir, coverage_dic, intergenic_buffer, intergenic_region_size, genes_per_chromosome, gene_infos):
	intergenic_window_coverages = []
	individual_pos_coverage = []
	nonzero_positions = []

	for contig in genes_per_chromosome:
		coverage_lookup = dict(coverage_dic[contig])
		gene_list = genes_per_chromosome[contig]
		for i, gene in enumerate(gene_list):
			if i + 1 >= len(gene_list):  #skip the last gene in the contig
				break
			next_gene = gene_list[i + 1]
			intergenic_start = gene_infos[gene]['end']
			intergenic_end = gene_infos[next_gene]['start']
			full_intergenic_coverage_slice = [(pos, coverage_lookup[pos]) for pos in range(intergenic_start, intergenic_end + 1) if pos in coverage_lookup]
			if len(full_intergenic_coverage_slice) > intergenic_buffer and len(full_intergenic_coverage_slice) > intergenic_region_size:
				sliced_intergenic_coverage = full_intergenic_coverage_slice[intergenic_buffer:(len(full_intergenic_coverage_slice)-intergenic_buffer)]#leave a buffer near the gene regions for both the genes
				chops = int(len(sliced_intergenic_coverage)/intergenic_region_size)#chop the intergenic region into bits of the same size so that a window of similar size can be compared in the sliding window strategy
				# Trim to only the evenly divisible portion before chunking to avoid inclusion of partial final chunks
				trimmed = sliced_intergenic_coverage[:chops * intergenic_region_size]
				chop_start = 0
				while chops > 0:
					chop_region = trimmed[chop_start:chop_start+intergenic_region_size]
					cumulative_coverage = sum(cov for _, cov in chop_region)
					for pos, cov in chop_region:
						individual_pos_coverage.append(cov)
						if cov >0:
							nonzero_positions.append(pos)
					avg_intergenic_window_coverage = cumulative_coverage / len(chop_region)
					if avg_intergenic_window_coverage != 0.0:#adding only non-zero average coverage values to the background noise
						intergenic_window_coverages.append(avg_intergenic_window_coverage)
					chop_start += intergenic_region_size
					chops -= 1
	# compute distances between consecutive nonzero positions
	gaps = [nonzero_positions[i + 1] - nonzero_positions[i]
			for i in range(len(nonzero_positions) - 1)]
	print(f"Mean gap:   {np.mean(gaps):.1f}bp")
	print(f"Median gap: {np.median(gaps):.1f}bp")
	print(f"Std gap:    {np.std(gaps):.1f}bp")
	modal_background_cov = mode(individual_pos_coverage)
	print(f'Modal background coverage is {modal_background_cov}')
	non_zero_cov_points = []
	zero_cov_points = []
	for each in individual_pos_coverage:
		if each == 0:
			zero_cov_points.append(each)
		else:
			non_zero_cov_points.append(each)
	percent_zero = (len(zero_cov_points) / len(individual_pos_coverage))*100
	percent_nonzero = (len(non_zero_cov_points) / len(individual_pos_coverage))*100
	print(f"Percentage of intergenic positions with zero coverage is {percent_zero}")
	print(f"Percentage of intergenic positions with nonzero coverage is {percent_nonzero}")
	cov_arr = np.array(intergenic_window_coverages)
	q1, q3 = np.percentile(cov_arr, [25, 75])
	iqr = q3 - q1
	print(f"IQR effect analysis with intergenic buffer of {intergenic_buffer}")
	intergenic_window_coverages_wo_outliers = cov_arr[cov_arr <= q3 + 1.5*iqr]
	print(f"Total intergenic windows:        {len(intergenic_window_coverages)}")
	print(f"Windows after IQR filter:        {len(intergenic_window_coverages_wo_outliers)}")
	print(f"Windows removed by IQR:          {len(intergenic_window_coverages) - len(intergenic_window_coverages_wo_outliers)}")
	print(f"Retention rate:                  {100 * len(intergenic_window_coverages_wo_outliers) / len(intergenic_window_coverages):.1f}%")
	intergenic_array = np.array(intergenic_window_coverages_wo_outliers)
	percentile_99 = np.percentile(intergenic_array, 99)

	#plotting the distribution of intergenic window coverages with outliers
	plt.figure(figsize=(8,5))
	sns.histplot(intergenic_window_coverages, kde=True)
	plt.xlim(0, np.percentile(intergenic_window_coverages, 99))
	plt.xlabel('Average_coverage')
	plt.ylabel('Frequency')
	plt.title('Distribution of average coverage values in intergenic region')
	plt.tight_layout()
	fig_save = os.path.join(outdir,'distribution_intergenic_with_outliers.png')
	plt.savefig(fig_save, dpi=600)

	# plotting the distribution of intergenic window coverages without outliers
	plt.figure(figsize=(8, 5))
	sns.histplot(intergenic_window_coverages_wo_outliers, kde=True)
	plt.xlim(0, np.percentile(intergenic_window_coverages_wo_outliers, 99))
	plt.xlabel('Average_coverage')
	plt.xlabel('Average_coverage')
	plt.ylabel('Frequency')
	plt.title('Distribution of average coverage values in intergenic region')
	plt.tight_layout()
	fig_save = os.path.join(outdir, 'distribution_intergenic_without_outliers.png')
	plt.savefig(fig_save, dpi=600)
	z_ig = percent_zero / 100  # fraction of zero cov positions in the intergenic region
	return intergenic_window_coverages_wo_outliers, z_ig, percent_zero, percent_nonzero

def generate_plot( other_tool, other_tool_tss_pos, values, svalues, fig_file, atg_pos, cov_walk_start, boundary_tss_for_plot, basal_tss_pos, elevated_tss_pos, accelerated_tss_pos, genomic_start, genomic_end, gene, orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot, introns_to_plot ):
	"""! @brief generate a coverage plot """
	fig, (ax1, ax_features) = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios':[4,1.5]}, sharex=True, layout = 'constrained')
	ax1.plot(values, color="black", linestyle="solid")  # coverage of aligned bases
	ax2 = ax1.twinx()
	ax2.plot(svalues, color="red", linestyle="dotted")  # coverage of spanning reads
	ax2.plot([atg_pos, atg_pos], [0, max(svalues + values)], color="green", linestyle="dotted", label="ATG")  # ATG position
	ax2.plot([cov_walk_start, cov_walk_start], [0, max(svalues + values)], color="orange", linestyle="dotted", label="Coverage walk origin")  # 5'UTR start or end or gene start or end position depending on strandedness and 5'UTR annotation being present for the gene's most upstream or downstream transcripts
	if basal_tss_pos is not None:
		ax2.plot([basal_tss_pos, basal_tss_pos], [0, max(svalues + values)], color="brown", linestyle="dotted", label="Basal TSS")  # basal TSS position
	if elevated_tss_pos is not None:
		ax2.plot([elevated_tss_pos, elevated_tss_pos], [0, max(svalues + values)], color="blue", linestyle="dotted", label="Elevated TSS")  # elevated TSS position
	if accelerated_tss_pos is not None:
		ax2.plot([accelerated_tss_pos, accelerated_tss_pos], [0, max(svalues + values)], color="pink", linestyle="dotted",label="Accelerated TSS")  # accelerated TSS position
	if other_tool_tss_pos is not None:
		ax2.plot([other_tool_tss_pos, other_tool_tss_pos], [0, max(svalues + values)], color="red", linestyle="dotted",label=f"{other_tool} TSS")  # other tool TSS position
	handles, labels = ax2.get_legend_handles_labels()
	fig.legend(handles, labels, loc = 'outside upper center', ncol = len(labels), fontsize = 10, frameon = False)
	ax1.set_title(gene + " (" + orientation +")",fontsize=12)
	# feature track plotting with features represented as arrows
	for mrna_start, mrna_end in mrnas_to_plot:
		feature_length = mrna_end - mrna_start
		if feature_length != 0:
			if orientation == '+':
				dx = feature_length
				x_origin = mrna_start
				xdot_origin = boundary_tss_for_plot
				predicted_tss_mrna_feature = mrna_end - boundary_tss_for_plot
			elif orientation == '-':
				dx = -feature_length
				x_origin = mrna_end
				xdot_origin = boundary_tss_for_plot
				predicted_tss_mrna_feature = -(boundary_tss_for_plot - mrna_start)
			head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
			arrow = FancyArrow(x=x_origin, y=0.2, dx=dx, dy=0, width=0.2, head_width=0.3, head_length=head_length, length_includes_head=True, facecolor='steelblue', alpha=0.5, edgecolor='black', linewidth=0.5)
			ax_features.add_patch(arrow)
			"""
			if predicted_tss_mrna_feature != 0:
				dotted_head_length = min(50, (abs(predicted_tss_mrna_feature)*0.2))
				dotted_arrow = FancyArrow(x=xdot_origin, y=0.2, dx=predicted_tss_mrna_feature, dy=0, width=0.2, head_width=0.3, head_length=dotted_head_length , length_includes_head=True, facecolor='none', alpha=0.5, edgecolor='steelblue', linewidth=0.5, linestyle='dotted')
				ax_features.add_patch(dotted_arrow)
			"""
	for cds_start, cds_end in cds_to_plot:
		feature_length = cds_end - cds_start
		if feature_length != 0:
			if orientation == '+':
				dx = feature_length
				x_origin = cds_start
			elif orientation == '-':
				dx = -feature_length
				x_origin = cds_end
			head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
			arrow = FancyArrow(x=x_origin, y=0.5, dx=dx, dy=0, width=0.2, head_width=0.3, head_length=head_length, length_includes_head=True, facecolor='lightgreen', alpha=0.5, edgecolor='black', linewidth=0.5)
			ax_features.add_patch(arrow)
	for five_utr_start, five_utr_end in five_utr_to_plot:
		feature_length = five_utr_end - five_utr_start
		if feature_length != 0:
			if orientation == '+':
				dx = feature_length
				x_origin = five_utr_start
				xdot_origin = boundary_tss_for_plot
				predicted_tss_five_utr_feature = five_utr_end - boundary_tss_for_plot
			elif orientation == '-':
				dx = -feature_length
				x_origin = five_utr_end
				xdot_origin = boundary_tss_for_plot
				predicted_tss_five_utr_feature = -(boundary_tss_for_plot - five_utr_start)
			head_length = min(50, feature_length * 0.2)  # cap arrowhead at 20% of feature width
			arrow = FancyArrow(x=x_origin, y=0.8, dx=dx, dy=0, width=0.2, head_width=0.3, head_length = head_length, length_includes_head=True, facecolor='salmon', alpha=0.5, edgecolor='black', linewidth=0.5)
			ax_features.add_patch(arrow)
			"""
			if predicted_tss_five_utr_feature != 0:
				dotted_head_length = min(50, (abs(predicted_tss_five_utr_feature) * 0.2))
				dotted_arrow = FancyArrow(x=xdot_origin, y=0.8, dx=predicted_tss_five_utr_feature, dy=0, width=0.2, head_width=0.3,head_length=dotted_head_length, length_includes_head=True, facecolor='none', alpha=0.5,edgecolor='salmon', linewidth=0.5, linestyle='dotted')
				ax_features.add_patch(dotted_arrow)
			"""
	if introns_to_plot:
		for intron_start, intron_end in introns_to_plot:
			ax_features.plot(
				[intron_start, intron_end], [1.1, 1.1],
				color='black', linewidth=2, solid_capstyle='butt', alpha=0.7, zorder=1
			)

	# extend vertical lines into feature axis
	ax_features.axvline(atg_pos, color="green", linestyle="dotted")
	ax_features.axvline(cov_walk_start, color="orange", linestyle="dotted")
	colour = 'black'
	if boundary_tss_for_plot == basal_tss_pos:
		colour = 'brown'
	elif boundary_tss_for_plot == elevated_tss_pos:
		colour = 'blue'
	elif boundary_tss_for_plot == accelerated_tss_pos:
		colour = 'pink'
	ax_features.axvline(boundary_tss_for_plot, color=colour, linestyle="dotted")


	if introns_to_plot:
		ax_features.set_yticks([0.2, 0.5, 0.8, 1.1])
		ax_features.set_yticklabels(['mRNA', 'CDS', "5'UTR", "Intron"], fontsize=8)
	else:
		ax_features.set_yticks([0.2, 0.5, 0.8])
		ax_features.set_yticklabels(['mRNA', 'CDS', "5'UTR"], fontsize=8)
	ax1.set_xlabel("Position in genomic region from " + str(genomic_start) + " to " + str(genomic_end))
	ax1.set_ylabel("Aligned RNA-seq coverage")
	ax1.yaxis.label.set_color('black')
	ax2.set_ylabel("Spanning RNA-seq coverage")
	ax2.yaxis.label.set_color('red')

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

#function to find the gene expression level of a gene of interest with respect to the background
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

#function to perform the kolmogorov smirnov test
def do_ks_test(normalized_cumulative_sum_array):
	cov_sum0 = normalized_cumulative_sum_array[0] #y0
	cov_sum_max = max(normalized_cumulative_sum_array) #ymax
	cov_len = normalized_cumulative_sum_array.size - 1 #xmax
	slope = (cov_sum_max - cov_sum0) / cov_len #((ymax-y0)/xmax)
	if slope == 0:
		pval_ks = 1
		diagonal_line_array = np.array([])
	else:
		diagonal_line_array = np.linspace(0, normalized_cumulative_sum_array.max(), normalized_cumulative_sum_array.size)#the uniform diagonal baseline with the ideal slope
		stat_ks, pval_ks = stats.ks_2samp(normalized_cumulative_sum_array, diagonal_line_array)
	return pval_ks, diagonal_line_array

#function to compute accelerated tss
def get_accelerated_tss(ks_pval_threshold, gene, coverage_slice_list_for_accelerated_tss_arr, lookahead, coverage_slice_list_for_accelerated_tss,strand):
	accelerated_tss = None
	if len(coverage_slice_list_for_accelerated_tss) == 0:
		accelerated_tss = None
	else:
		consecutive_cov_differences = np.insert(np.ediff1d(coverage_slice_list_for_accelerated_tss_arr), 0, 0)#impute the first value as 0 since the difference calculation reduces the length of the original array by 1
		cumulative_sum_of_cov_diff = np.cumsum(np.abs(consecutive_cov_differences))
		if cumulative_sum_of_cov_diff.max() == 0:
			crs = np.array([])
		else:
			crs = cumulative_sum_of_cov_diff/cumulative_sum_of_cov_diff.max()
			pval_ks, diagonal_line_array = do_ks_test(crs)
			if pval_ks < ks_pval_threshold:
				print(f'KS test passed for {gene}. Searching for accelerated TSS.')
				i=0
				index_mean_dic={}
				window_means=[]
				while i <= (len(consecutive_cov_differences) - lookahead):
					array_window = consecutive_cov_differences[i:i+lookahead]
					mean_diff_per_window = np.mean(array_window)
					index_mean_dic[i]=mean_diff_per_window
					window_means.append(mean_diff_per_window)
					i+=1

				window_means_arr=np.array(window_means)
				if window_means_arr.size > 0:#if (len(consecutive_cov_differences) - lookahead) < 0 then window_means_arr turns out empty. this could happen for cases wherein the elevated tss is already close to the ATG. So skip accelerated TSS detection for such cases.
					max_window_mean=np.max(window_means_arr)
					for i in index_mean_dic:
						if index_mean_dic[i] == max_window_mean:
							selected_index=i
							break
					target_window_for_steepest_rise_point = np.array(consecutive_cov_differences[selected_index:selected_index+lookahead])
					steepest_rise_point = np.max(target_window_for_steepest_rise_point)
					consecutive_cov_differences_list = consecutive_cov_differences.tolist()
					index_steepest_rise_point = consecutive_cov_differences_list.index(steepest_rise_point)
					positions = np.array([pos for pos, cov in coverage_slice_list_for_accelerated_tss])
					accelerated_tss = positions[index_steepest_rise_point]#the genomic position corresponding to the point of steepest increase
	return accelerated_tss

def run_fwd_analysis(other_tool, other_tss_dic,ks_pval, strength, lookahead, output_folder, pvalue,
					 intergenic_region_size, slide_step, intergenic_window_coverages, coverage_lookup, gene,
					 cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size,
					 hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig,
					 genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos,
					 gene_atg_dic, cds_infos,percentdiff_threshold):
	"""! @brief run analysis on forward strand """
	if gene in other_tss_dic:
		other_tool_tss = other_tss_dic[gene]
	else:
		other_tool_tss = None

	hard_cutoff_reached = False
	intron_boundary_marker = {}#dictionary to collect intron boundary positions where key is donor splice site and value is acceptor splice site
	most_upstream_pos = start
	final_pos_status = False
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
					pass
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
			while cov_per_contig[ current_position - 2 ] < mincov:#check if upstream position has low coverage
				current_position -= 1
				if current_position == hard_cutoff:
					hard_cutoff_reached = 'True;'+str(hard_cutoff)
					break
			avg_gap_coverage = sum( scov_per_contig[ current_position-1:most_upstream_pos-1 ] )/(most_upstream_pos-current_position)
			#average coverage in intron should be very low
			if current_position > min_exon_size and avg_gap_coverage > mincov:
			# --- check coverage gaps for (canonical) splice sites to continue across introns --- #
				donor_splice_site = seq_per_contig[current_position - 1:current_position + 1].upper()  # this should be GT for it to be canonical splice site
				acceptor_splice_site = seq_per_contig[most_upstream_pos - 2:most_upstream_pos].upper()  # this should be AG for it to be canonical splice site; changed the slicing from -3 and -1 respectively to -2 and nothing as only then the intron acceptor boundaries are correctly considered inclusively
				if splicesites == "off":  # ignore check for canonical splice sites #tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
					if most_upstream_pos - current_position < tolerated_gap:
						intron_boundary_marker[intron_donor_splice_site_position] = intron_acceptor_splice_site_position
						most_upstream_pos = current_position - 1
					else:#if splice sites off and tolerated gap check fails, make final pos status true; making the check stringent for the loop with cov less than mincov
						final_pos_status = True
				elif donor_splice_site == "GT" and acceptor_splice_site == "AG":
					print(f"Canonical donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Canonical acceptor splice site ending at {most_upstream_pos}: " + acceptor_splice_site + '\n')
					intron_boundary_marker[current_position] = most_upstream_pos
					most_upstream_pos = current_position - 1
				elif donor_splice_site != "GT" or acceptor_splice_site != "AG":#tolerate gap check is activated only when non-canonical sites appear or when splicesites is off
					print(f"Warning! Either one or both the donor and acceptor splice sites are non-canonical."+ '\n')
					print(f"Donor splice site starting at {current_position}: " + donor_splice_site+ '\n')
					print(f"Acceptor splice site ending at {most_upstream_pos}: " + acceptor_splice_site+ '\n')
					if most_upstream_pos - current_position < tolerated_gap:
						intron_boundary_marker[current_position] = most_upstream_pos
						most_upstream_pos = current_position - 1
					else:#if non canonical splice sites show up and tolerated gap check fails, make final pos status true; making the check stringent for the loop with cov less than mincov
						final_pos_status = True
			else:
				final_pos_status = True
		else:
			final_pos_status = True

	print("Walk TSS position of " + gene + ": " + str(most_upstream_pos))
	walk_tss = most_upstream_pos
	# collect the intron positions as a set
	intron_pos_set = set()
	for intron_start_key, intron_end_val in intron_boundary_marker.items():
		for p in range(intron_start_key, intron_end_val + 1):
			intron_pos_set.add(p)  # collect every intronic position

	# one-time prefix-sum setup, replacing per-window O(W) list rebuilds to speed up the coverage lookups and hence the overall TSS annotation efficiency
	region_lo = most_upstream_pos
	built_hi = region_lo  # tracks how far the arrays are currently built -- grows on demand, no longer a fixed ceiling

	value_prefix = [0]  # builds incrementally now, so value_arr/count_arr as separate intermediate lists are no longer needed; #initial value of cumulative sum is fixed as 0
	count_prefix = [0] # intial value of counts corresponding to cumulative sum value_prefix is fixed to be 0


	def extend_prefix_to(new_hi):  # NEW
		nonlocal built_hi
		if new_hi <= built_hi:
			return
		for pos in range(built_hi, new_hi):
			if pos in coverage_lookup and pos not in intron_pos_set:
				value_prefix.append(value_prefix[-1] + coverage_lookup[pos])
				count_prefix.append(count_prefix[-1] + 1)
			else:
				value_prefix.append(value_prefix[-1])
				count_prefix.append(count_prefix[-1])
		built_hi = new_hi

	initial_hi = start + strength * slide_step + intergenic_region_size  # CHANGED: now just a starting guess, not a hard limit
	extend_prefix_to(initial_hi)

	def window_avg(pos0):
		# avg_coverage for the window starting at pos0, or None if fully intronic/absent
		lo = pos0 - region_lo
		hi = lo + intergenic_region_size
		if hi >= len(count_prefix):  # NEW: extend past the guess if this window needs positions not yet built
			extend_prefix_to(region_lo + hi + 1)
		window_count = count_prefix[hi] - count_prefix[lo]
		if window_count == 0:
			return None
		return float((value_prefix[hi] - value_prefix[lo]) / window_count)

	def accelerated_tss_slice(lo_pos, hi_pos):
		# (pos, cov) list and coverage array for [lo_pos, hi_pos), excluding
		# intronic/absent positions -- exact equivalent of the old list comprehension
		lo = lo_pos - region_lo
		hi = hi_pos - region_lo
		if hi >= len(count_prefix):  # extend if this range reaches past what's built so far
			extend_prefix_to(region_lo + hi + 1)
		positions = np.arange(lo_pos, hi_pos)
		values = np.diff(np.array(value_prefix[lo:hi + 1]))
		valid = np.diff(np.array(count_prefix[lo:hi + 1])).astype(bool)
		filtered_positions = positions[valid]
		filtered_covs = values[valid]
		return list(zip(filtered_positions.tolist(), filtered_covs.tolist())), filtered_covs

	# window sliding strategy to compare against background intergenic noise
	window_pass = 0
	basal_tss = None
	elevated_tss = None
	accelerated_tss = None

	window_starts = []
	window_starts.append(most_upstream_pos)
	while most_upstream_pos <= start:
		avg_coverage = window_avg(most_upstream_pos)
		if avg_coverage is None:
			most_upstream_pos += slide_step
			window_starts.append(most_upstream_pos)
			continue
		B = len(intergenic_window_coverages)
		idx = bisect.bisect_left(intergenic_window_coverages, avg_coverage)
		hits = B - idx
		psig = float(hits / B)
		if psig < pvalue:
			print(f'psig is {psig}')
			if window_pass == 0:  # the window of the walk-based TSS itself is in the signal region
				print("No window sliding happened. TSS based on coverage walk only.")
				elevated_tss = most_upstream_pos
				coverage_slice_list_for_accelerated_tss, coverage_slice_list_for_accelerated_tss_arr = accelerated_tss_slice(elevated_tss, start)
				accelerated_tss = get_accelerated_tss(ks_pval, gene, coverage_slice_list_for_accelerated_tss_arr, lookahead,coverage_slice_list_for_accelerated_tss, '+')
				tsr = 'NA'
				break
			elif window_pass != 0:  # the window of the walk-based TSS is in the intergenic background noise region
				first_most_upstream_pos = most_upstream_pos
				window_pass_signal_strength = 0
				window_pass_check = 0
				while window_pass_check < strength:
					most_upstream_pos = most_upstream_pos + slide_step
					avg_coverage = window_avg(most_upstream_pos)
					if avg_coverage is None:
						most_upstream_pos += slide_step
						window_starts.append(most_upstream_pos)
						continue
					idx = bisect.bisect_left(intergenic_window_coverages, avg_coverage)
					hits = B - idx
					psig = float(hits / B)
					if psig < pvalue:
						window_pass_signal_strength += 1
					window_pass_check += 1
				if window_pass_signal_strength == strength:  # if a window with higher signal then the noise is found, check if the 3 (default; can be changed) consecutive windows next to it are also signal to determine confidently that elevated TSS
					coverage_slice_list_for_accelerated_tss, coverage_slice_list_for_accelerated_tss_arr = accelerated_tss_slice(first_most_upstream_pos, start)
					accelerated_tss = get_accelerated_tss(ks_pval, gene, coverage_slice_list_for_accelerated_tss_arr, lookahead,coverage_slice_list_for_accelerated_tss, '+')
					if slide_step == 1:
						elevated_tss = first_most_upstream_pos
						tsr = 'NA'
					else:
						if len(window_starts) > 1:
							elevated_tss = first_most_upstream_pos
							tsr = f'{window_starts[-2]} to {elevated_tss}'
					break

				else:  # if the three consecutive windows fail to show mean cov values, slide further to look at the next consecutive windows
					window_pass += 1
					# print(f'sliding window pass {window_pass}')
					most_upstream_pos = most_upstream_pos + slide_step
					window_starts.append(most_upstream_pos)
					if most_upstream_pos > start:
						print(f'Caution: No sustained signal region found upstream of {gene}. Returning only the basal TSS from coverage walk.')
						basal_tss = walk_tss
						tsr = 'NA'
						break

		else:  # the walk-based TSS window is in the basal region and sliding window needs to be adopted to capture the elevated and accelerated TSS points (if any)
			basal_tss = walk_tss
			window_pass += 1
			# print(f'sliding window pass {window_pass}')
			most_upstream_pos = most_upstream_pos + slide_step
			window_starts.append(most_upstream_pos)
			if most_upstream_pos > start:
				print(f'Caution: No background elevated signal region found upstream of {gene}. Returning only the basal TSS from coverage walk.')
				basal_tss = walk_tss
				tsr = 'NA'
				break
	# --- generate figures to visualize coverage around the TSS for manual inspection --- #
	tss_list = []
	if basal_tss:
		tss_list.append(basal_tss)
	if elevated_tss:
		tss_list.append(elevated_tss)
	if accelerated_tss:
		tss_list.append(accelerated_tss)
	if other_tool_tss:
		tss_list.append(other_tool_tss)
	starting_tss_for_plot = min(tss_list)

	transcript_list = transcripts_per_gene[gene]
	mrna_dic = {}
	cds_dic = {}
	five_utr_dic = {}
	for each in transcript_list:
		if each in mrna_infos:
			mrna_dic[each] = (mrna_infos[each]['start'], mrna_infos[each]['end'])
		if each in cds_infos:
			cds_dic[each] = cds_infos[each]
		if each in five_utr_infos:
			five_utr_dic[each] = (five_utr_infos[each])
	if starting_tss_for_plot > flank_region_for_plot:
		plot_start_region = starting_tss_for_plot - flank_region_for_plot
	else:
		plot_start_region = 0
	plot_end_region = max(start + flank_region_for_plot, atg_genomic_pos + flank_region_for_plot)
	mrnas_to_plot = []
	cds_to_plot = []
	five_utr_to_plot = []
	for transcript, (mrna_start, mrna_end) in mrna_dic.items():
		if mrna_end >= plot_start_region and mrna_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
			if mrna_start >= plot_start_region and mrna_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (mrna_end - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
				mrnas_to_plot.append((0, (mrna_end - plot_start_region)))
			elif mrna_start >= plot_start_region and mrna_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (plot_end_region - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				mrnas_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, cds_list in cds_dic.items():
		for (cds_start,cds_end) in cds_list:  # two level loop for cds_dic alone since cds_dic structure is a list of tuples per transcript similar to the cds_infos structure from which it is derived
			if cds_end >= plot_start_region and cds_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
				if cds_start >= plot_start_region and cds_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (cds_end - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
					cds_to_plot.append((0, (cds_end - plot_start_region)))
				elif cds_start >= plot_start_region and cds_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (plot_end_region - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					cds_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, five_utr_list in five_utr_dic.items():
		for (five_utr_start, five_utr_end) in five_utr_list: # two level loop for five_utr_dic alone since five_utr_dic structure is a list of tuples per transcript similar to the five_utr_infos structure from which it is derived
			if five_utr_end >= plot_start_region and five_utr_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
				if five_utr_start >= plot_start_region and five_utr_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
					five_utr_to_plot.append(((five_utr_start - plot_start_region), (five_utr_end - plot_start_region)))
				elif five_utr_start <= plot_start_region and five_utr_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
					five_utr_to_plot.append((0, (five_utr_end - plot_start_region)))
				elif five_utr_start >= plot_start_region and five_utr_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
					five_utr_to_plot.append(((five_utr_start - plot_start_region), (plot_end_region - plot_start_region)))
				elif five_utr_start <= plot_start_region and five_utr_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					five_utr_to_plot.append((0, (plot_end_region - plot_start_region)))

	introns_to_plot = []
	for intron_start_key, intron_end_val in intron_boundary_marker.items():
		if intron_end_val >= plot_start_region and intron_start_key <= plot_end_region:  # primary bounds check
			if intron_start_key >= plot_start_region and intron_end_val <= plot_end_region:  # case1: fully inside
				introns_to_plot.append(((intron_start_key - plot_start_region), (intron_end_val - plot_start_region)))
			elif intron_start_key <= plot_start_region and intron_end_val <= plot_end_region:  # case2: start out of bounds
				introns_to_plot.append((0, (intron_end_val - plot_start_region)))
			elif intron_start_key >= plot_start_region and intron_end_val >= plot_end_region:  # case3: end out of bounds
				introns_to_plot.append(((intron_start_key - plot_start_region), (plot_end_region - plot_start_region)))
			elif intron_start_key <= plot_start_region and intron_end_val >= plot_end_region:  # case4: spans whole plot
				introns_to_plot.append((0, (plot_end_region - plot_start_region)))
	boundary_tss_for_plot = starting_tss_for_plot - plot_start_region
	basal_tss_pos=None
	elevated_tss_pos=None
	accelerated_tss_pos=None
	other_tool_tss_pos = None
	values = cov_per_contig[plot_start_region:plot_end_region]
	svalues = scov_per_contig[plot_start_region:plot_end_region]
	atg_pos = atg_genomic_pos - plot_start_region
	if basal_tss:
		basal_tss_pos = basal_tss - plot_start_region
	if elevated_tss:
		elevated_tss_pos = elevated_tss - plot_start_region
	if accelerated_tss:
		accelerated_tss_pos = accelerated_tss - plot_start_region
	if other_tool_tss:
		other_tool_tss_pos = other_tool_tss - plot_start_region
	cov_walk_start = start - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "+"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end + 1]
	try:
		generate_plot(other_tool, other_tool_tss_pos, values, svalues, fig_file, atg_pos, cov_walk_start, boundary_tss_for_plot, basal_tss_pos, elevated_tss_pos, accelerated_tss_pos, genomic_start, genomic_end, gene,orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot, introns_to_plot)
	except:
		print("ERROR: plot failed" + gene)

	basal_tss_yr_compliant = False
	elevated_tss_yr_compliant = False
	accelerated_tss_yr_compliant = False
	other_tool_tss_yr_compliant = False

	if basal_tss is not None:  # seq_per_contig slicing is 0-based while tss position is genomic index based. So seq_per_contig[basal_tss-1] is the base at TSS and seq_per_contig[basal_tss-2] is the base just preceding the TSS
		if ((seq_per_contig[basal_tss - 2].upper() == 'T' or seq_per_contig[basal_tss - 2].upper() == 'C') and (seq_per_contig[basal_tss - 1].upper() == 'A' or seq_per_contig[basal_tss - 1].upper() == 'G')):
			basal_tss_yr_compliant = True
	else:
		basal_tss_yr_compliant = 'NA'
	if elevated_tss is not None:
		if ((seq_per_contig[elevated_tss - 2].upper() == 'T' or seq_per_contig[elevated_tss - 2].upper() == 'C') and (seq_per_contig[elevated_tss - 1].upper() == 'A' or seq_per_contig[elevated_tss - 1].upper() == 'G')):
			elevated_tss_yr_compliant = True
	else:
		elevated_tss_yr_compliant = 'NA'
	if accelerated_tss is not None:
		if ((seq_per_contig[accelerated_tss - 2].upper() == 'T' or seq_per_contig[accelerated_tss - 2].upper() == 'C') and (seq_per_contig[accelerated_tss - 1].upper() == 'A' or seq_per_contig[accelerated_tss - 1].upper() == 'G')):
			accelerated_tss_yr_compliant = True
	else:
		accelerated_tss_yr_compliant = 'NA'

	if other_tool_tss is not None:
		if ((seq_per_contig[other_tool_tss - 2].upper() == 'T' or seq_per_contig[other_tool_tss - 2].upper() == 'C') and (seq_per_contig[other_tool_tss - 1].upper() == 'A' or seq_per_contig[other_tool_tss - 1].upper() == 'G')):
			other_tool_tss_yr_compliant = True
	else:
		other_tool_tss_yr_compliant = 'NA'

	return {'TSS': walk_tss, 'start': start,'end': end}, walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant, accelerated_tss_yr_compliant, other_tool_tss_yr_compliant, hard_cutoff_reached


def run_rev_analysis(other_tool, other_tss_dic, ks_pval, strength, lookahead, output_folder, pvalue,
					 intergenic_region_size, slide_step, intergenic_window_coverages, coverage_lookup, gene,
					 cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size,
					 hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_genomic_pos, contig,
					 genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos,
					 gene_atg_dic, cds_infos, percentdiff_threshold):
	"""! @brief run analysis on reverse strand """
	if gene in other_tss_dic:
		other_tool_tss = other_tss_dic[gene]
	else:
		other_tool_tss = None

	hard_cutoff_reached = False
	intron_boundary_marker = {}
	most_downstream_pos = end
	final_pos_status = False
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
					pass
			if most_downstream_pos == hard_cutoff:  # stop if end of contig/pseudochromosome is reached
				hard_cutoff_reached = 'True;'+str(hard_cutoff)
				break
		print(f'most downstream pos is {most_downstream_pos}')
		# --- try to cross intron --- #
		current_position = most_downstream_pos + 1  # most_downstream_pos has coverage above cutoff (position, not index!)
		if current_position < hard_cutoff:
			while cov_per_contig[current_position] < mincov:  # check if downstream position has low coverage
				current_position += 1  # move one step downstream
				if current_position == hard_cutoff:
					print(f'current position is {current_position}')
					hard_cutoff_reached = 'True;'+str(hard_cutoff)
					break
			avg_gap_coverage = sum( scov_per_contig[ most_downstream_pos-1 : current_position-1 ] )/(current_position - most_downstream_pos)
			if current_position < (len(seq_per_contig) - min_exon_size) and avg_gap_coverage > mincov:
				donor_splice_site = seq_per_contig[current_position - 2:current_position].upper()
				acceptor_splice_site = seq_per_contig[most_downstream_pos - 1:most_downstream_pos + 1].upper()
				if splicesites == "off":
					if current_position - most_downstream_pos < tolerated_gap:
						intron_boundary_marker[most_downstream_pos] = current_position #the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
						most_downstream_pos = current_position + 1
					else:
						final_pos_status = True
				elif donor_splice_site == 'AC' and acceptor_splice_site == "CT":
					print(f"Low coverage regime: Canonical donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Low coverage regime: Canonical acceptor splice site ending at {most_downstream_pos}: " + acceptor_splice_site + '\n')
					intron_boundary_marker[most_downstream_pos] = current_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
					most_downstream_pos = current_position + 1
				elif donor_splice_site != 'AC' or acceptor_splice_site != "CT":
					print(f"Low coverage regime: Warning! Either one or both the donor and acceptor splice sites are non-canonical." + '\n')
					print(f"Low coverage regime: Donor splice site starting at {current_position}: " + donor_splice_site + '\n')
					print(f"Low coverage regime: Acceptor splice site ending at {most_downstream_pos}: " + acceptor_splice_site + '\n')
					if current_position - most_downstream_pos < tolerated_gap:
						intron_boundary_marker[most_downstream_pos] = current_position  # the key is acceptor site and value is donor site since donor is in high coords in rev strand and acceptor is in low coords in rev strand
						most_downstream_pos += 1
					else:
						final_pos_status = True
			else:
				final_pos_status = True
		else:
			final_pos_status = True

	print("Walk TSS position of " + gene + ": " + str(most_downstream_pos))
	walk_tss = most_downstream_pos

	# collect the intron positions as a set
	intron_pos_set = set()
	for intron_start_key, intron_end_val in intron_boundary_marker.items():
		for p in range(intron_start_key, intron_end_val + 1):
			intron_pos_set.add(p)  # collect every intronic position

	# one-time prefix-sum setup, replacing per-window O(W) list rebuilds to speed up the coverage lookups and hence the overall TSS annotation efficiency
	region_hi = most_downstream_pos
	built_lo = region_hi # tracks how far the arrays are currently built -- grows on demand, no longer a fixed ceiling

	value_prefix = [0]#initial value of cumulative sum is fixed as 0
	count_prefix = [0]# intial value of counts corresponding to cumualtive sum vlaue_prefix is fixed to be 0

	def extend_prefix_to(new_lo):
		nonlocal built_lo
		if new_lo >= built_lo:
			return
		for pos in range(built_lo, new_lo, -1):
			if pos in coverage_lookup and pos not in intron_pos_set:
				value_prefix.append(value_prefix[-1] + coverage_lookup[pos])
				count_prefix.append(count_prefix[-1] + 1)
			else:
				value_prefix.append(value_prefix[-1])
				count_prefix.append(count_prefix[-1])
		built_lo = new_lo

	initial_lo = end - strength*slide_step - intergenic_region_size
	extend_prefix_to(initial_lo)

	def window_avg(pos0):
		d = region_hi - pos0
		hi_idx = d + intergenic_region_size
		if hi_idx >= len(count_prefix):
			extend_prefix_to(pos0 - intergenic_region_size)
		window_count = count_prefix[hi_idx] - count_prefix[d]
		if window_count == 0:
			return None
		return float((value_prefix[hi_idx] - value_prefix[d]) / window_count)

	def accelerated_tss_slice(lo_pos, hi_pos):
		idx_start = region_hi - hi_pos
		idx_end = region_hi - lo_pos - 1
		if idx_end + 1 >= len(count_prefix):
			extend_prefix_to(lo_pos)
		idx_range = np.arange(idx_start, idx_end + 1)
		positions = region_hi - idx_range
		values = np.diff(np.array(value_prefix[idx_start:idx_end + 2]))
		valid = np.diff(np.array(count_prefix[idx_start:idx_end + 2])).astype(bool)
		filtered_positions = positions[valid]
		filtered_covs = values[valid]
		return list(zip(filtered_positions.tolist(), filtered_covs.tolist())), filtered_covs


	# window sliding strategy to compare against background intergenic noise
	window_pass = 0
	basal_tss = None
	elevated_tss = None
	accelerated_tss = None
	secondary_accelerated_tss = None

	window_starts = []
	window_starts.append(most_downstream_pos)
	while most_downstream_pos >= end:
		avg_coverage = window_avg(most_downstream_pos)
		if avg_coverage is None:
			most_downstream_pos -= slide_step
			window_starts.append(most_downstream_pos)
			continue
		B = len(intergenic_window_coverages)
		idx = bisect.bisect_left(intergenic_window_coverages, avg_coverage)
		hits = B - idx
		psig = float(hits / B)
		if psig < pvalue:
			print(f'psig is {psig}')
			if window_pass == 0:  # the window of the walk-based TSS itself is in the signal region
				print("No window sliding happened. TSS based on coverage walk only.")
				elevated_tss = most_downstream_pos
				coverage_slice_list_for_accelerated_tss, coverage_slice_list_for_accelerated_tss_arr = accelerated_tss_slice(end, elevated_tss)
				accelerated_tss = get_accelerated_tss(ks_pval, gene, coverage_slice_list_for_accelerated_tss_arr, lookahead,coverage_slice_list_for_accelerated_tss, '-')
				tsr = 'NA'
				break
			elif window_pass != 0:  # the window of the walk-based TSS is in the intergenic background noise region
				first_most_downstream_pos = most_downstream_pos
				window_pass_signal_strength = 0
				window_pass_check = 0
				while window_pass_check < strength:
					most_downstream_pos = most_downstream_pos - slide_step
					avg_coverage = window_avg(most_downstream_pos)
					if avg_coverage is None:
						most_downstream_pos -= slide_step
						window_starts.append(most_downstream_pos)
						continue
					idx = bisect.bisect_left(intergenic_window_coverages, avg_coverage)
					hits = B - idx
					psig = float(hits / B)
					if psig < pvalue:
						window_pass_signal_strength += 1
					window_pass_check += 1
				if window_pass_signal_strength == strength:  # if a window with higher signal then the noise is found, check if the 3 (default; can be changed) consecutive windows next to it are also signal to determine confidently that elevated TSS
					coverage_slice_list_for_accelerated_tss, coverage_slice_list_for_accelerated_tss_arr = accelerated_tss_slice(end, first_most_downstream_pos)
					accelerated_tss = get_accelerated_tss(ks_pval, gene, coverage_slice_list_for_accelerated_tss_arr, lookahead,coverage_slice_list_for_accelerated_tss, '-')
					if slide_step == 1:
						elevated_tss = first_most_downstream_pos
						tsr = 'NA'
					else:
						if len(window_starts) > 1:
							elevated_tss = first_most_downstream_pos
							tsr = f'{window_starts[-2]} to {elevated_tss}'
					break
				else:
					window_pass += 1
					# print(f'sliding window pass {window_pass}')
					most_downstream_pos = most_downstream_pos - slide_step
					window_starts.append(most_downstream_pos)
					if most_downstream_pos < end:
						print(f'Caution: No sustained signal region found upstream of {gene}. Returning only the basal TSS from coverage walk.')
						basal_tss = walk_tss
						tsr = 'NA'
						break
		else:
			basal_tss = walk_tss
			window_pass += 1
			# print(f'sliding window pass {window_pass}')
			most_downstream_pos = most_downstream_pos - slide_step
			window_starts.append(most_downstream_pos)
			if most_downstream_pos < end:
				print(f'Caution: No background elevated signal region found upstream of {gene}. Returning only the basal TSS from coverage walk.')
				basal_tss = walk_tss
				tsr = 'NA'
				break

	# --- generate figures to visualize coverage around the TSS for manual inspection --- #
	tss_list = []
	if basal_tss:
		tss_list.append(basal_tss)
	if elevated_tss:
		tss_list.append(elevated_tss)
	if accelerated_tss:
		tss_list.append(accelerated_tss)
	if other_tool_tss:
		tss_list.append(other_tool_tss)
	end_tss_for_plot = max(tss_list)
	transcript_list = transcripts_per_gene[gene]
	mrna_dic = {}
	cds_dic = {}
	five_utr_dic = {}
	for each in transcript_list:
		if each in mrna_infos:
			mrna_dic[each] = (mrna_infos[each]['start'], mrna_infos[each]['end'])
		if each in cds_infos:
			cds_dic[each] = cds_infos[each]
		if each in five_utr_infos:
			five_utr_dic[each] = (five_utr_infos[each])
	plot_start_region = min(end - flank_region_for_plot, atg_genomic_pos - flank_region_for_plot)
	if end_tss_for_plot < (len(seq_per_contig) - flank_region_for_plot):
		plot_end_region = end_tss_for_plot + flank_region_for_plot
	else:
		plot_end_region = len(seq_per_contig)

	mrnas_to_plot = []
	cds_to_plot = []
	five_utr_to_plot = []
	for transcript, (mrna_start, mrna_end) in mrna_dic.items():
		if mrna_end >= plot_start_region and mrna_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
			if mrna_start >= plot_start_region and mrna_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (mrna_end - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
				mrnas_to_plot.append((0, (mrna_end - plot_start_region)))
			elif mrna_start >= plot_start_region and mrna_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
				mrnas_to_plot.append(((mrna_start - plot_start_region), (plot_end_region - plot_start_region)))
			elif mrna_start <= plot_start_region and mrna_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				mrnas_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, cds_list in cds_dic.items():
		for (cds_start,
			 cds_end) in cds_list:  # two level loop for cds_dic alone since cds_dic structure is a list of tuples per transcript similar to the cds_infos structure from which it is derived
			if cds_end >= plot_start_region and cds_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
				if cds_start >= plot_start_region and cds_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (cds_end - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
					cds_to_plot.append((0, (cds_end - plot_start_region)))
				elif cds_start >= plot_start_region and cds_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
					cds_to_plot.append(((cds_start - plot_start_region), (plot_end_region - plot_start_region)))
				elif cds_start <= plot_start_region and cds_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					cds_to_plot.append((0, (plot_end_region - plot_start_region)))

	for transcript, five_utr_list in five_utr_dic.items():
		for (five_utr_start,five_utr_end) in five_utr_list:  # two level loop for five_utr_dic alone since five_utr_dic structure is a list of tuples per transcript similar to the five_utr_infos structure from which it is derived
			if five_utr_end >= plot_start_region and five_utr_start <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
				if five_utr_start >= plot_start_region and five_utr_end <= plot_end_region:  # case1 where feature is entirely within the plot bounds
					five_utr_to_plot.append(((five_utr_start - plot_start_region), (five_utr_end - plot_start_region)))
				elif five_utr_start <= plot_start_region and five_utr_end <= plot_end_region:  # case2 where feature start is out of the plot bounds
					five_utr_to_plot.append((0, (five_utr_end - plot_start_region)))
				elif five_utr_start >= plot_start_region and five_utr_end >= plot_end_region:  # case3 where feature end is out of the plot bounds
					five_utr_to_plot.append(((five_utr_start - plot_start_region), (plot_end_region - plot_start_region)))
				elif five_utr_start <= plot_start_region and five_utr_end >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
					five_utr_to_plot.append((0, (plot_end_region - plot_start_region)))

	introns_to_plot = []
	for intron_start_key, intron_end_val in intron_boundary_marker.items():
		if intron_end_val >= plot_start_region and intron_start_key <= plot_end_region:  # primary check to see if feature is within the bounds of the plot
			if intron_start_key >= plot_start_region and intron_end_val <= plot_end_region:  # case1 where feature is entirely within the plot bounds
				introns_to_plot.append(((intron_start_key - plot_start_region), (intron_end_val - plot_start_region)))
			elif intron_start_key <= plot_start_region and intron_end_val <= plot_end_region:  # case2 where feature start is out of the plot bounds
				introns_to_plot.append((0, (intron_end_val - plot_start_region)))
			elif intron_start_key >= plot_start_region and intron_end_val >= plot_end_region:  # case3 where feature end is out of the plot bounds
				introns_to_plot.append(((intron_start_key - plot_start_region), (plot_end_region - plot_start_region)))
			elif intron_start_key <= plot_start_region and intron_end_val >= plot_end_region:  # case4 where both feature start and end are out of the plot bounds making the feature span the entire plot boundary
				introns_to_plot.append((0, (plot_end_region - plot_start_region)))
	boundary_tss_for_plot = end_tss_for_plot - plot_start_region
	basal_tss_pos=None
	elevated_tss_pos=None
	accelerated_tss_pos=None
	other_tool_tss_pos = None
	values = cov_per_contig[plot_start_region:plot_end_region]
	svalues = scov_per_contig[plot_start_region:plot_end_region]
	atg_pos = atg_genomic_pos - plot_start_region
	if basal_tss:
		basal_tss_pos = basal_tss - plot_start_region
	if elevated_tss:
		elevated_tss_pos = elevated_tss - plot_start_region
	if accelerated_tss:
		accelerated_tss_pos = accelerated_tss - plot_start_region
	if other_tool_tss:
		other_tool_tss_pos = other_tool_tss - plot_start_region
	cov_walk_start = end - plot_start_region
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "-"
	dna_sequence_for_plot = genome_seq[contig][genomic_start:genomic_end + 1]
	try:
		generate_plot(other_tool, other_tool_tss_pos,values, svalues, fig_file, atg_pos, cov_walk_start, boundary_tss_for_plot, basal_tss_pos, elevated_tss_pos, accelerated_tss_pos, genomic_start, genomic_end, gene,orientation, dna_sequence_for_plot, mrnas_to_plot, cds_to_plot, five_utr_to_plot, introns_to_plot)
	except:
		print("ERROR: plot failed" + gene)

	basal_tss_yr_compliant = False
	elevated_tss_yr_compliant = False
	accelerated_tss_yr_compliant = False
	other_tss_yr_compliant = False

	if basal_tss is not None:  # seq_per_contig slicing is 0-based while tss position is genomic index based. Since this is the segment for reverse strand gene, the bases should follow reverse complementation. So seq_per_contig[basal_tss-1] is the base at TSS and seq_per_contig[basal_tss] is the base just succeeding the TSS
		if ((seq_per_contig[basal_tss].upper() == 'A' or seq_per_contig[basal_tss].upper() == 'G') and (
				seq_per_contig[basal_tss - 1].upper() == 'T' or seq_per_contig[basal_tss - 1].upper() == 'C')):
			basal_tss_yr_compliant = True
	else:
		basal_tss_yr_compliant = 'NA'
	if elevated_tss is not None:
		if ((seq_per_contig[elevated_tss].upper() == 'A' or seq_per_contig[elevated_tss].upper() == 'G') and (
				seq_per_contig[elevated_tss - 1].upper() == 'T' or seq_per_contig[elevated_tss - 1].upper() == 'C')):
			elevated_tss_yr_compliant = True
	else:
		elevated_tss_yr_compliant = 'NA'
	if accelerated_tss is not None:
		if ((seq_per_contig[accelerated_tss].upper() == 'A' or seq_per_contig[accelerated_tss].upper() == 'G') and (
				seq_per_contig[accelerated_tss - 1].upper() == 'T' or seq_per_contig[
			accelerated_tss - 1].upper() == 'C')):
			accelerated_tss_yr_compliant = True
	else:
		accelerated_tss_yr_compliant = 'NA'

	if other_tool_tss is not None:
		if ((seq_per_contig[other_tool_tss].upper() == 'A' or seq_per_contig[other_tool_tss].upper() == 'G') and (seq_per_contig[other_tool_tss - 1].upper() == 'T' or seq_per_contig[other_tool_tss - 1].upper() == 'C')):
			other_tss_yr_compliant = True
	else:
		other_tss_yr_compliant = 'NA'

	return {'TSS': walk_tss, 'start': start,'end': end}, walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant, accelerated_tss_yr_compliant, other_tss_yr_compliant, hard_cutoff_reached


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
	

def extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, tss, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size ):
	"""! @brief extract promoter region """

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
	if moods == 'moods-dna.py':
		cmd = moods + ' -m ' + pfm_folder + '/*.pfm ' + '-s ' + sequence + ' -p ' + str(pvalue) + ' -o ' + output_file
	else:
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
def promoter_motif_analysis (numcols, upstream_slice, downstream_slice, bg_scores, tss, tss_type, gene, orientation, promoter_seq, downstream_to_tss, moods, pvalue, pfm_config_dic, tmp_folder, output_folder):
	tss_neighbourhood = os.path.join(output_folder,f'{gene}_{tss_type}_tss_neighbourhood.png')
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
		tmp.write(f">{gene}_{tss_type}\n{full_seq}\n")
		tmp.flush()
		tmp.close()
		seq_len = len(full_seq)
		downstream_len = len(downstream_to_tss)

		sorted_hits_dic = {}
		cumulative_promoter_motif_score = 0
		for elements in pfm_config_dic:
			pfm_folder = pfm_config_dic[elements][0]
			output_file = os.path.join(tmp_folder, f"{gene}_{tss_type}_tss_moods_{elements}.txt")
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
		fig, axes = plt.subplots(nrows, ncols,figsize=(6 * ncols, 3 * nrows),squeeze=False, layout = 'constrained')
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

		# Hide any unused subplot panels
		for idx in range(n_motifs, nrows * ncols):
			row, col = divmod(idx, ncols)
			axes[row][col].set_visible(False)

		fig.legend(handles=legend_elements, loc = 'outside upper center', ncol=len(legend_elements))
		plt.savefig(tss_neighbourhood, dpi=600, bbox_inches = 'tight')
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

	if "--protein_encoding" in arguments:#yes or no to consider only protein encoding genes for TSS analysis
		protein_encoding = arguments[arguments.index("--protein_encoding")+1]
	else:
		protein_encoding = "no"

	if '--compare_tss' in arguments:
		other_tool = arguments[arguments.index('--compare_tss') + 1]  # specify epd tssfinder and so on
	else:
		other_tool = None

	other_tss_dic = {}
	if '--compare_tss_input' in arguments:
		other_tss_file = arguments[arguments.index('--compare_tss_input') + 1]  # full path to two column tab separated config file where column one is GOI and column two is the TSS for it from the source being compared against
		with open(other_tss_file, "r") as f:
			for line in f:
				parts = line.strip('\n').split()
				other_tss_dic[parts[0]] = int(parts[1])  # first column goi and second column tss position
	else:
		other_tss_file = None
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

	if '--parallel' in arguments: #full path to GNU parallel
		parallel = arguments[arguments.index('--parallel') + 1]
	else:
		parallel = 'parallel'

	if '--sample_support' in arguments: #yes or no for sample support-based tss scoring; default is no
		tss_scoring = arguments[arguments.index('--sample_support') + 1]
	else:
		tss_scoring = 'no'

	if '--m' in arguments:
		m = arguments[arguments.index('--m') + 1]
	else:
		m = "5000000"

	if '--threads' in arguments:
		t = int(arguments[arguments.index('--threads') + 1])
	else:
		t = 4

	if '--intron_percentile_cutoff' in arguments:
		percentile_cut = int(arguments[arguments.index('--intron_percentile_cutoff') + 1])
	else:
		percentile_cut = 99

	parallelize_cov_generation = False
	if '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]
		if os.path.isfile(bam_file):
			parallelize_cov_generation = False
		elif os.path.isdir(bam_file):
			parallelize_cov_generation = True

		if '--bam_is_sorted' in arguments:
			bam_sorted_status = True
		else:
			bam_sorted_status = False
			sorted_bam_folder = os.path.join(output_folder,'Sorted_BAM_files')
			if not os.path.exists(sorted_bam_folder):
				os.makedirs(sorted_bam_folder)
		
		if not bam_sorted_status:	#sorting the BAM file if it was not sorted already
			print ("sorting BAM file(s) ...")
			if parallelize_cov_generation:# to handle folder of BAM files sorting
				files = [str(p) for p in Path(bam_file).iterdir()]
				for f in files:
					file_name = Path(f).stem #stem removes the last extension from the file name
					sorted_bam_file = os.path.join(sorted_bam_folder,f'{file_name}_sorted.bam')
					cmd = samtools + " sort -m " + str(m) + " --threads " + str(t) + " " + f + " > " + sorted_bam_file
					p = subprocess.Popen(args=cmd, shell=True)
					p.communicate()
			else:#to handle single BAM file sorting
				file_name = Path(bam_file).stem  # stem removes the last extension from the file name
				sorted_bam_file = os.path.join(sorted_bam_folder, f'{file_name}_sorted.bam')
				cmd = samtools + " sort -m " + str(m) + " --threads " + str(t) + " " + bam_file + " > " + sorted_bam_file
				p = subprocess.Popen( args= cmd, shell=True )
				p.communicate()
		t_cov_start = time.perf_counter()
		bam_files = os.path.join(output_folder, 'BAM_files.txt')
		if os.path.isfile(bam_file):
			if bam_sorted_status:
				bam_list = [bam_file]
			else:
				bam_list = [sorted_bam_file]
		elif os.path.isdir(bam_file):
			if bam_sorted_status:
				bam_list = [
					os.path.join(bam_file, f)
					for f in os.listdir(bam_file)
					if f.endswith('.bam')
				]
			else:
				bam_list = [
					os.path.join(sorted_bam_folder, f)
					for f in os.listdir(sorted_bam_folder)
					if f.endswith('.bam')
				]
		with open(bam_files, 'w') as out:
			for f in bam_list:
				if f.endswith('bam'):
					out.write(f + '\n')
		aligned_read_coverage = os.path.join(output_folder,'Aligned_reads_cov')
		if not os.path.exists(aligned_read_coverage):
			os.makedirs(aligned_read_coverage)
			cov_file = construct_coverage_files (bam_list, bam_files, aligned_read_coverage, bedtools, parallelize_cov_generation, parallel,t, tss_scoring, 'aligned')
		spanning_read_coverage = os.path.join(output_folder, 'Spanning_reads_cov')
		if not os.path.exists(spanning_read_coverage):
			os.makedirs(spanning_read_coverage)
			scov_file = construct_coverage_files (bam_list, bam_files, spanning_read_coverage, bedtools, parallelize_cov_generation, parallel,t, tss_scoring, 'spanning')
		t_cov_end = time.perf_counter()
		time_cov_file_creation = t_cov_end - t_cov_start
		print(f'time taken for cov and scov file generation is {time_cov_file_creation} seconds')

	elif '--sra_folder' not in arguments:
		coverage_file = arguments[ arguments.index('--cov')+1 ]
		scoverage_file = arguments[arguments.index('--scov') + 1]
		if os.path.isfile(coverage_file) and os.path.isfile(scoverage_file):
			cov_file = coverage_file
			scov_file = scoverage_file
		elif os.path.isdir(coverage_file) and os.path.isdir(scoverage_file):
			input_aligned_read_coverage = coverage_file
			input_spanning_read_coverage = scoverage_file
			aligned_read_coverage = os.path.join(output_folder, 'Aligned_reads_cov')
			if not os.path.exists(aligned_read_coverage):
				os.makedirs(aligned_read_coverage)
			aligned_cov_files_list = [
				os.path.join(coverage_file, f)
				for f in os.listdir(coverage_file)
				if f.endswith('.cov')
			]
			spanning_read_coverage = os.path.join(output_folder, 'Spanning_reads_cov')
			if not os.path.exists(spanning_read_coverage):
				os.makedirs(spanning_read_coverage)
			spanning_cov_files_list = [
				os.path.join(scoverage_file, f)
				for f in os.listdir(scoverage_file)
				if f.endswith('.cov')
			]
			cov_file = sum_coverage_files(aligned_cov_files_list, input_aligned_read_coverage, aligned_read_coverage, parallel, t, tss_scoring, 'aligned')
			scov_file = sum_coverage_files(spanning_cov_files_list, input_spanning_read_coverage, spanning_read_coverage, parallel,t, tss_scoring, 'spanning')

	#coverage cutoff
	if '--mincov' in arguments:
		mincov = int( arguments[ arguments.index('--mincov')+1 ] )
	else:
		mincov = 1

	if '--coverage_difference' in arguments:#threshold to be used as percentage difference threshold between aligned and spanning read coverage; default is 10%
		percentdiff_threshold = float(arguments[arguments.index("--coverage_difference")+1])
	else:
		percentdiff_threshold = 10
	
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
	#flag for determining background strength for promoter motif random background seq retrieval
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
		pattern_names = arguments[arguments.index('--fastq_pattern')+1]#specify fastq read file name pattern for the paired end files separated by commas without spaces like - _pass_1,_pass_2_
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

	#background percentage threshold for basal vs elevated transcription regions identification
	if '--background_percentage' in arguments:
		background_percentage = float(arguments[arguments.index('--background_percentage')+1])
	else:
		background_percentage = 0.05 #5% is default

	#p-value for Kolmogorov Smirnov distribution fit test
	if '--ks_pval' in arguments:
		ks_pval = float(arguments[arguments.index('--ks_pval')+1])
	else:
		ks_pval = 0.01

	#p-value threshold for MOODS analysis
	if '--moods_pval' in arguments:
		pvalue = float(arguments[arguments.index('--moods_pval')+1])
	else:
		pvalue = 0.01

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

	#flag to fix the size of intergenic region, intron, exon tiles chosen as RNA-seq intergenic background units
	if '--background_unit' in arguments:
		intergenic_region_size = int(arguments[arguments.index('--background_unit')+1])
	else:
		intergenic_region_size = 10

	#flag to trim bases from intron ends to avoid splice junctions from showing up in intron coverage analysis
	if '--intron_trim' in arguments:
		trim = int(arguments[arguments.index('--intron_trim')+1])
	else:
		trim = 5

	#flag to fix the size of window sliding steps for comparing against intergenic background
	if '--slide' in arguments:
		slide_step = int(arguments[arguments.index('--slide')+1])
	else:
		slide_step = 1

	#flag to determine the number of consecutive windows that must have mean cov values above the background
	if '--signal_strength' in arguments:
		strength = int(arguments[arguments.index('--signal_strength')+1])
	else:
		strength = 3

	#flag to determine the number of bases or positions to be looked ahead for determining the minima point accompanying the steepest positive transition in coverage
	if '--lookahead' in arguments:
		lookahead = int(arguments[arguments.index('--lookahead')+1])
	else:
		lookahead = 20

	#flag to take a buffer region before slicing of the intergenic region to avoid read through signal interruptions from actual exonic regions
	if '--buffer' in arguments:
		intergenic_buffer = int(arguments[arguments.index('--buffer')+1])
	else:
		intergenic_buffer = 500

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
			cmd = star+' --runMode genomeGenerate --genomeDir '+star_indexing_folder+' --genomeFastaFiles '+fasta_file+' --sjdbGTFfile '+gff_file+' --sjdbGTFtagExonParentTranscript Parent --runThreadN '+str(t)+' --genomeSAindexNbases '+str(index_bases)
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
			print(f"[INFO] Using read patterns: R1='{pattern_names_list[0]}', R2='{pattern_names_list[1]}'")
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
						cmd = 'ulimit -n 4096 && '+star+' --runMode alignReads --genomeDir '+star_indexing_folder+' --outSAMtype BAM SortedByCoordinate --readFilesIn '+r1+' '+r2+' --runThreadN '+str(t)+' --outFileNamePrefix '+prefix+' --readFilesCommand zcat --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --alignIntronMax '+ str(intron_cutoff)
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
			with open(bam_files, 'w') as out:
				for f in bam_list:
					if f.endswith('bam'):
						out.write(f + '\n')

			aligned_read_coverage = os.path.join(output_folder, 'Aligned_reads_cov')
			if not os.path.exists(aligned_read_coverage):
				os.makedirs(aligned_read_coverage)
				cov_file = construct_coverage_files (bam_list, bam_files, aligned_read_coverage, bedtools, parallelize_cov_generation,parallel, t, tss_scoring, 'aligned')
			spanning_read_coverage = os.path.join(output_folder, 'Spanning_reads_cov')
			if not os.path.exists(spanning_read_coverage):
				os.makedirs(spanning_read_coverage)
				scov_file = construct_coverage_files (bam_list, bam_files, spanning_read_coverage, bedtools, parallelize_cov_generation,parallel, t, tss_scoring, 'spanning')

	# --- load data --- #

	coverage = load_coverage( cov_file)
	scoverage = load_coverage( scov_file)

	#obtain a dictionary with values being list of tuples of position and coverage per position per contig
	coverage_dic = {}
	with open(cov_file, 'r') as f:
		line = f.readline()
		while line:
			parts = line.strip().split()
			try:
				coverage_dic[parts[0]].append((int(parts[1]), int(parts[2])))
			except KeyError:
				coverage_dic.update({parts[0]: [(int(parts[1]), int(parts[2]))]})
			line = f.readline()

	#convert the position-coverage tuples to dictionaries per contig for easy coverage lookup later
	lookup_dic = {}
	for contig in coverage_dic:
		lookup_dic[contig] = dict(coverage_dic[contig])

	gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos, protein_coding_genes_per_chromosome = load_gene_infos( protein_encoding, gff_file, child_attribute, child_parent_linker, parent_attribute )
	genome_seq, seq_counter = load_sequences( fasta_file )
	t_intergenic_start = time.perf_counter()
	print(f'Retrieving intergenic background seqs')
	intergenic_window_coverages, z_ig, ig_percent_zero, ig_percent_nonzero = get_intergenic_background_seqs(output_folder, coverage_dic, intergenic_buffer, intergenic_region_size,genes_per_chromosome, gene_infos)
	intergenic_window_coverages = sorted(intergenic_window_coverages)#sorting the intergenic window coverages list for the downstream bisect operations
	print(f'Completed retrieving intergenic background seqs')
	t_intergenic_end = time.perf_counter()
	t_intergenic_time = t_intergenic_end - t_intergenic_start
	print(f'time taken to complete intergenic background seqs: {t_intergenic_time}')
	if promoter_analysis == 'yes':
		print("Retrieving background seqs.")
		background_seq_len = upstream_slice + downstream_slice
		background_seqs = get_random_background_seqs(genome_seq, background_seq_len, background_strength)
		print("Completed retrieving background seqs.")
		print("Starting MOODS scoring of background seqs.")
		bg_scores = compute_background_moods_scores(background_seqs, pfm_config_dic, tmp_folder, output_folder, moods,pvalue)
		print("Completed MOODS scoring of background seqs.")

	# run analysis per gene of interest
	gene_exp_status_dic = {}
	inflection_dic = {}
	candidate_tss_dic = {}
	tss_compare_dic = {}
	results = {}
	motifs = []
	motif_scores = {}
	pvalue_dic = {}
	full_seq_pos_strand_gene = {}
	full_seq_neg_strand_gene = {}
	promoter_status_dic = {}
	promoter_dic = {}
	# confidence_score_dic = {}
	tss_confidence_dic = {}
	isoforms_dic = {}
	distance_dic = {}
	avg_coverage_dic = {}
	five_utr_dic = {}
	percentile_dic = {}
	canonical_hits_dic = {}
	tflank = 0
	texp = 0
	toverlap = 0
	tanalysis = 0
	textract = 0
	t_tss_start_analysis = time.perf_counter()
	for gene in goi:
		if protein_encoding == "yes":
			if gene not in transcripts_per_gene:
				print(f'You chose only protein encoding gene analysis mode. But your GOI {gene} not found in protein-encoding gene group. Skipping TSS analysis for {gene}\n')
				continue
		try:
			cov_per_contig = coverage[gene_infos[gene]['chromosome']]  # get coverage of the sequence that harbours the gene of interest
			scov_per_contig = scoverage[gene_infos[gene]['chromosome']]  # get spanning read coverage of the sequence that harbours the gene of interest
			seq_per_contig = genome_seq[gene_infos[gene]['chromosome']]  # get the sequence of the contig/pseudochromosome that harbours the gene of interest
			# code block to check if the goi has annotated 5'UTR and if yes take the most upstream/ downstream 5'UTR start/ end as start or end according to + or - strand orientation
			# Initialize with gene coordinates as default with 5'UTR checks downstream
			gene_start, gene_end, orientation, contig = gene_infos[gene]['start'], gene_infos[gene]['end'], gene_infos[gene]['orientation'], gene_infos[gene]['chromosome']  # get information about gene of interest if it does not have 5'UTRs annotated
			coverage_lookup = lookup_dic[contig]
			# store ATG position as fixed reference before start may be modified by 5'UTR code block
			if gene in gene_atg_dic:
				atg_pos = gene_atg_dic[gene]
			else:
				atg_pos = gene_start if orientation == '+' else gene_end  # fallback to gene boundary if no CDS found
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
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No 5'UTR annotated. {gene} end used for TSS prediction."
			elif coverage_walk_origin == 'cds':
				if gene not in five_utr_dic and orientation == '+':
					five_utr_dic[gene] = f"No CDS annotated. {gene} start used for TSS prediction."
				if gene not in five_utr_dic and orientation == '-':
					five_utr_dic[gene] = f"No CDS annotated. {gene} end used for TSS prediction."
			tflank_start = time.perf_counter()
			if protein_encoding == "no":
				upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes( gene, gene_infos, genes_per_chromosome, window )
			else:
				upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes(gene, gene_infos,protein_coding_genes_per_chromosome,window)
			tflank_end = time.perf_counter()
			tflank += (tflank_end - tflank_start)
			avg_cov_gene = sum(cov_per_contig[gene_start - 1:gene_end]) / (gene_end - (gene_start - 1))  # get average coverage of the gene of interest for confidence thresholding
			texp_start = time.perf_counter()
			gene_exp_status = find_gene_exp_level(intergenic_window_coverages, genes_per_chromosome, coverage_lookup, gene_start, gene_end, gene, intergenic_region_size, background_percentage)
			texp_end = time.perf_counter()
			texp += (texp_end - texp_start)
			gene_exp_status_dic[gene] = gene_exp_status
			#check if goi is an overlapping gene
			#TSS-blocking overlap types
			SKIP_TSS_TYPES = {'head_head', 'head_into_neighbor', 'same_strand', 'nested'}#'tail_tail' and 'tail_into_neighbor' overlap types are safe for TSS analysis
			goi_strand = gene_infos[gene]['orientation']
			blocking_overlaps = 0
			toverlap_start = time.perf_counter()
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
			toverlap_end = time.perf_counter()
			toverlap += (toverlap_end - toverlap_start)
			if blocking_overlaps > 0:
				print(f"{gene} has {blocking_overlaps} overlap(s). TSS analysis skipped.")
				continue
			if gene_exp_status == None or gene_exp_status == 'low':
				print(f"{gene} shows lower expression level with respect to the background. TSS analysis skipped.")
				continue
			if orientation == "+":	#only works on forward strand
				fig_file = output_folder + gene + ".png"
				if upstream_gene:
					hard_cutoff = gene_infos[ upstream_gene ]['end']

				else:
					hard_cutoff = 1
				print(f"hardcutoff of {gene} is set to {hard_cutoff}")
				tanalysis_start = time.perf_counter()
				result, walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant, accelerated_tss_yr_compliant, other_tool_tss_yr_compliant, hard_cutoff_reached = run_fwd_analysis( other_tool, other_tss_dic, ks_pval, strength, lookahead, output_folder, background_percentage, intergenic_region_size, slide_step, intergenic_window_coverages, coverage_lookup, gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos, percentdiff_threshold )
				tanalysis_end = time.perf_counter()
				tanalysis += (tanalysis_end - tanalysis_start)
				if not other_tss_dic:
					tss_compare_dic[gene] = [walk_tss, basal_tss, elevated_tss, accelerated_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant,accelerated_tss_yr_compliant, hard_cutoff_reached]
				elif other_tss_dic:
					tss_compare_dic[gene] = [walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant,elevated_tss_yr_compliant, accelerated_tss_yr_compliant, other_tool_tss_yr_compliant, hard_cutoff_reached]
				tss_list = {}
				if basal_tss:
					tss_list['basal'] = basal_tss
				else:
					tss_list['basal'] = None
				if elevated_tss:
					tss_list['elevated'] = elevated_tss
				else:
					tss_list['elevated'] = None
				if accelerated_tss:
					tss_list['accelerated'] = accelerated_tss
				else:
					tss_list['accelerated'] = None
				if other_tool_tss:
					tss_list[other_tool] = other_tool_tss
				else:
					tss_list[other_tool] = None
				full_seq_pos_strand_gene[gene]={}
				promoter_status_dic[gene] = {}
				promoter_dic[gene] = {}
				motif_scores[gene] = {}
				percentile_dic[gene] = {}
				canonical_hits_dic[gene] = {}
				textract_start = time.perf_counter()
				for tss_type in tss_list:
					if tss_list[tss_type]:
						promoter_status, promoter, downstream_to_tss, full_seq = extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, tss_list[tss_type], orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size )
						full_seq_pos_strand_gene[gene][tss_type] = full_seq
						promoter_status_dic[gene][tss_type] = promoter_status
						promoter_dic[gene][tss_type] = promoter
						if promoter_analysis == 'yes':
							if os.path.exists(moods):
								cumulative_promoter_motif_score, percentile_of_promoter_score, canonical_hits = promoter_motif_analysis(numcols, upstream_slice, downstream_slice, bg_scores, tss_list[tss_type], tss_type, gene, orientation,promoter, downstream_to_tss, moods, pvalue, pfm_config_dic, tmp_folder,output_folder)
								motif_scores[gene][tss_type] = cumulative_promoter_motif_score
								percentile_dic[gene][tss_type] = percentile_of_promoter_score
								canonical_hits_dic[gene][tss_type] = canonical_hits
							else:
								print('MOODS not found. Promoter analysis not possible.')
					else:
						full_seq_pos_strand_gene[gene][tss_type] = None
						promoter_status_dic[gene][tss_type] = None
						promoter_dic[gene][tss_type] = None
						if promoter_analysis == 'yes':
							motif_scores[gene][tss_type] = None
							percentile_dic[gene][tss_type] = None
							canonical_hits_dic[gene][tss_type] = None
				results.update( { gene: result } )
				textract_end = time.perf_counter()
				textract += (textract_end - textract_start)
			else:	#solution for reverse strand genes
				fig_file = output_folder + gene + ".png"
				if downstream_gene:
					hard_cutoff = gene_infos[ downstream_gene ]['start']
				else:
					hard_cutoff = len( seq_per_contig )
				print(f"hardcutoff of {gene} is set to {hard_cutoff}")
				tanalysis_start = time.perf_counter()
				result, walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant, accelerated_tss_yr_compliant, other_tss_yr_compliant, hard_cutoff_reached = run_rev_analysis( other_tool, other_tss_dic, ks_pval, strength, lookahead, output_folder ,background_percentage, intergenic_region_size, slide_step, intergenic_window_coverages, coverage_lookup, gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites, atg_pos, contig, genome_seq, gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene, five_utr_infos, gene_atg_dic, cds_infos, percentdiff_threshold )
				tanalysis_end = time.perf_counter()
				tanalysis += (tanalysis_end - tanalysis_start)
				if not other_tss_dic:
					tss_compare_dic[gene] = [walk_tss, basal_tss, elevated_tss, accelerated_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant,accelerated_tss_yr_compliant, hard_cutoff_reached]
				elif other_tss_dic:
					tss_compare_dic[gene] = [walk_tss, basal_tss, elevated_tss, accelerated_tss, other_tool_tss, basal_tss_yr_compliant, elevated_tss_yr_compliant, accelerated_tss_yr_compliant,other_tss_yr_compliant, hard_cutoff_reached]
				tss_list = {}
				if basal_tss:
					tss_list['basal'] = basal_tss
				else:
					tss_list['basal'] = None
				if elevated_tss:
					tss_list['elevated'] = elevated_tss
				else:
					tss_list['elevated'] = None
				if accelerated_tss:
					tss_list['accelerated'] = accelerated_tss
				else:
					tss_list['accelerated'] = None
				if other_tool_tss:
					tss_list[other_tool] = other_tool_tss
				else:
					tss_list[other_tool] = None
				full_seq_neg_strand_gene[gene] = {}
				promoter_status_dic[gene] = {}
				promoter_dic[gene] = {}
				motif_scores[gene] = {}
				percentile_dic[gene] = {}
				canonical_hits_dic[gene] = {}
				textract_start = time.perf_counter()
				for tss_type in tss_list:
					if tss_list[tss_type]:
						promoter_status, promoter, downstream_to_tss, full_seq = extract_promoter_region( upstream_slice, downstream_slice, gene, start, end, result, tss_list[tss_type], orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size, downstream_size )
						full_seq_neg_strand_gene[gene][tss_type] = full_seq
						promoter_status_dic[gene][tss_type] = promoter_status
						promoter_dic[gene][tss_type] = promoter
						if promoter_analysis == 'yes':
							if os.path.exists(moods):
								cumulative_promoter_motif_score, percentile_of_promoter_score, canonical_hits = promoter_motif_analysis(numcols, upstream_slice, downstream_slice, bg_scores, tss_list[tss_type], tss_type, gene, orientation,promoter, downstream_to_tss, moods, pvalue, pfm_config_dic, tmp_folder,output_folder)
								motif_scores[gene][tss_type] = cumulative_promoter_motif_score
								percentile_dic[gene][tss_type] = percentile_of_promoter_score
								canonical_hits_dic[gene][tss_type] = canonical_hits
							else:
								print('MOODS not found. Promoter analysis not possible.')
					else:
						full_seq_neg_strand_gene[gene][tss_type] = None
						promoter_status_dic[gene][tss_type] = None
						promoter_dic[gene][tss_type] = None
						if promoter_analysis == 'yes':
							motif_scores[gene][tss_type] = None
							percentile_dic[gene][tss_type] = None
							canonical_hits_dic[gene][tss_type] = None
				results.update( { gene: result } )
				textract_end = time.perf_counter()
				textract += (textract_end - textract_start)
			# calculating confidence score
			isoforms = len(transcripts_per_gene[gene])  # no. of isoforms per gene

			isoforms_dic[gene]=isoforms
			coverage_dic[gene]=avg_cov_gene
		except KeyError as e:
			print(f"Missing gene error: {gene}, missing key: {e}")
			traceback.print_exc()
	t_tss_end_analysis = time.perf_counter()
	time_tss_analysis = t_tss_end_analysis - t_tss_start_analysis
	print(f'time taken for tss analysis of all gois: {time_tss_analysis} seconds')
	print(f'time taken for flanking gene analysis of all gois: {tflank} seconds')
	print(f'time taken for expression level analysis of all gois: {texp} seconds')
	print(f'time taken for overlap type assessment of all gois: {toverlap} seconds')
	print(f'time taken for tss analysis of all gois: {tanalysis} seconds')
	print(f'time taken for promoter extraction of all gois: {textract} seconds')

	# --- report TSS in output file --- #
	final_output_file = os.path.join(output_folder, "Results.tsv")
	final_promoter_analysis_file = os.path.join(output_folder, "Promoter_analysis_results.tsv")
	pos_strand_tss_neighbourhood_file = output_folder + "Positive_strand_gene_promoter_up_downstream_slice_seqs.fasta"
	neg_strand_tss_neighbourhood_file = output_folder + "Negative_strand_gene_promoter_up_downstream_slice_seqs.fasta"
	with open (pos_strand_tss_neighbourhood_file, 'w') as out:
		for gene, tss_type_dic in full_seq_pos_strand_gene.items():
			for tss_type, seq in tss_type_dic.items():
				if seq is not None:
					out.write(f">{gene}_{tss_type}\n{seq}\n")
	with open (neg_strand_tss_neighbourhood_file, 'w') as out:
		for gene, tss_type_dic in full_seq_neg_strand_gene.items():
			for tss_type, seq in tss_type_dic.items():
				if seq is not None:
					out.write(f">{gene}_{tss_type}\n{seq}\n")
	if not other_tss_dic:
		with open( final_output_file, "w" ) as out:
			out.write("\t".join(
				["GeneID", "Average gene coverage", "Gene expression level", "No.of isoforms",
				 "Walk TSS", "Basal TSS", "Elevated TSS","Accelerated TSS",
				 "Basal TSS YR compliance", "Elevated TSS YR compliance","Accelerated TSS YR compliance", "Hard cutoff reached"
				 "Basal promoter status", "Elevated promoter status", "Accelerated promoter status",
				 "Basal promoter sequence", "Elevated promoter sequence", "Accelerated promoter sequence",
				 "Additional comments"]) + "\n")
			for gene in list(results.keys()):
				final_results = []
				if gene in results:
					final_results.extend([gene, str(coverage_dic[gene]), str(gene_exp_status_dic[gene]), str(isoforms_dic[gene]),
					str(tss_compare_dic[gene][0]), str(tss_compare_dic[gene][1]), str(tss_compare_dic[gene][2]), str(tss_compare_dic[gene][3]),
					str(tss_compare_dic[gene][4]), str(tss_compare_dic[gene][5]), str(tss_compare_dic[gene][6]), str(tss_compare_dic[gene][7])])
					if gene in promoter_status_dic:
						for tss_type, status in promoter_status_dic[gene].items():
							final_results.append(str(status))

					if gene in promoter_dic:
						for tss_type, seq in promoter_dic[gene].items():
							final_results.append(str(seq))
					final_results.append(str(five_utr_dic[gene]))
					out.write("\t".join(str(x) for x in final_results) + "\n")#explicitly converting every result entry into string to avoid boolean value errors
		if promoter_analysis == 'yes':
			with open (final_promoter_analysis_file, 'w') as out:
				out.write( "\t".join( [ "Basal promoter motif score", "Elevated promoter motif score", "Accelerated promoter motif score",
										"Basal promoter motif score percentile", "Elevated promoter motif score percentile", "Accelerated promoter motif score percentile",
										"Basal promoter canonical hits", "Elevated promoter canonical hits", "Accelerated promoter canonical hits",])+ "\n")
				for main_gene in list( results.keys() ):
					final_results = []
					for gene, tss_type_dic in motif_scores.items():
						if gene == main_gene:
							for tss_type, score in tss_type_dic.items():
								final_results.extend([score])
							break
					for gene, tss_type_dic in percentile_dic.items():
						if gene == main_gene:
							for tss_type, percentile in tss_type_dic.items():
								final_results.extend([percentile])
							break
					for gene, tss_type_dic in canonical_hits_dic.items():
						if gene == main_gene:
							for tss_type, canonical_hits in tss_type_dic.items():
								final_results.extend([canonical_hits])
							break
					out.write("\t".join(str(x) for x in final_results) + "\n")

	elif other_tss_dic:
		with open( final_output_file, "w" ) as out:
			out.write("\t".join(
				["GeneID", "Average gene coverage", "Gene expression level", "No.of isoforms",
				 "Walk TSS", "Basal TSS", "Elevated TSS","Accelerated TSS", f"{other_tool} TSS"
				 "Basal TSS YR compliance", "Elevated TSS YR compliance","Accelerated TSS YR compliance", f"{other_tool} TSS YR compliance","Hard cutoff reached"
				 "Basal promoter status", "Elevated promoter status", "Accelerated promoter status",
				 "Basal promoter sequence", "Elevated promoter sequence", "Accelerated promoter sequence",
				 "Additional comments"]) + "\n")
			for gene in list(results.keys()):
				final_results = []
				if gene in results:
					final_results.extend([gene, str(coverage_dic[gene]), str(gene_exp_status_dic[gene]), str(isoforms_dic[gene]),
					str(tss_compare_dic[gene][0]), str(tss_compare_dic[gene][1]), str(tss_compare_dic[gene][2]), str(tss_compare_dic[gene][3]),
					str(tss_compare_dic[gene][4]), str(tss_compare_dic[gene][5]), str(tss_compare_dic[gene][6]), str(tss_compare_dic[gene][7]),
					str(tss_compare_dic[gene][8]), str(tss_compare_dic[gene][9])])
					if gene in promoter_status_dic:
						for tss_type, status in promoter_status_dic[gene].items():
							final_results.append(str(status))

					if gene in promoter_dic:
						for tss_type, seq in promoter_dic[gene].items():
							final_results.append(str(seq))
					final_results.append(str(five_utr_dic[gene]))
					out.write("\t".join(str(x) for x in final_results) + "\n")#explicitly converting every result entry into string to avoid boolean value errors
		if promoter_analysis == 'yes':
			with open (final_promoter_analysis_file, 'w') as out:
				out.write( "\t".join( [ "Basal promoter motif score", "Elevated promoter motif score", "Accelerated promoter motif score",f"{other_tool} TSS promoter motif score",
										"Basal promoter motif score percentile", "Elevated promoter motif score percentile", "Accelerated promoter motif score percentile", f"{other_tool} TSS promoter motif score percentile",
										"Basal promoter canonical hits", "Elevated promoter canonical hits", "Accelerated promoter canonical hits",])+ "\n")
				for main_gene in list( results.keys() ):
					final_results = []
					for gene, tss_type_dic in motif_scores.items():
						if gene == main_gene:
							for tss_type, score in tss_type_dic.items():
								final_results.extend([score])
							break
					for gene, tss_type_dic in percentile_dic.items():
						if gene == main_gene:
							for tss_type, percentile in tss_type_dic.items():
								final_results.extend([percentile])
							break
					for gene, tss_type_dic in canonical_hits_dic.items():
						if gene == main_gene:
							for tss_type, canonical_hits in tss_type_dic.items():
								final_results.extend([canonical_hits])
							break
					out.write("\t".join(str(x) for x in final_results) + "\n")

if '--bam' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--scov' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--run_mode' in sys.argv and '--sra_folder' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main(sys.argv)
else:
	sys.exit( __usage__ )
