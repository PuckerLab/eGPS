### Boas Pucker ###
### pucker@uni-bonn.de ###
__version__ = "v0.0231"

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


import os, sys, subprocess, gzip
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
	
	#add support for gzip-compressed files
	
	coverage_per_seq = {}
	with open( cov_file, "r" ) as f:
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
	genes_per_chromosome = {}
	transcripts_per_gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[2].lower() == "gene":	#could be extended to other feature types
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
			line = f.readline()
	for chromosome in genes_per_chromosome: #sort the genes in each contig/ chromosome in the ascending order of start positions
		genes_per_chromosome[chromosome].sort(key=lambda gene: gene_infos[gene]['start'])
	return gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene

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


def generate_plot( values, svalues, fig_file, atg_pos, genomic_start, genomic_end, gene, orientation ):
	"""! @brief generate a coverage plot """
	
	fig, ax1 = plt.subplots()
	ax1.plot( values, color="black", linestyle="solid" )	#coverage of aligned bases
	ax2 = ax1.twinx()
	ax2.plot( svalues, color="red", linestyle="dotted" )	#coverage of spanning reads
	ax2.plot( [ atg_pos, atg_pos ], [ 0, max( svalues+values ) ], color="green", linestyle="dotted" )	#ATG position
	
	ax1.set_title( gene + "   (" + orientation + ")" )
	ax1.set_xlabel( "position in genomic region from " + str( genomic_start ) + " to " + str( genomic_end ) )
	ax1.set_ylabel( "aligned RNA-seq coverage" )
	ax1.yaxis.label.set_color('black')
	ax2.set_ylabel( "spanning RNA-seq coverage" )
	ax2.yaxis.label.set_color('red')
	
	fig.savefig( fig_file, dpi=300 )


def run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites ):
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
			print( "donor splice site: " + donor_splice_site )
			print( "acceptor splice site: " + acceptor_splice_site )
			print("current position is " + str( current_position ) )
			print ("most upstream position is " + str( most_upstream_pos ) )
			if donor_splice_site == "GT" and acceptor_splice_site == "AG":
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
	plot_end_region = start + flank_region_for_plot
	
	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = len( svalues )-flank_region_for_plot
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "+"
	
	try:
		generate_plot( values, svalues, fig_file, atg_pos, genomic_start, genomic_end, gene, orientation )
	except:
		print( "ERROR: plot failed" + gene )
		
	return { 'TSS': most_upstream_pos, 'start': start, 'end': end }


def run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites ):
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
			print( "donor splice site: " + donor_splice_site )
			print( "acceptor splice site: " + acceptor_splice_site )
			if donor_splice_site == "AC" and acceptor_splice_site == "CT":	#reverse sequences of GT-AG
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
	plot_start_region = end - flank_region_for_plot
	if most_downstream_pos < ( len( seq_per_contig ) - flank_region_for_plot ):
		plot_end_region = most_downstream_pos + flank_region_for_plot
	else:
		plot_end_region = len( seq_per_contig )
		
	values = cov_per_contig[ plot_start_region:plot_end_region ]
	svalues = scov_per_contig[ plot_start_region:plot_end_region ]
	atg_pos = flank_region_for_plot + 0
	genomic_start, genomic_end = plot_start_region, plot_end_region
	orientation = "-"
	
	try:
		generate_plot( values, svalues, fig_file, atg_pos, genomic_start, genomic_end, gene, orientation )
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
	

def extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size=50, max_promoter_size=1000 ):
	"""! @brief extract promoter region """
	
	tss = result['TSS']
	gene_start = result['start']
	gene_end = result['end']
	if orientation == "+":	#forward strand	
		if hard_cutoff - tss > min_promoter_size:
			if hard_cutoff - tss > max_promoter_size:
				promoter = seq_per_contig[ tss-max_promoter_size:tss ]
				promoter_status = True
			else:
				promoter = seq_per_contig[ hard_cutoff:tss ]
				promoter_status = True
		else:
			#no promoter detected (returning everything upstream of start codon
			promoter = seq_per_contig[ hard_cutoff:gene_start ]
			promoter_status = False
	else:	#reverse strand
		if tss - hard_cutoff > min_promoter_size:
			if tss - hard_cutoff > max_promoter_size:
				promoter = seq_per_contig[ tss:tss+max_promoter_size ]
				promoter_status = True
			else:
				promoter = seq_per_contig[ tss:hard_cutoff ]
				promoter_status = True
		else:
			#no promoter detected (returning everything upstream of start codon
			promoter = seq_per_contig[ gene_end:hard_cutoff ]
			promoter_status = False
	return promoter_status, promoter


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
	
	if '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]
				
		if '--bam_is_sorted' in arguments:
			bam_sorted_status = True
		else:
			bam_sorted_status = False
		
		if '--samtools' in arguments:
			samtools = arguments[ arguments.index( '--samtools' )+1 ]
		else:
			samtools = "samtools"
	
		if '--bedtools' in arguments:
			bedtools = arguments[ arguments.index( '--bedtools' )+1 ]
		else:
			bedtools = "genomeCoverageBed"
	
		if '--m' in arguments:
			m = arguments[ arguments.index( '--m' )+1 ]
		else:
			m = "5000000000"
		
		if '--threads' in arguments:
			t = arguments[ arguments.index( '--threads' )+1 ]
		else:
			t = "4"
		
		if not bam_sorted_status:	#sorting the BAM file if it was not sorted already
			print ("sorting BAM file ...")
			sorted_bam_file = output_folder + "sorted.bam"
			cmd = samtools + " sort -m " + m + " --threads " + t + " " + bam_file + " > " + sorted_bam_file
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
		
	else:
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

	# --- load data --- #
	coverage = load_coverage( cov_file, input_mode )
	scoverage = load_coverage( scov_file, input_mode )
	
	gene_infos, genes_per_chromosome, mrna_infos, transcripts_per_gene = load_gene_infos( gff_file )
	
	genome_seq = load_sequences( fasta_file )
	
	#run analysis per gene of interest
	results = {}
	confidence_score_dic = {}
	for gene in goi:
		try:
			cov_per_contig = coverage[ gene_infos[ gene ]['chromosome'] ]	#get coverage of the sequence that harbours the gene of interest
			scov_per_contig = scoverage[ gene_infos[ gene ]['chromosome'] ]	#get spanning read coverage of the sequence that harbours the gene of interest
			seq_per_contig = genome_seq[ gene_infos[ gene ]['chromosome'] ]	#get the sequence of the contig/pseudochromosome that harbours the gene of interest
			start, end, orientation = gene_infos[ gene ]['start'], gene_infos[ gene ]['end'], gene_infos[ gene ]['orientation']	#get information about gene of interest
			upstream_gene, downstream_gene, upstream_gene_list, downstream_gene_list = find_flanking_genes( gene, gene_infos, genes_per_chromosome, window )
			avg_cov_gene = sum(cov_per_contig[start-1:end])/(end-(start-1)) #get average coverage of the gene of interest for confidence thresholding
			# check if goi is an overlapping gene
			overlaps = 0
			for ugene in upstream_gene_list:
				up_end = gene_infos[ ugene ]['end']
				if start <= up_end:
					overlaps += 1
			#this code section is commented to disable downstream overlap check as upstream overlap check is sufficient condition
			"""
			for dgene in downstream_gene_list:
				down_start = gene_infos[ dgene ]['start']
				if end >= down_start:
					overlaps += 1
			"""
			if overlaps > 0: # skip the TSS analysis for the overlapping goi and continue with the for loop iteration of the enxt goi
				print(f"{gene} is overlapping. TSS analysis skipped for this gene.")
				continue
			if orientation == "+":	#only works on forward strand
				fig_file = output_folder + gene + ".png"
				if upstream_gene:
					hard_cutoff = gene_infos[ upstream_gene ]['end']

				else:
					hard_cutoff = 1
				result = run_fwd_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites )
				promoter_status, promoter = extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size )
				result.update( { 'promoter_status': promoter_status } )
				result.update( { 'promoter': promoter } )
				results.update( { gene: result } )
			else:	#solution for reverse strand genes
				fig_file = output_folder + gene + ".png"
				if downstream_gene:
					hard_cutoff = gene_infos[ downstream_gene ]['start']
				else:
					hard_cutoff = len( seq_per_contig )
				result = run_rev_analysis( gene, cov_per_contig, scov_per_contig, seq_per_contig, start, end, fig_file, mincov, min_exon_size, hard_cutoff, flank_region_for_plot, tolerated_gap, splicesites )
				promoter_status, promoter = extract_promoter_region( result, orientation, hard_cutoff, seq_per_contig, min_promoter_size, max_promoter_size )
				result.update( { 'promoter_status': promoter_status } )
				result.update( { 'promoter': promoter } )
				results.update( { gene: result } )
			# calculating confidence score
			isoforms = len(transcripts_per_gene[gene])  # no. of isoforms per gene
			distance = abs((results[gene]['TSS']) - start)  # distance between predicted TSS and gene start
			total_contribution = avg_cov_gene + distance + isoforms
			coverage_score = (avg_cov_gene / total_contribution)  # contribution of coverage to total score
			distance_score = 1 - (distance / total_contribution)  # contribution of distance to total score; 1 - is used since distance is inversely related to total score
			isoform_score = 1 - (isoforms / total_contribution)  # contribution of no. of isoforms to total score; 1 - is used since no. of isoforms is inversely related to total score
			TSS_confidence_score = (coverage_score * distance_score * isoform_score) ** (1 / 3) # confidence score is a gemotric mean of the individual factor scores; each factor is independent of one another and must contribute well for the overal confidence making gemoetric mean and the multiplicative approach preferred over arithmetic mean and the additive approach
			confidence_score_dic[gene]=TSS_confidence_score
		except KeyError:
			print( "Missing gene error: " + gene )
			# If gene was added to results but confidence calculation failed, add fall back
			if gene in results and gene not in confidence_score_dic:
				confidence_score_dic[gene] = "NA"
	# --- report TSS in output file --- #
	final_output_file = output_folder + "results.txt"
	with open( final_output_file, "w" ) as out:
		out.write( "\t".join( [ "GeneID", "TSS", "TSS_confidence_Score", "Start", "End", "PromoterStatus", "Promoter" ] ) + "\n" )
		for gene in list( results.keys() ):
			out.write( "\t".join( [ 	gene,
												str( results[ gene ]['TSS'] ),
									   			str(confidence_score_dic[gene]),
												str( results[ gene ]['start'] ),
												str( results[ gene ]['end'] ),
												str( results[ gene ]['promoter_status'] ),
												str( results[ gene ]['promoter'] )
										] ) + "\n" )


if '--bam' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--scov' in sys.argv and '--out' in sys.argv and '--goi' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
