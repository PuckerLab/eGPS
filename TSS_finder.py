
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
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[2] == "gene":	#could be extended to other feature types
					ID = parts[-1].split('ID=')[-1]
					if ";" in ID:
						ID = ID.split(';')[0]
					gene_infos.update( { ID: { 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } } )
			line = f.readline()
	return gene_infos


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



def main( arguments ):
	"""! @brief run everything """
	
	#input: BAM/COV file via --in with correct file name extension
	#or BAM file via --bam
	#or COV file via --cov
	
	if '--in' in arguments:
		input_file = arguments[ arguments.index('--in')+1 ]
		try:
			if input_file.split('.')[-1].lower() == "bam":
				bam_file = input_file
				cov_file = ""
				input_mode = "bam"
			elif input_file.split('.')[-1].lower() == "cov":
				cov_file = input_file
				cov_file = ""
				input_mode = "cov"
			elif input_file.split('.')[-1].lower() == "gz" and input_file.split('.')[-2].lower() == "bam":
				cov_file = input_file
				bam_file = ""
				input_mode = "bam_gz"
			elif input_file.split('.')[-1].lower() == "gz" and input_file.split('.')[-2].lower() == "cov":
				cov_file = input_file
				bam_file = ""
				input_mode = "cov_gz"
			
			else:
				sys.exit( "ERROR: unrecognized input file type. If you file is a BAM or COV file, please ensure that it has the proper file name extension (.bam or .cov). Gziped versions are supported (.bam.gz or .cov.gz). #else" )
		except:
			sys.exit( "ERROR: unrecognized input file type. If you file is a BAM or COV file, please ensure that it has the proper file name extension (.bam or .cov). Gziped versions are supported (.bam.gz or .cov.gz). #except" )
	elif '--bam' in arguments:
		bam_file = arguments[ arguments.index('--bam')+1 ]
		if bam_file.split('.')[-1].lower() == "gz":
			input_mode = "bam_gz"
		else:
			input_mode = "bam"
		cov_file = ""
	elif '--cov' in arguments:
		cov_file = arguments[ arguments.index('--cov')+1 ]
		if cov_file.split('.')[-1].lower() == "gz":
			input_mode = "cov_gz"
		else:
			input_mode = "cov"
		bam_file = ""
	
	
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
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	#coverage cutoff
	if '--mincov' in arguments:
		mincov = int( arguments[ arguments.index('--mincov')+1 ] )
	else:
		mincov = 1
	
	#ADD TOOL FOR BAM TO COV CONVERSION
	
	#CONVERT BAM TO COV IF NECESSARY
	
	# --- load data --- #
	coverage = load_coverage( cov_file, input_mode )
	
	gene_infos = load_gene_infos( gff_file )
	
	genome_seq = load_sequences( fasta_file )
	
	#run analysis per gene of interest
	
	#walk coverage upstream of transcription start
	
	#check coverage gaps for (canonical) splice sites to continue across introns
	
	#identify most 5' position with coverage above coverage cutoff
	
	#(invert everything for genes on reverse strand)
	
	#report TSS in output file
	
	#generate figures to visualize coverage around the TSS for manual inspection (coverage plot for whole gene + zoomed version for TSS)


