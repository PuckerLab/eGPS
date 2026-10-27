### v0.1 ###
### Shakunthala Natarajan ###

import argparse
import sys, os, glob
import pandas as pd
import numpy as np

def find_sample_name(filepath):
	"""! @brief find sample name """
	basefilename = os.path.basename(filepath)
	sample = basefilename.replace('_reads_aligned.cov','')
	return sample

def get_data_frame(filepath):
	sample = find_sample_name(filepath)
	"""! @brief get data frame """
	df = pd.read_csv(filepath, sep='\t', header=None, names=['col1','col2',sample],dtype={'col1':str,'col2':str,sample:float})
	return df

def merge_files(file_list,output_file):
	"""! @brief merge files """
	anchor = get_data_frame(file_list[0])
	merged = anchor.set_index(['col1','col2'])

	for filepath in file_list[1:]:
		df = get_data_frame(filepath)
		merged = merged.join(df.set_index(['col1','col2']),how='left')#to merge all the rows in each file even when some rows don't have corresponding values that will be handled with NaN; But this won't be the case for legitimate coverage files
	merged.reset_index().to_csv(output_file,sep='\t',index=False)

def compute_avg_gene_cov_cv(file_list, contig_boundary_dic):
	across_samples_cv = []
	for file_path in file_list:
		sum_cov = 0
		num_pos = 0
		with open(file_path, 'r') as file:
			for line in file:
				parts = line.strip().split()
				contig = parts[0]
				pos = int(parts[1])
				cov = float(parts[2])
				if contig in contig_boundary_dic:
					for (start, end) in contig_boundary_dic[contig]:
						if start <= pos <= end:
							sum_cov += cov
							num_pos += 1
		if num_pos > 0:
			avg_gene_cov_per_sample = sum_cov / num_pos
		else:
			avg_gene_cov_per_sample = np.nan
		across_samples_cv.append(avg_gene_cov_per_sample)

	across_samples_cv_arr = np.array(across_samples_cv)
	#compute of CV of avg gene cov across samples as percentage
	CV = float(((np.std(across_samples_cv_arr, ddof=1))/np.mean(across_samples_cv_arr)))*100 if np.mean(across_samples_cv_arr)!=0 else np.nan
	print(across_samples_cv_arr)
	print("Mean =", np.mean(across_samples_cv_arr))
	print("Std  =", np.std(across_samples_cv_arr, ddof=1))
	print(f'CV(%) of avg gene cov across samples is {CV}')
	if CV > 10:
		print('gene is perturbable across the given samples')
	else:
		print('gene is constant across the given samples')
	return CV
def in_line_file_read_write_support_calculation(file_list,contig_pos_dic):
	buffer = []
	buffer_size = 10000
	filehandles = [open(path) for path in file_list]
	support = 0
	cov_value = []
	while True:
		current_rows = []
		for f in filehandles:
			line = f.readline()
			if line == "":
				break#break out of for loop if end of file is reached
			current_rows.append(line)
		if len(current_rows) < len(filehandles):
			break#break out of while loop if end of file is reached
		for row in current_rows:
			if int(row.split('\t')[1]) in contig_pos_dic and row.split('\t')[0]==contig_pos_dic[int(row.split('\t')[1])]:
				cov=float(row.split('\t')[2].rstrip('\n'))
				if cov!=0:
					support+=1
				cov_value.append(cov)
	max_cov_value = max(cov_value)
	sum_norm_cov=0
	for each in cov_value:
		norm_cov = float((each)/max_cov_value)
		sum_norm_cov+=norm_cov
	print(f'TSS sample support is {support}')
	print(f'TSS normalized sample support is {sum_norm_cov}')



def main(arguments):
	inputdir = arguments[arguments.index('--inputdir')+1]
	#output_file = arguments[arguments.index('--out')+1]
	file_list = sorted(glob.glob(os.path.join(inputdir,'*_reads_aligned.cov')))
	#CV = compute_avg_gene_cov_cv(file_list, {'Chr5': [(2803833, 2805606)]} )
	in_line_file_read_write_support_calculation(file_list,{2803971:'Chr5'})

if '--inputdir' in sys.argv:
	main(sys.argv)


