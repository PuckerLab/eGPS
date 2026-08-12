#!/usr/bin/env python3
"""
Parallelizes coverage file generation using bedtools genomeCoverageBed
"""
import subprocess
from pathlib import Path


def parallelize_coverage_file_generation ( bam_files, cores, read_coverage_folder, parallel,bedtools, coverage_type):
	Path(read_coverage_folder).mkdir(parents=True, exist_ok=True)
	if coverage_type == 'aligned':
		flags, suffix = "-d -split -ibam", "aligned"
	elif coverage_type == 'spanning':
		flags, suffix = "-d -ibam", "spanning"
	else:
		raise ValueError(f"Coverage type must be 'aligned' or 'spanning', got {coverage_type}")

	cmd = f'{parallel} -j {cores} "{bedtools} {flags} {{}} > {read_coverage_folder}/{{/.}}_{suffix}.cov" :::: {bam_files}'
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()