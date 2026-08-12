#!/usr/bin/env python3
"""
Concatenates all chunk_XX.cov files (in numeric order) into one final
cumulative.cov
"""
import glob
import os
import shutil
import sys


def concatenate_chunks(chunk_sums_dir, out_path):
	chunk_files = sorted(glob.glob(os.path.join(chunk_sums_dir, "chunk_*.cov")))
	with open(out_path, "wb") as out:
		for chunk_file in chunk_files:
			with open(chunk_file, "rb") as f:
				shutil.copyfileobj(f, out)  # streams in chunks, doesn't load whole file into memory
	return out_path


if __name__ == "__main__":
	chunk_sums_dir, out_path = sys.argv[1:3]
	concatenate_chunks(chunk_sums_dir, out_path)
