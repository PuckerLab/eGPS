#!/usr/bin/env python3
"""
Sums one chunk across all samples by seeking directly into each original
.cov file at its byte offset for that chunk.
"""
import sys
import glob
import os
import subprocess
from pathlib import Path


def sum_chunk_seek(cov_dir, offsets_dir, idx, chunk_size, out_dir):
	idx = int(idx)
	chunk_size = int(chunk_size)

	cov_files = sorted(glob.glob(os.path.join(cov_dir, "*.cov")))
	handles = []
	for cov in cov_files:
		base = os.path.splitext(os.path.basename(cov))[0]
		with open(os.path.join(offsets_dir, base + ".offsets")) as f:
			offsets = [int(line) for line in f]
		fh = open(cov, "rb")
		fh.seek(offsets[idx])
		handles.append(fh)

	out_path = os.path.join(out_dir, f"chunk_{idx:02d}.cov")
	with open(out_path, "w") as out:
		for _ in range(chunk_size):
			lines = [h.readline() for h in handles]
			if not lines[0]:
				break  # last chunk can be shorter than chunk_size
			fields0 = lines[0].split()
			chrom, pos = fields0[0].decode(), fields0[1].decode()
			total = sum(int(line.split()[2]) for line in lines)
			out.write(f"{chrom}\t{pos}\t{total}\n")

	for h in handles:
		h.close()


def run_sum_chunk_seek_parallel(cov_dir, offsets_dir, cores, chunk_size, out_dir, parallel):
	"""
	Runs sum_chunk_seek() across ALL <cores> chunk indices via GNU parallel
	"""
	Path(out_dir).mkdir(parents=True, exist_ok=True)
	script_path = Path(__file__).resolve()
	subprocess.run(
		[
			parallel, "-j", str(cores),
			"python3", str(script_path), cov_dir, offsets_dir, "{}", str(chunk_size), out_dir,
			":::", *[str(i) for i in range(int(cores))],
		],
		check=True,
	)


if __name__ == "__main__":
	cov_dir, offsets_dir, idx, chunk_size, out_dir = sys.argv[1:6]
	sum_chunk_seek(cov_dir, offsets_dir, idx, chunk_size, out_dir)
