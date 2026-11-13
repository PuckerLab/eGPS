#!/usr/bin/env python3
"""
Builds per-sample byte-offset indexes in parallel (phase 1 of the chunked,
seek-based summation pipeline)

Uses awk to do the actual per-file scan
(fast, purpose-built for this) and GNU parallel to distribute the work across
cores

"""
import subprocess
import sys
from pathlib import Path


AWK_PROGRAM = """
BEGIN { offset = 0; printed = 0 }
{
	if ((NR - 1) % cs == 0 && printed < cores) { print offset; printed++ }
	offset += length($0) + 1
}
END {
	while (printed < cores) { print offset; printed++ }
	print offset
}
"""


def build_offsets(cov_file, chunk_size, cores, out_path):
	"""Run the awk offset-builder on ONE file (plain or .gz), writing the result to out_path."""
	with open(out_path, "w") as out:
		cat = subprocess.Popen(["gzip", "-cdf", cov_file], stdout=subprocess.PIPE)
		subprocess.run(
			["awk", "-v", f"cs={chunk_size}", "-v", f"cores={cores}", AWK_PROGRAM],
			stdin=cat.stdout, stdout=out, check=True,
		)
		cat.stdout.close()
		if cat.wait() != 0:
			raise subprocess.CalledProcessError(cat.returncode, "gzip -cdf " + cov_file)


def run_build_offsets_parallel(cov_files, chunk_size, cores, offsets_dir, parallel):
	#Runs build_offsets() across ALL cov_files via GNU parallel,

	Path(offsets_dir).mkdir(parents=True, exist_ok=True)
	script_path = Path(__file__).resolve()

	subprocess.run(
		[
			parallel, "-j", str(cores),
			"python3", str(script_path), "{}", str(chunk_size), str(cores),
			f"{offsets_dir}/{{/.}}.offsets",
			":::", *[str(f) for f in cov_files],
		],
		check=True,
	)


if __name__ == "__main__":
	cov_file, chunk_size, cores, out_path = sys.argv[1:5]
	build_offsets(cov_file, int(chunk_size), int(cores), out_path)
