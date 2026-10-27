"""
Second-pass repositioning of the walk TSS using read-pair identity.

Rationale
---------
When the coverage walk stops at the hard cutoff it did not terminate on evidence --
it ran out of room, because the neighbouring gene's coverage never fell below mincov.
Coverage alone cannot say which gene a read in the intergenic space came from. This is
a major problem for eukaryotes with compact genoems taht have small intergenic regions
like Saccharomyces cerevisiae. The method discussed below is hence applicable only
to compact genome cases and only with paired-end RNA-seq reads.

Matching the read pair origin with the gene it maps to, can help with this. Both mates of a 
fragment come from the same transcript, so a mate that lands unambiguously inside 
the GOI labels its partner as the GOI's, and a mate inside the neighbour labels 
its partner as the neighbour's. Walking back from the cutoff toward the gene, the 
position where the GOI's reads start to outnumber the neighbour's is the honest 
edge of the GOI's own transcript.

Limitation: a mate can only anchor a partner within about one fragment length of the
CDS. Beyond that both mates lie in the gap and nothing is assignable. This suits
compact genomes (short intergenic regions, short 5'UTRs).
"""

import numpy as np
import pysam


def build_mate_identity_counts(bam_file, contig, region_lo, region_hi,
                               goi_start, goi_end, nbr_start, nbr_end,
                               min_mapq=20):
	"""
	ONE pass over the BAM for this region.

	For every genomic position in [region_lo, region_hi] (1-based, inclusive) count:
	  goi_counts[i] = reads covering that position whose MATE lies inside the GOI
	  nbr_counts[i] = reads covering that position whose MATE lies inside the neighbour

	Index i = position - region_lo.

	Reads whose mate lies in neither gene (both mates in the gap, unpaired, mate on a
	different contig) are deliberately NOT counted: they carry no evidence either way.
	"""
	n = region_hi - region_lo + 1
	# difference arrays: +1 at the read's first base, -1 just past its last base,
	# then one cumsum turns them into per-position counts. Avoids a loop per base.
	goi_diff = np.zeros(n + 1, dtype=np.int32)
	nbr_diff = np.zeros(n + 1, dtype=np.int32)

	with pysam.AlignmentFile(bam_file, 'rb') as bam:
		# pysam fetch is 0-based half-open; genomic coordinates are 1-based inclusive
		for read in bam.fetch(contig, max(0, region_lo - 1), region_hi):
			if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
				continue
			if read.mapping_quality < min_mapq:
				continue
			if not read.is_paired or read.mate_is_unmapped:
				continue
			if read.next_reference_id != read.reference_id:
				continue  # mate on another contig -> chimeric/mismapped, not usable

			mate_pos = read.next_reference_start + 1  # 0-based -> 1-based
			if goi_start <= mate_pos <= goi_end:
				diff = goi_diff
			elif nbr_start <= mate_pos <= nbr_end:
				diff = nbr_diff
			else:
				continue  # unassignable

			# get_blocks() gives the ALIGNED blocks only, so a spliced read does not
			# count across the intron it skips. reference_start..reference_end would.
			for blk_start, blk_end in read.get_blocks():
				lo = max(blk_start + 1, region_lo)   # 0-based -> 1-based first base
				hi = min(blk_end, region_hi)         # blk_end is 0-based exclusive == 1-based inclusive
				if hi < lo:
					continue
				diff[lo - region_lo] += 1
				diff[hi - region_lo + 1] -= 1

	return np.cumsum(goi_diff)[:n], np.cumsum(nbr_diff)[:n]


def reposition_by_mate_identity(goi_counts, nbr_counts, region_lo, scan_positions,slide_step,
                                window_size, intron_pos_set=None,
                                gene='', verbose=True):
	"""
	Walk OUTWARD from the gene toward the hard cutoff and stop where the neighbour
	takes over the reads. The last position still owned by the GOI is the boundary.

	Direction matters and is the reason this needs no read-count threshold. Scanning
	inward (cutoff -> gene) returns the FIRST window the GOI wins, so a single stray
	GOI-anchored read sitting in an otherwise empty stretch wins before the real
	boundary is ever reached. Scanning outward never evaluates that stray read,
	because the scan has already stopped at the crossover.

	scan_positions must be ordered FROM the gene TOWARD the cutoff:
	    + strand:  range(start, walk_tss - 1, -1)     (descending)
	    - strand:  range(end,   walk_tss + 1)         (ascending)
	The winning window's LAST element is returned -- its outer edge, the furthest
	point from the gene that the GOI still owns.

	Per window, over the non-intronic positions it contains:
	    sum(goi) >  sum(nbr)  -> the GOI owns this window; keep going
	    sum(goi) <= sum(nbr)  -> the neighbour has taken over; stop
	Note sum(goi) == sum(nbr) == 0 also stops the scan, and that is deliberate: an
	empty stretch means the GOI's fragments no longer reach here, which is exactly
	the boundary being looked for.

	A window that is mostly intronic is SKIPPED rather than decided, so a 5'UTR intron
	does not end the scan prematurely. Positions stay in genomic order, so a window is
	always window_size consecutive bases and never jumps across an intron.
	"""
	n = goi_counts.size
	introns = intron_pos_set or set()

	if len(scan_positions) < window_size:
		if verbose:
			print(f'  [mate identity] scan range shorter than one window; not repositioned')
		return None, None, 0

	last_pass = None
	last_share = None
	n_windows_passed = 0
	n_windows_skipped = 0

	w = 0
	limit = len(scan_positions) - window_size
	while w <= limit:
		idx = [p - region_lo for p in scan_positions[w:w + window_size]
			   if 0 <= p - region_lo < n and p not in introns]
		if len(idx) < window_size // 2:
			n_windows_skipped += 1          # mostly intronic: no verdict, keep going
			w += slide_step
			continue
		g = int(goi_counts[idx].sum())
		b = int(nbr_counts[idx].sum())
		if g > b:
			last_pass = scan_positions[w + window_size - 1]   # outer edge of this window
			last_share = g / (g + b)
			n_windows_passed += 1
			w += slide_step
		else:
			break                            # the neighbour owns it from here outward

	if verbose:
		if last_pass is None:
			print(f'  [mate identity] the GOI does not own even the window at the gene; '
				  f'not repositioned')
		else:
			print(f'  [mate identity] GOI owns {n_windows_passed} window(s) outward from the gene '
				  f'({n_windows_skipped} skipped as intronic); boundary at {last_pass}')

	return last_pass, last_share, n_windows_passed
