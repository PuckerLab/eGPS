# TSS_finder
tool for transcription start site (TSS) detection in plants based on RNA-seq data



```
Usage:
  python3 TSS_finder.py [--bam <BAM_FILE> | --cov <COV_FILE> --scov <SCOV_FILE> ] --fasta <FASTA_FILE> --gff <GFF_FILE> --out <DIR>

Mandatory:
  --cov              STR       Aligned reads coverage file (COV)
  --scov             STR       Spanning read coverage file (COV)
  --bam              STR       BAM file to automatically create the coverage file]
  --fasta            STR       FASTA assembly file for read mapping
  --gff              STR       GFF file with gene information
  --out              STR       Output directory

Optional:
  --mincov           STR       Minimal coverage [1]
  --samtools         STR       Full path to samtools (if not in your $PATH)
  --bedtools         STR       Full path to bedtools (if not in your $PATH)
  --bam_is_sorted    STR       Do not sort BAM file
  --m                INT       Memory for sorting via samtools[5000000000]
	--threads <NUMBER_THREADS_FOR_SAMTOOLS_SORTING>[4]
	--minexon <MINIMAL_EXON_SIZE>[10]
	--flanksize <FLANKING_REGION_SIZE>[50]
	--gapsize <COVERAGE_GAP_SIZE>[5]
```

# References
Cite this repository.

