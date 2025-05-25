# TSS_finder
This tool was developed to identify transcription start sites (TSS) in plants based on RNA-seq data coverage. The approach looks for the most 5' position upstream of a start codon (ATG) that is covered by RNA-seq reads. The coverage of aligned reads and the coverage of spanning reads are considered in this analysis. Input for the analysis is an RNA-seq read mapping (BAM file). The genome sequence and positions of genes are required for the TSS identification.


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
  --m                INT       Memory for sorting via samtools [5000000000]
  --threads          INT       Number of threads for samtools [4]
  --minexon          INT       Minimal exon size [10]
  --flanksize        INT       Flanking region size in plot [50]
  --gapsize          INT       Coverage gap size [5]
```





# References
Cite this repository.

