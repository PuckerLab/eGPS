# eGPS: eukaryotic Gene Promoter Seeker

## Background

This tool was developed to identify transcription start sites (TSS) in plants based on RNA-seq data coverage. The coverage of aligned reads and the coverage of spanning reads are considered in this analysis. The approach looks for the most 5' position upstream of a start codon (ATG) that is covered by RNA-seq reads. Input sample-specific background baseline coverage determination is achieved via a sliding-window approach . Input for the analysis is an RNA-seq read mapping (BAM file). The genome sequence and positions of genes are required for the TSS identification. 

Depending on the background basal coverage baseline, three possible TSS positions can be reported per gene - basal, elevated and accelerated TSS. Accordingly, the respective promoters are extracted and reported. Coverage plots of the reported TSS regions are also provided as outputs, encompassing sub-plots of gene feature boundaries of mRNA, 5'UTR, CDS retrived from the user-provided GFF3 file.

**Basal TSS position** - First position with coverage in the basal coverage region i.e. this region shows no significant elevation in coverage with respect to the intergenic region average coverage values that are used as a baseline

**Elevated TSS position** - First position with coverage in the elevated coverage region i.e. this region shows significant elevation in coverage with respect to the intergenic region average coverage values that are used as a baseline

**Accelerated TSS position** - In case an elevated TSS is detected, then the elevated TSS and coverage walk origin are taken as boundaries and this region is first investigated to see if the average coverage values per window in this region deviate from uniform distribution via a Kolmogorov-Smirnov distribution fit test (adopted from 1. Cass, A. A. & Xiao, X. mountainClimber Identifies Alternative Transcription Start and Polyadenylation Sites in RNA-Seq. Cell Systems 9, 393-400.e6 (2019).). In case, the KS-test p-value is statistically significant (ks_pval turns out less than the default threshold of 0.01), then this region is investigated to find the accelerated TSS, the minima point accompanying the steepest increase in coverage.

The rationale behind the multi-regime TSS analysis stems from the heterogeneity in the nature and noise of the input RNA-seq samples and previous reports of multiple TSS sites for eukaryotic genes. Thus, the approach adopted by eGPS helps obtain putative TSS positions for a gene at different coverage levels.

An optional promoter analysis is facilitated by integrating MOODS, that helps determine motif hits in a given sequence, along with rich motif density plots. In case, motifs for certain transcription factor binding sites are already known to exist in the promoter sequences, this can very well be used as a secondary test to determine the most confident promoter sequence from the different TSS position associated promoter sequences, and thereby the relatively higher confidence TSS values.

## Usage


```
Usage:
  python3 eGPS.py [--sra_folder <READ_FILES> | --bam <BAM_FILE> | --cov <COV_FILE> --scov <SCOV_FILE> ] --fasta <FASTA_FILE> --gff <GFF_FILE> --out <DIR>

Mandatory:
  --sra_folder       STR       Folder encompassing sub-folders of SRA files for RNA-seq mapping
  --cov              STR       Aligned bases coverage file (COV)
  --scov             STR       Spanning read coverage file (COV)
  --bam              STR       BAM file to automatically create the coverage file
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
  --splicesites      STR       Handling of splice sites[strict](strict|off)
  --intron_percentile_cutoff  STR  Intron size percentile cutoff for RNA-seq mapping
  --neighbourhood    STR      Number of neighbourhood genes to be considered for overlapping gene check analysis
  --min_promoter_size  INT    Minimum length of promoter to be extracted [50]
  --max_promoter_size  INT    Maximum length of promoter to be extracted [1000]
  --background    INT    Number of random background sequences to be considered for motif scoring in promoter analysis [1000]
  --upstream_slice  INT  Length of promoter region to be considered for motif analysis [200]
  --downstream_slice  INT  Length of region downstream to the identified promoter to be considered for motif analysis [50]
  --aligner    STR    Option to choose between aligners HISAT2 and STAR [STAR]
  --STAR    STR    Full path to STAR
  --index_bases  
```

`--cov` full path to the file containing the coverage based on aligned bases.

`--scov` full path to the file containing the coverage based on spanning reads.

`--bam` full path to the BAM file.

`--fasta` full path to the FASTA file containing the genome sequence.

`--gff` full path to the GFF file containing information about all genes.

`--out` full path to output folder. This folder will be created if it does not exist already.

`--mincov` minimal coverage for 5'-UTR identification. Default: 1.

`--samtools` full path to samtools. Default: samtools.

`--bedtools` full path to bedtools genomeCoverageBed. Default: genomeCoverageBed.

`--m` specifies the memory for BAM sorting. Default: 5000000000.

`--threads` specifies the number of threads for sorting of a BAM file via samtools. Default: 4.

`--minexon` specifies the minimal exon size [bp]. Default: 10.

`--flanksize` specifies the flanking size for the plot. Default: 50.

`--gapsize` specifies the gap size in coverage due to sequence variants between RNA-seq reads and reference. Default: 5.

`--splicesites` specifies the handling of putative introns. Modes: strict, off. strict enforces a check for canonical splice sites at the ends of a putative intron. off enables the consideration of introns without canonical splice sites. Default: strict.


## References
Cite this repository.

