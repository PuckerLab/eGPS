# eGPS: eukaryotic Gene Promoter Seeker

## Background

This tool was developed to identify transcription start sites (TSS) in plants based on RNA-seq data coverage. The coverage of aligned reads and the coverage of spanning reads are considered in this analysis. The approach looks for the most 5' position upstream of a start codon (ATG) that is covered by RNA-seq reads. Input sample-specific background baseline coverage determination is achieved via a sliding-window approach . Input for the analysis is an RNA-seq read mapping (BAM file). The genome sequence and positions of genes are required for the TSS identification. 

Depending on the background basal coverage baseline, three possible TSS positions can be reported per gene - basal, elevated and accelerated TSS. Accordingly, the respective promoters are extracted and reported. Coverage plots of the reported TSS regions are also provided as outputs, encompassing sub-plots of gene feature boundaries of mRNA, 5'UTR, CDS retrived from the user-provided GFF3 file.

**Basal TSS position** - First position with coverage in the basal coverage region i.e. this region shows no significant elevation in coverage with respect to the intergenic region average coverage values that are used as a baseline

**Elevated TSS position** - First position with coverage in the elevated coverage region i.e. this region shows significant elevation in coverage with respect to the intergenic region average coverage values that are used as a baseline

**Accelerated TSS position** - In case an elevated TSS is detected, then the elevated TSS and coverage walk origin are taken as boundaries and this region is first investigated to see if the average coverage values per window in this region deviate from uniform distribution via a Kolmogorov-Smirnov (KS) distribution fit test (adopted from 1. Cass, A. A. & Xiao, X. mountainClimber Identifies Alternative Transcription Start and Polyadenylation Sites in RNA-Seq. Cell Systems 9, 393-400.e6 (2019).). In case, the KS-test p-value is statistically significant (ks_pval turns out less than the default threshold of 0.01), then this region is investigated to find the accelerated TSS, the minima point accompanying the steepest increase in coverage.

The rationale behind the multi-regime TSS analysis stems from the heterogeneity in the nature and noise of the input RNA-seq samples and previous reports of multiple TSS sites for eukaryotic genes. Thus, the approach adopted by eGPS helps obtain putative TSS positions for a gene at different coverage levels.

An optional promoter analysis is facilitated by integrating MOODS, that helps determine motif hits in a given sequence, along with rich motif density plots. In case, motifs for certain transcription factor binding sites are already known to exist in the promoter sequences, this can very well be used as a secondary test to determine the most confident promoter sequence from the different TSS position associated promoter sequences, and thereby the relatively higher confidence TSS positions.

## Workflow

<img width="5748" height="5973" alt="tssfinder_workflow drawio" src="https://github.com/user-attachments/assets/41a09b2a-0756-44fb-b00c-5ef5803504f9" />

## Installation

**(1) Manual installation**

```
git clone https://github.com/bpucker/eGPS
```
**Mandatory dependencies**
-> Tools - samtools, bedtools
-> Python libraries - pandas (v2.3.1 or greater), matplotlib (v3.10.5 or greater)

**Optional dependencies**
-> Tools - STAR/ HISAT2, MOODS

**(2) Installation in a conda environment**

This method of installation installs all the dependencies in a conda environment using the environment.yml file in this repository

```
git clone https://github.com/bpucker/eGPS

cd eGPS

conda env create -f environment.yml

conda activate egps

```

## Usage


```
Usage:

  python3 eGPS.py [--sra_folder <READ_FILES> | --bam <BAM_FILE> | --cov <COV_FILE> --scov <SCOV_FILE> ] \
                  --fasta <FASTA_FILE> --gff <GFF_FILE> --out <DIR>

Mandatory:

  --sra_folder                STR       Folder encompassing sub-folders of SRA files
                                        for RNA-seq mapping

  --cov                       STR       Aligned bases coverage file (COV)

  --scov                      STR       Spanning read coverage file (COV)

  --bam                       STR       BAM file to automatically create the
                                        coverage file

  --fasta                     STR       FASTA assembly file for read mapping

  --gff                       STR       GFF file with gene information

  --out                       STR       Output directory

Optional:

  --run_mode                  STR       Mode option for running the script [find_tss]

  --mincov                    STR       Minimal coverage [1]

  --samtools                  STR       Full path to samtools (if not in your $PATH)

  --bedtools                  STR       Full path to bedtools (if not in your $PATH)

  --bam_is_sorted             STR       Do not sort BAM file

  --m                         INT       Memory for sorting via samtools [5000000000]

  --threads                   INT       Number of threads for samtools [4]

  --minexon                   INT       Minimal exon size [10]

  --flanksize                 INT       Flanking region size in plot [50]

  --gapsize                   INT       Coverage gap size [5]

  --splicesites               STR       Handling of splice sites[strict](strict|off)

  --intron_percentile_cutoff  STR       Intron size percentile cutoff for RNA-seq mapping

  --neighbourhood             STR       Number of neighbourhood genes to be considered
                                        for overlapping gene check analysis [5]

  --min_promoter_size         INT       Minimum length of promoter to be extracted [50]

  --max_promoter_size         INT       Maximum length of promoter to be extracted [1000]

  --background                INT       Number of random background sequences to be
                                        considered for motif scoring in
                                        promoter analysis [1000]

  --upstream_slice            INT       Length of promoter region to be considered for
                                        motif analysis [200]

  --downstream_slice          INT       Length of region downstream to the identified
                                        promoter to be considered for
                                        motif analysis [50]

  --aligner                   STR       Option to choose between aligners -
                                        STAR and HISAT2 [STAR]

  --STAR                      STR       Full path to STAR

  --HISAT2                    STR       Full path to HISAT2

  --index_bases               STR       Parameter for genome index generation for
                                        RNA-seq mapping [12]

  --fastq_pattern             STR       SRA FASTQ file's naming pattern
                                        [_pass_1, _pass_2]

  --analyse_promoter          STR       Option to activate promoter analysis [no]

  --moods                     STR       Full path to the MOODS script

  --PFM                       STR       Full path to config file for promoter
                                        motif analysis

  --moods_pval                FLOAT     Moods threshold for discovering
                                        motif hits [0.01]

  --background_percentage     FLOAT     Maximum fraction of intergenic background
                                        fragments that have average coverage values
                                        that meet or exceed the candidate window's
                                        average coverage [0.05]

  --ks_pval                   FLOAT     p-value cut-off for the Kolmogorov
                                        Smirnov test [0.01]

  --background_unit           INT       Length of windows to be considered for sliding
                                        window approach used to infer gene expression
                                        levels and distinguish the coverage region into
                                        basal, elevated and accelerated regions [10]

  --slide                     INT       Progression interval of a sliding window [1]

  --signal_strength           INT       Number of consecutive windows needed to cross
                                        the baseline for marking the
                                        elevated TSS position [3]

  --lookahead                 INT       Number of positions to be looked ahead
                                        for determining the accelerated TSS [20]

  --buffer                    INT       Number of bases to be ignored at the
                                        intergenic region ends [500]
```

## More details

`--cov` full path to the file containing the coverage based on aligned bases.

`--scov` full path to the file containing the coverage based on spanning reads.

`--bam` full path to the BAM file.

`--fasta` full path to the FASTA file containing the genome sequence.

`--gff` full path to the GFF file containing information about all genes.

`--out` full path to output folder. This folder will be created if it does not exist already.

`--run_mode` Two modes are available -> find_tss | make_bam 
             find_tss is the default; make_bam can be activated if you want to perform the RNA-seq mapping and then proceed to TSS analysis

`--mincov` minimal coverage for 5'-UTR identification. Default: 1.

`--samtools` full path to samtools. Default: samtools.

`--bedtools` full path to bedtools genomeCoverageBed. Default: genomeCoverageBed.

`--m` specifies the memory for BAM sorting. Default: 5000000000.

`--threads` specifies the number of threads for sorting of a BAM file via samtools. Default: 4.

`--minexon` specifies the minimal exon size [bp]. Default: 10.

`--flanksize` specifies the flanking size for the plot. Default: 50.

`--gapsize` specifies the gap size in coverage due to sequence variants between RNA-seq reads and reference. Default: 5.

`--splicesites` specifies the handling of putative introns. Modes: strict, off. strict enforces a check for canonical splice sites at the ends of a putative intron. off enables the consideration of introns without canonical splice sites. Default: strict.

`--intron_percentile_cutoff' specifies the percentile cutoff to be used for determining the maximum intron size for RNA-seq mapping. By default 99th percentile of intron sizes is used.

`--neighbourhood` is the number of genes to be considered as neighbours for the gene of interest (GOI) both up and downstream. This facilitates overlap checking analysis in the specified neighbourhood window. By default overlap analysis is done across 5 genes up and downstream, and if one of the following overlap types is found, the TSS analysis is skipped - 

  -> head-head overlap: There is an overlap between the GOI's start and one of the neighbour's start

  -> head into neighbour: The gene start is within the bounds of one of its neighbours

  -> same strand: Overlap between the start and end of two genes in the neighbourhood with one being the GOI and both the genes belonging to the same strand

  -> nested: One of the two genes in a neighbourhood window are nested within the other gene with either of them being the GOI 

`--background`: MOODS provides a score for every motif it finds in the given DNA sequence based on the position frequency matrix (PFM) supplied. But this motif can also be found randomly in any other random DNA sequence and not necessarily in the promoter sequence alone. Hence for every cumulative motif scores computed for a fixed length of bases in and downstream of a promoter sequence, the MOODS analysis is also repeated for n random fragments of the same size from across the genome to compute the percentile of the promoter's reported cumulative motif score against a set of random bacgkround scores. n is 1000 by default.

`--upstream_slice` and `--downstream_slice` are the number of bases to be sliced upstream of the TSS and downstream of it, for MOODS analysis. The cumulative sum of these two values will be used as the length of the random sequence fragments that will be retrieved from across the genome for teh background MOODS analysis. These flags' default values are 200 bp and 50 bp by default.

`--PFM` is the full path to a config TXT file for promoter motif analysis. The structure of this config file is as follows - It is a TAB separated TXT file with the following columns:

   -> **Motif_name  Path_to_PFM_file	Upstream_boundary	Downstream_boundary  Direction_sensitivity**
   
   -> Motif_name is the name of your motif of interest like TATA
   
   -> MOODS accepts PFM files in the JASPAR raw .pfm file format
   
   -> Upstream and downstream boundaries are respectively the canonical or previously reported position boundaries for occurrence of the motif of your interest upstream of and downstream to the determined TSS of the extracted promoter at hand. Both these columns should have integer values. If you want a motif element to be analysed only in the upstream end then put in 0 in the downstream boundary column and vice versa if you want only downstream analysis.

   -> Direction sensitivity means that some motif's are meaningful only if they are discovered in the same strand as the GOI. Specify yes if the motif is direction sensitive and no if it is not.

`--moods_pval` is the p-value cut-off for MOODS to report significant motif hits. The lower this value, stricter is the motif hit analysis. MOODS has a default of 0.001 of this p-value. But in some experimentally verified plant promoter sequences known to have TATA box motifs, this p-value could not retrieve the hits. Hence the p-value threshold for the MOODS analysis within eGPS has been increased.

`--background_percentage` is the maximum fraction of fixed size fragements in the intergenic region that can have the same or greater average coverage value as the candidate window from a gene of interest being analysed for TSS detection. 

`--signal_strength` is the number of consecutive windows that must show an elevated coverage value against the intergenic baseline for a position to be deemed as the start of the elevated coverage region or the elevated TSS

`--background_unit` is the length of a window for the sliding window analysis. Since coverage values are being dealt with in eGPS, taking the coverage of a single position or the average coverage of a large number of positions will be extremely prone to outliers in the highly random coverage value dataset. Hence to leverage the power of average and to also counter against the randomness associated at single position resolution, a window of size 10 bp is chosen as a default while also keeping the computational time in perspective. 

`--slide` determines the number of bases a sliding window progresses across a region. It is set to 1 by default so that the TSS reporting can happen at a single base resolution.

`--lookahead` is the number of bases to be looked ahead to adopt a cumulative coverage sum analysis to determine the minima point accompanying the steepest jump in coverage as the accelerated TSS position.

`--buffer` is the number of bases to be skipped before and after the beginning and end of a gene to obtain a typical intergenic region free of read-through signal interruption from actual exonic regions.

  
## Third party tool references

-> Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352

-> Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033

-> Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635

-> Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4

-> Janne Korhonen, Petri Martinmäki, Cinzia Pizzi, Pasi Rastas, Esko Ukkonen, MOODS: fast search for position weight matrix matches in DNA sequences, Bioinformatics, Volume 25, Issue 23, December 2009, Pages 3181–3182, https://doi.org/10.1093/bioinformatics/btp554

## References
Cite this repository.

