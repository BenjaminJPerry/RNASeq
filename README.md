# RNASeq Analysis using Rsubread and DESeq2 or ALDEx2

# Description
Template repo for running RNASeq Analysis. The /example directory contains the code used for RNASeq analysis.

## Data Type
Strand specific paired-end RNASeq data from dUTP library preparation method
* R2 maps to transcript forward strand.
* R1 maps to transcript reverse strand.

Reference genome in .fna format.
Genome annotation in .gtf format, with features of interest flagged as 'exon'. You can use RSEM to make a gtf file from a gff file, make sure to use the "--make-genes-as-transcript" option when you make it.
Experiment file in .csv format indicating: UniqueSampleID, condition, path/R1.fastq.gz, path/R2.fastq.gz.  
## Overview
### Read Alignment and Feature Counts
    1. Rsubread builds index
    2. Aligns paired reads
    3. Feature counts generated with fragments allowed to map across multiple features

### Differential Expression Analysis
    DESeq2 - comparison of two conditions indicated in metadata file.

### KO Functional Informations
    KO numbers are appended to the DE analysis output file.  
    A KO mapping file is created with sig. DE genes colour coded.
