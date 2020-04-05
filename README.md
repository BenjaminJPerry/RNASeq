# RNASeq Analysis using Rsubread and DESeq2 or ALDEx2

# Description
Template repo for running RNASeq Analysis. The /example directory contains the code used for R7AStar2020 RNASeq analysis.

## Data Type
Strand specific paired-end RNASeq data from dUTP librariy preparation method:  
    * R2 maps to transcript forward strand.
    * R1 maps to transcript reverse strand.  
Reference genome in .fasta format.  
Genome annotation in .gtf format, with features of interest flagged as 'exon'.  
Experiment file in .csv format indicating: SampleID, condition, path/R1.fastq.gz, path/R2.fastq.gz.  
## Overview
### Read Alignment and Feature Counts
    1. Rsubread builds index
    2. Align paired reads
    3. Feature counts generated with fragments allowed to map across multiple features
### Differential Expression Analysis
    DESeq2 - comparison of two conditions indicated in metadata file.
    ALDEx2 - comparison of two conditions indicated in metadata file.
    ALDEx2 - comparison of multiple conditions indicated in metadata file.
