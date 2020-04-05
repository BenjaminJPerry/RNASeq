# RNASeq Analysis using Rsubread and DESeq2 or ALDEx2

# Description
Templea repo for running RNASeq Analysis. The /exmple directory contains the code used for R7AStar2020 RNASeq analysis.

## Data Type
Strand specific paired-end RNASeq data from dUTP method librariy preparation method. Meaning, R2 maps to transcript forward strands and R1 maps to transcript reverse strand.
Reference genome sequence in .fasta format.
Genome annotation prepared in .gtf format, with features of interest flagged as 'exon'.
A .csv experiment file indicating SampleID, condition, path/R1.fastq.gz, path/R2.fastq.gz.

## Overview
### Read Alignment and Feature Counts
    1. Rsubread builds index
    2. Align paired reads
    3. Feature counts generated with fragments allowed to map across multiple features
### Differential Expression Analysis
    DESeq2 - comparison of two conditions indicated in metadata file.
    ALDEx2 - comparison of two conditions indicated in metadata file.
    ALDEx2 - comparison of multiple conditions indicated in metadata file.
