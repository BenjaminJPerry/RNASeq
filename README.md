# Template for RNASeq Analysis using Rsubread and DESeq2 or ALDEx2

# Description

# Data Type
Strand specific paired-end RNASeq data from dUTP method librariy preparation method. Meaning, R2 maps to transcript forward strands and R1 maps to transcript reverse strand.  
Reference genome sequence in .fasta format.  
Genome annotation prepared in .gtf format, with features of interest flagged as 'exon'.  
Experiment metadata file which indicates SampleID, condition, path/R1.fastq.gz, path/R2.fastq.gz.  

#Overview  
1. Read Alignment and Feature Counts
    + Rsubread builds index
    + Align paired reads
    + Feature counts generated with fragments allowed to map across multiple features
2. Differential Expression Analysis
    + DESeq2 - comparison of two conditions indicated in metadata file.
    + ALDEx2 - comparison of two conditions indicated in metadata file.
    + ALDEx2 - comparison of multiple conditions indicated in metadata file.
