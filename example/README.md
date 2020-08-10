## To run RNASeq analysis:

1. Install tidyverse, DESeq2, and Rsubread R packages.
2. Move reads files into /[analysis] directory (./); make sure reads files match those listed in experiment.metadata.csv and have unique sample names.
3. Open terminal and navigate to /[analysis] directory.
4. Run command: "Rscript --verbose RNASeqAnalysis.R | tee RNASeqAnalysis.R.log".

## Example Analysis directory
### For RNASeq analysis the following dir structure is used:
```bash
 [analysis]/  
 ├── sample1-A_S1_L001_R1_001.fastq.gz  
 ├── sample1-A_S1_L001_R2_001.fastq.gz  
 ├── sample1-B_S2_L001_R1_001.fastq.gz  
 ├── sample1-B_S2_L001_R2_001.fastq.gz  
 ├── sample2-A_S3_L001_R1_001.fastq.gz  
 ├── sample2-A_S3_L001_R2_001.fastq.gz  
 ├── sample2-B_S4_L001_R1_001.fastq.gz  
 ├── sample2-B_S4_L001_R2_001.fastq.gz  
 ├── experiment.metadata.csv  
 ├── README.txt  
 ├── ref/  
 │    ├── Reference.fna  
 │    ├── Reference.gtf  
 │    ├── KO.Detailed.txt  
 └── RNASeqAnalysis.R  
```

## Output Files:

### Rsubread alignments and alignment summary files:
* sample*.subjunc.sort.bam
* sample*.subjunc.sort.bam.summary

### Rsubread featureCount summary files:
* featureCounts.summary.txt

### DE analysis output tables:
* DE.results.rsubread.txt
* DE.results.rsubread.KEGGMAP.txt
