## To run RNASeq analysis:

1. Install tidyverse, DESeq2, and Rsubread R packages.
2. Move read files into /analysis directory (./); make sure read files match those listed in R7AStar.experiment.metadata.csv.
3. Open Terminal and navigate to /analysis directory.
4. Run command: "Rscript --verbose RNASeqAnalysis.R | tee RNASeqAnalysis.R.log".

## Output Files:

### Rsubread alignments and alignment summary files:
* R7A*.subjunc.sort.bam
* R7A.*.subjunc.sort.bam.summary

### Rsubread featureCount summary files:
* featureCounts.summary.txt
* featureCounts.summary.antisense.txt

### DESeq2 analysis summary files:
* R7A.DE.results.rsubread.IndFiltF.cooksF.alpha-0.1.lfcT-0.txt
* R7A.DE.results.rsubread.antisense.IndFiltF.cooksF.alpha-0.1.lfcT-0.txt

### KEGG Mapper files of signifcantly DE genes for mapping:
* R7A.DE.results.rsubread.IndFiltF.cooksF.alpha-0.1.lfcT-0.KEGGMAP.txt
* R7A.DE.results.rsubread.antisense.IndFiltF.cooksF.alpha-0.1.lfcT-0.KEGGMAP.txt

### DESeq2 analysis files filtered for transfer related genes listed in loci.rsubread.txt:
* R7AStar.DE.rsubread.out.summary.txt
* R7AStar.DE.rsubread.antisense.out.summary.txt
