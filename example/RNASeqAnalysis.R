list.files()
list.files()
# install.packages("tidyverse")
# install.packages("BiocManager")
# BiocManager::install(pkgs = c("Rsubread", "DESeq2"))

library(tidyverse)
library(Rsubread)
library(DESeq2)

makeCounts <- function(refIndex,read1,read2,bamOut,RefGTF,antisense = F) {
    if (antisense == T) {
      strandedness = 1
    } else {
      strandedness = 2 # NEB dUTP Strandedness
    }
    
    align.stat <- subjunc(
      index = refIndex,
      readfile1 = read1,
      readfile2 = read2,
      output_file = bamOut,
      input_format = "gzFASTQ",
      output_format = "BAM",
      nBestLocations = 1,
      maxFragLength = 800,
      nthreads = 12,
      sortReadsByCoordinates = T
    )
    
    counts.stat <- featureCounts(
      files = bamOut,
      annot.ext = RefGTF,
      isGTFAnnotationFile = T,
      minMQS = 10,
      strandSpecific = strandedness,
      isPairedEnd = T,
      requireBothEndsMapped = T,
      checkFragLength = T,
      maxFragLength = 800,
      allowMultiOverlap = T,
      nthreads = 12
    )
    
    return(counts.stat)
  }

# For R7AStar 2020 RNASeq Analysis The Followig dir structure was used:
# ./
# ├── loci.rsubread.txt
# ├── R7A-392-A_S1_L001_R1_001.fastq.gz
# ├── R7A-392-A_S1_L001_R2_001.fastq.gz
# ├── R7A-392-B_S2_L001_R1_001.fastq.gz
# ├── R7A-392-B_S2_L001_R2_001.fastq.gz
# ├── R7A-Star-4510-A_S3_L001_R1_001.fastq.gz
# ├── R7A-Star-4510-A_S3_L001_R2_001.fastq.gz
# ├── R7A-Star-4510-B_S4_L001_R1_001.fastq.gz
# ├── R7A-Star-4510-B_S4_L001_R2_001.fastq.gz
# ├── R7AStar.experiment.metadata.csv
# ├── README.txt
# ├── ref
# │   ├── R7A2020.fna
# │   ├── R7A2020.fseA.gtf
# │   ├── R7A2020.KO.Detailed.txt
# └── RNASeqAnalysis.R

# Build Reference for subjunc alignment
R7A2020Index <- "ref/R7A2020"
buildindex(basename = "ref/R7A2020",
           reference = "ref/R7A2020.fna")
R7AGTF <- "ref/R7A2020.fseA.gtf"

# Prepare the Experimental Metadata File
RNASEQ <- read_csv(
  file = "R7AStar.experiment.metadata.csv",
  col_names = T,
  col_types = "cfcc",
  trim_ws = T
)
# Out put file paths based on the samples in the metadata file
RNASEQ$BAM <- paste(RNASEQ$sample, "subjunc", "sort", "bam", sep = ".")

# Compute feature counts object for  samples
countsObj <- makeCounts(
  refIndex = R7A2020Index,
  read1 = RNASEQ$R1,
  read2 = RNASEQ$R2,
  bamOut = RNASEQ$BAM,
  RefGTF = R7AGTF,
  antisense = F
)

# Print our featureCounts stats
StatsTable <- as.data.frame(x = countsObj$stat)
write_tsv(StatsTable, path = "featureCounts.summary.txt")

### Prepare DESeqDataSet object
# Tidy the count data
countData <- as.data.frame(countsObj$counts)
colnames(countData) <- RNASEQ$sample
countData$ID <- rownames(countData)
countData <- countData %>% select(ID, R7A392A, R7A392B, R7A4510A, R7A4510B)
condition <- factor(RNASEQ$condition)

# Merge the dds Object
ddsR7A <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition, tidy = T)
head(ddsR7A@assays@data@listData[["counts"]])

### Start DESeq2 DE Analysis
# Filter out features with no transcripts
keep <- rowSums(counts(ddsR7A)) > 1
ddsR7AFilter <- ddsR7A[keep,]
nrow(ddsR7A)
nrow(ddsR7AFilter)

# Make a regular log normalized counts table for cursor comparison
rlogR7A <- rlog(ddsR7AFilter, blind = F)
head(assay(rlogR7A), 10)
rlogR7AsampleDist <- dist(t(assay(rlogR7A)))
# Show the euclidean distance between samples at all genes after normlizaion
rlogR7AsampleDist
# The distances between samples is logical given the samples

# Go ahead with the DE analysis
# DESeq 2 DE computation on the raw counts
ddsR7AFilter <- DESeq(ddsR7AFilter)
# Compute DE expression results object
resR7AFilter <- results(
  ddsR7AFilter,
  contrast = c("condition", "R7AStar", "R7A"),
  cooksCutoff = F,
  test = "Wald",
  independentFiltering = F,
  pAdjustMethod = 'BH'
)

### Prepare Files to Export
resOut <- as.data.frame(resR7AFilter)
resOut$ID <- row.names(resOut)
resOut <- resOut %>% select(ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
normCounts <- as.data.frame(counts(ddsR7AFilter, normalized=T))
normCounts$ID <- row.names(normCounts)
normCounts <- normCounts %>% select(ID,
                                    normalize.counts.R7A392A = R7A392A,
                                    normalize.counts.R7A392B = R7A392B,
                                    normalize.counts.R7A4510A = R7A4510A,
                                    normalize.counts.R7A4510B = R7A4510B)

DEout <- left_join(x = countData, y = normCounts, by = "ID")
DEout <- left_join(x = DEout, y = resOut, by = "ID")

# Read in KO annotations
KO <- read_delim(file = "ref/R7A2020.KO.Detailed.txt", col_names = F, delim = "\t")
KO <- KO %>% select(ID = X1, KO = X2, Details = X3, KEGG.Score = X4, KO.Alt = X5, Alt.Score = X6)

# Upating ID strings in DEout
NewIDS <- DEout$ID
NewIDS <- as.data.frame(str_split(NewIDS, pattern = "-", simplify = T))
DEout$ID <- NewIDS$V2

# Merge annotations with DE out
DEKO <- left_join(DEout, KO, by = "ID")
DEKO$CC <- "#00cc00"

# Assign color codes for criteria here
DEKO$CC <- replace(x = DEKO$CC, DEKO$log2FoldChange > 0, "#ff0000")
DEKO$CC <- replace(x = DEKO$CC, DEKO$log2FoldChange < 0, "#0000ff")
DEKO$CC <- replace(x = DEKO$CC, DEKO$padj > 0.05, "#00cc00")

# Print out master table
write_delim(DEKO, path = "R7A.DE.results.rsubread.lfcT-0.CooksF.IndFiltF.txt", col_names = T, delim = "\t")

# Filter genes with no KO codes and print KEGG Mapper Input File
KEGGMAP <- DEKO %>% filter(!is.na(KO)) %>% select(KO, CC)
write_delim(KEGGMAP, path = "R7A.DE.results.rsubread.lfcT-0.CooksF.IndFiltF.KEGGMAP.txt", col_names = F, delim = "\t")


### Compute DESeq2 Analysis for Antisense Transcription ###

# Out put file paths based on the samples in the metadata file
RNASEQ$BAM <- paste(RNASEQ$sample, "subjunc", "antisense", "sort", "bam", sep = ".")

# Compute feature counts object for  samples
countsObj <- makeCounts(
  refIndex = R7A2020Index,
  read1 = RNASEQ$R1,
  read2 = RNASEQ$R2,
  bamOut = RNASEQ$BAM,
  RefGTF = R7AGTF,
  antisense = T
)

# Print our featureCounts stats
StatsTable <- as.data.frame(x = countsObj$stat)
write_tsv(StatsTable, path = "featureCounts.summary.antisense.txt")

### Prepare DESeqDataSet object
# Tidy the count data
countData <- as.data.frame(countsObj$counts)
colnames(countData) <- RNASEQ$sample
countData$ID <- rownames(countData)
countData <-
  countData %>% select(ID, R7A392A, R7A392B, R7A4510A, R7A4510B)
condition <- factor(RNASEQ$condition)

# Merge the dds Object
ddsR7A <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition, tidy = T)
head(ddsR7A@assays@data@listData[["counts"]])

### Start DESeq2 DE Analysis
# Filter out features with no transcripts
keep <- rowSums(counts(ddsR7A)) > 1
ddsR7AFilter <- ddsR7A[keep,]
nrow(ddsR7A)
nrow(ddsR7AFilter)

# Make a regular log normalized counts table for cursor comparison
rlogR7A <- rlog(ddsR7AFilter, blind = F)
head(assay(rlogR7A), 10)
rlogR7AsampleDist <- dist(t(assay(rlogR7A)))
# Show the euclidean distance between samples at all genes after normlizaion
rlogR7AsampleDist
# The distances between samples is logical given the samples

# Go ahead with the DE analysis
# DESeq 2 DE computation on the raw counts
ddsR7AFilter <- DESeq(ddsR7AFilter)
# Compute DE expression results object
resR7AFilter <- results(
  ddsR7AFilter,
  contrast = c("condition", "R7AStar", "R7A"),
  cooksCutoff = F,
  test = "Wald",
  independentFiltering = F,
  pAdjustMethod = 'BH'
)

### Prepare Files to Export
resOut <- as.data.frame(resR7AFilter)
resOut$ID <- row.names(resOut)
resOut <- resOut %>% select(ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
normCounts <- as.data.frame(counts(ddsR7AFilter, normalized=T))
normCounts$ID <- row.names(normCounts)
normCounts <- normCounts %>% select(ID, normalize.counts.R7A392A = R7A392A,
                                    normalize.counts.R7A392B = R7A392B,
                                    normalize.counts.R7A4510A = R7A4510A,
                                    normalize.counts.R7A4510B = R7A4510B)
DEout <- left_join(x = countData, y = normCounts, by = "ID")
DEout <- left_join(x = DEout, y = resOut, by = "ID")

# Read in KO annotations
KO <- read_delim(file = "ref/R7A2020.KO.Detailed.txt", col_names = F, delim = "\t")
KO <- KO %>% select(ID = X1, KO = X2, Details = X3, KEGG.Score = X4, KO.Alt = X5, Alt.Score = X6)

# Upating ID strings in DEout
NewIDS <- DEout$ID
NewIDS <- as.data.frame(str_split(NewIDS, pattern = "-", simplify = T))
DEout$ID <- NewIDS$V2

# Merge annotations with DE out
DEKO <- left_join(DEout, KO, by = "ID")
DEKO$CC <- "#00cc00"

# Assign color codes for criteria here
DEKO$CC <- replace(x = DEKO$CC, DEKO$log2FoldChange > 0, "#ff0000")
DEKO$CC <- replace(x = DEKO$CC, DEKO$log2FoldChange < 0, "#0000ff")
DEKO$CC <- replace(x = DEKO$CC, DEKO$padj > 0.05, "#00cc00")

# Print out master table
write_delim(DEKO, path = "R7A.DE.results.rsubread.antisense.lfcT-0.CooksF.IndFiltF.txt", col_names = T, delim = "\t")

# Filter genes with no KO codes and print KEGG Mapper Input File
KEGGMAP <- DEKO %>% filter(!is.na(KO)) %>% select(KO, CC)
write_delim(KEGGMAP, path = "R7A.DE.results.rsubread.antisense.lfcT-0.CooksF.IndFiltF.KEGGMAP.txt", col_names = F, delim = "\t")


### filter outputs for making summary table ###

system("grep -f loci.rsubread.txt R7A.DE.results.rsubread.lfcT-0.CooksF.IndFiltF.txt > R7AStar.DE.rsubread.out.summary.txt")
system("grep -f loci.rsubread.txt R7A.DE.results.rsubread.antisense.lfcT-0.CooksF.IndFiltF.txt > R7AStar.DE.rsubread.antisense.out.summary.txt")
system("rm *antisense.sort.bam*")
