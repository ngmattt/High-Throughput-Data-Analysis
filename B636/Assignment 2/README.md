The first part of the assignment which included data acquisition, fastqc and trimming, indexing and alignment, sam/bam conversion and featureCounts was performed in Quartz HPC. The second part of the assignment which included differential gene expression analysis and GO enrichment was done on local system via RStudio. 

Modules that were loaded into Quartz HPC included: 
```bash
module load sra-toolkit
module load fastqc
module load trimgalore
module load python
module load hisat/2.2.1
module load samtools
module load subread
```

This is the beginning of the first part of the pipeline. 

Data acquisition via sra-toolkit:
```bash
cd /N/project/NGS-JangaLab/Matthew/assignment2

mkdir -p fastq

for srr in SRR22269883 SRR22269882 SRR22269881 SRR22269880 SRR22269879 SRR22269878 SRR22269877 SRR22269876 SRR22269875 SRR22269874 SRR22269873 SRR22269872; do
  fasterq-dump $srr -O fastq --split-files --threads $SLURM_CPUS_PER_TASK
done

```

Fastqc 
```bash
mkdir -p fastqc_reports
fastqc -t 8 -o fastqc_reports/ *.fastq
```

Trimming
```bash
RAW_DIR=/N/project/NGS-JangaLab/Matthew/assignment2/fastq
OUT_DIR=/N/project/NGS-JangaLab/Matthew/assignment2/trimmed_fastq

mkdir -p $OUT_DIR

cd $RAW_DIR

for sample in *_1.fastq
do
  base=$(basename $sample _1.fastq)
  echo "Processing $base ..."
  trim_galore --fastqc --cores $SLURM_CPUS_PER_TASK \
              --output_dir $OUT_DIR \
              $sample
done
```


Next, the reference genome and gtf files were downloaded via wget.
Reference genome was indexed via hisat/2.2.1

```bash
mkdir -p genome
cd genome

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
gunzip hg38.ensGene.gtf.gz


hisat2-build hg38.fa hg38_index
```

Then, HISAT2 was used for aligning SARS-CoV-2 RNA-seq samples
```bash
# Paths
BASE_DIR=/N/project/NGS-JangaLab/Matthew/assignment2
GENOME_INDEX=/N/project/NGS-JangaLab/Matthew/assignment2/genome/hg38_index
FASTQ_DIR=${BASE_DIR}/trimmed_fastq
ALIGN_DIR=${BASE_DIR}/alignments

mkdir -p $ALIGN_DIR

# Alignment
for fq in ${FASTQ_DIR}/*_1_trimmed.fq; do
    base=$(basename $fq _1_trimmed.fq)

    hisat2 -p $SLURM_CPUS_PER_TASK \
        -x $GENOME_INDEX \
        -U $fq \
        -S ${ALIGN_DIR}/${base}.sam

done
```

Taking the sam files from the hisat2 alignment and converting to bam files

```bash
ALIGN_DIR=/N/project/NGS-JangaLab/Matthew/assignment2/alignments
BAM_DIR=/N/project/NGS-JangaLab/Matthew/assignment2/bam

mkdir -p $BAM_DIR

for sam in ${ALIGN_DIR}/*.sam; do
    base=$(basename $sam .sam)
    samtools sort -@ $SLURM_CPUS_PER_TASK -o ${BAM_DIR}/${base}.bam $sam
    samtools index ${BAM_DIR}/${base}.bam
done
```
After conversion to bam files, featureCounts was used for gene-level quantification
```bash
BASE_DIR=/N/project/NGS-JangaLab/Matthew/assignment2
GENOME_DIR=/N/project/NGS-JangaLab/Matthew/assignment2/genome
BAM_DIR=${BASE_DIR}/bam
OUT_DIR=${BASE_DIR}/counts

mkdir -p $OUT_DIR

featureCounts \
    -a ${GENOME_DIR}/hg38.ensGene.gtf \
    -o ${OUT_DIR}/gene_counts.txt \
    -T $SLURM_CPUS_PER_TASK \
    -g gene_id \
    -t exon \
    -s 0 \
    ${BAM_DIR}/*.bam
```
The featureCounts process generated two files: "gene_counts.txt" and "gene_counts.txt.summary" which were then transfered from Quartz to local machine for DESeq2 analysis and functional enrichment.

```bash
scp ngmat@quartz.uits.iu.edu:/N/project/NGS-JangaLab/Matthew/assignment2/counts/gene_counts.txt /N/project/NGS-JangaLab/Matthew/assignment2/counts/gene_counts.txt.summary ./Downloads
```

The DESeq2 analysis was done in RStudio with R version 4.5.1 

These were the libraries installed for this analysis in R
```r
# Load libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
```

Reading and cleaning the data
```r
# Read in count data
countData <- read.delim("/Users/matthewng/Downloads/gene_counts.txt", comment.char = "#")
rownames(countData) <- countData$Geneid
countData <- countData[,7:ncol(countData)]  # keep only sample columns

# Read sample metadata
colData <- read.csv("/Users/matthewng/Downloads/assignment_2_info.csv", stringsAsFactors = TRUE)

# Match sample order
colnames(countData) <- colData$SRA.Accession

# Clean up Condition and Time.Point columns
colData$Condition <- trimws(as.character(colData$Condition))  # remove trailing spaces
colData$Condition <- factor(colData$Condition, levels = c("Mock", "SARS-CoV-2"))

# Convert Time.Point numeric → factor with "H" suffix
colData$Time.Point <- paste0(colData$Time.Point, "H")
colData$Time.Point <- factor(colData$Time.Point, levels = c("24H", "72H"))

# Verify the structure again
str(colData)
table(colData$Condition)
table(colData$Time.Point)
```

DESeq2 analysis
```r
# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = colData,
                              design = ~ Condition + Time.Point)

# Run DESeq normalization & modeling
dds <- DESeq(dds)
resultsNames(dds)

# SARS-CoV-2 vs Mock
res1 <- results(dds, name = "Condition_SARS.CoV.2_vs_Mock", alpha = 0.05)

# Sort by adjusted p-value
res1 <- res1[order(res1$padj), ]

# Filter for significant DEGs: padj < 0.05 and |log2FC| ≥ 1
sig_res1 <- subset(res1, padj < 0.05 & abs(log2FoldChange) >= 1)

# Save results
write.csv(as.data.frame(res1),
          "DEG_Mock_vs_SARS-CoV-2_all.csv",
          row.names = TRUE)

write.csv(as.data.frame(sig_res1),
          "DEG_Mock_vs_SARS-CoV-2_sig_padj0.05_LFC1.csv",
          row.names = TRUE)

# Print summary
cat("Summary: Mock vs SARS-CoV-2\n")
summary(res1)
cat("\nNumber of significant DEGs:", nrow(sig_res1), "\n\n")

# SARS-CoV-2 24H vs 72H

if ("Time.Point_72H_vs_24H" %in% resultsNames(dds)) {
  res2 <- results(dds, name = "Time.Point_72H_vs_24H", alpha = 0.05)
} else {
  res2 <- results(dds, contrast = c("Time.Point", "72H", "24H"), alpha = 0.05)
}

# Sort and filter
res2 <- res2[order(res2$padj), ]
sig_res2 <- subset(res2, padj < 0.05 & abs(log2FoldChange) >= 1)

# Save results
write.csv(as.data.frame(res2),
          "DEG_SARS-CoV-2_24h_vs_72h_all.csv",
          row.names = TRUE)

write.csv(as.data.frame(sig_res2),
          "DEG_SARS-CoV-2_24h_vs_72h_sig_padj0.05_LFC1.csv",
          row.names = TRUE)

# Print summary
cat("Summary: 24H vs 72H\n")
summary(res2)
cat("\nNumber of significant DEGs:", nrow(sig_res2), "\n")

```

GO enrichment
```r
# Prepare significant gene lists

sig_res1 <- res1[which(res1$padj < 0.05 & !is.na(res1$padj)), ]
sig_res2 <- res2[which(res2$padj < 0.05 & !is.na(res2$padj)), ]

genes1 <- gsub("\\..*", "", rownames(sig_res1)) 
genes2 <- gsub("\\..*", "", rownames(sig_res2))

cat("Significant DEGs:\n")
cat("Mock vs SARS-CoV-2:", length(genes1), "\n")
cat("24H vs 72H:", length(genes2), "\n")

# Map Ensembl to Entrez IDs

genes1_map <- bitr(genes1, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes2_map <- bitr(genes2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("\nMapped successfully:\n")
cat("SARS-CoV-2 vs Mock:", nrow(genes1_map), "mapped out of", length(genes1), "\n")
cat("SARS-CoV-2 24H vs 72H:", nrow(genes2_map), "mapped out of", length(genes2), "\n")


# Run GO enrichment

ego1 <- enrichGO(
  gene = genes1_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego2 <- enrichGO(
  gene = genes2_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Save results and visualize

write.csv(as.data.frame(ego1), "GO_Mock_vs_SARS-CoV-2.csv", row.names = FALSE)
write.csv(as.data.frame(ego2), "GO_24H_vs_72H.csv", row.names = FALSE)

# Top categories
library(enrichplot)
dotplot(ego1, showCategory = 15, title = "GO Enrichment: Mock vs SARS-CoV-2")
dotplot(ego2, showCategory = 15, title = "GO Enrichment: 24H vs 72H")
```
 



