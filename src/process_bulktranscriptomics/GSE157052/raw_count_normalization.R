library(DESeq2)

rm(list = ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

setwd('/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/BT_OMV/THP_1_analysis/GEO157052')

#read raw counts
raw_count <- read.csv('GSE157052_RNAseq_countdata.txt', header = TRUE, sep = ",")
conditions <- factor(c("THP1_0H", "THP1_0H_2", "THP1_0H_3"))

dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = DataFrame(conditions) , design = ~ conditions)