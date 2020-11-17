library(DESeq2)

rm(list = ls())


setwd('/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/BT_OMV/THP_1_analysis/GSE157052')

#read raw counts
raw_count <- read.csv('GSE157052_RNAseq_countdata.txt', header = TRUE, sep = ",", row.names = 1)
boxplot(log2(raw_count))

# no metadata table -> coldata = describe a mock tabel with 2 columns (1 - THP_0H---rows matching with conditions 2 - anything, WT / 16H (may more type))
meta <- read.csv('sample_condition.csv', header = TRUE)

#design matrix: melyikhez kepest vizsgalja melyiket
design_list <- as.matrix(c("Samples", "Conditions"))


dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = DataFrame(meta) , design = ~ Conditions)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

boxplot(log2(normalized_counts))

write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)



