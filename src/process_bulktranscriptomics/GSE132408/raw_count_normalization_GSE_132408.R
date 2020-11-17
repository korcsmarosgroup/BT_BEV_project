library(DESeq2)

rm(list = ls())

setwd('path_for_source_file')

#read raw counts
raw_count <- read.csv('GSE132408_rawcounts.csv', header = TRUE, sep = ",", row.names = 1)
boxplot(log2(raw_count))

# no metadata table -> coldata = describe a mock tabel with 2 columns (1 - samples, 2 - conditions)
meta <- read.csv('sample_conditions.csv', header = TRUE)

dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = DataFrame(meta) , design = ~ Conditions)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

boxplot(log2(normalized_counts))

write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)
