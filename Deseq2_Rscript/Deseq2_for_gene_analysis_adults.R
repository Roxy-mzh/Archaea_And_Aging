#How to analyze your gene catalog output from ATLAS 
#Normalization output exporting
# Method used for normalization in Deseq: relative log expression (RLE)


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")


#Read the tsv file (gene_annotation) and import as as matrix
gene_data <- as.matrix(read.csv("unnormalized_annotated_genes_adults.csv", header = TRUE, sep = ",", dec = ".", row.names="gene_id"))
View(gene_data)

#Import the metadata file 

metadata <- read.csv("metadata_adults.csv", header = TRUE, sep = ",", dec = ".", row.names="Sample_ID")
View(metadata)

#Define the columns and metadata (condition or type or age, or study, etc)
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

metadata$Methane_production <- factor(metadata$Methane_production)

levels(metadata$Methane_production)


#With the count matrix, gene_data, and the sample information, metadata, we can construct a DESeqDataSet:
#count data is the gene_data (gene count)
#colData is the metadata file
#design is the grouping
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = gene_data,
                              colData = metadata,
                              design = ~ Methane_production)
dds
#To get the normalized counts use one of the followings. Both outputs will be the same. 

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts_adults.csv")

#OR
dds <- DESeq(dds)
sizeFactors(dds)
normalized_counts_DESeq <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_DESeq, file="normalized_counts_adults_v2.csv")

#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

#Pre-filtering
#present in at least N samples (A recommendation for the minimal number of samples is to specify the smallest group size: in my case number of high methane
#The count of 10 is a reasonable choice for bulk RNA-seq (min number of reads) #number of methane producers in adults=32 
smallestGroupSize <- 32 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
#choose a reference level 
dds$Methane_production <- relevel(dds$Methane_production, ref = "yes")

#------------------------------------------------------------------

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Exploring and exporting results
plotMA(res, ylim=c(-2,2))


#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]

#Exporting results to CSV files
write.csv(as.data.frame(resOrdered), file="methane_production_Deseq2_adults_results.csv")
