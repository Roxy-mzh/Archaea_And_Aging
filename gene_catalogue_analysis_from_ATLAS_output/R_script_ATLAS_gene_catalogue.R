# library(plotly) # don't load due to import conflicts
library(heatmaply)
library(dplyr) # dpyr masks select from plotly
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggbeeswarm)
library(pheatmap)
library(grid)
# library(vegan)
library(useful)
library(kableExtra)

library(rhdf5)
library(dplyr)
library(tibble)

library(arrow)
library(yaml)

##### Download all the necessary files from Atlas output and assign to new names:
coverage_stats <- "gene_coverage_stats.parquet"
abundance_file <- "median_coverage.h5"
sample_stats <- read.table ("sample_coverage_stats.tsv", sep='\t', header=T, row.names=1)
kegg <- "kegg.parquet"
readcounts <- "read_counts.tsv" #is located in stats folder
gene_stats <- arrow::read_parquet("gene_coverage_stats.parquet")

#counts: Genecatalog/counts/Nmapped_reads.h5
#geneinfo: Genecatalog/clustering/orf_info.parquet
#eggnog: Genecatalog/annotations/eggNOG.parquet
#cazy: Genecatalog/annotations/dram/cazy.parquet 
#pfam: Genecatalog/annotations/dram/pfam.parquet

# get dimension of data

h5overview <- rhdf5::h5ls(abundance_file)
dim <- h5overview[1, "dim"] %>%
  stringr::str_split(" x ", simplify = T) %>%
  as.numeric()

Ngenes <- dim[1]
Nsamples <- dim[2]

cat("The genecatalog contains", Ngenes, "genes and", Nsamples, "samples.\n")

#Because the dimensions of the genecatalog are huge (even more so with more samples) but
#many genes are detected only in a subset of samples,
#I optimized the file format to allow for fast loading of a subset of the data.

#However we still want information from all the genes.
#The file `r genecatalog_files$sample_stats` contains stats per sample of the genecatalog.
#Especially the number of genes that are detected in each sample and the total coverage which we will use for normalization.

#Similarly the file `r gene_catalog_files$gene_stats` contains stats per gene, e.g. the number of samples in which the gene is detected.

#Let's first look at the stats per sample.

### Gene stats per sample

# check that we have expected column names
assertthat::assert_that("Sum_coverage" %in% colnames(sample_stats))

# get copies per million for normalization (divides the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM))
total_covarage <- sample_stats[, "Sum_coverage"]
names(total_covarage) <- rownames(sample_stats)

head(sample_stats)

# Plot the two numeric columns side by side using facet_wrap
# !!!Not working!!!!
df_melt <- sample_stats %>%
  select(all_of(c("Sum_coverage", "Genes_nz_coverage"))) %>%
  print() %>%
  rename(c(
    Sum_coverage = "Total coverage",
    Genes_nz_coverage = "N detected genes"
  )) %>%
  reshape2::melt()

ggplot(df_melt, aes(x = variable, y = value)) +
  geom_beeswarm(size = 3, alpha = 0.7, cex = 5) +
  facet_wrap(vars(variable), ncol = 2, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme_minimal()

### Stats per gene

#Samples_nz_coverage: Number of samples in which the gene has a non-zero coverage
#Samples_nz_counts: Number of samples in which the gene has a non-zero counts
#Sum_coverage: Sum of the coverage of the gene in all samples
#The values for `Samples_nz_coverage` and `Samples_nz_counts` are not the same
#because if there are only a view reads mapped to a gene but less than halve of the gene is covered the median coverage is zero.
# Assuming your tibble is named 'gene_stats', and you want to create histograms for the specified columns

head(gene_stats)

ggplot(gene_stats, aes(x = log10(Sum_coverage))) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of log(Sum_coverage)",
    x = "log(Sum_coverage)",
    y = "Frequency"
  )

ggplot(gene_stats, aes(x = Samples_nz_coverage)) +
  geom_histogram(binwidth = 10, fill = "green", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Samples_nz_coverage",
    x = "Samples_nz_coverage",
    y = "Frequency"
  )

# Assuming your tibble is named 'gene_stats', and you want to create histograms for the specified columns
ggplot(gene_stats, aes(x = log10(Sum_coverage))) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of log(Sum_coverage)",
    x = "log(Sum_coverage)",
    y = "Frequency"
  )

ggplot(gene_stats, aes(x = Samples_nz_coverage)) +
  geom_histogram(binwidth = 10, fill = "green", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Samples_nz_coverage",
    x = "Samples_nz_coverage",
    y = "Frequency"
  )

# Load the ggplot2 library if not already loaded
library(ggplot2)

# Create a 2D histogram for log10(Sum_coverage) vs. GC
ggplot(gene_stats, aes(x = GC, y = log10(Sum_coverage))) +
  geom_bin2d(bins = 30, aes(fill = ..density..)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "2D Histogram: log10(Sum_coverage) vs. GC",
    x = "GC",
    y = "log10(Sum_coverage)"
  )

# Create a 2D histogram for log10(Sum_coverage) vs. Length
ggplot(gene_stats, aes(x = Length, y = log10(Sum_coverage))) +
  geom_bin2d(bins = 30, aes(fill = ..density..)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "2D Histogram: log10(Sum_coverage) vs. Length",
    x = "Length",
    y = "log10(Sum_coverage)"
  )

# Create a 2D histogram for log10(Sum_coverage) vs. Samples_nz_coverage
ggplot(gene_stats, aes(x = Samples_nz_coverage, y = log10(Sum_coverage))) +
  geom_bin2d(bins = 30, aes(fill = ..density..)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "2D Histogram: log10(Sum_coverage) vs. Samples_nz_coverage",
    x = "Samples_nz_coverage",
    y = "log10(Sum_coverage)"
  )

## Mapping rate

#Let's check the fraction of reads that could be mapped to the genecatalog for each sample.

#```{r, mapping_rate, fig.width=2.5, fig.height=4}

read_stats <- read_tsv("read_counts.tsv", show_col_types = FALSE) %>%
  filter(Step == "QC") %>%
  mutate(Total_reads = Reads_pe * 2 + Reads_se) %>%
  .[, c("Sample", "Total_reads", "Reads_pe", "Reads_se")]

read_stats <- read_stats %>%
  left_join(rownames_to_column(sample_stats, "Sample")) %>%
  mutate(Mapping_rate = Total_counts / Total_reads * 100)



# add metadata to read_stats
read_stats <- read_stats
read_stats[, "x"] <- "Genecatalog"



plt <- ggplot(read_stats, aes(
  y = Mapping_rate,
  x = x, text = paste("Sample:", Sample)
)) +
  ylim(c(50, 100)) +
  # xlim(c(-0.1,0.1)) +
  labs(x = NULL) +
  geom_beeswarm(cex = 2, size = 2) +
  theme_minimal()

print(plt)

# ggplotly(plt,tooltip = c('text','y' ))

# Load abundance of a subset of genes 

#One could load the whole gene catalog like this:

# Load full genecatalog matrix, 

data <- h5read(abundance_file, "data")

attributes= h5readAttributes(abundance_file, "data")

colnames(data) <- attributes$sample_names

#but usually only a subset of genes is of interest.
#The following code allows to load only subset of genes.

#I see two ways to subset the genes:

#1. Load only genes based on gene stats, e.g. genes that are found in at least 10% of the samples
#2. Load only genes that have annotations, e.g. KEGG annotations


# Some helper functions to convert between gene names and gene numbers

GeneNr2GeneName <- function(GeneNumbers, MaxGeneNumber = Ngenes) {
  paste0("Gene", formatC(format = "d", GeneNumbers, flag = "0", width = ceiling(log10(MaxGeneNumber))))
}

GeneName2GeneNr <- function(GeneNames) {
  as.integer(str_sub(GeneNames, start = 5))
}


# helper function to load subset of genes from hdf5 file
load_subset_of_genes <- function(abundance_file, indexes_of_genes_to_load) {
  cat("Load ", length(indexes_of_genes_to_load), " genes\n")


  data <- rhdf5::h5read(
    file = abundance_file, name = "data",
    index = list(indexes_of_genes_to_load, NULL)
  )

  # add sample names
  attributes <- rhdf5::h5readAttributes(abundance_file, "data")
  colnames(data) <- attributes$sample_names
  rownames(data) <- GeneNr2GeneName(indexes_of_genes_to_load)

  return(data)
}


## 1. Load only genes based on gene stats
#here we load only genes that are found in at least 2 samples.


# renmove "Gene" from the gene name to get only the GeneNr
gene_stats["GeneNr"] <- GeneName2GeneNr(gene_stats$GeneName)


# use filtering to make your dataset load-able!
# filter genes that are present in at least 90% of the samples
filtered_genes <- gene_stats %>%
  filter(Samples_nz_coverage > 0.9 * Nsamples) %>%
  dplyr::pull(GeneNr)


# take a subset to speed up caluclation !!!! remove this line
# filtered_genes <- filtered_genes[1:1000]

cat("Select", length(filtered_genes), " genes.\n")
gene_coverage <- load_subset_of_genes(abundance_file, filtered_genes)
gene_coverage[1:5, 1:5]


#Lets calculate the Shannon diversity of the genes. 

shannon_diversity <- vegan::diversity(t(gene_coverage), index = "shannon")
hist(shannon_diversity)

## Normalisation
#Normalize gene counts to **gene counts per million**, which is analogous to transcripts per million.


# This is a complex, but efficient way of deviding each column by the total coverage of the sample and multiplying with 1Million
gene_gcpm <- gene_coverage %*% diag(1 / total_covarage[colnames(gene_coverage)]) * 1e6
colnames(gene_gcpm) <- colnames(gene_coverage)
gene_gcpm[1:5, 1:5]

# Save the normalized data as tsv.gz file
#(if you dont need annotations at this step, you can just go ahead with this table)
#write_delim (as.data.frame(gene_gcpm),
#file.path ("~/R/2023_Roxy_test", "Filtered_genes_cpm.tsv.gz"), delim = "\t")

write.csv(gene_gcpm,"gene_gcpm.csv")

rm(gene_coverage, gene_gcpm, gene_stats)
gc()

## 2. Load only genes with annotations

#Here we load KEGG annotations.
#As alternative one can load the cazy annotation. 

#However, CAZy has multiple annotations per gene.
#I need to write a function how to deal with this.

annotations <- read_parquet(kegg)
annotations["GeneName"] <- GeneNr2GeneName(annotations$GeneNr)

# take a subset to speed up caluclation !!!! remove this line
# annotations <- annotations[1:1000, ]

gene_nrs_with_annotations <- annotations %>%
  pull("GeneNr") %>%
  unique()

annotated_genes <- load_subset_of_genes(abundance_file, gene_nrs_with_annotations)

annotated_genes[1:10, 1:5]
write.csv(annotated_genes,"unnormalized_annotated_genes.csv") #Take this if you want to merge with other datasets, etc)




# Normalize

annotated_genes_gcpm <- annotated_genes %*% diag(1 / total_covarage[colnames(annotated_genes)]) * 1e6
colnames(annotated_genes_gcpm) <- colnames(annotated_genes)

annotated_genes_gcpm[1:5, 1:5]

#In theory, genes can have multiple annotations we need to create a matrix with genes as columns and annotations as rows
#KEGG annotations seem to be unique.

# Create a binary annotation matrix
gene_names_with_annotations <- unique(annotations$GeneName)
ko_ids <- unique(annotations$ko_id)

annotation_matrix <- Matrix::sparseMatrix(
  i = match(annotations$ko_id, ko_ids),
  j = match(annotations$GeneName, gene_names_with_annotations),
  x = 1,
  dims = c(length(ko_ids), length(gene_names_with_annotations))
)

colnames(annotation_matrix) <- gene_names_with_annotations
rownames(annotation_matrix) <- ko_ids

annotation_matrix[1:5, 1:5]

#cat(" Genes have max ", max(colSums(annotation_matrix)), " annotations.\n")

#Now let's multiply the gene coverage with the annotation matrix to get the abundance per annotation.
#**If a gene has multiple annotations the abundance is counted for *each* annotation.**


# make shure the data aligns
assertthat::assert_that(all(colnames(annotation_matrix) == rownames(annotated_genes_gcpm)))

# multiply gene coverage with annotation matrix,
annotation_cpm <- annotation_matrix %*% annotated_genes_gcpm

annotation_cpm <- as.matrix(annotation_cpm)

annotation_cpm[1:5, 1:5]


str(annotation_cpm)
head(annotation_cpm)

annotation_cpm$first_column <- as.character(annotation_cpm$first_column)


# Save the normalized data as csv file

write.csv(annotation_cpm,"annotation_cpm.csv")


# take the log10 of the data
log_annotation_cpm <- log10(annotation_cpm + 1)

# create a heatmap
heatmap(log_annotation_cpm,
        col = colorRampPalette(c("white", "blue"))(100),
        margins = c(10, 10),
        main = "Data in annotations [log10(cmp+1)]"
)

# you need metadata for the next step!
# Differencial abundance with Deseq2


#Load the metadata

metadata <- read.csv(file.path(data_dir, "metadata.csv")) %>%
  column_to_rownames("sample_accession")

head(metadata)

group_variable <- "group"

#Deseq2 requires the unnormalized annotations counts. let's calculate them.

assertthat::assert_that(all(colnames(annotation_matrix) == rownames(annotated_genes)))

# multiply gene coverage with annotation matrix,
annotation_genecounts <- annotation_matrix %*% annotated_genes

annotation_genecounts <- as.matrix(annotation_genecounts)


# Save the normalized data as tsv.gz file
write_delim(as.data.frame(annotation_genecounts),
  file.path(getwd(), "Kegg_counts.tsv.gz"),
  delim = "\t"
)

annotation_genecounts[1:5, 1:5]


#To perform a differential abundance testing in R with abundance data.
#Deseq2 required the unnormalized data.
#You will need the DESeq2 package for performing differential abundance analysis.
#If you haven't already, install and load the package:
  

# BiocManager::install("DESeq2")
library(DESeq2)

#Prepare Data:
  

Description <- annotations %>%
  select(ko_id, kegg_hit) %>%
  distinct(ko_id,.keep_all = TRUE) %>%
  column_to_rownames("ko_id") %>%
  as.data.frame()

# Assuming annotation_cpm is your count data (CPM)
dds <- DESeqDataSetFromMatrix(countData = annotation_genecounts,
                              colData = metadata[colnames(annotation_genecounts),],
                              rowData = Description,
                              design = ~ group)




#Normalize the data using the DESeq2 normalization methods, and estimate differencial abundance.


dds <- DESeq(dds)

res <- results(dds,contrast=c(group_variable,"AD","HC")) %>% as.data.frame()

res %>% merge( Description, by = "row.names", all = TRUE) %>% 
  column_to_rownames("Row.names") -> res

res %>% arrange(padj) %>% select(all_of(c("baseMean", "log2FoldChange", "padj","kegg_hit"))) %>% head(n=10)


ma_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "MA Plot", x = "Mean of normalized counts", y = "log2(Fold Change)")

# Display the MA plot
print(ma_plot)



## Alternative ways to normalize

#An other interesting way to normalize would be to normalize to single copy genes. 
#Example the one from [MuSiCC](https://github.com/borenstein-lab/MUSiCC/blob/master/musicc/data/uscg_76_kegg_min_2313.lst)

#MUSiCC: A marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome. Ohad Manor and Elhanan Borenstein. Genome Biology