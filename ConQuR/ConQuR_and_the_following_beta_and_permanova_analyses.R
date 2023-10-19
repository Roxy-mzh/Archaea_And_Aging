
.libPaths()
install.packages(c("quantreg", "cqrReg", "glmnet", "dplyr", "doParallel", "gplots", "vegan", "ade4", "compositions", "randomForest", "ROCR", "ape", "GUniFrac", "fastDummies"), lib="/Users/roxy/Library/R/arm64/4.3/library")
install.packages("coda.base") 
install.packages("doParallel")
devtools::install_github("wdl2459/ConQuR")
devtools::install_github("wdl2459/ConQuR", build_vignettes = TRUE, force=TRUE)                
#Due to technical issues, always library doParallel together with ConQuR, from an R session:
some_packages <- c('ConQuR', 'doParallel')
lapply(some_packages, library, character.only=TRUE)


library(readr)
library(tibble)



species_data <- read_tsv("OTU_table_modified_for_ConQuR.txt", show_col_types = FALSE)%>%
  column_to_rownames(var = "...1") %>%
  as.data.frame()
view(species_data)
#Define taxa by writing the number of columns (how many taxa do you have from 1 to XXX)

taxa = species_data[, 1:4253]

#Show some of the taxa to see if everything is working
taxa[146:150, 1:5]
#Define the batchid=study column; At this step, always do factor() for discrete variables, and do droplevels() to drop unused factor levels.
species_data$batchid <- as.factor(species_data$batchid)
species_data$batchid <- droplevels(species_data$batchid)
batchid = species_data[, 'batchid']
summary(batchid)

#All other metadata are considered as covariates. Specify them. Numerical metadata should be changed to binary data. 
#always do factor() for discrete variables, and do droplevels() to drop unused factor levels. (for each metadata)
species_data$sex <- as.factor(species_data$sex)
species_data$sex <- droplevels(species_data$sex)

species_data$age_group <- as.factor(species_data$age_group)
species_data$age_group <- droplevels(species_data$age_group)

species_data$high_methane_phenotype <- as.factor(species_data$high_methane_phenotype)
species_data$high_methane_phenotype <- droplevels(species_data$high_methane_phenotype)


covar = species_data[, c('sex', 'age_group', 'high_methane_phenotype')]
summary(covar)

#We first use ConQuR (default fitting strategy) to remove batch effects from the taxa read count table, with the randomly picked reference batch: Batch 0. 
#For the default, standard logistic regression and standard quantile regression are used, interpolation of the piece-wise estimation is not used when estimating the conditional quantile functions (both original and batch-free) of the taxonmic read count for each sample and each taxon.
#option(warn = -1) is used in the vignette to suppress a benign warning from quantile regression indicating the solution is not unique.
#takes a while depending on number of samples and taxa
#Choose batch_ref based on prior knowledge or try several options, there is no default
options(warn=-1)
taxa_corrected1 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="2")

#Show some of the taxa
taxa_corrected1[146:150, 1:5]


#We then apply ConQuR with the penalized fitting strategy
#The penalized fitting strategy tackles the high-dimensional covariate problem, 
#helps to stablize the estimates, and prevents over-correction
#takes long
options(warn=-1)
taxa_corrected2 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="2",
                         logistic_lasso=T, quantile_type="lasso", interplt=T) 

#logistic LASSO regression  e.g., c(T, F)
#types of quantile regression, e.g., c("standard", "lasso", "composite").
#whether use interpolation in the piece-wise estimation, e.g., c(T, F).

#Show some of the taxa
taxa_corrected2[146:150, 1:5]

write.csv(taxa_corrected1, file ="otu_table_ConQuR_default_batch_ref_2.csv")
write.csv(taxa_corrected2, file ="otu_table_ConQuR_penalized_batch_ref_2.csv")


#Check how ConQuR works by several ways
#----------------------------------------------------------------------------
#1

#we compare the corrected taxa read count table to the original one, 
#checking whether the batch effects are eliminated and the key variable's effects are preserved
#Whether the PCoA plot is more unified with reference to batchid

par(mfrow=c(2, 3)) #The pcoA plots will be shown in three columns and two rows 

#Plot PcoA with Bray-curtis distance. The batch numbers should be close in the ConQuR default and Penalized
Plot_PCoA(TAX=taxa, factor=batchid, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, main="ConQuR (Penalized), Bray-Curtis")

#Plot PcoA with the Euclidean distance between clr-transformed compositions. The batch numbers should be close in the ConQuR default and Penalized. 
#Penalized should be the closest
Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")

#--------------------------------------------------------------------------
#2

#Whether the variability explained by batchid (quantified by PERMANOVA R2) is decreased and that explained by the key variable is preserved. (in my case age_group)
#In the original taxa read count table, with the key variable sbp be the "n"th variable in covar: 
#(in my case age group is the 3rd variable in covar)
#The value of R2 for batch should be less than key index in both default and penalized ConQuR as compared with  the original

PERMANOVA_R2(TAX=taxa, batchid=batchid, covariates=covar, key_index=3)

#In the corrected taxa read count table, by ConQuR (default):

PERMANOVA_R2(TAX=taxa_corrected1, batchid=batchid, covariates=covar, key_index=3)
#In the corrected taxa read count table, by ConQuR (penalized):

PERMANOVA_R2(TAX=taxa_corrected2, batchid=batchid, covariates=covar, key_index=3)
#----------------------------------------------------------------------------

