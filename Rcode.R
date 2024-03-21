# Loading packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("synapter", type = "binary")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

library(AnnotationDbi)
library(org.Hs.eg.db)

library(readxl)
library(dplyr)
library(DESeq2)
library(vsn)
library(readxl)
library(ggplot2)
library(tidyverse)
library(grid)
library(ComplexHeatmap)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(Rtsne)

# Loading Count Matrix (data frame)
count_matrix_orig <- read.csv("BMMF_PBMC_breast_CountMatrix.csv")
class(count_matrix_orig)
dim(count_matrix_orig)
str(count_matrix_orig)

count_mat_df <- count_matrix_orig

# Extract sample names from the first row of each column and set them as column names
colnames(count_mat_df) <- sapply(count_mat_df[1, ], 
                                      function(x) sub("\\./(.*)_.*", "\\1", x))

# Remove the first row after assigning column names
count_mat_df <- count_mat_df[-1, ]

# Set Geneid column as row names
rownames(count_mat_df) <- count_mat_df$Geneid

# Delete the first column
count_mat_df <- count_mat_df[, -1]

# Convert character values to numerical values
count_mat_df <- count_mat_df %>% mutate_at(1:length(count_mat_df),
                                           as.numeric)

str(count_mat_df)
head(count_mat_df)
dim(count_mat_df)


# Loading meta data (originally converted .txt file to .xlsx file)
meta_data <- read_xlsx('SraRunTable.xlsx')
class(meta_data)
metadata_df <- as.data.frame(meta_data)
class(metadata_df)
str(metadata_df)

# Set Samples column as row names
rownames(metadata_df) <- metadata_df$Run
dim(metadata_df)
head(metadata_df)

# Basic Quality Check of Feature Count Matrix and Meta Data
# Check if row names in metadata_df matches col names in count data
all(colnames(count_mat_df) %in% rownames(metadata_df))
# Checking order of row names and column names
all(rownames(metadata_df) == colnames(count_mat_df))


# Building DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = count_mat_df, 
                              colData = metadata_df,
                              design = ~ disease_state)
dds

# Removing Low Counts Reads Genes (leaky expression of genes)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds


#Differential Analysis of Gene (compare cancer with normal)
# Set reference value for DEG analysis 
dds$disease_state <- relevel(dds$disease_state, ref = "normal")

# Perform Differential Analysis of Gene
deg <- DESeq(dds)

# Getting results and sorting as data frame
DESeqRes <- results(deg, alpha = 0.05)
DESeqRes <- DESeqRes[!is.na(DESeqRes$padj),]
DESeqRes <- DESeqRes[order(DESeqRes$padj),]

summary(DESeqRes)
DESeqRes_df <- as.data.frame(DESeqRes)
class(DESeqRes_df)

# Converting Gene ID to Gene Name
# Create new column Symbol for Gene Name in data frame
DESeqRes_df_Genenames <- DESeqRes_df
DESeqRes_df_Genenames$Symbol <- mapIds(org.Hs.eg.db, rownames(DESeqRes_df_Genenames),
                                  keytype = "ENSEMBL",
                                  column = "SYMBOL")

str(DESeqRes_df_Genenames)

#Getting idea about best Genes
best_genes_names <- DESeqRes_df_Genenames %>%
  arrange(padj) %>%
  head(20)

str(best_genes_names)

# Writing to .csv file
write.csv(DESeqRes_df_Genenames, "output_files/DESeqRes_df_Genenames.csv")
write.csv(best_genes_names, "output_files/best_genes_names.csv")



