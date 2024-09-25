# Following JMMP (LabBCES)

# Set working directory
#setwd("/home/santosrac/Repositories/maize_microbiome_transcriptomics/PCA")
setwd("/home/rsantos/Repositories/maize_microbiome_transcriptomics/poster_presentation_notebooks/plant_center_retreat_2024/")

# Importing required libraries
library(plotly)
library(DESeq2)
library(ggplot2)
library(htmlwidgets)

kremling_expression_data <- read.csv("/home/rsantos/Repositories/maize_microbiome_transcriptomics/poster_presentation_notebooks/plant_center_retreat_2024/expression_counts_day_night_pca.tsv",
                                     sep="\t", header=TRUE, row.names=1)
wallace_otu_data <- read.csv("/home/rsantos/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/summed_day_night_otu_counts_conekt_pca.tsv",
                             sep="\t", header=TRUE, row.names=1)
sample_metadata <- read.csv("/home/rsantos/Repositories/maize_microbiome_transcriptomics/poster_presentation_notebooks/plant_center_retreat_2024/sample_metadata_subpopulations_pca.txt",
                            sep="\t", header=TRUE, row.names=1)

library(dplyr)

sample_metadata_filtered <- filter(sample_metadata, rownames(sample_metadata) %in% colnames(kremling_expression_data))
head(sample_metadata_filtered)
dim(sample_metadata_filtered)

expression_data_sorted <- kremling_expression_data[ , rownames(sample_metadata_filtered)]
wallace_otu_data_sorted <- wallace_otu_data[ , rownames(sample_metadata_filtered)]

head(expression_data_sorted)
dim(expression_data_sorted)
head(wallace_otu_data_sorted)
dim(wallace_otu_data_sorted)

expression_raw <- DESeqDataSetFromMatrix(countData = expression_data_sorted,
                                         colData = sample_metadata_filtered,
                                         design = ~ subpopulation)

expression_data <- estimateSizeFactors(expression_raw)
expression_vst <- varianceStabilizingTransformation(expression_data)
df_vst <- assay(expression_vst)

# Make PCA with DESEQ2 function
pca_data <- plotPCA(expression_vst, intgroup = c("subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(subpopulation)))
p <- p + geom_point(size=3)
p <- p + labs(title = "Maize genotype subpopulations (Gene Expression)", color="subpopulation")
p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
p <- p + theme_bw()
p <- p + theme( text = element_text(family = "Times New Roman", size=22),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) 

saveWidget(ggplotly(p   ), file = "expression_pca.html")

otu_raw <- DESeqDataSetFromMatrix(countData = wallace_otu_data_sorted,
                                  colData = sample_metadata_filtered,
                                  design = ~ subpopulation)

otu_data <- estimateSizeFactors(otu_raw)
otu_vst <- varianceStabilizingTransformation(otu_data)
df_vst <- assay(otu_vst)

# Make PCA with DESEQ2 function
pca_data <- plotPCA(otu_vst, intgroup = c("subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(subpopulation)))
p <- p + geom_point(size=3)
p <- p + labs(title = "Maize genotype subpopulations (OTU Count)", color="subpopulation")
p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
p <- p + theme_bw()
p <- p + theme( text = element_text(family = "Times New Roman", size=22),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) 

saveWidget(ggplotly(p   ), file = "otu_pca.html")
