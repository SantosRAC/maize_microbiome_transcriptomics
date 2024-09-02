# Following JMMP (LabBCES)

# Set working directory
setwd("/home/santosrac/Projects/UGA_RACS/PCA_RNAseq_16S")

# Importing required libraries
library(DESeq2)
library(wesanderson)
library(ggplot2)

# Expression data (raw counts)
expression_data_day = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/correlations_rnaseq_metataxonomics/kremling_expression_v5_day.tsv",
			  sep="\t", header=TRUE, row.names=1)

# OTU data (raw abundance counts)
wallace_otu_data_day = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/otu_table_merged_counts_day_filtered.tsv",
                          sep="\t", header=TRUE, row.names=1)

# Sample metadata (all samples)
sample_metadata = read.csv("/home/santosrac/Projects/UGA_RACS/PCA_RNAseq_16S/sample_groups_flint_garcia_romay.txt",
			  sep="\t", header=TRUE, row.names=1)

# First need to filter sample metadata to have only the subset of day samples


sample_metadata_filtered <- sample_metadata[rownames(sample_metadata) %in% colnames(expression_data_day), ]

expression_data_day_sorted <- expression_data_day[ , rownames(sample_metadata_filtered)]

wallace_otu_data_day_sorted <- wallace_otu_data_day[ , rownames(sample_metadata_filtered)]

# Visualizing first rows of both dataframes
head(sample_metadata_filtered)
head(expression_data_day_sorted)
head(wallace_otu_data_day)

# Checking that dimensions are correct
print(dim(sample_metadata_filtered))
print(dim(expression_data_day_sorted))
print(dim(wallace_otu_data_day))

raw <- DESeqDataSetFromMatrix(countData = expression_data_day_sorted,
                              colData = sample_metadata_filtered,
                              design = ~ Subpopulation)

data <- estimateSizeFactors(raw)
vst <- varianceStabilizingTransformation(data)
df_vst <- assay(vst)
write.table(df_vst, file = paste0("kremling_vst_counts_day.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

# Make PCA with DESEQ2 function
pca_data <- plotPCA(vst, intgroup = c("Subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(Subpopulation)))
p <- p + geom_point(size=3)
p <- p + labs(title = "Maize genotype subpopulations (Gene Expression)", color="Subpopulation")
p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
p <- p + theme_bw()
p <- p + theme( text = element_text(family = "Times New Roman", size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))

ggsave(p, filename = "pca_maize_subpopulations_expression.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

raw1 <- DESeqDataSetFromMatrix(countData = wallace_otu_data_day_sorted,
                              colData = sample_metadata_filtered,
                              design = ~ Subpopulation)

data1 <- estimateSizeFactors(raw1)
vst1 <- varianceStabilizingTransformation(data1)
df_vst1 <- assay(vst1)
write.table(df_vst1, file = paste0("otus_vst_counts_day.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

# Make PCA with DESEQ2 function
pca_data1 <- plotPCA(vst1, intgroup = c("Subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar1 <- round(100 * attr(pca_data1, "percentVar"))
p1 <- ggplot(pca_data1, aes(x = PC1, y = PC2, color = as.factor(Subpopulation)))
p1 <- p1 + geom_point(size=3)
p1 <- p1 + labs(title = "Maize genotype subpopulations (OTUs)", color="Subpopulation")
p1 <- p1 + xlab((paste0("PC1 : ", percentVar1[1],"%")))
p1 <- p1 + ylab((paste0("PC2 : ", percentVar1[2],"%")))
p1 <- p1 + theme_bw()
p1 <- p1 + theme( text = element_text(family = "Times New Roman", size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))

ggsave(p1, filename = "pca_maize_subpopulations_otu.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

