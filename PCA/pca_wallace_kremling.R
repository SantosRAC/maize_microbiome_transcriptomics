# Following JMMP (LabBCES)

# Set working directory
setwd("/home/santosrac/Repositories/maize_microbiome_transcriptomics/PCA")

# Importing required libraries
library(DESeq2)
#library(wesanderson)
library(ggplot2)

# Expression data (raw counts)
expression_data_day = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/correlations_rnaseq_metataxonomics/kremling_expression_v5_day.tsv",
			  sep="\t", header=TRUE, row.names=1)

# OTU data (raw abundance counts)
wallace_otu_data_day = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/otu_table_merged_counts_day_filtered.tsv",
                          sep="\t", header=TRUE, row.names=1)

# OTU data (genus abundance counts)
wallace_otu_data_all_genus = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/collapse_counts_gtdb/genus_counts.tsv", sep="\t", header=TRUE, row.names=1)

# OTU data (family abundance counts)
wallace_otu_data_all_family = read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/collapse_counts_gtdb/family_counts.tsv", sep="\t", header=TRUE, row.names=1)


# Sample metadata (all samples)
sample_metadata = read.csv("/home/santosrac/Projects/UGA_RACS/PCA_RNAseq_16S/sample_groups_flint_garcia_romay.txt",
			  sep="\t", header=TRUE, row.names=1)

# First need to filter sample metadata to have only the subset of day samples


sample_metadata_filtered <- sample_metadata[rownames(sample_metadata) %in% colnames(expression_data_day), ]

expression_data_day_sorted <- expression_data_day[ , rownames(sample_metadata_filtered)]

wallace_otu_data_day_sorted <- wallace_otu_data_day[ , rownames(sample_metadata_filtered)]

wallace_otu_data_genus_filtered_day <- wallace_otu_data_all_genus[ , colnames(wallace_otu_data_all_genus) %in% rownames(sample_metadata_filtered)]
wallace_otu_data_genus_filtered_day_sorted <- wallace_otu_data_genus_filtered_day[ , rownames(sample_metadata_filtered)]

wallace_otu_data_family_filtered_day <- wallace_otu_data_all_family[ , colnames(wallace_otu_data_all_family) %in% rownames(sample_metadata_filtered)]
wallace_otu_data_family_filtered_day_sorted <- wallace_otu_data_family_filtered_day[ , rownames(sample_metadata_filtered)]

# Visualizing first rows of both dataframes
head(expression_data_day_sorted)
head(wallace_otu_data_day)
head(wallace_otu_data_family_filtered_day_sorted)
head(sample_metadata_filtered)
head(wallace_otu_data_genus_filtered_day_sorted)

# Checking that dimensions are correct
print(dim(expression_data_day_sorted))
print(dim(wallace_otu_data_day))
print(dim(wallace_otu_data_family_filtered_day_sorted))
print(dim(sample_metadata_filtered))
print(dim(wallace_otu_data_genus_filtered_day_sorted))

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

raw2 <- DESeqDataSetFromMatrix(countData = wallace_otu_data_family_filtered_day_sorted,
                              colData = sample_metadata_filtered,
                              design = ~ Subpopulation)

data2 <- estimateSizeFactors(raw2)
vst2 <- varianceStabilizingTransformation(data2)
df_vst2 <- assay(vst2)
write.table(df_vst2, file = paste0("otus_vst_counts_day_family.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

# Make PCA with DESEQ2 function
pca_data2 <- plotPCA(vst2, intgroup = c("Subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar2 <- round(100 * attr(pca_data2, "percentVar"))
p2 <- ggplot(pca_data2, aes(x = PC1, y = PC2, color = as.factor(Subpopulation)))
p2 <- p2 + geom_point(size=3)
p2 <- p2 + labs(title = "Maize genotype subpopulations (Family)", color="Subpopulation")
p2 <- p2 + xlab((paste0("PC1 : ", percentVar2[1],"%")))
p2 <- p2 + ylab((paste0("PC2 : ", percentVar2[2],"%")))
p2 <- p2 + theme_bw()
p2 <- p2 + theme( text = element_text(family = "Times New Roman", size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))

ggsave(p2, filename = "pca_maize_subpopulations_family.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

raw3 <- DESeqDataSetFromMatrix(countData = wallace_otu_data_genus_filtered_day_sorted,
                              colData = sample_metadata_filtered,
                              design = ~ Subpopulation)

data3 <- estimateSizeFactors(raw3)
vst3 <- varianceStabilizingTransformation(data3)
df_vst3 <- assay(vst3)
write.table(df_vst3, file = paste0("otus_vst_counts_day_genus.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

# Make PCA with DESEQ2 function
pca_data3 <- plotPCA(vst3, intgroup = c("Subpopulation"), returnData = TRUE)

# Plot customized PCA
percentVar3 <- round(100 * attr(pca_data3, "percentVar"))
p3 <- ggplot(pca_data3, aes(x = PC1, y = PC2, color = as.factor(Subpopulation)))
p3 <- p3 + geom_point(size=3)
p3 <- p3 + labs(title = "Maize genotype subpopulations (Genus)", color="Subpopulation")
p3 <- p3 + xlab((paste0("PC1 : ", percentVar3[1],"%")))
p3 <- p3 + ylab((paste0("PC2 : ", percentVar3[2],"%")))
p3 <- p3 + theme_bw()
p3 <- p3 + theme( text = element_text(family = "Times New Roman", size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))

ggsave(p3, filename = "pca_maize_subpopulations_genus.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)
