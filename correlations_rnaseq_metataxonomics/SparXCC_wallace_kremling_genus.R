library(CompoCor)

wallace_genus_counts <- read.csv("/home/santosrac/Repositories/maize_microbiome_transcriptomics/16S_wallace2018/combine_day_night_samples/collapse_counts_gtdb/genus_counts.tsv",
				        sep="\t", header=TRUE, row.names=1)
kremling_expression_counts_day <- read.csv("/home/santosrac/tmp/kremling_expression_v5_day.tsv",
					sep="\t", header=TRUE, row.names=1)

wallace_genus_counts_day_filtered <- wallace_genus_counts[ , colnames(wallace_genus_counts) %in% colnames(kremling_expression_counts_day)]
wallace_genus_counts_day_sorted <- wallace_genus_counts_day_filtered[ , colnames(kremling_expression_counts_day)]

dim(wallace_genus_counts_day_sorted)
dim(kremling_expression_counts_day)

print(head(wallace_genus_counts_day_sorted))
print(head(kremling_expression_counts_day))

wallace_genus_counts_day_filtered_transposed <- t(wallace_genus_counts_day_sorted)
kremling_expression_counts_day_transposed <- t(kremling_expression_counts_day)

if(identical(row.names(wallace_genus_counts_day_filtered_transposed), row.names(kremling_expression_counts_day_transposed))) {
 SparXCC_output <- SparXCC_base(wallace_genus_counts_day_filtered_transposed,
                               kremling_expression_counts_day_transposed,
                               pseudo_count=1,
                               var_min=1e-05,
                               Find_m=FALSE)
 write.table(SparXCC_output$cor, file="SparXCC_genus_output.txt", sep="\t")
}
