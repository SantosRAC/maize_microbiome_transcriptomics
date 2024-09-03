library(CompoCor)

wallace_otu_counts_day_filtered <- read.csv("/home/santosrac/tmp/otu_table_merged_counts_day_filtered.tsv",
					sep="\t", header=TRUE, row.names=1)
kremling_expression_counts_day <- read.csv("/home/santosrac/tmp/kremling_expression_v5_day.tsv",
					sep="\t", header=TRUE, row.names=1)

wallace_otu_counts_day_filtered_transposed <- t(wallace_otu_counts_day_filtered)
kremling_expression_counts_day_transposed <- t(kremling_expression_counts_day)

if(identical(row.names(wallace_otu_counts_day_filtered_transposed), row.names(kremling_expression_counts_day_transposed))) {
 SparXCC_output <- SparXCC_base(wallace_otu_counts_day_filtered_transposed,
                               kremling_expression_counts_day_transposed,
                               pseudo_count=1,
                               var_min=1e-05,
                               Find_m=FALSE)
 write.table(SparXCC_output$cor, file="SparXCC_output.txt", sep="\t")
}
