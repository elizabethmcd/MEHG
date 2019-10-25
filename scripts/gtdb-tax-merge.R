library(tidyverse)

metadata = read.csv("files/stats/medium-mehg-metadata-stats.csv")
tax = read_delim("files/stats/combined_gtdbk_results.tsv", delim="\t")
colnames(tax)[1] = c("genome_name")
tax_names = tax[,c(1,2)]

merged = left_join(tax_names, metadata)
colnames(merged)[2] = c("gtdb_taxonomy")
df = merged[,c(1,5,3,4,6,2,7:25)]
colnames(df)[5] = c("ncbi_lineage")

write.csv(df, "files/stats/combined-metadata-gtdbtk-taxonomy.csv", quote = FALSE, row.names=FALSE)
