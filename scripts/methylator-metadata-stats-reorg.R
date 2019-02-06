library(tidyverse)

# Create master metadata table of methylators with taxonomy and quality statistics

checkm = read.table("files/all-meth-checkm-stats.txt", sep=" ", header=FALSE)
metadata = read.table("files/all-metadata-gtdb-names.tsv", sep="\t", header=TRUE)
colnames(checkm) = c("genome_name", "classification", "genome_size", "completeness", "redundancy")

full = left_join(checkm, metadata)
meta = full[, c("genome_name", "original_name", "genome_size", "completeness", "redundancy", "taxonomy", "gtdb_taxonomy")]
meta$genome_size = meta$genome_size / 1000000
meta$genome_size = formatC(meta$genome_size, digits=3, format="fg")

# Medium quality bins for metabolic summaries
medium = meta %>% filter(meta$completeness > 50 & meta$redundancy < 10)

medium$filename = paste(medium$genome_name, ".fna", sep="")

write_delim(medium, "files/medium-qual-methylator-metadata.csv", delim=",")

# High quality bins 
high = meta %>% filter(meta$completeness > 90 & meta$redundancy < 5)

tax.clean = separate(meta, gtdb_taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax.clean = tax.clean[,-13]
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="d__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="p__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern="c__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean,gsub, pattern="o__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern = "f__", replacement=""))
tax.clean = as.data.frame(sapply(tax.clean, gsub, pattern = "g__", replacement=""))

meta.clean = tax.clean[,-c(7,9,10,11,12)]
colnames(meta.clean) = c("genome_name", "original_name", "genome_size", "completeness", "redundancy", "NCBI_taxonomy", "GTDB_phylum")

meta.full = left_join(meta.clean, metadata)
meta.full = meta.full[,-8]
meta.final = meta.full[, c("genome_name", "original_name", "genome_size", "completeness", "redundancy", "NCBI_taxonomy", "gtdb_taxonomy", "GTDB_phylum")]
# Get rid of underscores

# Output tables 
write_delim(meta.final, "files/all-methylator-metadata.tsv", delim="\t")
write_delim(high, "files/high-qual-methylator-metadata.tsv", delim="\t")