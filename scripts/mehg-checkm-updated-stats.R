library(tidyverse)

# merge metadata with qual information
metadata = read.csv("files/methylator-updated-hgcA-metadata.csv")
stats = read.delim("files/stats/meth-checkm-stats.tsv", sep="\t", header=FALSE)
colnames(stats) <- c("genome_name", "Completeness", "Contamination")

# merge
fullstats = left_join(stats, metadata)
medstats = fullstats %>% filter(Completeness > 50) %>% filter(Contamination < 15)
highqual = fullstats %>% filter(Completeness > 90) %>% filter(Contamination < 5)

# write out different subsets
write_delim(fullstats, "files/stats/complete-mehg-metadata-stats.csv", delim=",", quote_escape="none")
write_delim(medstats, "files/stats/medium-mehg-metadata-stats.csv", delim=",", quote_escape = "none")
write_delim(highqual, "files/stats/high-qual-mehg-metadata-stats.csv", delim=",", quote_escape = "none")
