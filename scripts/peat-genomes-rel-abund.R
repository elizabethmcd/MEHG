# 2018-10-11 RStudyGroup Example

library(tidyverse)
library(reshape2)

# Taking three metadata files, merge together, and create an organized relative abundance table

# Note, removed quotes around relabund and stats files downloaded from Nature supplementary material with sed
# Read in supplementary files and metadata
genomestats <- read.delim("files/Tyson2018-peat-metadata/Tyson2018-fixed-stats.txt", header=TRUE, sep="\t")
relabund <- read.delim("files/Tyson2018-peat-metadata/Tyson2018-fixed-rel-abund.txt", header=TRUE, sep="\t")
ncbi <- read.delim("files/Tyson2018-peat-metadata/PRJNA386568_AssemblyDetails.txt", header=TRUE, sep="\t", row.names=NULL)
# Methylators accession numbers 
methacc <- read.delim("files/Tyson2018-peat-metadata/peat-meth-accessions.txt", header=FALSE, sep="\t")

# rename colnames
colnames(ncbi) <- colnames(ncbi)[-1]
ncbi <- ncbi[,-7]
colnames(genomestats)[10] <- c("BioSample")
colnames(ncbi)[1] <- c("assembly_accession")

# Duplicated methylator lines 
dedupmeth <- as.data.frame(unique(methacc[,1]))
colnames(dedupmeth) <- c("assembly_accession")

# Join stats and assembly files 
ncbistats <- inner_join(ncbi, genomestats, by="BioSample")
colnames(ncbistats)[7] <- c("genome")

# Final merged dataset 
final <- inner_join(ncbistats, relabund, by="genome")

# Drop unecessary columns
subset <- final[,c(7,16,18)]

# Relative abundance table of all genomes 
fintable <- as.data.frame(spread(subset, key="genome", value="relative_abundance", fill=0))
rownames(fintable) <- fintable[,1]

# Select only methylating genomes identified from assembly accession numbers 
methncbi <- inner_join(ncbi, dedupmeth, by="assembly_accession")
methstats <- inner_join(methncbi, genomestats, by="BioSample")
colnames(methstats)[7] <- c("genome")
methgenomes <- as.data.frame(methstats[,7])
colnames(methgenomes) <- c("genome")
methgenomes$genome <- as.character(methgenomes$genome)
methrelabund <- inner_join(methgenomes, relabund, by="genome")
methrelabund$phylum <- gsub("_.*$", "", methrelabund$genome)

# Drop unecessary columns
colnames(methrelabund)
methsub <- methrelabund[,c(1,2,4)]
methreltable <- as.data.frame(spread(methsub, key="genome", value="relative_abundance", fill=0))

methrelabund %>% ggplot(aes(x=sample, y=relative_abundance, group=genome, colour=phylum)) + geom_point()


                                                                                                                        