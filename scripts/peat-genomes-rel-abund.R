# 2018-10-11 RStudyGroup Example

library(tidyverse)
library(reshape2)

# Taking three metadata files, merge together, and create an organized relative abundance table

# Read in data
genomestats <- read.delim("Tyson2018-fixed-stats.txt", header=TRUE, sep="\t")
relabund <- read.delim("Tyson2018-fixed-rel-abund.txt", header=TRUE, sep="\t")
ncbi <- read.delim("PRJNA386568_AssemblyDetails.txt", header=TRUE, sep="\t", row.names=NULL)

# rename colnames
colnames(ncbi) <- colnames(ncbi)[-1]
ncbi <- ncbi[,-7]
colnames(genomestats)[10] <- c("BioSample")

# Join stats and assembly files 
ncbistats <- inner_join(ncbi, genomestats, by="BioSample")
test <- left_join(ncbi, genomestats, by="BioSample")
colnames(ncbistats)[7] <- c("genome")

# Drop unecessary columns
subset <- final[,c(7,16,18)]

# Final merged dataset 
final <- inner_join(ncbistats, relabund, by="genome")

# Relative abundance table
fintable <- as.data.frame(spread(subset, key="genome", value="relative_abundance", fill=0))
rownames(fintable) <- fintable[,1]



