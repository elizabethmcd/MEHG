# Exploration of hgcA hits in UBA dataset
library(dplyr)
library(ggplot2)
library(data.table)

# Methylating bacterial bins detected in Parks 2017 UBA dataset
meth_full<- read.csv("genomes-hgcA-table.csv", sep=",", header=TRUE)
# Parks full metagenomic dataset supplement
env_supp <- read.csv("Parks2017-dataset-metadata.csv", sep=",", header=TRUE)
# Column name of SRA same in both
colnames(meth_full)[3] <- "SRA"
colnames(env_supp)[1] <- "SRA"
master <- left_join(meth_full, env_supp, by="SRA")
# Archaeal methylating bins in UBA dataset
arch_meths <- read.csv("ARCHgenomes-hgcA.csv", sep=",", header=TRUE)
colnames(arch_meths)[3] <- "SRA"
# Archaeal and bacterial bins
all_meths <- bind_rows(meth_full, arch_meths)
# Master file of all meths with metadata
master_meths <- left_join(all_meths, env_supp, by="SRA")

# Genome Bins and Environmental Data
meth_summ <- master_meths %>% select(UBA.Genome.ID,NCBI.Organism.Name,CheckM.Completeness,CheckM.Contamination,Genome.Size..bp.,Experiment.Title,Study.Title)
colnames(meth_summ) <- c("UBA-Bin", "Classification", "Completeness", "Contamination", "Mbp", "Experiment", "Title")
meth_summ$Mbp <- meth_summ$Mbp/1000000
meth_summ$Mbp <- format(meth_summ$Mbp, digits=3, format='fg')
ordered <- meth_summ %>% setorder(-Completeness)
above90 <- meth_summ %>% filter(Completeness > 90)
highqual <- meth_summ %>% filter(Completeness > 90 & Contamination < 5)

# Output files
write.table(bin_meta, file="UBA-methylating-bins-metadata.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(highqual, file="UBA-methylating-high-qual-bins.csv", sep=",", row.names=FALSE, quote=FALSE)
