# Figures for assessing hgcA primers

# data files
delta = read.csv("results/primer-results/hgcA-delta-qPCR-primers-hits-metadata.txt", header=FALSE)
arch = read.csv("results/primer-results/arch-hgcA-qPCR-hits-metadata.txt", header=FALSE)
firm = read.csv("results/primer-results/hgcA-firmicutes-qPCR-hits-metadata.txt", header=FALSE)

names = c("file", "name", "phyla")
colnames(delta) = names
colnames(arch) = names
colnames(firm) = names

# tables of phyla names
deltaPhy = table(delta$phyla)
archPhy = table(arch$phyla)
firmPhy = table(firm$phyla)

# labels
deltaLbls = paste(names(deltaPhy), " ", deltaPhy)
archLbls = paste(names(archPhy), " ", archPhy)
firmLbls = paste(names(firmPhy))

# pie chart figures
pie(deltaPhy, labels=deltaLbls)
pie(archPhy, labels=FALSE)
