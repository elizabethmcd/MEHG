library(tidyverse)
library(reshape2)
library(plyr)

metadata = read.csv("~/Desktop/methylator-metadata.csv")
environment = count(metadata, "Code")
phyla = count(metadata, "Phylum")

# Environments plot
phaColfunc = colorRampPalette(c("palegreen3", "slateblue"))
bar_env = environment %>% ggplot(aes(x=reorder(Code, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(11))) + coord_flip() + scale_y_continuous(limits=c(0,125), expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_env


# phyla plot
phylaColfunc = colorRampPalette(c("khaki1", "orchid4"))
bar_phyla = phyla %>% ggplot(aes(x=reorder(Phylum, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(28))) + coord_flip() + scale_y_continuous(limits=c(0,200), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_phyla

# reference isolate phyla 
references = metadata %>% filter(Code=="Isolate")
reference_phyla = count(references, "Phylum")
bar_ref = reference_phyla %>% ggplot(aes(x=reorder(Phylum, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(4))) + coord_flip() + scale_y_continuous(limits=c(0,30), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_ref

# Save plots
ggsave(file="~/Desktop/environment-barplot.png", bar_env, width=15, height=10, units=c("cm"))
ggsave(file="~/Desktop/phyla-barplot.png", bar_phyla, width=15, height=10, units=c("cm"))
ggsave(file="~/Desktop/reference-phya-barplot.png", bar_ref, width=7, height=4, units=c("cm"))
