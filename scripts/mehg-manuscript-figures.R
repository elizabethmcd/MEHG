library(tidyverse)
library(reshape2)
library(plyr)
library(ggplot2)
library()

metadata = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/methylator-metadata.csv")
environment = count(metadata, "Code")
group = count(metadata, "Group")
study = count(metadata, "Study")

# Environments plot
phaColfunc = colorRampPalette(c("palegreen3", "slateblue"))
bar_env = environment %>% ggplot(aes(x=reorder(Code, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(11))) + coord_flip() + scale_y_continuous(limits=c(0,125), expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_env 


# phyla plot
phylaColfunc = colorRampPalette(c("khaki1", "orchid4"))
bar_phyla = group %>% ggplot(aes(x=reorder(Group, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(15))) + coord_flip() + scale_y_continuous(limits=c(0,200), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_phyla

# reference isolate phyla 
references = metadata %>% filter(Code=="Isolate")
reference_phyla = count(references, "Phylum")
bar_ref = reference_phyla %>% ggplot(aes(x=reorder(Phylum, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(4))) + coord_flip() + scale_y_continuous(limits=c(0,30), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_ref

# CP phyla
cp = metadata %>% filter(Group=="CP")
cp_phyla = count(cp, "Phylum")
bar_cp= cp_phyla %>% ggplot(aes(x=reorder(Phylum, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(9))) + coord_flip() + scale_y_continuous(breaks=c(0:10), expand= c(0,0))
bar_cp 

# Save plots
ggsave(file="/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/figs/2019-06-17-figs/environment-barplot.png", bar_env, width=20, height=10, units=c("cm"))
ggsave(file="/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/figs/2019-06-17-figs/phyla-barplot.png", bar_phyla, width=20, height=10, units=c("cm"))
ggsave(file="/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/figs/2019-06-17-figs/reference-phya-barplot.png", bar_ref, width=7, height=4, units=c("cm"))
