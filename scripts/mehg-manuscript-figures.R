library(tidyverse)
library(reshape2)
library(plyr)

metadata = read.csv("methylator-metadata.csv")
environment = count(metadata, "Code")
phyla = count(metadata, "Phylum")

phaColfunc = colorRampPalette(c("palegreen3", "slateblue"))
bar_env = environment %>% ggplot(aes(x=reorder(Code, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(11))) + coord_flip() + scale_y_continuous(limits=c(0,125), expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
ggsave(file="environment-barplot.png", bar_env, width=30, height=20, units=c("cm"))

bar_phyla = phyla %>% ggplot(aes(x=reorder(Phylum, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(28))) + coord_flip() + scale_y_continuous(limits=c(0,200), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
bar_phyla
