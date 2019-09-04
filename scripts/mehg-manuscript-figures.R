library(tidyverse)
library(reshape2)
library(plyr)

metadata = read.csv("files/methylator-updated-hgcA-metadata.csv")
environment = count(metadata, "code")
group = count(metadata, "Group")
study = count(metadata, "study")

# Environments plot
phaColfunc = colorRampPalette(c("palegreen3", "slateblue"))
bar_env = environment %>% ggplot(aes(x=reorder(code, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(13))) + coord_flip() + scale_y_continuous(limits=c(0,215), expand = c(0, 0)) + theme_classic()
bar_env 

# phyla plot
phylaColfunc = colorRampPalette(c("khaki1", "orchid4"))
bar_phyla = group %>% ggplot(aes(x=reorder(Group, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(18))) + coord_flip() + scale_y_continuous(limits=c(0,400), expand= c(0,0)) + theme_classic()
bar_phyla

# reference isolate phyla 
references = metadata %>% filter(code=="Isolate")
reference_phyla = count(references, "Phyla")
bar_ref = reference_phyla %>% ggplot(aes(x=reorder(Phyla, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(8))) + coord_flip() + scale_y_continuous(limits=c(0,100), expand= c(0,0), breaks=seq(0,100,10)) + theme_classic()
bar_ref

# CP phyla
cp = metadata %>% filter(Group=="Candidate Phyla")
cp_phyla = count(cp, "Phyla")
bar_cp= cp_phyla %>% ggplot(aes(x=reorder(Phyla, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(21))) + coord_flip() + scale_y_continuous(limits=c(0,6), expand= c(0,0)) + theme_classic()
bar_cp 

# genome statistics
# all
# genome size
genome_size = metadata %>% select(Group, size)
genome_size_avgs = aggregate(genome_size[2], list(genome_size$Group), mean)
all_avgs = genome_size %>% ggplot(aes(x=Group, y=size)) + geom_boxplot(fill="gray")
all_avgs_theme = all_avgs + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1))
all_avgs_theme
# gc content
gc_all = metadata %>% select(Group, gc)
gc_meth = gc_all %>% ggplot(aes(x=Group, y=gc)) + geom_boxplot(fill="gray")
gc_meth_theme = gc_meth + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1))
gc_meth_theme
# just isolates
# genome size
isolate_sizes = metadata %>% filter(code=="Isolate") %>% select(Group, size)
isolate_avgs = isolate_sizes %>% ggplot(aes(x=Group, y=size)) + geom_boxplot(fill="gray")
isolate_avgs_theme = isolate_avgs + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1))
isolate_avgs_theme
# gc content 
gc_isolates = metadata %>% filter(code=="Isolate") %>% select(Group, gc)
gc_meth_isolates = gc_isolates %>% ggplot(aes(x=Group, y=gc)) + geom_boxplot(fill="gray")
gc_isolates_theme = gc_meth_isolates + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1))
gc_isolates_theme
# Save plots
# bar plots
ggsave(file="figs/2019-09-02-environment-barplot.png", bar_env, width=20, height=10, units=c("cm"))
ggsave(file="figs/2019-09-02-phyla-barplot.png", bar_phyla, width=20, height=10, units=c("cm"))
ggsave(file="figs/2019-09-02-reference-phyla-barplot.png", bar_ref, width=7, height=4, units=c("cm"))
ggsave(file="figs/2019-09-02-cp-phyla-barplot.png", bar_cp, width=7, height=10, units=c("cm"))
# stats
ggsave(file="figs/2019-09-02-all-genome-sizes.png", all_avgs_theme, width=15, height=10, units=c("cm"))
ggsave(file="figs/2019-09-02-all-gc.png", gc_meth_theme, width=15, height=10, units=c("cm"))
ggsave(file="figs/2019-09-02-isolates-sizes.png", isolate_avgs_theme, width=10, height=7, units=c("cm"))
ggsave(file="figs/2019-09-02-isolates-gc.png", gc_isolates_theme, width=10, height=7, units=c("cm"))