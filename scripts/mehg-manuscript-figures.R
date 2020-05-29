library(tidyverse)
library(reshape2)
library(plyr)

metadata = read.csv("files_and_results//stats/mehg-final-dataset-metadata.csv")
environment = count(metadata, "code")
group = count(metadata, "Group")
phyla = count(metadata, "Phyla")
study = count(metadata, "study")
author = count(metadata, "reference")

# Environments plot
phaColfunc = colorRampPalette(c("palegreen3", "slateblue"))
bar_env = environment %>% ggplot(aes(x=reorder(code, freq), y=freq)) + geom_bar(stat="identity", fill=(phaColfunc(12))) + coord_flip() + scale_y_continuous(limits=c(0,215), expand = c(0, 0)) + theme_classic()
bar_env 

# phyla plot
phylaColfunc = colorRampPalette(c("khaki1", "orchid4"))
bar_phyla = group %>% ggplot(aes(x=reorder(Group, freq), y=freq)) + geom_bar(stat="identity", fill=(phylaColfunc(14))) + coord_flip() + scale_y_continuous(limits=c(0,400), expand= c(0,0)) + theme_classic()
bar_phyla

# reference isolate phyla 
references = metadata %>% filter(code=="Isolate")
reference_phyla = count(references, "Phyla")
bar_ref = reference_phyla %>% ggplot(aes(x=reorder(Phyla, freq), y=freq)) + geom_bar(stat="identity", fill="azure4") + coord_flip() + scale_y_continuous(limits=c(0,100), expand= c(0,0), breaks=seq(0,100,10)) + theme_classic()
bar_ref

# other with low amounts
other = metadata %>% filter(Group=="Other")
other_phyla = count(other, "Phyla")
bar_other= other_phyla %>% ggplot(aes(x=reorder(Phyla, freq), y=freq)) + geom_bar(stat="identity", fill="azure4") + coord_flip() + scale_y_continuous(limits=c(0,10), expand= c(0,0), breaks=seq(0,10,1)) + theme_classic()
bar_other

# genome statistics
# all
# genome size
genome_size = metadata %>% select(Group, size) %>% filter(Group != "Other")
genome_size_avgs = aggregate(genome_size[2], list(genome_size$Group), mean)
all_avgs = genome_size %>% ggplot(aes(x=Group, y=size)) + geom_boxplot(fill="gray")
all_avgs_theme = all_avgs + scale_y_continuous(breaks=seq(0,15,3)) + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1))
all_avgs_theme
# gc content
gc_all = metadata %>% select(Group, gc) %>% filter(Group != "Other")
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

# stats for "other" 
#size
other_sizes = metadata %>% filter(Group=="Other") %>% select(Phyla, size)
other_avgs = other_sizes %>% ggplot(aes(x=Phyla, y=size)) + geom_boxplot(fill="gray")
other_avgs_theme = other_avgs + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1)) 
other_avgs_theme
#gc
other_gc = metadata %>% filter(Group=="Other") %>% select(Phyla, gc)
other_gc_avgs = other_gc %>% ggplot(aes(x=Phyla,y=gc)) + geom_boxplot(fill="gray")
other_gc_theme = other_gc_avgs + theme_classic() + theme(axis.text.x = element_text(angle = 85, hjust = 1)) 
other_gc_theme

# quality stats of MAGs
mags = metadata %>% filter(code != "Isolate")
mag_qual = mags %>% ggplot(aes(x=Completeness, y=Contamination, color=Group)) + geom_point(size=2) + scale_x_continuous(breaks=seq(0,100,10)) + theme_classic()
mag_qual

# Save plots
# bar plots
ggsave(file="figs/2019-10-21-environment-barplot.png", bar_env, width=20, height=10, units=c("cm"))
ggsave(file="figs/2019-10-21-phyla-barplot.png", bar_phyla, width=20, height=10, units=c("cm"))
ggsave(file="figs/2019-10-21-reference-phyla-barplot.png", bar_ref, width=7, height=4, units=c("cm"))
ggsave(file="figs/2019-10-21-other-phyla-barplot.png", bar_other, width=10, height=10, units=c("cm"))

# stats
ggsave(file="figs/2019-10-21-all-genome-sizes.png", all_avgs_theme, width=15, height=10, units=c("cm"))
ggsave(file="figs/2019-10-21-all-gc.png", gc_meth_theme, width=15, height=10, units=c("cm"))
ggsave(file="figs/2019-10-21-isolates-sizes.png", isolate_avgs_theme, width=10, height=7, units=c("cm"))
ggsave(file="figs/2019-10-21-isolates-gc.png", gc_isolates_theme, width=10, height=7, units=c("cm"))
ggsave(file="figs/2019-10-21-other-phyla-sizes.png", other_avgs_theme, width=7,height=10,units=c("cm"))
ggsave(file="figs/2019-10-21-other-phyla-gc.png", other_gc_theme,width=7,height=10,units=c("cm"))
ggsave(file="figs/2019-10-21-mags-qual-stats.png", mag_qual, width=15,height=10,units=c("cm"))
