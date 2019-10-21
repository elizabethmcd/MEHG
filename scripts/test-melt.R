library(tidyverse)
library(reshape2)

groups = readLines("/Users/emcdaniel/Desktop/McMahon-Lab/Pipelines/metabolisHMM/bin/groups.txt")
markers = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/Pipelines/metabolisHMM/bin/TEST1/results/custom-markers-results.csv", header=TRUE)
metadata = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/Pipelines/metabolisHMM/data/bac-genomes.csv", header=FALSE)
colnames(metadata) = c("genome", "group")
colnames(markers)[1] = c("genome")

merged_table = left_join(metadata, markers)
presence_absence = as.data.frame(lapply(merged_table[3:ncol(merged_table)], function(x) ifelse(x>1, 1, x)))
presence_absence$group <- metadata$group
exclude <- ncol(presence_absence) - 1
agg <- aggregate(presence_absence[1:exclude], list(presence_absence$group), mean)
colnames(agg)[1] <- "group"
table_melted <- melt(agg, id.vars="group")
table_melted$group <- factor(table_melted$group, levels=as.factor(groups))

table_melted <- table_melted %>% mutate(group = factor(group),group = factor(group, levels = rev(as.factor(groups))))

plot <- table_melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black') + scale_fill_viridis_c(alpha=1,begin=0,end=1,direction=-1) + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))      
plot_formatted <- plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=85, hjust=0), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + labs(x=NULL,y=NULL) + guides(fill = guide_colorbar(nbin = 10)) + scale_y_discrete(expand=c(0,0), order(rev(as.factor(groups))))
ggsave(file=figureOut, plot_formatted, height=15, width=55, units=c("cm"))
plot
plot_formatted
