library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)

markers = read.csv("~/Desktop/custom-markers-results.csv")
metadata = read.csv("files/stats/quals/mehg-markers-metadata.csv", header=FALSE, stringsAsFactors = FALSE)
colnames(markers)[1] = c("genome")
colnames(metadata) = c("genome", "group")
merged = left_join(metadata, markers)
presence_absence = data.frame(lapply(merged[3:ncol(merged)], function(x) ifelse(x>1, 1, x)), stringsAsFactors = FALSE)
presence_absence$group = metadata$group
pa = data.frame(presence_absence, stringsAsFactors = FALSE)
group_means = colMeans(presence_absence[1:27])
test = data.frame(rbind(pa, group_means), stringsAsFactors = FALSE)
test[525,28] = "total"
exclude <- ncol(presence_absence) - 1
agg <- aggregate(test[1:exclude], list(test$group), mean)

# individual melts for heatmap colors
  # carbon & wlj pathways
carbon = cbind(agg$Group.1, agg[,c("codhC_TIGR00316", "codhD_TIGR00381", "codh_catalytic_TIGR01702", "folD", "metF", "coxA_TIGR02891", "coxB_TIGR02866")])
colnames(carbon) = c("group", "codhC", "codhD", "codh catalytic", "folD", "metF", "coxA", "coxB")
carbon.melted = melt(carbon, id.vars= "group") %>%  mutate(group=factor(group), group = factor(group, levels = rev(levels(group))))
carbon.plot = carbon.melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() +  scale_fill_gradient(low="gray85", high="brown4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
carbon.formatted = carbon.plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0))
carbon.clean = carbon.formatted + theme(legend.position="none")
  # sulfur
sulfur = cbind(agg$Group.1, agg[,c("dsrA_TIGR02064", "dsrB_TIGR02066", "dsrD_PF08679", "sat_TIGR00339")])
colnames(sulfur) = c("group", "dsrA", "dsrB", "dsrD", "sat")
sulfur.melted = melt(sulfur, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group))))
sulfur.plot = sulfur.melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() +  scale_fill_gradient(low="grey85", high="orchid4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
sulfur.formatted = sulfur.plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Sulfur", y=NULL) +  scale_y_discrete(expand=c(0,0))
sulfur.clean = sulfur.formatted + theme(legend.position="none")
  # nitrogen
nitrogen = cbind(agg$Group.1, agg[,c("nifA_Mo_TIGR01282", "nifB_Mo_TIGR01286", "nifH_TIGR01287", "nirB_TIGR02374", "nirD_TIGR02378")])
colnames(nitrogen) = c("group", "nifA", "nifB", "nifH", "nirB", "nirD")
nitrogen.melted <- melt(nitrogen, id.vars="group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group))))
nitrogen.plot <- nitrogen.melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() +  scale_fill_gradient(low="grey85", high="steelblue4") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
nitrogen.formatted = nitrogen.plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + labs(x="Nitrogen", y=NULL) +  scale_y_discrete(expand=c(0,0))
nitrogen.clean = nitrogen.formatted + theme(legend.position="none")
  # metal resistance
metals = cbind(agg$Group.1, agg[,c("hgcA", "arsC_thioredoxin", "sulfur_dioxygenase_sdo")])
colnames(metals) = c("group", "hgcA", "arsC", "sdo")  
metals.melted = melt(metals, id.vars = "group") %>% mutate(group=factor(group), group = factor(group, levels = rev(levels(group))))
metals.plot = metals.melted %>% ggplot(aes(x=variable, y=(group), fill=value)) + geom_tile(color='black',size=0.5,aes(width=1, height=1)) + coord_fixed() +  scale_fill_gradient(low="grey85", high="midnightblue") + theme(panel.grid = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
metals.formatted = metals.plot + scale_x_discrete(position="top", expand=c(0,0)) + theme(axis.text.x.top=element_text(angle=85, hjust=0, face="italic"), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + labs(x=NULL, y=NULL) +  scale_y_discrete(expand=c(0,0))
metals.clean = metals.formatted + theme(legend.position="none")

# stitch
carbon.tmp = ggplot_build(carbon.clean)
nitrogen.tmp = ggplot_build(nitrogen.clean)
sulfur.tmp = ggplot_build(sulfur.clean)
metals.tmp = ggplot_build(metals.clean)
g1 = ggplot_gtable(metals.tmp) ; g2 = ggplot_gtable(carbon.tmp) ; g3 = ggplot_gtable(nitrogen.tmp) ; g4 = ggplot_gtable(sulfur.tmp)
n1 <- length(metals.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n2 <- length(carbon.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n3 <- length(nitrogen.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])
n4 <- length(sulfur.tmp[["panel"]][["ranges"]][[1]][["x.major_source"]])


g = cbind(g1, g2, g3, g4, size="min")
ggsave(file="~/Desktop/test.png", plot=g, width=20, height=10, units=c("cm"))
