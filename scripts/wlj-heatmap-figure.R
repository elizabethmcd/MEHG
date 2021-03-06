library(tidyverse)
library(reshape2)
library(viridis)

# Analysis for WLJ markers

# Load marker counts for WLJ and all markers and metadata 
wlj = read_csv("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/results/metabolic-results/2019-06-20-wlj-results.csv")
metadata = read_csv("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/methylator-metadata.csv")
all_markers = read_csv("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/results/metabolic-results/2019-06-07-metabolic-results.csv")
colnames(wlj)[1] = "genome"
colnames(all_markers)[1] = "genome"
colnames(metadata)[1] = "genome"

# Merge into one for wlj
phyla = metadata %>% select(genome, Phylum)
colnames(phyla) = c("genome", "phyla")
marker_table = inner_join(phyla, wlj)

# Merge for all markers to get hgcA column
full_table = inner_join(phyla, all_markers)

####### WLJ statistics
# Investigate percentages for each marker by presence/absence, ignore copy #s, WLJ
presence_absence = as.data.frame(lapply(marker_table[3:24], function(x) ifelse(x>1, 1, x)))
sort(colMeans(presence_absence))
# aggregate by phyla for percentages
sumphyla_pres_ab = presence_absence
sumphyla_pres_ab$phyla = marker_table$phyla
sumphyla = aggregate(sumphyla_pres_ab[1:22], list(sumphyla_pres_ab$phyla), sum)

######## all markers statistics
all_pres_ab = as.data.frame(lapply(full_table[3:132], function(x) ifelse(x>1, 1, x)))
sumphyla_all = all_pres_ab
sumphyla_all$phyla = full_table$phyla
sumAll = aggregate(sumphyla_all[1:130], list(sumphyla_all$phyla), sum)

# percentage by total genomes in phyla
sumphyla$hgcA = sumAll$hgcA
pcts_phyla = as.data.frame(lapply(sumphyla[,-1], function(x) {
  (x / sumphyla$hgcA) * 100
}))
pcts_phyla$hgcA_count = sumphyla$hgcA
pcts_phyla_table = as.data.frame(lapply(pcts_phyla, round, 2))
pcts_phyla_table$phyla = sumphyla$Group.1

# percentage of marker by total genomes
genome_sums = as.data.frame(colSums(presence_absence))
genome_sums$average = lapply(genome_sums$`colSums(presence_absence)`, function(x) x/518 * 100)
genome_sums$average = lapply(genome_sums$average, round, 2)

# percentage of marker at phyla level
pres_ab_phyla = as.data.frame(lapply(pcts_phyla_table[1:22], function(x) ifelse(x>0, 1, x)))
pres_ab_phyla["percent_phyla", ] = (colSums(pres_ab_phyla[1:22]) / 23) * 100
pres_ab_phyla["percent_genomes", ] = genome_sums$average

# put together in table
pcts_phyla_table["percent_phyla", ] = pres_ab_phyla["percent_phyla", ]
pcts_phyla_table["percent_genomes", ] = pres_ab_phyla["percent_genomes", ]

# export raw table to clean up
write_csv(pcts_phyla_table, "~/Desktop/mehg-wlj-percentages-raw.csv")

# read in cleaned table
master_table = read_csv("~/Desktop/mehg-wlj-percentages-cleaned.csv")
master_table_v1 = master_table %>% column_to_rownames("phyla")
master_table_v2 = master_table_v1[,order(-master_table_v1[nrow(master_table_v1),])]
master_table_ordered = master_table_v2
master_table_ordered$phyla = master_table$phyla

write_csv(master_table_ordered, "~/Desktop/mehg-wlj-stats-ordered.csv")

# order of pathways
wlj_ordered = master_table_ordered[,c("phyla", "fdhA", "fdhB", "fhs", "folD", "metF", "fwdA", "fwdB", "fwdC", "fwdD", "fwdF", "ftr", "mch", "mtd", "mer", "acsA", "acsB", "cdhA", "cdhC", "cdhD", "cdhE", "cdhB", "acsE")]

# melt for heatmap
table_melted = melt(wlj_ordered, id.vars="phyla")
table_melted$phyla = factor(wlj_ordered$phyla, level = c('Acidobacteria', 'Actinobacteria', 'Aminicenantes', 'Bacteroidetes', 'Chlorobi', 'Chloroflexi', 'Deltaproteobacteria', 'Eisenbacteria', 'Elusimicrobia', 'Euryarchaeota', 'FCPU426', 'Fibrobacteres', 'Firestonebacteria', 'Firmicutes', 'KSB1', 'Margulisbacteria', 'Nitrospirae', 'PVC', 'Raymondbacteria', 'Spirochaetes', 'Synergistaceae', 'Taylorbacteria', 'Unclassified', 'Woesearchaeota', 'WOR1', 'Phyla', 'Genomes'))
marker_plot = ggplot(table_melted, aes(x=variable, y=fct_rev(phyla), fill=value)) + geom_tile(color="black") + scale_fill_viridis(option="viridis", alpha=1, begin=0, end=1, direction=-1) + theme_bw()
marker_plot2 = marker_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
marker_plot2

marker_plot_no_legends = marker_plot2 +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position="none")
marker_plot_no_legends

ggsave(marker_plot2, file="~/Desktop/wlj-heatmap.png", height=20, width=40, units=c("cm"))

ggsave(marker_plot_no_legends, file="~/Desktop/wlj-no-legends-heatmap.png", height=20, width=40, units=c("cm"))

# weird actinos
actinos = c("bacteria00190",
            "bacteria12801",
            "bacteria12873",
            "bacteria13742",
            "bacteria17823",
            "bacteria20068",
            "bacteria20100",
            "bacteria20693",
            "bacteria20701",
            "bacteria20733",
            "bacteria20819",
            "bacteria23279",
            "bacteria30001",
            "bacteria30002")
actino_wlj = wlj %>% filter(genome %in% actinos)
thermo = wlj %>% filter(genome == c("bacteria30002","bacteria20068"))
