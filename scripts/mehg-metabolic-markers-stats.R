library(tidyverse)
library(reshape2)
library(viridis)

# Load marker counts and metadata 
metabolic_markers = read_csv("results/mebolic-results/2019-03-11-metabolic-marker-results.csv")
metadata = read_csv("~/Desktop/methylator-metadata.csv")
colnames(metabolic_markers)[1] = "genome"
colnames(metadata)[1] = "genome"
wlj = read_csv("results/mebolic-results/2019-03-11-WLJ-markers-fixed.csv")

# Merge into one
phyla = metadata %>% select(genome, Phylum)
colnames(phyla) = c("genome", "phyla")
marker_table = inner_join(phyla, metabolic_markers)

# Counts with genome IDs
count_data = metabolic_markers %>% column_to_rownames("genome")

# Investigate percentages for each marker by presence/absence, ignore copy #s
presence_absence = as.data.frame(lapply(count_data, function(x) ifelse(x>1, 1, x)))
sort(colMeans(presence_absence))
    # aggregate by phyla for percentages
sumphyla_pres_ab = presence_absence
sumphyla_pres_ab$phyla = marker_table$phyla
sumphyla = aggregate(sumphyla_pres_ab[1:130], list(sumphyla_pres_ab$phyla), sum)

# get rid of ones with all zeros
presence_absence_no_all_zeros = as.data.frame(presence_absence[,-(which(colSums(presence_absence) == 0))])
presence_absence_no_all_zeros$phyla = marker_table$phyla
no_zeros = aggregate(presence_absence_no_all_zeros[1:106], list(presence_absence_no_all_zeros$phyla), sum)

# percentage by total genomes in phyla
pcts_phyla = as.data.frame(lapply(no_zeros[,-1], function(x) {
  (x / no_zeros$hgcA) * 100
}))
pcts_phyla$hgcA_count = no_zeros$hgcA
pcts_phyla_table = as.data.frame(lapply(pcts_phyla, round, 2))
pcts_phyla_table$phyla = no_zeros$Group.1

# percentage of marker by phyla and genomes
phyla_percentages = pcts_phyla[,-c(107,108)]
avg = as.data.frame(lapply(no_zeros[2:107], function(x) x/512 * 100))
avg_sigf = as.data.frame(lapply(avg, round, 2))
avg_sigf$phyla = no_zeros$Group.1
pres_ab_phyla = as.data.frame(lapply(no_zeros[2:107], function(x) ifelse(x>1, 1, x)))
pres_ab_phyla$phyla = no_zeros$Group.1
phyla_percentages["percent_phyla", ] = (colSums(pres_ab_phyla[1:106]) / 28) * 100
phyla_percentages["percent_genomes", ] = (colSums(no_zeros[2:107]) / 512 ) * 100

# raw of both to stitch together
write_csv(phyla_percentages, "~/Desktop/mehg-metabolism-percentages-phyla-genomes.csv")
write_csv(pcts_phyla_table, "~/Desktop/mehg-metabolism-percentages-each-phyla.csv")

# reimport to sort the markers by highest > lowest percentage
master_table = read_csv("~/Desktop/mehg-metabolism-stats-raw.csv")
master_table_v1 = master_table %>% column_to_rownames("phyla")
master_table_v2 = master_table_v1[,order(-master_table_v1[nrow(master_table_v1),])]
master_table_ordered = master_table_v2
master_table_ordered$phyla = master_table$phyla

write_csv(master_table_ordered, "~/Desktop/mehg-metabolism-stats-ordered.csv")

# heatmap of ordered marker percentages
table_melted = melt(master_table_ordered, id.vars="phyla")
table_melted$phyla = factor(table_melted$phyla, level = c('Acidobacteria', 'Actinobacteria', 'Aminicenantes', 'Bacteroidetes', 'Chlorobi', 'Chloroflexi', 'Deltaproteobacteria', 'Eisenbacteria', 'Elusimicrobia', 'Euryarchaeota', 'FCPU426', 'Fibrobacteres', 'Firestonebacteria', 'Firmicutes', 'KSB1', 'Lentisphaerae', 'Margulisbacteria', 'Nitrospirae', 'Planctomycetes', 'PVC', 'Raymondbacteria', 'Spirochaetes', 'Synergistaceae', 'Taylorbacteria', 'Unclassified', 'Verrucomicrobia', 'Woesearchaeota', 'WOR1', 'Phyla', 'Genomes'))
marker_plot = ggplot(table_melted, aes(x=variable, y=fct_rev(phyla), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="plasma", alpha=1, begin=0, end=1, direction=-1) + theme_bw()
marker_plot2 = marker_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
marker_plot2

ggsave(marker_plot2, file="~/Desktop/metabolic-markers-heatmap.png", height=20, width=40, units=c("cm"))

# WLJ pathway characterization
wlj_counts = wlj %>% column_to_rownames("genome")
wlj_presence_absence = as.data.frame(lapply(wlj_counts, function(x) ifelse(x>1, 1, x)))
as.data.frame(sort(colMeans(wlj_presence_absence)))
colMeans(wlj_presence_absence)

# folD seems to be important, which ones are missing it 
genomes_no_folD = wlj %>%  filter(folD==0) %>% select(genome)
no_folD_taxonomy = inner_join(genomes_no_folD, metadata)
