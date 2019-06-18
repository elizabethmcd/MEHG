library(tidyverse)
library(reshape2)
library(viridis)

# Load marker counts and metadata 
metabolic_markers = read_csv("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/results/metabolic-results/2019-06-07-metabolic-results.csv")
metadata = read_csv("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/methylator-metadata.csv")
colnames(metabolic_markers)[1] = "genome"
colnames(metadata)[1] = "genome"

# Merge into one
phyla = metadata %>% select(genome, Phylum)
colnames(phyla) = c("genome", "phyla")
marker_table = inner_join(phyla, metabolic_markers)

# Counts with genome IDs
count_data = metabolic_markers %>% column_to_rownames("genome")

# Investigate percentages for each marker by presence/absence, ignore copy #s
presence_absence = as.data.frame(lapply(count_data, function(x) ifelse(x>1, 1, x)))
sort(colMeans(presence_absence))

# Get rid of columns that don't really add to the analysis
greater20 = presence_absence[,colMeans(presence_absence) > 0.2]

    # aggregate by phyla for percentages
sumphyla_pres_ab = greater20
sumphyla_pres_ab$phyla = marker_table$phyla
sumphyla = aggregate(sumphyla_pres_ab[1:32], list(sumphyla_pres_ab$phyla), sum)

# get rid of ones with all zeros
## already done when set filter for markers have to be set above a total percentage
# presence_absence_no_all_zeros = as.data.frame(presence_absence[,-(which(colSums(presence_absence) == 0))])
# presence_absence_no_all_zeros$phyla = marker_table$phyla
# no_zeros = aggregate(presence_absence_no_all_zeros[1:106], list(presence_absence_no_all_zeros$phyla), sum)

# percentage by total genomes in phyla
pcts_phyla = as.data.frame(lapply(sumphyla[,-1], function(x) {
  (x / sumphyla$hgcA) * 100
}))
pcts_phyla$hgcA_count = sumphyla$hgcA
pcts_phyla_table = as.data.frame(lapply(pcts_phyla, round, 2))
pcts_phyla_table$phyla = sumphyla$Group.1

# percentage of marker by total genomes
genome_sums = as.data.frame(colSums(greater20))
genome_sums$average = lapply(genome_sums$`colSums(greater20)`, function(x) x/518 * 100)
genome_sums$average = lapply(genome_sums$average, round, 2)

# percentage of marker at phyla level
pres_ab_phyla = as.data.frame(lapply(pcts_phyla_table[1:32], function(x) ifelse(x>1, 1, x)))
pres_ab_phyla["percent_phyla", ] = (colSums(pres_ab_phyla[1:32]) / 23) * 100
pres_ab_phyla["percent_genomes", ] = genome_sums$average

# put together in table
pcts_phyla_table["percent_phyla", ] = pres_ab_phyla["percent_phyla", ]
pcts_phyla_table["percent_genomes", ] = pres_ab_phyla["percent_genomes", ]

# export raw table to clean up
write_csv(pcts_phyla_table, "~/Desktop/mehg-metabolism-percentages-raw.csv")

# reimport to sort the markers by highest > lowest percentage
master_table = read_csv("~/Desktop/mehg-metabolism-percentages-cleaned.csv")
master_table_v1 = master_table %>% column_to_rownames("phyla")
master_table_v2 = master_table_v1[,order(-master_table_v1[nrow(master_table_v1),])]
master_table_ordered = master_table_v2
master_table_ordered$phyla = master_table$phyla

write_csv(master_table_ordered, "~/Desktop/mehg-metabolism-stats-ordered.csv")

# heatmap of ordered marker percentages
no_totals = master_table_ordered[,c(-1)]
throwout_no_cutoffs = no_totals[,c(-2, -3, -4)]
table_melted = melt(throwout_no_cutoffs, id.vars="phyla")
table_melted$phyla = factor(throwout_no_cutoffs$phyla, level = c('Acidobacteria', 'Actinobacteria', 'Aminicenantes', 'Bacteroidetes', 'Chlorobi', 'Chloroflexi', 'Deltaproteobacteria', 'Eisenbacteria', 'Elusimicrobia', 'Euryarchaeota', 'FCPU426', 'Fibrobacteres', 'Firestonebacteria', 'Firmicutes', 'KSB1', 'Margulisbacteria', 'Nitrospirae', 'PVC', 'Raymondbacteria', 'Spirochaetes', 'Synergistaceae', 'Taylorbacteria', 'Unclassified', 'Woesearchaeota', 'WOR1', 'Phyla', 'Genomes'))
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
