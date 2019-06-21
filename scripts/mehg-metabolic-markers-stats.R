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

# Investigate percentages for each marker by presence/absence, ignore copy #s
presence_absence = as.data.frame(lapply(marker_table[3:132], function(x) ifelse(x>1, 1, x)))
sort(colMeans(presence_absence))

# Get rid of columns with WLJ from full metabolic analysis, including curated WLJ proteins for more closer look
# also remove markers that don't have a cutoff and throw things off
no_wlj = presence_absence %>% select(c(-codhC_TIGR00316, -codhD_TIGR00381, -codh_catalytic_TIGR01702, -fdhA_TIGR01591, -fdh_thiol_id_TIGR02819, -fdhB_TIGR01582, -fdhC_TIGR01583, -carbon_monoxide_dehydrogenase_coxL_TIGR02416, -carbon_monoxide_dehydrogenase_coxM, -carbon_monoxide_dehydrogenase_coxS))
# Get rid of columns that don't really add to the analysis
greater20 = no_wlj[,colMeans(no_wlj) > 0.2]

    # aggregate by phyla for percentages
sumphyla_pres_ab = greater20
sumphyla_pres_ab$phyla = marker_table$phyla
sumphyla = aggregate(sumphyla_pres_ab[1:24], list(sumphyla_pres_ab$phyla), sum)

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
pres_ab_phyla = as.data.frame(lapply(pcts_phyla_table[1:24], function(x) ifelse(x>0, 1, x)))
pres_ab_phyla["percent_phyla", ] = (colSums(pres_ab_phyla[1:24]) / 23) * 100
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

# heatmap of marker percentages 
throwout_no_cutoffs = master_table_ordered[,c(-2, -3, -4)]

# split by function
# hgcA > cytochromes > hydrogenases > nitrogen > sulfur > followed by wlj results parsed in a separate script and heatmaps stitched together

ordered = throwout_no_cutoffs[,c("phyla", "hgcA", "cydA_PF01654", "cydB_TIGR00203", "coxB_TIGR02866", "nrfA_PF02335", "cyoE_TIGR01473","Hydrogenase_Group_1", "Hydrogenase_Group_3b", "Hydrogenase_Group_3c", "Hydrogenase_Group_3d", "Hydrogenase_Group_4", "nifA_Mo_TIGR01282", "nifB_Mo_TIGR01286", "nifH_TIGR01287", "sulfur_dioxygenase_sdo", "ars_thioredoxin_TIGR02691", "cysC_TIGR00455", "cysN_TIGR02034")]

# melt for heatmap
table_melted = melt(ordered, id.vars="phyla")
table_melted$phyla = factor(ordered$phyla, level = c('Acidobacteria', 'Actinobacteria', 'Aminicenantes', 'Bacteroidetes', 'Chlorobi', 'Chloroflexi', 'Deltaproteobacteria', 'Eisenbacteria', 'Elusimicrobia', 'Euryarchaeota', 'FCPU426', 'Fibrobacteres', 'Firestonebacteria', 'Firmicutes', 'KSB1', 'Margulisbacteria', 'Nitrospirae', 'PVC', 'Raymondbacteria', 'Spirochaetes', 'Synergistaceae', 'Taylorbacteria', 'Unclassified', 'Woesearchaeota', 'WOR1', 'Phyla', 'Genomes'))
marker_plot = ggplot(table_melted, aes(x=variable, y=fct_rev(phyla), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="plasma", alpha=1, begin=0, end=1, direction=-1) + theme_bw()
marker_plot2 = marker_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
marker_plot2

metabolism_no_legends = marker_plot2 +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position="none")

onecolor =ggplot(table_melted, aes(x=variable, y=fct_rev(phyla), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="grey", high="steelblue") + theme_bw()
onecolor2 = onecolor +  theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
onecolor2

onecolorNolegends = onecolor2 +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position="none")
onecolorNolegends

ggsave(marker_plot2, file="~/Desktop/metabolic-markers-heatmap.png", height=20, width=40, units=c("cm"))

ggsave(metabolism_no_legends, file="~/Desktop/metabolic-markers-no-legends.png", height=20, width=40, units=c("cm"))

ggsave(onecolorNolegends, file="~/Desktop/metabolic-markers-heatmap-one-color.png", height=20, width=40, units=c("cm"))

ggsave(onecolor2, file="~/Desktop/metabolic-markers-one-color-legends.png", height=20, width=40, units=c("cm"))
