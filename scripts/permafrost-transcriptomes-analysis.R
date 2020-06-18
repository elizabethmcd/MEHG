library(tximport)
library(readr)
library(tibble)
library(dplyr)
library(reshape)
library(gtable)
library(grid)
library(gridExtra)
library(viridis)
library(tidyverse)

# kallisto directory for mapped transcriptomes to methylators
meth_dir <- "kallisto_mapping/"
meth_samples <- read.table(file.path(meth_dir, "permafrost-metadata.txt"), header=TRUE)
meth_files <- file.path(meth_dir, meth_samples$sample, "abundance.h5")
rownames(meth_samples) <- meth_samples$sample
names(meth_files) <- meth_samples$sample
# check if in the cloud backup
meth.kallisto <- tximport(meth_files, type="kallisto", txOut = TRUE)

# counts file
meth.counts <- as.data.frame(meth.kallisto)
final.meth.counts <- rownames_to_column(meth.counts, var="genome_name")
meth.counttable <- final.meth.counts[, c(1:27)]

# ids to pick certain genomes
ids <- read.delim("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/kallisto_mapping/perma-meth-ids.txt", sep="\t", header=TRUE)
counttable.ids <- left_join(meth.counttable, ids)

## Merge with annotations to get gene lengths
meth_prokka = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files_and_results/permafrostTranscription/permafrost-methylator-prokka-annotations.txt", sep="\t", header=FALSE)
colnames(meth_prokka) <- c("genome_name", "prokka_annotation", "size_bp", "accession")

# permafrost results from archived set of metadata
# merge with metadata to get genome sizes to get relative expression by phyla
meth_metadata = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files_and_results/archived/archived-metadata/archived-methylator-metadata.csv")
peat_phyla = meth_metadata %>% select("genome_name", "Phylum")
colnames(peat_phyla) <- c("genome", "Phylum")
meth_counts_totals_table = left_join(counttable.ids, peat_phyla)
meth_expression_average = aggregate(meth_counts_totals_table[2:27], list(meth_counts_totals_table$Phylum), sum)
genome_expression_average = aggregate(meth_counts_totals_table[2:27], list(meth_counts_totals_table$genome), sum)
meth_expression_average$avg_expression = rowSums(meth_expression_average[2:27]) / 26
genome_expression_average$avg_expression = rowSums(genome_expression_average[2:27]) / 26

# hgcA locus tags 
hgcA = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files_and_results/archived/tree-tax-files/hgcA-locus-tags-phyla.txt", sep="\t", header=FALSE)
colnames(hgcA) = c("genome_name", "locus_tag", "Phylum")
hgcA = hgcA %>% select(genome_name, locus_tag)
mehg_peat = meth_metadata %>% filter(Study=="Woodcroft2018") 
peat_locus_tags = left_join(mehg_peat, hgcA)
peat_hgcA = peat_locus_tags %>% select(c(Phylum, locus_tag))

# counts of hgcA
colnames(counttable.ids)[1] = c("locus_tag")
hgcA_counts = left_join(peat_hgcA, counttable.ids)
write.csv(hgcA_counts, "/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files_and_results/permafrostTranscription/meth-counts-hgcA-norm.csv", quote=FALSE, row.names=FALSE)
hgcA_non_agg_zeros = hgcA_counts %>% select(c(Phylum, locus_tag, 4,6,9,10,12,14,16,17,18,20:26,28))
hgcA_phyla_total = aggregate(hgcA_counts[3:28], list(hgcA_counts$Phylum),sum)
hgcA_no_zeros = hgcA_phyla_total %>% select(c(Group.1, 3,5,8,9,11,13,15,16,17,19:25,27))
hgcA_counts.m = melt(hgcA_phyla_total, id.vars="Group.1")
hgcA_no_zeros.m = melt(hgcA_no_zeros, id.vars="Group.1")
sample_list = c("abundance.20120600_S2M",
                "abundance.20120800_S1M",
                "abundance.20120600_S2D",
                "abundance.20120700_S1D",
                "abundance.20120800_S1X",
                "abundance.20100900_E2S",
                "abundance.20110600_E1M",
                "abundance.20120600_E1M",
                "abundance.20120700_E3M",
                "abundance.20120800_E2M",
                "abundance.20100900_E1D",
                "abundance.20110600_E1D",
                "abundance.20110700_E1D",
                "abundance.20120600_E1D",
                "abundance.20120700_E3D",
                "abundance.20120800_E2D",
                "abundance.20120800_E3D")
hgcA_no_zeros.m$variable = factor(hgcA_no_zeros.m$variable, levels = c(sample_list))

# total hgcA by sample
sample_totals = as.data.frame(colSums(hgcA_no_zeros[,c(2:15)]))
sample_totals$sample = rownames(sample_totals)
colnames(sample_totals) = c("total", "sample")
sample_totals$sample = factor(sample_totals$sample, levels=c(sample_list))

# heatmap of hgcA total expression
hgcA_plot = ggplot(hgcA_counts.m, aes(fct_rev(Group.1), y=variable, fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis",alpha=1, begin=0, end=0.95, direction=-1) + theme_bw()
hgcA_plot2 = hgcA_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
hgcA_plot2

# with no zeros of total expression by phyla
no_zeros_plot = ggplot(hgcA_no_zeros.m, aes(x=(Group.1), y=fct_rev(variable), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis",alpha=1, begin=0, end=1, direction=-1, limits=c(0,165), breaks=c(0, 25, 50, 75, 100, 165), rescaler = function(x, to = c(0, 0.8), from = NULL) {
  ifelse(x<100, 
         scales::rescale(x,
                         to = to,
                         from = c(min(x, na.rm = TRUE), 100)),
         1)
  }) + theme_bw()
no_zeros_plot
no_zeros_plot2 = no_zeros_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10)) 
no_zeros_plot3 = no_zeros_plot2 + theme_gray() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 7))
no_zeros_plot3

# total hgcA by sample
sample_totals = as.data.frame(colSums(hgcA_no_zeros[,c(2:18)]))
sample_totals$sample = rownames(sample_totals)
colnames(sample_totals) = c("total", "sample")
sample_totals$sample = factor(sample_totals$sample, levels=c(sample_list))

# plot bar graph of total hgcA per sample
colors = c("bog", "bog", "fen", "fen", "bog", "bog", "bog", "fen", "fen", "fen","fen","fen","fen","fen", "fen", "fen", "fen")
hgcA_sample = sample_totals %>% ggplot(aes(x=fct_rev(sample), y=total, fill=colors)) + geom_bar(stat="identity", width=.70) + coord_flip() + scale_y_continuous(limits=c(0,200), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + scale_fill_manual("legend", values=c("bog"="darkgreen","fen"="midnightblue"))
hgcA_sample

# average expression of each methylator
avg_expr = meth_expression_average %>% ggplot(aes(x=Group.1, y=avg_expression)) + geom_bar(stat="identity", fill="purple4") + theme(axis.text.x= element_text(angle=85, hjust=1))
avg_expr

# save plots indvidually
ggsave(file="~/Desktop/hgcA-TPM-totals.png", hgcA_sample, width=15, height=15, units=c("cm"))
ggsave(file="~/Desktop/hgcA-expression-heatmap.png", no_zeros_plot2, width=15, height=10, units=c("cm"))
ggsave(file="~/Desktop/average-expression-methylators.png", avg_expr, width=15, height=15, units=c("cm"))

# combine plots
# cleanup heatmap
tmp <- ggplot_gtable(ggplot_build(no_zeros_plot3))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
hm.clean <- no_zeros_plot3 +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")
hm.clean
# x axis clean
avg_clean = avg_expr + theme(axis.line=element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks=element_blank(),
                             axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                             legend.position="none",
                             panel.background=element_blank(),
                             panel.border=element_blank(),
                             panel.grid.major=element_blank(),
                             panel.grid.minor=element_blank(),
                             plot.background=element_blank())
# y axis clean
hgcA_sample_clean = hgcA_sample + theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.title.y=element_blank(),
                                      legend.position="none",
                                      panel.background=element_blank(),
                                      panel.border=element_blank(),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      plot.background=element_blank())
hgcA_sample_clean
# together
grid = grid.arrange(avg_clean, legend, hm.clean, hgcA_sample_clean, nrow=2, ncol=2, widths=c(30,40), heights=c(40,60))
# save grid 
ggsave(file="figs/2020-05-22-permafrost-transcription-grid.png", grid, height=10, width=15, units=c("cm"))


### comparison to housekeeping gene expression
rpoB = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files_and_results/permafrostTranscription/rpoB-locus-tags.tsv", header=TRUE, sep="\t")
peat_rpoB = rpoB %>% select(c(Phyla, locus_tag))
rpoB_counts = left_join(peat_rpoB, counttable.ids)
rpoB_no_zeros = rpoB_counts %>% select(c(Phyla, locus_tag, 4,6,9,10,12,14,16,17,18,20:26,28))
rpoB_agg = aggregate(rpoB_no_zeros[3:19], list(rpoB_no_zeros$Phyla), sum)
rpoB.m = melt(rpoB_agg, id.vars="Group.1")
rpoB.m$variable = factor(rpoB.m$variable, levels=c(sample_list))
rpo.sub = rpoB.m %>% filter(value < 1000)

hgcA.rpo = hgcA_no_zeros.m %>% filter(Group.1 %in% c("Acidobacteria", "Actinobacteria", "Aminicenantes", "Bacteroidetes", "Chlorobi", "Chloroflexi", "Deltaproteobacteria", "Elusimicrobia", "Euryarchaeota", "FCPU426", "Fibrobacteres", "Nitrospirae"))
hgcA.test = hgcA.rpo %>% select(variable, value)
write.csv(hgcA.rpo, "~/Desktop/hgcA-counts.csv", quote=FALSE, row.names=FALSE)
write.csv(rpo.sub, "~/Desktop/rpo-counts.csv", quote=FALSE, row.names=FALSE)
colnames(hgcA.test) = c("variable", "hgcA.counts")
merged.counts = left_join(rpo.sub, hgcA.test)

hgcA.rpo$marker = c("hgcA")
colnames(hgcA.rpo)[1] = c("phyla")
rpo.sub$marker = c("rpo")
colnames(rpo.sub)[1] = c("phyla")

# separate totals
  #rpoB
rpoB.plt = ggplot(rpo.sub, aes(x=Group.1, y=value)) + geom_boxplot()
rpoB.plt + geom_jitter(shape=15, position=position_jitter(0.2)) + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))
  # hgcA
hgcA.rpo.plt = ggplot(hgcA.rpo, aes(x=Group.1, y=value)) + geom_boxplot()
hgcA.rpo.plt + geom_jitter(shape=15, position=position_jitter(0.2)) + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))

# together with variable as marker instead of sampling depth
combined = rbind(hgcA.rpo, rpo.sub)
comb.plt = ggplot(combined, aes(x=phyla, y=value, fill=marker)) + geom_boxplot() + geom_jitter()
comb.plt
comb.form = comb.plt + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))
comb.form
ggsave(comb.form, filename="/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/permafrostTranscription/hgcA-vs-rpoB-expression.png", width=30, height=20, units=c("cm"))

# divide hgcA by rpoB
hgcA.counts <- hgcA.rpo %>% select(phyla, variable, value)
rpo.counts <- rpo.sub %>% select(phyla, variable, value)
colnames(hgcA.counts)[3] <- c("hgcA.val")
colnames(rpo.counts)[3] <- c("rpo.val")
hgcA_rpo_combined <- cbind(hgcA.counts, rpo.counts$rpo.val)
colnames(hgcA_rpo_combined) <- c("phyla", "sample", "hgcA.val", "rpo.val")
hgcA_rpo_combined$hgcA_over_rpo <- hgcA_rpo_combined$hgcA.val / hgcA_rpo_combined$rpo.val
hgcA_rpo_combined[is.na(hgcA_rpo_combined)] <- 0
hgcA_rpo_combined$sample = factor(hgcA_rpo_combined$sample, levels=c(sample_list))
hgcA_rpo_combined_cleaned <- hgcA_rpo_combined %>% mutate(sample=str_replace_all(sample, "abundance.", ""))

# Plot of hgcA normalized by rpoB facet by sample and color by phyla

hgcA_rpoB_facet <- hgcA_rpo_combined_cleaned %>% ggplot(aes(x=hgcA.val, y=rpo.val, color=phyla)) + geom_point() + facet_wrap( ~ sample, scales="free") + 
  scale_x_continuous(limits=c(0,175)) + scale_y_continuous(limits=c(0,450))

hgcA_vs_rpoB_all <- hgcA_rpo_combined %>% ggplot(aes(x=hgcA.val, y=rpo.val, color=phyla)) + geom_point(size=3) + theme_classic() + theme(legend.position="bottom")

ggsave(filename="figs/hgcA_vs_rpoB_facet_sample.png", plot=hgcA_rpoB_facet, width=20, height=15, units=c("cm"))
ggsave(filename="figs/hgcA_vs_rpoB_all.png", plot=hgcA_vs_rpoB_all, width=15, height=8, units=c("cm"))
ggsave(filename="figs/hgcA_vs_rpoB_all_configured.png", plot=hgcA_vs_rpoB_all, width=10, height=15, units=c("cm"))

# hgcB locus tags 

old_metadata <- read.csv("files_and_results/archived/archived-metadata/medium-qual-methylator-metadata.csv") %>% select(genome_name, original_name)
peat_split <- separate(data= peat_hgcA, col=locus_tag, into=c("genome_name", "hgcA_locus_tag"), sep="_", extra="merge")
old_metadata_peat <- left_join(peat_split, old_metadata)

hgcAB_loci <- read.csv("files_and_results/annotations/hgcAB_loci.csv")
colnames(hgcAB_loci)[1] = c("original_name")

peat_hgcAB_loci <- left_join(old_metadata_peat, hgcAB_loci)
hgcB_locus_tags <- peat_hgcAB_loci %>% select(Phylum, genome_name, hgcb_locus_tag)
hgcB_locus_tags$hgcB_tag <- paste(hgcB_locus_tags$genome_name, hgcB_locus_tags$hgcb_locus_tag, sep="_")

peat_hgcB_tags <- hgcB_locus_tags %>% select(Phylum, hgcB_tag)
colnames(peat_hgcB_tags) <- c("Phylum", "locus_tag")

hgcB_counts <- left_join(peat_hgcB_tags, counttable.ids) %>% drop_na()

hgcB_phyla_total = aggregate(hgcB_counts[3:28], list(hgcB_counts$Phylum), sum)

hgcB_cross <- hgcB_counts[-2]
hgcB.m <- melt(hgcB_phyla_total, id.vars="Group.1", var="sample")
hgcB.m$sample = factor(hgcB.m$sample, levels=c(sample_list))
hgcB_table = hgcB.m %>% drop_na()

hgcB_expression <- hgcB_table %>%  ggplot(aes(x=Group.1, y=fct_rev(sample), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis",alpha=1, begin=0, end=1, direction=-1, limits=c(0,220), breaks=c(0, 25, 50, 75, 100, 120, 220), rescaler = function(x, to = c(0, 0.8), from = NULL) {
  ifelse(x<120, 
         scales::rescale(x,
                         to = to,
                         from = c(min(x, na.rm = TRUE), 100)),
         1)
}) + theme_bw() + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))

# regulator locus tags

regulator_metadata = read.csv("files_and_results/annotations/flanking-regulator-metadata.csv") %>% select(genome, reg_locus_tag)
colnames(regulator_metadata) = c("original_name", "reg_locus_tag")
peat_regulator = left_join(old_metadata_peat, regulator_metadata) %>% drop_na() %>% select(Phylum, genome_name, reg_locus_tag)
peat_regulator$locus_tag = paste(peat_regulator$genome_name, peat_regulator$reg_locus_tag, sep="_")
reg_tags = peat_regulator %>% select(Phylum, locus_tag)
reg_counts = left_join(reg_tags, counttable.ids) %>% drop_na()
reg_phyla_total = aggregate(reg_counts[3:28], list(reg_counts$Phylum), sum)
reg.m = melt(reg_phyla_total, id.vars="Group.1", var="sample")
reg.m$sample = factor(reg.m$sample, levels=c(sample_list))
reg_table = reg.m %>% drop_na()

reg_expression <- reg_table %>%  ggplot(aes(x=Group.1, y=(fct_rev(sample)), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis",alpha=1, begin=0, end=1, direction=-1, limits=c(0,100)) + theme_bw() + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))

ggsave(file="figs/hgcB-expression.png", plot=hgcB_expression, width=15, height=10, units=c("cm"))
ggsave(file="figs/reg-expression.png", plot=reg_expression, width=15, height=10, units=c("cm"))
