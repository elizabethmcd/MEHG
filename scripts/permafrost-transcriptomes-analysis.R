library(tximport)
library(readr)
library(tibble)
library(dplyr)
library(reshape2)
library(gtable)
library(grid)
library(gridExtra)
library(viridis)
library(ggplot2)

# kallisto directory for mapped transcriptomes to methylators
meth_dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/kallisto_mapping"
meth_samples <- read.table(file.path(meth_dir, "permafrost-metadata.txt"), header=TRUE)
meth_files <- file.path(meth_dir, meth_samples$sample, "abundance.h5")
rownames(meth_samples) <- meth_samples$sample
names(meth_files) <- meth_samples$sample
# check if in the cloud backup
meth.kallisto <- tximport(meth_files, type="kallisto", txOut = TRUE)

# counts file
meth.counts <- as.data.frame(meth.kallisto)
final.meth.counts <- rownames_to_column(meth.counts, var="genome_name")
meth.counttable <- final.meth.counts[, c(1, 28:53)]

# ids to pick certain genomes
ids <- read.delim("~/Desktop/McMahon-Lab/MeHg-Projects/MEHG/kallisto_mapping/perma-meth-ids.txt", sep="\t", header=TRUE)
counttable.ids <- left_join(meth.counttable, ids)

# TPM normalization
## Merge with annotations to divide by gene length for counts in all samples
# also check if this is in the cloud
meth_prokka = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/permafrostTranscription/permafrost-methylator-prokka-annotations.txt", sep="\t", header=FALSE)
colnames(meth_prokka) <- c("genome_name", "prokka_annotation", "size_bp", "accession")
meth_counts_annots <- left_join(meth_prokka, counttable.ids)
meth_counts_annots$size_kbp = meth_counts_annots$size_bp / 1000

# Divide by gene lengths
meth_counts_rpk = as.data.frame(lapply(meth_counts_annots[,c(-1,-2,-3, -4, -31, -32)], function(x) {
  (x / meth_counts_annots$size_kbp)
}))
# per million factor
meth_counts_rpk["PM", ] = (colSums(meth_counts_rpk[1:26]) / 1000000)
# divide by per million factor
meth_counts_tpm = as.data.frame(lapply(meth_counts_rpk, function(x) x/tail(x,1) ))
# get metadata back to aggregate by genome_id
meth_counts_tpm = meth_counts_tpm[-387690, ]
meth_counts_tpm$locus_tag = meth_counts_annots$genome_name
meth_counts_tpm$genome_name = meth_counts_annots$genome
meth_counts_totals = aggregate(meth_counts_tpm[1:26], list(meth_counts_tpm$genome_name), sum)

# permafrost results from archived set of metadata
# merge with metadata to get genome sizes to get relative expression by phyla
meth_metadata = read.csv("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/archived/archived-metadata/archived-methylator-metadata.csv")
peat_phyla = meth_metadata %>% select("genome_name", "Phylum")
colnames(meth_counts_totals)[1] = "genome_name"
meth_counts_totals_table = left_join(meth_counts_totals, peat_phyla)
meth_expression_average = aggregate(meth_counts_totals_table[2:27], list(meth_counts_totals_table$Phylum), sum)
meth_expression_average$avg_expression = rowSums(meth_expression_average[2:27]) / 26

# hgcA locus tags 
hgcA = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/archived/tree-tax-files/hgcA-locus-tags-phyla.txt", sep="\t", header=FALSE)
colnames(hgcA) = c("genome_name", "locus_tag", "Phylum")
hgcA = hgcA %>% select(genome_name, locus_tag)
mehg_peat = meth_metadata %>% filter(Study=="Woodcroft2018") 
peat_locus_tags = left_join(mehg_peat, hgcA)
peat_hgcA = peat_locus_tags %>% select(c(Phylum, locus_tag))

# counts of hgcA
meth_counts_norm = meth_counts_tpm
hgcA_counts = left_join(peat_hgcA, meth_counts_norm)
hgcA_non_agg_zeros = hgcA_counts %>% select(c(Phylum, locus_tag, 4,6,9,10,12,14,16,17,18,20:26,28))
hgcA_phyla_total = aggregate(hgcA_counts[3:28], list(hgcA_counts$Phylum),sum)
hgcA_no_zeros = hgcA_phyla_total %>% select(c(Group.1, 3,5,8,9,11,13,15,16,17,19:25,27))
hgcA_counts.m = melt(hgcA_phyla_total, id.vars="Group.1")
hgcA_no_zeros.m = melt(hgcA_no_zeros, id.vars="Group.1")
sample_list = c("counts.20120600_S2M",
                "counts.20120800_S1M",
                "counts.20120600_S2D",
                "counts.20120700_S1D",
                "counts.20120800_S1X",
                "counts.20100900_E2S",
                "counts.20110600_E1M",
                "counts.20120600_E1M",
                "counts.20120700_E3M",
                "counts.20120800_E2M",
                "counts.20100900_E1D",
                "counts.20110600_E1D",
                "counts.20110700_E1D",
                "counts.20120600_E1D",
                "counts.20120700_E3D",
                "counts.20120800_E2D",
                "counts.20120800_E3D")
hgcA_no_zeros.m$variable = factor(hgcA_no_zeros.m$variable, levels = c(sample_list))

# total hgcA by sample
sample_totals = as.data.frame(colSums(hgcA_no_zeros[,c(2:15)]))
sample_totals$sample = rownames(sample_totals)
colnames(sample_totals) = c("total", "sample")
sample_totals$sample = factor(sample_totals$sample, levels=c(sample_list))

# heatmap of hgcA total expression
hgcA_plot = ggplot(hgcA_counts.m, aes(x=fct_rev(Group.1), y=variable, fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="magma",alpha=1, begin=0, end=0.95, direction=-1) + theme_bw()
hgcA_plot2 = hgcA_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10))
hgcA_plot2

# with no zeros of total expression by phyla
no_zeros_plot = ggplot(hgcA_no_zeros.m, aes(x=Group.1, y=fct_rev(variable), fill=value)) + geom_tile(color="white") + scale_fill_viridis(option="viridis",alpha=1, begin=0, end=1, direction=-1) + theme_bw()
no_zeros_plot2 = no_zeros_plot + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 10)) 
no_zeros_plot3 = no_zeros_plot2 + theme_gray() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 7))
no_zeros_plot2

# total hgcA by sample
sample_totals = as.data.frame(colSums(hgcA_no_zeros[,c(2:18)]))
sample_totals$sample = rownames(sample_totals)
colnames(sample_totals) = c("total", "sample")
sample_totals$sample = factor(sample_totals$sample, levels=c(sample_list))

# plot bar graph of total hgcA per sample
colors = c("bog", "bog", "fen", "fen", "bog", "bog", "bog", "fen", "fen", "fen","fen","fen","fen","fen", "fen", "fen", "fen")
hgcA_sample = sample_totals %>% ggplot(aes(x=fct_rev(sample), y=total, fill=colors)) + geom_bar(stat="identity", width=.70) + coord_flip() + scale_y_continuous(limits=c(0,350), expand= c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + scale_fill_manual("legend", values=c("bog"="darkgreen","fen"="midnightblue"))
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
ggsave(file="~/Desktop/mehg-grid-axes.png", grid, height=10, width=15, units=c("cm"))


### comparison to housekeeping gene expression
rpoB = read.delim("/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/permafrostTranscription/rpoB-locus-tags.tsv", header=TRUE, sep="\t")
peat_rpoB = rpoB %>% select(c(Phyla, locus_tag))
rpoB_counts = left_join(peat_rpoB, meth_counts_norm)
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

# separate tottals
  #rpoB
rpoB.plt = ggplot(rpo.sub, aes(x=Group.1, y=value)) + geom_boxplot()
rpoB.plt + geom_jitter(shape=15, position=position_jitter(0.2)) + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))
  # hgcA
hgcA.rpo.plt = ggplot(hgcA.rpo, aes(x=Group.1, y=value)) + geom_boxplot()
hgcA.rpo.plt + geom_jitter(shape=15, position=position_jitter(0.2)) + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))

# together with variable as marker instead of sampling depth
hgcA.var = read.csv("~/Desktop/hgcA-counts-var.csv")
rpoB.var = read.csv("~/Desktop/rpo-counts-var.csv")
combined = rbind(hgcA.var, rpoB.var)
comb.plt = ggplot(combined, aes(x=phyla, y=value, fill=variable)) + geom_boxplot()
comb.form = comb.plt + theme_classic() + theme(axis.text.x=element_text(angle=85, hjust=1))
ggsave(comb.form, filename="/Users/emcdaniel/Desktop/McMahon-Lab/MeHg-Projects/MEHG/files/permafrostTranscription/hgcA-vs-rpoB-expression.png", width=30, height=20, units=c("cm"))
