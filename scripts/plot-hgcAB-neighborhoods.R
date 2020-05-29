library(ggplot2)
library(gggenes)

genes <- read.csv("files_and_results/annotations/ars/hgcAB-ars-neighborhoods.csv")

hgcA_alignment <- make_alignment_dummies(genes, aes(xmin=start, xmax=end, y=genome, id=gene), on = "hgcA")

ars <- ggplot(genes, aes(xmin=start, xmax=end, y=genome, fill=gene, label=gene)) + geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + geom_blank(data=hgcA_alignment) + facet_wrap(~ genome, scales="free", ncol = 1) + scale_fill_manual(values=c("mediumorchid3", "purple4", "plum4", "slateblue3", "cyan4", "lightseagreen", "snow3", "beige", "turquoise")) + theme_genes()

ggsave(filename="figs/hgcAB-ars-neighborhoods.png", plot=ars, width=17, height=15, units=c("cm"))
