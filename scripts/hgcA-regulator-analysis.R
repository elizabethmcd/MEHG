# merge regulator and hgca loci hits
library(tidyverse)
library(plyr)

reg = read.csv("files/annotations/reg_locus_tags.txt", sep="\t", header=FALSE)
hgca = read.csv("files/annotations/hgcA_list_locus_tags.txt", sep="\t", header=FALSE)
reg1 = reg %>% select(V1,V3)
hgca1 = hgca %>% select(V1,V3)
colnames(reg1) = c("genome", "reg_locus_tag")
colnames(hgca1) = c("genome", "hgca_locus_tag")
merged = left_join(reg1, hgca1)
# write out results to check if flanking
write.csv(merged, file="files/annotations/flanking-hgcA-regulator-results.csv", quote=FALSE, row.names=FALSE)
# merge with metadata
metadata = read.csv("files/stats/mehg-final-dataset-metadata.csv")
names = metadata %>% select(genome_name, Phyla, code)
colnames(names) = c("genome", "Group", "origin")
merge_names = left_join(merged, names)
# write out results with metadata
write.csv(merge_names, file="files/annotations/flanking-regulator-metadata.csv", quote=FALSE, row.names=FALSE)

# read in the counts
counts = read.csv("files/annotations/flanking-regulator-metadata-flanking-info.csv")
flanking = count(counts, "flanking")
# flanking true 112 total 283, 12 with 1 or 2 genes between regulator and hgcA

# merge with regulator annotations
annotations = read.csv("files/annotations/regulator_annotations.txt", sep="\t", header=FALSE)
colnames(annotations) = c("genome", "reg_locus_tag", "annotation")
annot_merged = left_join(merge_names, annotations)
annot_counts = left_join(counts, annotations)
flanking_true = annot_counts %>% filter(flanking == 't')

# hgcB results with hgcA results
hgcb = read.csv("files/annotations/hgcB_locus_tags.txt", sep="\t", header=FALSE)
colnames(hgcb) = c("genome", "hgcb_locus_tag")
hgcab = left_join(hgca1, hgcb)
write.csv(hgcab, file="files/annotations/hgcAB_loci.csv", quote=FALSE, row.names=FALSE)

