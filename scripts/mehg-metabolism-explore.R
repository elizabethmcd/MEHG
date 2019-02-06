# Characterizing Metabolic Capabilities of Methylating Organisms 

library(tidyverse)
library(reshape2)

metabolic_matrix = read.table("mehg-metabolic-marker-results.csv", sep=",", header=TRUE)
methylating_metadata = read.table("methylating-genomes-phyla.csv", sep=",", header=TRUE)

# Merge datasets 

metabolic_table = left_join(methylating_metadata, metabolic_matrix)
metabolic_table_phyla = metabolic_table %>% arrange(phylum)
without_phyla = subset(metabolic_table_phyla, select=-c(taxonomy, phylum))
metabolic_melted = melt(without_phyla, id.vars="genome")
melt_test = melt(metabolic_matrix, id.vars='genome')

# Large heatmap

p1 <- ggplot(metabolic_melted, aes(x=variable, y=genome, fill=value)) + geom_tile(color="white") + scale_fill_gradient(low="lightblue", high="midnightblue")
p2 <- p1 + theme(axis.text.x= element_text(angle=85, hjust=1)) + guides(fill = guide_colorbar(nbin = 4))
p2
ggsave("test.png", p2, width=60, height=90, units="cm")

# Carbon metabolism, monoxide, and fixation dataframes
carbon_metabolism = read.table("central-carbon-markers.csv", sep=",", header=TRUE)
carbon_monoxide = read.table("carbon-monoxide-markers.csv", sep=",", header=TRUE)
carbon_fixation = read.table("carbon-fixation-markers.csv", sep=",", header=TRUE)

central_carbon = left_join(methylating_metadata, carbon_metabolism)
CMO = left_join(methylating_metadata, carbon_monoxide)
CO2 = left_join(methylating_metadata, carbon_fixation)

# Deltaproteobacteria carbon dataframes 
central_delta = central_carbon %>% filter(phylum=="Deltaproteobacteria") %>% subset(select=-c(taxonomy, phylum))
cmo_delta = CMO %>% filter(phylum=="Deltaproteobacteria") %>% subset(select=-c(taxonomy, phylum))
co2_delta = CO2 %>% filter(phylum=="Deltaproteobacteria") %>% subset(select=-c(taxonomy, phylum))

central_delta_m = melt(central_delta, id.vars="genome")
cmo_delta_m = melt(cmo_delta, id.vars="genome")
co2_delta_m = melt(co2_delta, id.vars="genome")

central_delta_plot = ggplot(central_delta_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="darkcyan") + theme(axis.text.x= element_text(angle=85, hjust=1))
cmo_delta_plot = ggplot(cmo_delta_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="chartreuse1") + theme(axis.text.x= element_text(angle=85, hjust=1))
co2_delta_plot = ggplot(co2_delta_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="violetred4") + theme(axis.text.x= element_text(angle=85, hjust=1))

ggsave("delta-central-carbon.png", central_delta_plot, width=10, height=60, units="cm")
ggsave("delta-CMO.png", cmo_delta_plot, width=7, height=60, units="cm")
ggsave("delta-CO2.png", co2_delta_plot, width=10, height=60, units="cm")

# Methane and focus on Euryarchaeota 
methane_markers = read.table("methane-markers.csv", sep=",", header=TRUE)
methane_meta = left_join(methylating_metadata, methane_markers)

methane_arch = methane_meta %>% filter(phylum=="Euryarchaeota") %>% subset(select=-c(taxonomy, phylum))
central_arch = central_carbon %>% filter(phylum=="Euryarchaeota") %>% subset(select=-c(taxonomy, phylum))
cmo_arch = CMO %>% filter(phylum=="Euryarchaeota") %>% subset(select=-c(taxonomy, phylum))
co2_arch = CO2 %>% filter(phylum=="Euryarchaeota") %>% subset(select=-c(taxonomy, phylum))

methane_arch_m = melt(methane_arch, id.vars="genome")
central_arch_m = melt(central_arch, id.vars="genome")
cmo_arch_m = melt(cmo_arch, id.vars="genome")
co2_arch_m = melt(co2_arch, id.vars="genome")

methane_arch_plot = ggplot(methane_arch_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="gray20") + theme(axis.text.x= element_text(angle=85, hjust=1))
central_arch_plot = ggplot(central_arch_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="deepskyblue4") + theme(axis.text.x= element_text(angle=85, hjust=1))
co2_arch_plot = ggplot(co2_arch_m, aes(x=variable, y=genome, fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="deeppink4") + theme(axis.text.x= element_text(angle=85, hjust=1))

ggsave("methane-arch.png", methane_arch_plot, width=10, height=30, units="cm")
ggsave("central-arch.png", central_arch_plot, width=10, height=30, units="cm")
ggsave("co2-arch.png", co2_arch_plot, width=10, height=30, units="cm")

# Test guilds Deltaproteobacteria and Nitrospirae hgcA most similar
central_guild = central_carbon %>% filter(genome %in% c("bacteria08141",
                                                        "bacteria07777",
                                                        "bacteria07753",
                                                        "bacteria07729",
                                                        "bacteria07725",
                                                        "bacteria07754",
                                                        "bacteria07713",
                                                        "bacteria07719",
                                                        "bacteria08100",
                                                        "bacteria08090",
                                                        "bacteria08033",
                                                        "bacteria08032",
                                                        "bacteria08027",
                                                        "bacteria08020",
                                                        "bacteria08017",
                                                        "bacteria08015",
                                                        "bacteria08687",
                                                        "bacteria08685",
                                                        "bacteria15437"))
cmo_guild = CMO %>% filter(genome %in% c("bacteria08141",
                                         "bacteria07777",
                                         "bacteria07753",
                                         "bacteria07729",
                                         "bacteria07725",
                                         "bacteria07754",
                                         "bacteria07713",
                                         "bacteria07719",
                                         "bacteria08100",
                                         "bacteria08090",
                                         "bacteria08033",
                                         "bacteria08032",
                                         "bacteria08027",
                                         "bacteria08020",
                                         "bacteria08017",
                                         "bacteria08015",
                                         "bacteria08687",
                                         "bacteria08685",
                                         "bacteria15437"))
co2_guild = CO2 %>% filter(genome %in% c("bacteria08141",
                                         "bacteria07777",
                                         "bacteria07753",
                                         "bacteria07729",
                                         "bacteria07725",
                                         "bacteria07754",
                                         "bacteria07713",
                                         "bacteria07719",
                                         "bacteria08100",
                                         "bacteria08090",
                                         "bacteria08033",
                                         "bacteria08032",
                                         "bacteria08027",
                                         "bacteria08020",
                                         "bacteria08017",
                                         "bacteria08015",
                                         "bacteria08687",
                                         "bacteria08685",
                                         "bacteria15437"))                                      
sulfur = read.table("mehg-sulfur-markers.csv", sep=",", header=TRUE)
nitrogen = read.table("mehg-nitrogen-markers.csv", sep=",", header=TRUE)
S = left_join(methylating_metadata, sulfur)
N = left_join(methylating_metadata, nitrogen)
S_guild = S %>% filter(genome %in% c("bacteria08141",
                                     "bacteria07777",
                                     "bacteria07753",
                                     "bacteria07729",
                                     "bacteria07725",
                                     "bacteria07754",
                                     "bacteria07713",
                                     "bacteria07719",
                                     "bacteria08100",
                                     "bacteria08090",
                                     "bacteria08033",
                                     "bacteria08032",
                                     "bacteria08027",
                                     "bacteria08020",
                                     "bacteria08017",
                                     "bacteria08015",
                                     "bacteria08687",
                                     "bacteria08685",
                                     "bacteria15437"))
N_guild = N %>% filter(genome %in% c("bacteria08141",
                                     "bacteria07777",
                                     "bacteria07753",
                                     "bacteria07729",
                                     "bacteria07725",
                                     "bacteria07754",
                                     "bacteria07713",
                                     "bacteria07719",
                                     "bacteria08100",
                                     "bacteria08090",
                                     "bacteria08033",
                                     "bacteria08032",
                                     "bacteria08027",
                                     "bacteria08020",
                                     "bacteria08017",
                                     "bacteria08015",
                                     "bacteria08687",
                                     "bacteria08685",
                                     "bacteria15437"))
N_counts = subset(N_guild, select=-c(phylum, taxonomy))
S_counts = subset(S_guild, select=-c(phylum,taxonomy))
n_m = melt(N_counts, id.vars="genome")
s_m = melt(S_counts, id.vars="genome")
n_m$genome = factor(n_m$genome, levels = c("bacteria08141",
                                           "bacteria07777",
                                           "bacteria07753",
                                           "bacteria07729",
                                           "bacteria07725",
                                           "bacteria07754",
                                           "bacteria07713",
                                           "bacteria07719",
                                           "bacteria08100",
                                           "bacteria08090",
                                           "bacteria08033",
                                           "bacteria08032",
                                           "bacteria08027",
                                           "bacteria08020",
                                           "bacteria08017",
                                           "bacteria08015",
                                           "bacteria08687",
                                           "bacteria08685",
                                           "bacteria15437"))

s_m$genome = factor(s_m$genome, levels = c("bacteria08141",
                                           "bacteria07777",
                                           "bacteria07753",
                                           "bacteria07729",
                                           "bacteria07725",
                                           "bacteria07754",
                                           "bacteria07713",
                                           "bacteria07719",
                                           "bacteria08100",
                                           "bacteria08090",
                                           "bacteria08033",
                                           "bacteria08032",
                                           "bacteria08027",
                                           "bacteria08020",
                                           "bacteria08017",
                                           "bacteria08015",
                                           "bacteria08687",
                                           "bacteria08685",
                                           "bacteria15437"))

nitrogen_guild_plot = ggplot(n_m, aes(x=variable, y=fct_rev(genome), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="blue4") + theme(axis.text.x= element_text(angle=85, hjust=1))
nitrogen_guild_plot

sulfur_guild_plot = ggplot(s_m, aes(x=variable, y=fct_rev(genome), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="purple4") + theme(axis.text.x= element_text(angle=85, hjust=1))
sulfur_guild_plot

actino_core_guild = central_carbon %>% filter(genome %in% c("bacteria20678",
                                                            "bacteria20901",
                                                            "bacteria20275",
                                                            "bacteria20842",
                                                            "bacteria20480",
                                                            "bacteria20819",
                                                            "bacteria20689",
                                                            "bacteria20830",
                                                            "bacteria12873",
                                                            "bacteria13742",
                                                            "bacteria12801",
                                                            "bacteria07624",
                                                            "bacteria08617",
                                                            "bacteria12786"))
actino_s_guild = S %>% filter(genome %in% c("bacteria20678",
                                             "bacteria20901",
                                             "bacteria20275",
                                             "bacteria20842",
                                             "bacteria20480",
                                             "bacteria20819",
                                             "bacteria20689",
                                             "bacteria20830",
                                             "bacteria12873",
                                             "bacteria13742",
                                             "bacteria12801",
                                             "bacteria07624",
                                             "bacteria08617",
                                             "bacteria12786"))
actino_n_guild = N %>% filter(genome %in% c("bacteria20678",
                                            "bacteria20901",
                                            "bacteria20275",
                                            "bacteria20842",
                                            "bacteria20480",
                                            "bacteria20819",
                                            "bacteria20689",
                                            "bacteria20830",
                                            "bacteria12873",
                                            "bacteria13742",
                                            "bacteria12801",
                                            "bacteria07624",
                                            "bacteria08617",
                                            "bacteria12786"))
actino_core_counts = subset(actino_core_guild, select=-c(taxonomy,phylum))
actino_s_counts = subset(actino_s_guild, select=-c(taxonomy, phylum))
actino_n_counts = subset(actino_n_guild, select=-c(taxonomy, phylum))
actino_core_m = melt(actino_core_counts, id.vars="genome")
actino_sulfur_m = melt(actino_s_counts, id.vars="genome")
actino_n_m = melt(actino_n_counts, id.vars="genome")
actino_core_m$genome = factor(actino_core_m$genome, levels = c("bacteria20678",
                                                               "bacteria20901",
                                                               "bacteria20275",
                                                               "bacteria20842",
                                                               "bacteria20480",
                                                               "bacteria20819",
                                                               "bacteria20689",
                                                               "bacteria20830",
                                                               "bacteria12873",
                                                               "bacteria13742",
                                                               "bacteria12801",
                                                               "bacteria07624",
                                                               "bacteria08617",
                                                               "bacteria12786"))
actino_sulfur_m$genome = factor(actino_sulfur_m$genome, levels = c("bacteria20678",
                                                               "bacteria20901",
                                                               "bacteria20275",
                                                               "bacteria20842",
                                                               "bacteria20480",
                                                               "bacteria20819",
                                                               "bacteria20689",
                                                               "bacteria20830",
                                                               "bacteria12873",
                                                               "bacteria13742",
                                                               "bacteria12801",
                                                               "bacteria07624",
                                                               "bacteria08617",
                                                               "bacteria12786"))
actino_n_m$genome = factor(actino_n_m$genome, levels = c("bacteria20678",
                                                               "bacteria20901",
                                                               "bacteria20275",
                                                               "bacteria20842",
                                                               "bacteria20480",
                                                               "bacteria20819",
                                                               "bacteria20689",
                                                               "bacteria20830",
                                                               "bacteria12873",
                                                               "bacteria13742",
                                                               "bacteria12801",
                                                               "bacteria07624",
                                                               "bacteria08617",
                                                               "bacteria12786"))
actino_core_plot = ggplot(actino_core_m, aes(x=variable, y=fct_rev(genome), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="red4") + theme(axis.text.x= element_text(angle=85, hjust=1))
actino_s_plot = ggplot(actino_sulfur_m, aes(x=variable, y=fct_rev(genome), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="purple4") + theme(axis.text.x= element_text(angle=85, hjust=1))
actino_n_m = ggplot(actino_n_m, aes(x=variable, y=fct_rev(genome), fill=value)) + geom_tile(color="black") + scale_fill_gradient(low="white", high="navyblue") + theme(axis.text.x= element_text(angle=85, hjust=1))
actino_core_plot
actino_s_plot
actino_n_m

# All Actinobacteria summaries
methane_actino = methane_meta %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))
central_actino = central_carbon %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))
cmo_actino = CMO %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))
co2_actino = CO2 %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))
sulfur_meta = left_join(methylating_metadata, sulfur)
sulfur_actino = sulfur_meta %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))
nitrogen_meta = left_join(methylating_metadata, nitrogen)
nitrogen_actino = nitrogen_meta %>% filter(phylum=="Actinobacteria") %>% subset(select=-c(taxonomy, phylum))

