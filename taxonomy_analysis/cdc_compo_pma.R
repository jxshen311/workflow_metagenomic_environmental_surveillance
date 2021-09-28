# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 
# purpose: 



# Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(openxlsx) # handle xlsx files
library(tidyr)
library(rstatix)
library(ggpubr)
library(vegan)
library(tibble)
library(stringr)
library(reshape2)
library("RColorBrewer")
library(Rfast)
library(cowplot)
library(ape)  # pcoa
library(scales)  # extarct default colors of R
library(ggtext)
library(ggh4x)
library(Cairo)


# Functions ----




# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/compo_analysis")



# Import data ----
level7 <- read.xlsx("input/merged_SILVA132_level7.xlsx",
                    sheet = "merged_SILVA132_level7_bacteria", startRow = 1, colNames = TRUE,
                    rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                    skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                    check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

# remove all singletons
for (ii in 8:ncol(level7)) {level7[which(level7[,ii]<2),ii] <- 0}

# * to the finest available level (treat "unclassified..." & "... sp" the same)
level7 <- level7

# make "unclassified..." & "... sp" to NA
level7$species[which(grepl("sp", level7$species, fixed = TRUE))] <- NA

for (ii in 1:7){
  level7[which(grepl("Unclassified", level7[ ,ii], fixed = TRUE)), ii] <- NA
}

tag <- c("k","p",	"c",	"o",	"f",	"g",	"s")

for (ii in 1:7){
  level7[,ii] <- paste(tag[ii],level7[,ii],sep="_")
  level7[which(grepl("_NA", level7[ ,ii], fixed = TRUE)), ii] <- NA
}

level7 <- level7 %>% 
  group_by(domain, phylum,	class,	order,	family,	genus,	species) %>% 
  summarise(across(everything(), sum))


# calculate relative abundance
for(i in 8:ncol(level7)) {level7[,i] <- level7[,i]/sum(level7[,i])}

level7 <- level7 %>%
  unite(taxa, colnames(level7[,1:7]), sep = ";", remove = FALSE)

level7$taxa <- gsub(";NA","",level7$taxa)

level7$lowest_level <- gsub(".*;", "", level7$taxa)

level7 <- level7 %>%
  select(taxa, lowest_level, everything())  # re-order the columns

level7_fi <- level7 %>%
  ungroup() %>%
  select(-c("taxa", "domain","phylum","class","order","family", "genus","species"))

level7_fi <- column_to_rownames(level7_fi ,var = "lowest_level")

level7_fi <- as.data.frame(t(level7_fi))

level7_fi <- level7_fi[rowSums(is.na(level7_fi)) != ncol(level7_fi), ]  # Remove rows/columns where all data is NA

# Alpha diversity ----------------------------------------------------------------------------------------------------------------
# * Shannon ----
level7_fi_alpha <- as.data.frame(diversity(level7_fi, index="shannon"))
colnames(level7_fi_alpha) <- "shannon"

level7_fi_alpha <- rownames_to_column(level7_fi_alpha, var = "sample")

level7_fi_alpha$sample_2 <- sapply(str_split(level7_fi_alpha$sample, "[.]",  n = 2), `[`, 1)

mean(level7_fi_alpha[which(level7_fi_alpha$sample_2 == "SN"), 2])  # output: 2.764197; shannon of SP is not applicable

level7_fi_alpha_pma <- level7_fi_alpha %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
  mutate(category = substring(sample_2,1,1),
         pma = substring(sample_2,2,2))

level7_fi_alpha_pma <- level7_fi_alpha_pma %>% mutate_if(is.character, as.factor)

p_alpha <- ggboxplot(level7_fi_alpha_pma, x = "pma", y = "shannon", fill = "category", palette = "npg", legend = "none", xlab = "Sample", ylab = "Shannon diversity index", ggtheme = theme_pubr(border = TRUE), facet.by = c("category"), panel.labs = list(category = c("Spike+", "Spike-")), panel.labs.font = list(size = 12)) + scale_x_discrete(labels=c("N" = "PMA-", "P" = "PMA+")) + theme(panel.spacing=unit(0, "lines"))



p_alpha

ggpar(p_alpha,
      font.x = c(14),
      font.y = c(14),
      font.xtickslab = c(12),
      font.ytickslab = c(12))

#// ggsave("out_fig_pma/shannon.pdf", width =  4, height = 4, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/shannon.pdf", width =  5, height = 5, unit = "in", dpi = 300)

# * Richness ----
level7_fi_richness <- as.data.frame(apply(level7_fi[ , ]>0, 1, sum))  

colnames(level7_fi_richness) <- "richness"

level7_fi_richness <- rownames_to_column(level7_fi_richness, var = "sample")

level7_fi_richness$sample_2 <- sapply(str_split(level7_fi_richness$sample, "[.]",  n = 2), `[`, 1)

mean(level7_fi_richness[which(level7_fi_richness$sample_2 == "SN"), 2])  # 8

level7_fi_richness_pma <- level7_fi_richness %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
  mutate(fill = substring(sample_2,1,1))

mean_level7_fi_richness <- aggregate(richness ~  sample_2, level7_fi_richness_pma, mean)
mean_level7_fi_richness$richness <- round(mean_level7_fi_richness$richness, digits = 2)


(p_richness <- ggplot(level7_fi_richness_pma, aes(x=sample_2, y=richness)) +
    geom_boxplot(aes(fill = fill)) +
    geom_text(data = mean_level7_fi_richness, aes(label = richness, y = richness + 2)) +
    theme_bw() +
    labs(x="Sample", y="Richness") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          panel.grid = element_blank(), 
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none"))

#// ggsave("out_fig_pma/richness.pdf", width =  4, height = 4, unit = "in", dpi = 300)




# Beta diversity -----------------------------------------------------------------------------------------------------------------
# * Jaccard -----------------------------------------------------------------------------------------------------------------
dist_jaccard <- vegdist(level7_fi,  method = "jaccard")


pcoa_jaccard <- pcoa(dist_jaccard, correction="cailliez")
pcoa_jaccard$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_jaccard_axes <- as.data.frame(pcoa_jaccard$vectors[,c(1,2)])
pcoa_jaccard_axes$sample <- rownames(pcoa_jaccard_axes)

pcoa_jaccard_axes$sample_2 <- sapply(str_split(pcoa_jaccard_axes$sample, "[.]",  n = 2), `[`, 1)

pcoa_jaccard_axes_pma <- pcoa_jaccard_axes %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
  mutate(category = substring(sample_2,1,1), pma = substring(sample_2,2,2))



(p_beta_jaccard <- ggplot(pcoa_jaccard_axes_pma, aes(x=Axis.1, y=Axis.2, color=category, shape = pma)) +
    geom_point(size = 5, alpha = 0.8)+
    scale_color_discrete(name="Internal standard addition", breaks = c("I", "N"), labels = c("Y", "N") ) +
    scale_shape_discrete(name="PMA treatment", breaks = c("P", "N"), labels = c("Y", "N")) +
    theme_bw() +
    labs(x="Axis.1 (41.16%)", y="Axis.2 (7.85%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "bottom")) 

#// ggsave("out_fig_pma/beta_jaccard.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# box plot of distance
dist_jaccard_pma <- as.data.frame(as.matrix(dist_jaccard))

dist_jaccard_pma <- dist_jaccard_pma[c("IN.1", "IN.2", "IN.3", "IP.1", "IP.2", "IP.3", "NN.1", "NN.2", "NN.3", "NP.1", "NP.2", "NP.3"),c("IN.1", "IN.2", "IN.3", "IP.1", "IP.2", "IP.3", "NN.1", "NN.2", "NN.3", "NP.1", "NP.2", "NP.3")]

dist_jaccard_pma <- as.matrix(dist_jaccard_pma)
dist_jaccard_pma_melt <- melt(dist_jaccard_pma, varnames = c("row", "col"))


#removing self distances and duplicated distances
dist_jaccard_pma_melt  <- dist_jaccard_pma_melt[-which(duplicated(dist_jaccard_pma_melt$value) | dist_jaccard_pma_melt[,1] == dist_jaccard_pma_melt[,2]), ] 

dist_jaccard_pma_melt  <- dist_jaccard_pma_melt[-which(substring(dist_jaccard_pma_melt[,1], 1, 1) != substring(dist_jaccard_pma_melt[,2], 1, 1)), ]

dist_jaccard_pma_melt <- dist_jaccard_pma_melt %>%
  mutate(category = substring(row,1,1), 
         pma = ifelse(substring(row,2,2) == substring(col,2,2), 
                      substring(row,2,2),
                      "Between"))

dist_jaccard_pma_melt$pma <- factor(dist_jaccard_pma_melt$pma, levels = c("N", "P", "Between"))

dist_jaccard_pma_melt$category <- as.factor(dist_jaccard_pma_melt$category)


p_dist_jaccard <- ggboxplot(dist_jaccard_pma_melt, x = "pma", y = "value", fill = "category", palette = "npg", legend = "none", xlab = "Sample", ylab = "Jaccard distance", ggtheme = theme_pubr(border = TRUE), facet.by = c("category"), panel.labs = list(category = c("Spike+", "Spike-")), panel.labs.font = list(size = 12)) + scale_x_discrete(labels=c("N" = "PMA-", "P" = "PMA+", "Between" = "Between")) + theme(panel.spacing=unit(0, "lines"))

ggpar(p_dist_jaccard,
      font.x = c(14),
      font.y = c(14),
      font.xtickslab = c(12),
      font.ytickslab = c(12))

#// ggsave("out_fig_pma/distance_jaccard.pdf", width =  5, height = 5, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/b_distance_jaccard.pdf", width =  5, height = 5, unit = "in", dpi = 300)


# * Bray-Curtis -----------------------------------------------------------------------------------------------------------------
dist_bray <- vegdist(level7_fi,  method = "bray")

pcoa_bray <- pcoa(dist_bray, correction="cailliez")
pcoa_bray$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_bray_axes <- as.data.frame(pcoa_bray$vectors[,c(1,2)])
pcoa_bray_axes$sample <- rownames(pcoa_bray_axes)

pcoa_bray_axes$sample_2 <- sapply(str_split(pcoa_bray_axes$sample, "[.]",  n = 2), `[`, 1)

pcoa_bray_axes_pma <- pcoa_bray_axes %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
  mutate(category = substring(sample_2,1,1), pma = substring(sample_2,2,2))



(p_beta_bray <- ggplot(pcoa_bray_axes_pma, aes(x=Axis.1, y=Axis.2, color=category, shape = pma)) +
    geom_point(size = 5, alpha = 0.8)+
    scale_color_discrete(name="Internal standard addition", breaks = c("I", "N"), labels = c("Y", "N") ) +
    scale_shape_discrete(name="PMA treatment", breaks = c("P", "N"), labels = c("Y", "N")) +
    theme_bw() +
    labs(x="Axis.1 (47.01%)", y="Axis.2 (7.17%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "bottom")) 

#// ggsave("out_fig_pma/beta_bray.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# box plot of distance
dist_bray_pma <- as.data.frame(as.matrix(dist_bray))

dist_bray_pma <- dist_bray_pma[c("IN.1", "IN.2", "IN.3", "IP.1", "IP.2", "IP.3", "NN.1", "NN.2", "NN.3", "NP.1", "NP.2", "NP.3"),c("IN.1", "IN.2", "IN.3", "IP.1", "IP.2", "IP.3", "NN.1", "NN.2", "NN.3", "NP.1", "NP.2", "NP.3")]

dist_bray_pma <- as.matrix(dist_bray_pma)
dist_bray_pma_melt <- melt(dist_bray_pma, varnames = c("row", "col"))


#removing self distances and duplicated distances
dist_bray_pma_melt  <- dist_bray_pma_melt[-which(duplicated(dist_bray_pma_melt$value) | dist_bray_pma_melt[,1] == dist_bray_pma_melt[,2]), ] 

dist_bray_pma_melt  <- dist_bray_pma_melt[-which(substring(dist_bray_pma_melt[,1], 1, 1) != substring(dist_bray_pma_melt[,2], 1, 1)), ]

dist_bray_pma_melt <- dist_bray_pma_melt %>%
  mutate(category = substring(row,1,1), 
         pma = ifelse(substring(row,2,2) == substring(col,2,2), 
                      substring(row,2,2),
                      "Between"))

dist_bray_pma_melt$pma <- factor(dist_bray_pma_melt$pma, levels = c("N", "P", "Between"))



p_dist_bray <- ggboxplot(dist_bray_pma_melt, x = "pma", y = "value", fill = "category", palette = "npg", legend = "none", xlab = "Sample", ylab = "Bray–Curtis distance", ggtheme = theme_pubr(border = TRUE)) + 
  facet_wrap(
    nrow = 1,
    ncol = 2,
    ~category)

ggpar(p_dist_bray,
      font.x = c(14),
      font.y = c(14),
      font.xtickslab = c(12),
      font.ytickslab = c(12))

#// ggsave("out_fig_pma/distance_bray.pdf", width =  5, height = 5, unit = "in", dpi = 300)


# Biomass ---------------------------------------------------------------------------------------------------------------------------
qpcr <- read.xlsx("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/qpcr/cdc_qpcr_sample.xlsx",
                    sheet = 1, startRow = 1, colNames = TRUE,
                    rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                    skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                    check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)[,c(1,9,10)]

qpcr <- dplyr::rename(qpcr, "copy_number_per_250uL_sample" = "copy_number_per_250uL_sample.(account.for.concentration)")

# change CDC to AG
qpcr$sample[qpcr$sample == "CDC_1"] <- "AG_1"
qpcr$sample[qpcr$sample == "CDC_2"] <- "AG_2"
qpcr$sample[qpcr$sample == "CDC_3"] <- "AG_3"

qpcr$sample_2 <- sapply(strsplit(qpcr$sample, "_"), `[`, 1)

qpcr_pma <- qpcr %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN", "SN", "SP")) %>%
  mutate(category = substring(sample_2, 1,1), pma = substring(sample_2, 2,2))

qpcr_pma$pma[which(qpcr_pma$pma == "P")] <- "Y"


# * DNA ----
ttest.dna <- qpcr_pma %>%
  group_by(category) %>%
  t_test(ng_dna_per_250uL_sample ~ pma,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

p_dna <- ggbarplot(qpcr_pma, x = "category", y = "ng_dna_per_250uL_sample", fill = "pma", position = position_dodge(0.8), add = c("mean_se"), error.plot = "errorbar", palette = "jco", legend = "bottom", legend.title = "PMA treatment", xlab = "Sample type", ylab = "DNA amount per 250 \u03BCL sample (ng)", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("I" = "Sample + internal\nstandard", "N" = "Sample", "S" = "External standard"))

ttest.dna <- ttest.dna %>% add_xy_position(x = "category")

p_dna <- p_dna + stat_pvalue_manual(ttest.dna, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain") +
  font("legend.title", size = 12, face = "plain") +
  font("legend.text", size = 12, face = "plain")

#// ggsave("out_fig_pma/biomass_dna.pdf", width =  5, height = 5, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/biomass_dna.pdf", plot=p_dna, device=cairo_pdf, width = 5, height = 5, units = "in", dpi = 300) 

# standard reduction
mean(qpcr_pma$ng_dna_per_250uL_sample[which(qpcr_pma$sample_2 == "SP")])  # 12.00237
(mean(qpcr_pma$ng_dna_per_250uL_sample[which(qpcr_pma$sample_2 == "SN")]) - mean(qpcr_pma$ng_dna_per_250uL_sample[which(qpcr_pma$sample_2 == "SP")])) / mean(qpcr_pma$ng_dna_per_250uL_sample[which(qpcr_pma$sample_2 == "SN")])  # 0.954094

# * qpcr copy number ----
ttest.qpcr <- qpcr_pma %>%
  group_by(category) %>%
  t_test(copy_number_per_250uL_sample ~ pma,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

p_qpcr <- ggbarplot(qpcr_pma, x = "category", y = "copy_number_per_250uL_sample", fill = "pma", position = position_dodge(0.8), add = c("mean_se"), error.plot = "errorbar", palette = "jco", legend = "bottom", legend.title = "PMA treatment", xlab = "Sample type", ylab = "16S rRNA gene copy number per 250 \u03BCL sample", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("I" = "Sample + internal\nstandard", "N" = "Sample", "S" = "External standard"))

ttest.qpcr <- ttest.qpcr %>% add_xy_position(x = "category")

p_qpcr <- p_qpcr +  stat_pvalue_manual(ttest.qpcr, label = "p.adj.signif") +
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  font("xy.text", size = 12, face = "plain") +
  font("legend.title", size = 12, face = "plain") +
  font("legend.text", size = 12, face = "plain")

#// ggsave("out_fig_pma/biomass_qpcr.pdf", width =  5, height = 5, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/biomass_qpcr.pdf", plot=p_qpcr, device=cairo_pdf, width = 5, height = 5, units = "in", dpi = 300) 



# standard reduction
mean(qpcr_pma$copy_number_per_250uL_sample[which(qpcr_pma$sample_2 == "SP")])  # 629254.2
(mean(qpcr_pma$copy_number_per_250uL_sample[which(qpcr_pma$sample_2 == "SN")]) - mean(qpcr_pma$copy_number_per_250uL_sample[which(qpcr_pma$sample_2 == "SP")])) / mean(qpcr_pma$copy_number_per_250uL_sample[which(qpcr_pma$sample_2 == "SN")])  # 0.9998381


# Abundance change ------------------------------------------------------------------------------------------------------------------
# * Relative abundance --------------------------------------------------------------------------------------------------------------
level7_fi_pma <- level7_fi

level7_fi_pma <- rownames_to_column(level7_fi_pma, var = "sample")
level7_fi_pma$sample_2 <- sapply(str_split(level7_fi_pma$sample, "[.]",  n = 2), `[`, 1)

level7_fi_pma <- level7_fi_pma %>%
  filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
  select(-sample_2) %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  as.data.frame()

# order the taxa according to average abundance
level7_fi_pma$avg <- apply(level7_fi_pma, 1, mean)

level7_fi_pma <- level7_fi_pma[order(level7_fi_pma$avg, decreasing = T),]

apply(level7_fi_pma, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_fi_pma$avg != 0))  # 30
length(which(level7_fi_pma$avg == 0))  # 51
length(which(level7_fi_pma$avg > 0.001))  # 23
length(which(level7_fi_pma$avg > 0.01))  # 6


# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_fi_pma_mod <- level7_fi_pma[which(level7_fi_pma$avg > 0.001), ]

level7_fi_pma_mod = rbind(level7_fi_pma_mod,
                          1-apply(level7_fi_pma_mod, 2, sum))
rownames(level7_fi_pma_mod)[nrow(level7_fi_pma_mod)] = c("Others")

# delete avg column
level7_fi_pma_mod <- select(level7_fi_pma_mod, -avg)

# melt
level7_fi_pma_melt <- melt(t(level7_fi_pma_mod))

# add category and pma columns for facet_wrap()
level7_fi_pma_melt$category <- as.factor(substring(level7_fi_pma_melt$Var1,1,1))
level7_fi_pma_melt$pma <- as.factor(substring(level7_fi_pma_melt$Var1,2,4))
level7_fi_pma_melt$pma_group <- as.factor(substring(level7_fi_pma_melt$pma,1,1))

# barplot, with expanded color palette
colourCount = length(unique(level7_fi_pma_melt$Var2))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

# generate a number of most distinctive colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, colourCount)
show_col(col)

# specify facet labeller
labs.category <- c("Spike+", "Spike-")
names(labs.category) <- c("I", "N")

labs.pma_group <- c("PMA-", "PMA+")
names(labs.pma_group) <- c("N", "P")


p_rel_abun_bar <- ggplot(level7_fi_pma_melt, aes(x = pma, y = value, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.95) +
  labs(title = "", x = "", y = "Relative Abundance") +
  theme_bw() +  
  scale_fill_manual(values = col) +
  # scale_fill_manual(values = getPalette(colourCount)) +
  scale_x_discrete(expand = c(0.18,0.18)) + 
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "bottom",
        panel.spacing=unit(0, "lines"),
        strip.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow =6,byrow=TRUE)) + 
  facet_nested(~category + pma_group, scales = "free_x", labeller = labeller(category = labs.category, pma_group = labs.pma_group))
  #facet_grid(~category + pma_group, scales = "free_x")

p_rel_abun_bar

#// ggsave("out_fig_pma/rel_abun_bar.pdf", width =  9, height = 9, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abs_abun_bar.pdf", plot=p_abs_abun_bar, device=cairo_pdf, width = 10, height = 9, units = "in", dpi = 300) 





# * Rare taxa distribution ----------------------------------------------------------------------------------------------------------
taxa_no <- as.data.frame(colnames(level7_fi_pma)[1:ncol(level7_fi_pma)-1])
taxa_no <- rename(taxa_no, "sample" = "colnames(level7_fi_pma)[1:ncol(level7_fi_pma) - 1]")

for (ii in 1:(ncol(level7_fi_pma)-1)){
  taxa_no$abun_no[ii] <- sum(level7_fi_pma[ ,ii] >= 0.01 , na.rm=TRUE)
  taxa_no$moderate_no[ii] <- sum(level7_fi_pma[ ,ii] < 0.01 & level7_fi_pma[ ,ii] >= 0.0001, na.rm=TRUE)
  taxa_no$rare_1_no[ii] <- sum(level7_fi_pma[ ,ii] < 0.001 & level7_fi_pma[ ,ii] > 0, na.rm=TRUE)
  taxa_no$rare_2_no[ii] <- sum(level7_fi_pma[ ,ii] < 0.0001 & level7_fi_pma[ ,ii] > 0, na.rm=TRUE)
  taxa_no$no_zero_no[ii] <- sum(level7_fi_pma[ ,ii] != 0, na.rm=TRUE)
}
# no rare taxa, only moderate and abundant taxa




# * Absolute abundance --------------------------------------------------------------------------------------------------------------
qpcr$sample_3 <- gsub("_", ".", qpcr$sample)

level7_fi_pma_mod_abs <- level7_fi_pma_mod

for(ii in 1:ncol(level7_fi_pma_mod_abs)) {level7_fi_pma_mod_abs[ ,ii] <- level7_fi_pma_mod_abs[ ,ii] * qpcr[which(qpcr$sample_3 ==  colnames(level7_fi_pma_mod_abs)[ii] ) ,2] }

apply(level7_fi_pma_mod_abs , 2, sum)  # check

# level7_fi_pma_abs <- level7_fi
# 
# level7_fi_pma_abs <- rownames_to_column(level7_fi_pma_abs, var = "sample")
# level7_fi_pma_abs$sample_2 <- sapply(str_split(level7_fi_pma_abs$sample, "[.]",  n = 2), `[`, 1)
# 
# level7_fi_pma_abs <- level7_fi_pma_abs %>%
#   filter(sample_2 %in% c("IP", "IN", "NP", "NN")) %>%
#   select(-sample_2) %>%
#   column_to_rownames(var = "sample") 

# for(ii in 1:nrow(level7_fi_pma_abs)) {level7_fi_pma_abs[ii, ] <- level7_fi_pma_abs[ii, ] * qpcr[which(qpcr$sample_3 ==  rownames(level7_fi_pma_abs)[ii] ) ,2] }
# 
# apply(level7_fi_pma_abs , 1, sum)  # check
# 
# level7_fi_pma_abs <- as.data.frame(t(level7_fi_pma_abs))

level7_fi_pma_mod_abs_melt <- melt(t(level7_fi_pma_mod_abs))

# add category and pma columns for facet_wrap()
level7_fi_pma_mod_abs_melt$category <- as.factor(substring(level7_fi_pma_mod_abs_melt$Var1,1,1))
level7_fi_pma_mod_abs_melt$pma <- as.factor(substring(level7_fi_pma_mod_abs_melt$Var1,2,4))
level7_fi_pma_mod_abs_melt$pma_group <- as.factor(substring(level7_fi_pma_mod_abs_melt$pma,1,1))

# barplot, with expanded color palette
colourCount = length(unique(level7_fi_pma_mod_abs_melt$Var2))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

# specify facet labeller
labs.category <- c("Spike+", "Spike-")
names(labs.category) <- c("I", "N")

labs.pma_group <- c("PMA-", "PMA+")
names(labs.pma_group) <- c("N", "P")

p_abs_abun_bar <- ggplot(level7_fi_pma_mod_abs_melt, aes(x = pma, y = value, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.95) +
  labs(title = "", x = "", y = "Absolute abundance\n(16S rRNA gene copy number per 250 \u03BCL sample)") +
  theme_bw() +  scale_fill_manual(values = col) +
  scale_x_discrete(expand = c(0.18,0.18)) + 
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "bottom",
        panel.spacing=unit(0, "lines"),
        strip.text = element_text(size=12)) + 
  guides(fill=guide_legend(nrow =6,byrow=TRUE)) + 
  facet_nested(~category + pma_group, scales = "free_x", labeller = labeller(category = labs.category, pma_group = labs.pma_group))
 # facet_grid(~category + pma_group, scales = "free_x")

p_abs_abun_bar 

#// ggsave("out_fig_pma/abs_abun_bar.pdf", width =  10, height = 9, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abs_abun_bar.pdf", plot=p_abs_abun_bar, device=cairo_pdf, width = 10, height = 9, units = "in", dpi = 300) 



# * Arrange relative abundance + absolute abundance ----
p_abs_rel <- ggarrange(p_abs_abun_bar, p_rel_abun_bar, nrow = 1,
          common.legend = TRUE, legend="bottom")

#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abs_rel_abun_bar.pdf", plot=p_abs_rel, device=cairo_pdf, width = 10, height = 8, units = "in", dpi = 300) 

# heatmap
p_abs_abun_heatmap <- ggplot(level7_fi_pma_mod_abs_melt, aes(x = pma, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "", x = "Sample", y = "Taxa") +
  scale_fill_gradient(name = "",
                      trans = "log",
                      low = "black",
                      high = "red",
                      labels = scientific) +
  scale_y_discrete(limits = rev(levels(level7_fi_pma_mod_abs_melt$Var2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, face = "plain", angle = 0, hjust = 0),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + 
  facet_wrap(
    nrow = 1,
    ncol = 2,
    ~category)

p_abs_abun_heatmap

#// ggsave("out_fig_pma/abs_abun_heat.pdf", width =  11, height = 9, unit = "in", dpi = 300)




# IS taxa efficacy ------------------------------------------------------------------------------------------------------------------
# level7_fi_pma_is <- filter(level7_fi_pma, apply(level7_fi_pma[ ,7:12], 1, sum)  == 0)
level7_fi_pma_is <- select(level7_fi_pma, -avg)
                           
for(ii in 1:ncol(level7_fi_pma_is)) {level7_fi_pma_is[ ,ii] <- level7_fi_pma_is[ ,ii] * qpcr[which(qpcr$sample_3 ==  colnames(level7_fi_pma_is)[ii] ) ,2] }

apply(level7_fi_pma_is, 2, sum)  # check

#// write.xlsx(rownames_to_column(level7_fi_pma_is, var = "sample"), file = "out_metaxa_level7_fi_pma_abs.xlsx")




# Absolute abundance (all taxa, not group to "Others")-------------------------------------------------------------------------------
level7_fi_pma_abs <- level7_fi_pma_is

level7_fi_pma_abs$avg = apply(level7_fi_pma_abs, 1, mean)

# check the number of taxa
length(which(level7_fi_pma_abs$avg != 0))  # 30
length(which(level7_fi_pma_abs$avg == 0))  # 51
# test <- list()
# for(ii in 1:ncol(level7_fi_pma_is)){test[[ii]] <- sum(level7_fi_pma_abs[ ,ii] != 0)}

level7_fi_pma_abs <- level7_fi_pma_abs[order(level7_fi_pma_abs$avg, decreasing = T), ]
level7_fi_pma_abs <- level7_fi_pma_abs[which(level7_fi_pma_abs$avg != 0), ]

level7_fi_pma_abs <- select(level7_fi_pma_abs,-avg)

level7_fi_pma_abs_melt <- melt(t(level7_fi_pma_abs))

# add category and pma columns for facet_wrap()
level7_fi_pma_abs_melt$category <- as.factor(substring(level7_fi_pma_abs_melt$Var1,1,1))
level7_fi_pma_abs_melt$pma <- as.factor(substring(level7_fi_pma_abs_melt$Var1,2,4))

# heatmap
p_abs_abun_heatmap_2 <- ggplot(level7_fi_pma_abs_melt, aes(x = pma, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "", x = "Sample", y = "Taxa") +
  scale_fill_gradient(name = "",
                      trans = "log",
                      low = "black",
                      high = "red",
                      labels = scientific) +
  scale_y_discrete(limits = rev(levels(level7_fi_pma_abs_melt$Var2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_markdown(size = 12, face = "italic", colour = c('black',	'black',	'black',	'#6600CC',	'black',	'black',	'#6600CC',	'#009900',	'#6600CC',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#009900',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'black',	'black',	'black',	'#009900',	'black',	'black',	'black')),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, face = "plain", angle = 0, hjust = 0),
        legend.position = "right") + 
  facet_wrap(
    nrow = 1,
    ncol = 2,
    ~category)

#// ggsave("out_fig_pma/abs_abun_heat_not_group_others.pdf", width =  9, height = 7, unit = "in", dpi = 300)



# Abs abun change + all 30 taxa -----------------------------------------------------------------------------------------------
level7_fi_pma_abs_melt$pma_2 <- as.factor(substring(level7_fi_pma_abs_melt$pma,1,1))

level7_fi_pma_abs_melt_stat <- level7_fi_pma_abs_melt %>%
  group_by(Var2, category, pma_2) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(length(value)))
 
level7_fi_pma_abs_melt_stat$pma_2 <- factor(level7_fi_pma_abs_melt_stat$pma_2, levels = c("P", "N"))

level7_fi_pma_abs_melt_stat <- level7_fi_pma_abs_melt_stat %>%
  group_by(Var2, category) %>%
  filter(sum(value_mean)>0)

level7_fi_pma_abs_melt_stat$zero <- as.factor(ifelse(level7_fi_pma_abs_melt_stat$value_mean == 0, "Y", "N"))


  
# "Escherichia coli" is recovered in cultivation but not included here. Instead, it is considered as IS taxon. Because it is not detected in NP & NN.

p_abs_abun_change <- ggplot(level7_fi_pma_abs_melt_stat, aes(x = value_mean, y = Var2, color = pma_2, shape = zero)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(y=Var2, xmin=value_mean-value_se, xmax=value_mean+value_se), width=.2, position=position_dodge(0.05)) +
  geom_path(aes(group = Var2), arrow=arrow(length = unit(0.04, "inches"), type = "open"), color = "black") + 
  labs(title = "", x = "Absolute abundance", y = "Taxa") +
  scale_y_discrete(limits = rev(levels(level7_fi_pma_abs_melt$Var2))) +
  scale_color_discrete(name="PMA", breaks = c("P", "N"), labels = c("Yes", "No")) +
  scale_shape_discrete(name="Detection", breaks = c("N", "Y"), labels = c("Yes", "No")) +
 # scale_shape_discrete(guide = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_markdown(size = 12, face = "italic", colour = c('black',	'black',	'black',	'#6600CC',	'black',	'black',	'#6600CC',	'#009900',	'#6600CC',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#009900',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'black',	'black',	'black',	'#009900',	'black',	'black',	'black')),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12, face = "plain"),
        legend.text = element_text(size = 12, face = "plain", angle = 0, hjust = 0),
        legend.position = "right",
        panel.spacing=unit(0, "lines"),
        strip.text = element_text(size=12)) + 
  facet_grid(~category, labeller = labeller(category = c(I = "Spike+", N = "Spike-")))

# \n(16S rRNA gene copy number per 250 \u03BCL sample)
p_abs_abun_change


#// ggsave("out_fig_pma/abs_abun_change.pdf", width =  9, height = 7, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abs_abun_change.pdf", plot=p_abs_abun_change, device=cairo_pdf, width = 9, height = 7, units = "in", dpi = 300) 


## check where absolute abundance increases
abs_abun_change <- select(level7_fi_pma_abs_melt_stat, 1:4)

abs_abun_change  <- abs_abun_change  %>%
  group_by(Var2, category) %>%
  pivot_wider(names_from = pma_2, values_from = value_mean)

abs_abun_change <- mutate(abs_abun_change, PN = P-N, PN_percent = P/N)

# which(abs_abun_change$PN_percent <= 0.1 & abs_abun_change$PN_percent > 0)

#// write.xlsx(abs_abun_change, file = "out_pma_abs_abun_change.xlsx")
## check done



# Rel abun change + all 30 taxa -----------------------------------------------------------------------------------------------
level7_fi_pma <- level7_fi_pma[which(level7_fi_pma$avg != 0), ]

level7_fi_pma_all_taxa_melt <- melt(t(select(level7_fi_pma, -avg)))

# add category and pma columns for facet_wrap()
level7_fi_pma_all_taxa_melt$category <- as.factor(substring(level7_fi_pma_all_taxa_melt$Var1,1,1))
level7_fi_pma_all_taxa_melt$pma_2 <- as.factor(substring(level7_fi_pma_all_taxa_melt$Var1,2,2))

level7_fi_pma_all_taxa_melt_stat <- level7_fi_pma_all_taxa_melt %>%
  group_by(Var2, category, pma_2) %>%
  summarise(value_mean = mean(value), value_se = sd(value)/sqrt(length(value)))

level7_fi_pma_all_taxa_melt_stat <- level7_fi_pma_all_taxa_melt_stat %>%
  group_by(Var2, category) %>%
  filter(sum(value_mean)>0)

level7_fi_pma_all_taxa_melt_stat$pma_2 <- factor(level7_fi_pma_all_taxa_melt_stat$pma_2, levels = c("P", "N"))
level7_fi_pma_all_taxa_melt_stat$zero <- as.factor(ifelse(level7_fi_pma_all_taxa_melt_stat$value_mean == 0, "Y", "N"))


# "Escherichia coli" is recovered in cultivation but not included here. Instead, it is considered as IS taxon. Because it is not detected in NP & NN.

p_rel_abun_change <- ggplot(level7_fi_pma_all_taxa_melt_stat, aes(x = value_mean, y = Var2, color = pma_2, shape = zero)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(y=Var2, xmin=value_mean-value_se, xmax=value_mean+value_se), width=.2, position=position_dodge(0.05)) +
  geom_path(aes(group = Var2), arrow=arrow(length = unit(0.04, "inches"), type = "open"), color = "black") + 
  labs(title = "", x = "Relative abundance", y = "Taxa") +
  xlim(0,0.5) +
  scale_y_discrete(limits = rev(levels(level7_fi_pma_abs_melt$Var2))) +
  scale_color_discrete(name="PMA", breaks = c("P", "N"), labels = c("Yes", "No")) +
  scale_shape_discrete(name="Detection", breaks = c("N", "Y"), labels = c("Yes", "No")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_markdown(size = 12, face = "italic", colour = c('black',	'black',	'black',	'#6600CC',	'black',	'black',	'#6600CC',	'#009900',	'#6600CC',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#009900',	'#6600CC',	'black',	'#6600CC',	'#6600CC',	'#6600CC',	'black',	'#6600CC',	'black',	'black',	'black',	'#009900',	'black',	'black',	'black')),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12, face = "plain"),
        legend.text = element_text(size = 12, face = "plain", angle = 0, hjust = 0),
        legend.position = "right",
        panel.spacing=unit(0, "lines"),
        strip.text = element_text(size=12)) + 
  facet_grid(~category, labeller = labeller(category = c(I = "Spike+", N = "Spike-")))

p_rel_abun_change

#// ggsave("out_fig_pma/rel_abun_change.pdf", width =  9, height = 7, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/rel_abun_change.pdf", plot=p_rel_abun_change, device=cairo_pdf, width = 9, height = 7, units = "in", dpi = 300) 

# Arrange change (relative abundance + absolute abundance) ----
# p_abs_rel_change <- ggarrange(p_abs_abun_change, p_rel_abun_change + rremove("ylab") + rremove("y.text") + rremove("y.ticks"), nrow = 1, common.legend = TRUE, legend="right", align = "v")

p1 <- p_abs_abun_change + 
  theme(plot.margin=unit(c(5.5,-30,5.5,5.5), "pt"))

p2 <- p_rel_abun_change + 
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.margin=unit(c(5.5,5.5,5.5,-30), "pt"))


p_abs_rel_change <- grid.arrange(p1, p2, ncol=2, nrow =1, widths=c(7.1, 5), heights= 7 )

p_abs_rel_change

#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abs_rel_abun_change_point.pdf", plot=p_abs_rel_change, device=cairo_pdf, width = 14, height = 7, units = "in", dpi = 300)



## check where relative abundance increases
rel_abun_change <- select(level7_fi_pma_all_taxa_melt_stat, 1:4)

rel_abun_change  <- rel_abun_change  %>%
  group_by(Var2, category) %>%
  pivot_wider(names_from = pma_2, values_from = value_mean)

rel_abun_change <- mutate(rel_abun_change, PN = P-N, PN_percent = P/N)

#// write.xlsx(rel_abun_change, file = "out_pma_rel_abun_change.xlsx")
## check done


# Test ------------------------------------------------------------------------------------------------------------------------------
test <- level7_fi_pma[,1:6]
test <- as.data.frame(t(test))

dist_bray_test <- vegdist(test,  method = "bray")

