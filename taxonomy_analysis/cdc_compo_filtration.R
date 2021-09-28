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
library(dichromat)
library("RColorBrewer")
library(Rfast)
library(cowplot)
library(ape)  # pcoa
library(scales)  # extarct default colors of R
library("ggsci")
library(Cairo)
library(ggh4x)


# Functions ----

# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/compo_analysis")



# Import data ----
# check the data type of columns of interest by str()
# import data
level7 <- read.xlsx("input/merged_SILVA132_level7.xlsx",
                    sheet = "merged_SILVA132_level7_bacteria", startRow = 1, colNames = TRUE,
                    rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                    skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                    check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

# remove all singletons
for (ii in 8:ncol(level7)) {level7[which(level7[,ii]<2),ii] <- 0}

# * metaxa; not group ----
level7_id <- level7

# calculate relative abundance
for(i in 8:ncol(level7_id)) {level7_id[,i] <- level7_id[,i]/sum(level7_id[,i])}

level7_id <- level7_id[,-c(1:7)]


# * to the finest available level (treat "unclassified..." & "... sp" the same) ------------------------------------------------
level7_fi <- level7

# make "unclassified..." & "... sp" to NA
level7_fi$species[which(grepl("sp", level7_fi$species, fixed = TRUE))] <- NA

for (ii in 1:7){
    level7_fi[which(grepl("Unclassified", level7_fi[ ,ii], fixed = TRUE)), ii] <- NA
}
               
tag <- c("k","p",	"c",	"o",	"f",	"g",	"s")

for (ii in 1:7){
  level7_fi[,ii] <- paste(tag[ii],level7_fi[,ii],sep="_")
  level7_fi[which(grepl("_NA", level7_fi[ ,ii], fixed = TRUE)), ii] <- NA
}

level7_fi <- level7_fi %>% 
  group_by(domain, phylum,	class,	order,	family,	genus,	species) %>% 
  summarise(across(everything(), sum))

# # name NA as "unclassified to this level"
# level7$species[which(is.na(level7$species))] <- c("Unclassified to genus level")

# calculate relative abundance
for(i in 8:ncol(level7_fi)) {level7_fi[,i] <- level7_fi[,i]/sum(level7_fi[,i])}

# group the same entry together for "species" column
# level7 <- level7 %>%
#   group_by(species) %>%
#   summarise_at(colnames(level7)[2:ncol(level7)], sum)

level7_fi <- level7_fi %>%
  unite(taxa, colnames(level7_fi[,1:7]), sep = ";", remove = FALSE)

level7_fi$taxa <- gsub(";NA","",level7_fi$taxa)

level7_fi$lowest_level <- gsub(".*;", "", level7_fi$taxa)

level7_fi <- level7_fi %>%
  select(taxa, lowest_level, everything())  # re-oredr the columns

#// write.csv(level7_fi, file = "out_metaxa_level7_fi.csv")

# Filtration ----
# * metaxa - level7_id - effect on overall composition ----

level7_id_fil <- as.data.frame(t(level7_id))
level7_id_fil <- rownames_to_column(level7_id_fil, var = "sample")

level7_id_fil$sample_2 <- sapply(str_split(level7_id_fil$sample, "[.]", n = 2), `[`, 1)

level7_id_fil <- level7_id_fil %>%
  filter(sample_2 %in% c('NNF',
                         'NNN',
                         'NPF',
                         'NPN',
                         'INF',
                         'INN',
                         'IPF',
                         'IPN',
                         'SNF',
                         'SNN',
                         'SPF',
                         'SPN'))

level7_id_fil <- select(level7_id_fil, -sample_2)

level7_id_fil <- column_to_rownames(level7_id_fil,var = "sample")

level7_id_fil <- as.data.frame(t(level7_id_fil))

# order the taxa according to average abundance
level7_id_fil <- Filter(function(x)!all(is.na(x)), level7_id_fil)  # remove the two columns with all NA

level7_id_fil$avg <- apply(level7_id_fil, 1, mean)

level7_id_fil <- level7_id_fil[order(level7_id_fil$avg, decreasing = T),]

apply(level7_id_fil, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_id_fil$avg != 0))  # 43
length(which(level7_id_fil$avg == 0))  # 47
length(which(level7_id_fil$avg > 0.001))  # 32
length(which(level7_id_fil$avg > 0.01))  # 12

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_id_fil_mod <- level7_id_fil[which(level7_id_fil$avg>0.001), ]

level7_id_fil_mod = rbind(level7_id_fil_mod,
                          1-apply(level7_id_fil_mod, 2, sum))
rownames(level7_id_fil_mod)[nrow(level7_id_fil_mod)] = c("Others")

# delete avg column
level7_id_fil_mod <- select(level7_id_fil_mod, -avg)

# melt
level7_id_fil_melt <- melt(t(level7_id_fil_mod))

# add category and filtration columns for facet_wrap()
level7_id_fil_melt$category <- substring(level7_id_fil_melt$Var1,1,2)
level7_id_fil_melt$filtration <- substring(level7_id_fil_melt$Var1,3,nchar(as.character(level7_id_fil_melt$Var1)))

# plot, with expanded color palette
colourCount = length(unique(level7_id_fil_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_filtration_id <- ggplot(level7_id_fil_melt, aes(x = filtration, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Filtration", y = "Relative Abundance") +
  theme_bw() +  
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE)) +
  facet_wrap(
    nrow = 3,
    ncol = 2,
    ~category)

p_filtration_id

#// ggsave("out_fig/species_abundance_filtration_id_no_singleton.pdf", width =  9, height = 7, unit = "in", dpi = 300)

# ** PCoA plot ----
# *** level7_id_fil----
# **** level7_id_fil + jaccard ----
dist_id_fil <- vegdist(t(select(level7_id_fil,-avg)),  method = "jaccard")


pcoa_id_fil <- pcoa(dist_id_fil, correction="cailliez")
pcoa_id_fil$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_fil_axes <- as.data.frame(pcoa_id_fil$vectors[,c(1,2)])
pcoa_id_fil_axes$sample_id <- rownames(pcoa_id_fil_axes)

pcoa_id_fil_axes$category <- substring(pcoa_id_fil_axes$sample_id,1,2)
pcoa_id_fil_axes$category <- factor(pcoa_id_fil_axes$category, levels = c("IN", "IP", "NN", "NP", "SN", "SP"))

pcoa_id_fil_axes$filtration <- as.factor(substring(pcoa_id_fil_axes$sample_id,3,3))

(p_pcoa_id_fil <- ggplot(pcoa_id_fil_axes, aes(x=Axis.1, y=Axis.2, color=category, shape = filtration)) +
    geom_point(size = 5, alpha = 0.7)+
    scale_color_discrete(name="Category") +
    scale_shape_discrete(name="Filtration") + 
    theme_bw() +
    labs(x="Axis.1 (69.57%)", y="Axis.2 (12.89%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) 

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_1.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)

# **** level7_id_fil_sample + jaccard ----
dist_id_fil_sample <- vegdist(t(select(level7_id_fil,-c("SNF.1","SNF.2", "SNF.3", "SNN.1", "SNN.2", "SNN.3" ,"SPF.3" ,"SPN.3" ,"avg"))),  method = "jaccard")


pcoa_id_fil_sample <- pcoa(dist_id_fil_sample, correction="cailliez")
pcoa_id_fil_sample$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_fil_sample_axes <- as.data.frame(pcoa_id_fil_sample$vectors[,c(1,2)])
pcoa_id_fil_sample_axes$sample_id <- rownames(pcoa_id_fil_sample_axes)

pcoa_id_fil_sample_axes$category <- substring(pcoa_id_fil_sample_axes$sample_id,1,2)
pcoa_id_fil_sample_axes$category <- factor(pcoa_id_fil_sample_axes$category, levels = c("IN", "IP", "NN", "NP"))

pcoa_id_fil_sample_axes$filtration <- as.factor(substring(pcoa_id_fil_sample_axes$sample_id,3,3))

(p_pcoa_id_fil_sample <- ggplot(pcoa_id_fil_sample_axes, aes(x=Axis.1, y=Axis.2, color=category, shape = filtration)) +
    geom_point(size = 5, alpha = 0.7)+
    scale_shape_discrete(name="Filtration", breaks=c("F", "N"),labels=c("Yes", "No")) + 
    scale_colour_manual(name="Category", values = c("#F8766D","#B79F00", "#00BA38", "#00BFC4"),
                        breaks=c("IN", "IP", "NN", "NP"),
                        labels=c("S+,P-", "S+,P+", "S-,P-", "S-,P+")) +   
    theme_bw() +
    labs(x="Axis.1 (39.24%)", y="Axis.2 (16.13%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) 

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_jaccard_sample.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/c_filtration_id_no_singleton_pcoa_jaccard_sample.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# **** level7_id_fil_sample + Bray ----
dist_id_fil_sample_bray <- vegdist(t(select(level7_id_fil,-c("SNF.1","SNF.2", "SNF.3", "SNN.1", "SNN.2", "SNN.3" ,"SPF.3" ,"SPN.3" ,"avg"))),  method = "bray")


pcoa_id_fil_sample_bray <- pcoa(dist_id_fil_sample_bray, correction="cailliez")
pcoa_id_fil_sample_bray$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_fil_sample_bray_axes <- as.data.frame(pcoa_id_fil_sample_bray$vectors[,c(1,2)])
pcoa_id_fil_sample_bray_axes$sample_id <- rownames(pcoa_id_fil_sample_bray_axes)

pcoa_id_fil_sample_bray_axes$category <- substring(pcoa_id_fil_sample_bray_axes$sample_id,1,2)
pcoa_id_fil_sample_bray_axes$category <- factor(pcoa_id_fil_sample_bray_axes$category, levels = c("IN", "IP", "NN", "NP"))

pcoa_id_fil_sample_bray_axes$filtration <- as.factor(substring(pcoa_id_fil_sample_bray_axes$sample_id,3,3))

(p_pcoa_id_fil_sample_bray <- ggplot(pcoa_id_fil_sample_bray_axes, aes(x=Axis.1, y=Axis.2, color=category, shape = filtration)) +
    geom_point(size = 5, alpha = 0.7)+
    scale_shape_discrete(name="Filtration") + 
    scale_colour_manual(name="Category", values = c("#F8766D","#B79F00", "#00BA38", "#00BFC4")) +   
    theme_bw() +
    labs(x="Axis.1 (39.30%)", y="Axis.2 (15.88%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) 

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_bray_sample.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# **** level7_id_fil_sample + canberra ----
dist_id_fil_sample_canberra <- vegdist(t(select(level7_id_fil,-c("SNF.1","SNF.2", "SNF.3", "SNN.1", "SNN.2", "SNN.3" ,"SPF.3" ,"SPN.3" ,"avg"))),  method = "canberra")


pcoa_id_fil_sample_canberra <- pcoa(dist_id_fil_sample_canberra, correction="cailliez")
pcoa_id_fil_sample_canberra$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_fil_sample_canberra_axes <- as.data.frame(pcoa_id_fil_sample_canberra$vectors[,c(1,2)])
pcoa_id_fil_sample_canberra_axes$sample_id <- rownames(pcoa_id_fil_sample_canberra_axes)

pcoa_id_fil_sample_canberra_axes$category <- substring(pcoa_id_fil_sample_canberra_axes$sample_id,1,2)
pcoa_id_fil_sample_canberra_axes$category <- factor(pcoa_id_fil_sample_canberra_axes$category, levels = c("IN", "IP", "NN", "NP"))
pcoa_id_fil_sample_canberra_axes$filtration <- as.factor(substring(pcoa_id_fil_sample_canberra_axes$sample_id,3,3))


(p_pcoa_id_fil_sample_canberra <- ggplot(pcoa_id_fil_sample_canberra_axes, aes(x=Axis.1, y=Axis.2, color=category, shape = filtration)) +
    geom_point(size = 5, alpha = 0.7)+
    scale_shape_discrete(name="Filtration") + 
    scale_colour_manual(name="Category", values = c("#F8766D","#B79F00", "#00BA38", "#00BFC4")) +   
    theme_bw() +
    labs(x="Axis.1 (43.49%)", y="Axis.2 (12.20%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) 

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_canberra_sample.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# *** level7_id + jaccard ----
dist_id <- vegdist(t(Filter(function(x)!all(is.na(x)), level7_id)),  method = "jaccard")


pcoa_id <- pcoa(dist_id, correction="cailliez")

pcoa_id$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_axes <- as.data.frame(pcoa_id$vectors[,c(1,2)])
pcoa_id_axes$sample_id <- rownames(pcoa_id_axes)

pcoa_id_axes$sample_id <- sapply(str_split(pcoa_id_axes$sample_id, "[.]", n = 2), `[`, 1)

pcoa_id_axes$filtration <- ifelse(grepl("\\d", pcoa_id_axes$sample_id),
                                  "Filter",
                                  ifelse(nchar(pcoa_id_axes$sample_id) == 2, "N", substring(pcoa_id_axes$sample_id,3,3)))
pcoa_id_axes$filtration <- factor(pcoa_id_axes$filtration, levels = c("F", "N", "Filter"))



pcoa_id_axes$filtration_2 <- ifelse(nchar(pcoa_id_axes$sample_id) == 2, 
                                    "N", 
                                    substring(pcoa_id_axes$sample_id,3,nchar(pcoa_id_axes$sample_id)))
pcoa_id_axes$filtration_2 <- factor(pcoa_id_axes$filtration_2, levels = c("F", "N", "F100", "F80", "F41", "F5"))


pcoa_id_axes$category <- substring(pcoa_id_axes$sample_id,1,2)
pcoa_id_axes$category <- factor(pcoa_id_axes$category, levels = c( "IN", "IP", "NN", "NP", "SN", "SP", "AG"))


(p_pcoa_id <- ggplot(pcoa_id_axes, aes(x=Axis.1, y=Axis.2, color = category, shape = filtration)) +
    geom_point(size = 4, alpha = 0.7)+
    scale_colour_manual(name="Category", values = c("#F8766D", "#B79F00" ,"#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#671bd1")) +  
    scale_shape_discrete(name="Filtration") +
    theme_bw() +
    labs(x="Axis.1 (40.04%)", y="Axis.2 (8.57%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) # exported // text size was changed for e-poster

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_2.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


(p_pcoa_id_2 <- ggplot(pcoa_id_axes, aes(x=Axis.1, y=Axis.2, shape = filtration_2)) +
    geom_point(size = 4)+
    scale_shape_discrete(name="Filtration", solid=F) +
    theme_bw() +
    labs(x="Axis.1 (40.04%)", y="Axis.2 (8.57%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right")) # exported // text size was changed for e-poster

#// ggsave("out_fig/filtration_id_no_singleton_pcoa_3.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)

# *** level7_id_sample + jaccard ----
level7_id_sample <- level7_id[,1:93]
level7_id_sample <- Filter(function(x)!all(is.na(x)), level7_id_sample)


dist_id_sample <- vegdist(t(level7_id_sample),  method = "jaccard")


pcoa_id_sample <- pcoa(dist_id_sample, correction="cailliez")

pcoa_id_sample$values[c(1,2),3]  # get Rel_corr_eig

# Extract the plot scores from first two PCoA axes:
pcoa_id_sample_axes <- as.data.frame(pcoa_id_sample$vectors[,c(1,2)])
pcoa_id_sample_axes$sample_id <- rownames(pcoa_id_sample_axes)

pcoa_id_sample_axes$sample_id <- sapply(str_split(pcoa_id_sample_axes$sample_id, "[.]", n = 2), `[`, 1)

pcoa_id_sample_axes$category <- substring(pcoa_id_sample_axes$sample_id,1,2)
pcoa_id_sample_axes$category <- factor(pcoa_id_sample_axes$category, levels = c( "IN", "IP", "NN", "NP", "AG"))

pcoa_id_sample_axes$filtration <- ifelse(grepl("\\d", pcoa_id_sample_axes$sample_id),
                                  "Filter",
                                  ifelse(nchar(pcoa_id_sample_axes$sample_id) == 2, "N", substring(pcoa_id_sample_axes$sample_id,3,3)))

pcoa_id_sample_axes$filtration <- factor(pcoa_id_sample_axes$filtration, levels = c("F", "N", "Filter"))



pcoa_id_sample_axes$filtration_2 <- ifelse(nchar(pcoa_id_sample_axes$sample_id) == 2, 
                                    "N", 
                                    substring(pcoa_id_sample_axes$sample_id,3,nchar(pcoa_id_sample_axes$sample_id)))
pcoa_id_sample_axes$filtration_2 <- factor(pcoa_id_sample_axes$filtration_2, levels = c("F", "N", "F100", "F80", "F41", "F5"))

(p_pcoa_id_sample <- ggplot(pcoa_id_sample_axes, aes(x=Axis.1, y=Axis.2, color = category, shape = filtration)) +
    geom_point(size = 4, alpha = 0.7)+
    scale_colour_manual(name="Category", values = c("#F8766D", "#B79F00" ,"#00BA38", "#00BFC4", "#671bd1")) +  
    scale_shape_discrete(name="Filtration") +
    theme_bw() +
    labs(x="Axis.1 (22.43%)", y="Axis.2 (13.80%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right") +
    guides(colour = guide_legend(order = 1), 
           shape = guide_legend(order = 2)))

#// ggsave("out_fig/filtration_id_sample_no_singleton_pcoa_2.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)


# get color
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
# or change "div" to "qual"
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, 4)
show_col(col)

(p_pcoa_id_sample_2 <- ggplot(pcoa_id_sample_axes, aes(x=Axis.1, y=Axis.2, shape = filtration_2)) +
    geom_point(size = 4, aes(color = filtration_2), stroke = 0.8)+
    scale_shape_discrete(name="", solid = F, 
                         breaks=c("F", "N"  ,  "F100" ,"F80",  "F41" , "F5" ),
                         labels=c("Filtration+", "Filtration-"  , "Filter: 100 \u03BCm","Filter: 80 \u03BCm",  "Filter: 41 \u03BCm" , "Filter: 5 \u03BCm" )) +
    scale_color_manual(name="",
                       breaks=c("F", "N"  ,  "F100" ,"F80",  "F41" , "F5" ),
                       labels=c("Filtration+", "Filtration-"  , "Filter: 100 \u03BCm","Filter: 80 \u03BCm",  "Filter: 41 \u03BCm" , "Filter: 5 \u03BCm" ),
                      # values = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF")) + 
                      values = c("#0073C2FF", "#EFC000FF", "#A50026", "#8073AC", "#4D9221" ,"#C51B7D")) + 
    theme_bw() +
    labs(x="Axis.1 (22.43%)", y="Axis.2 (13.80%)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "right"))

#// ggsave("out_fig/filtration_id_sample_no_singleton_pcoa_3.pdf", width =  6.5, height = 5, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/d_filtration_id_sample_no_singleton_pcoa_3.pdf", plot=p_pcoa_id_sample_2, device=cairo_pdf, width = 7, height = 5, units = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/d_filtration_id_sample_no_singleton_pcoa_3_color.pdf", plot=p_pcoa_id_sample_2, device=cairo_pdf, width = 7, height = 5, units = "in", dpi = 300)

# * metaxa - level7 - effect on specific taxa ----------------------------------------------------------------------------------
# ** relative abundance (bar + heatmap) ----------------------------------------------------------------------------------------
level7_fi_fil <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_fi_fil <- column_to_rownames(level7_fi_fil, var = "lowest_level")

level7_fi_fil <- as.data.frame(t(level7_fi_fil))
level7_fi_fil <- rownames_to_column(level7_fi_fil, var = "sample")

level7_fi_fil$sample_2 <- sapply(str_split(level7_fi_fil$sample, "[.]", n = 2), `[`, 1)

level7_fi_fil <- level7_fi_fil %>%
  filter(sample_2 %in% c('NNF',
                         'NNN',
                         'NPF',
                         'NPN',
                         'INF',
                         'INN',
                         'IPF',
                         'IPN',
                         'SNF',
                         'SNN',
                         'SPF',
                         'SPN'))

level7_fi_fil <- select(level7_fi_fil, -sample_2)

level7_fi_fil <- column_to_rownames(level7_fi_fil,var = "sample")

level7_fi_fil <- as.data.frame(t(level7_fi_fil))

level7_fi_fil <- Filter(function(x)!all(is.na(x)), level7_fi_fil)  # delete the columns with all as NA

# order the taxa according to average abundance
level7_fi_fil$avg <- apply(level7_fi_fil, 1, mean)

level7_fi_fil <- level7_fi_fil[order(level7_fi_fil$avg, decreasing = T),]

apply(level7_fi_fil, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_fi_fil$avg != 0))  # 40
length(which(level7_fi_fil$avg == 0))  # 41
length(which(level7_fi_fil$avg > 0.001))  # 29
length(which(level7_fi_fil$avg > 0.01))  # 11

# keep only the taxa with relative abundance larger than 1%; the rest grouped to "Others"
level7_fi_fil_mod <- level7_fi_fil[which(level7_fi_fil$avg>0.01), ]

level7_fi_fil_mod = rbind(level7_fi_fil_mod,
              1-apply(level7_fi_fil_mod, 2, sum))
rownames(level7_fi_fil_mod)[nrow(level7_fi_fil_mod)] = c("Others")

# delete avg column
level7_fi_fil_mod <- select(level7_fi_fil_mod, -avg)

# melt
level7_fi_fil_melt <- melt(t(level7_fi_fil_mod))

# add category and filtration columns for facet_wrap()
level7_fi_fil_melt$category <- as.factor(substring(level7_fi_fil_melt$Var1,1,2))
level7_fi_fil_melt$filtration <-  as.factor(substring(level7_fi_fil_melt$Var1,3,nchar(as.character(level7_fi_fil_melt$Var1))))
level7_fi_fil_melt$filtration_2 <-  as.factor(substring(level7_fi_fil_melt$Var1,3,3))

level7_fi_fil_melt$is_cat <- as.factor(substring(level7_fi_fil_melt$category,1,1))
level7_fi_fil_melt$pma <- as.factor(substring(level7_fi_fil_melt$category,2,2))




# barplot, with expanded color palette
colourCount = length(unique(level7_fi_fil_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_filtration <- ggplot(level7_fi_fil_melt, aes(x = filtration, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Filtration", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE)) +
  facet_grid(labeller = labeller(is_cat = c(I = "Sample + internal standard", N = "Sample", S = "External standard"), pma = c(N = "PMA-", P = "PMA+")), rows = vars(is_cat), cols = vars(pma))

p_filtration

#// ggsave("out_fig/species_abundance_filtration_no_singleton.pdf", width =  9, height = 9, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/species_abundance_filtration_no_singleton.pdf", width =  8, height = 8, unit = "in", dpi = 300)


# heatmap
p_filtration_heatmap <- ggplot(level7_fi_fil_melt, aes(x = filtration, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "", x = "Filtration", y = "Taxa") +
  scale_fill_gradient(name = "",
                      low = "white",
                      high = "#542788") +
  scale_y_discrete(limits = rev(levels(level7_fi_fil_melt$Var2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, face = "plain"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + 
  # guides(fill=guide_legend(nrow =3,byrow=TRUE)) +
  facet_wrap(
    nrow = 3,
    ncol = 2,
    ~category)

p_filtration_heatmap
#// ggsave("out_fig/heatmap_abundance_filtration_no_singleton.pdf", width =  9, height = 9, unit = "in", dpi = 300)


# ** t test
# ttest.filtration <- level7_fi_fil_melt %>%
#   group_by(Var2, category) %>%
#   t_test(value ~ filtration_2,  paired = TRUE) %>%
#   adjust_pvalue(method = "BH") %>%
#   add_significance()
#// ERROR: not enough 'x' observations

# ** absolute abundance (heatmap) ----------------------------------------------------------------------------------------
qpcr_2 <- read.xlsx("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/qpcr/cdc_qpcr_sample.xlsx",
                    sheet = 1, startRow = 1, colNames = TRUE,
                    rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                    skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                    check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

# rename col name
qpcr_2 <-dplyr::rename(qpcr_2, "copy_number_per_250uL_sample" = "copy_number_per_250uL_sample.(account.for.concentration)")
# change CDC to AG
qpcr_2$sample[qpcr_2$sample == "CDC_1"] <- "AG_1"
qpcr_2$sample[qpcr_2$sample == "CDC_2"] <- "AG_2"
qpcr_2$sample[qpcr_2$sample == "CDC_3"] <- "AG_3"

qpcr_2$sample_2 <- sapply(strsplit(qpcr_2$sample, "_"), `[`, 1)

qpcr_2 <- qpcr_2 %>%
  select(sample, copy_number_per_250uL_sample, sample_2) %>%
  filter(sample_2 %in% c('NNF',
                         'NNN',
                         'NPF',
                         'NPN',
                         'INF',
                         'INN',
                         'IPF',
                         'IPN',
                         'SNF',
                         'SNN',
                         'SPF',
                         'SPN'))

qpcr_2$sample <- gsub("_",".",qpcr_2$sample)

# calculate absolute abundance
level7_fi_fil_mod_abs <- level7_fi_fil_mod

for(ii in 1:ncol(level7_fi_fil_mod_abs)) {level7_fi_fil_mod_abs[ ,ii] <- level7_fi_fil_mod_abs[ ,ii] * qpcr_2[which(qpcr_2$sample ==  colnames(level7_fi_fil_mod_abs)[ii] ) ,2] }

apply(level7_fi_fil_mod_abs , 2, sum)  # check


# melt
level7_fi_fil_melt_abs <- melt(t(level7_fi_fil_mod_abs))

# add category and filtration columns for facet_wrap()
level7_fi_fil_melt_abs$category <- as.factor(substring(level7_fi_fil_melt_abs$Var1,1,2))
level7_fi_fil_melt_abs$filtration <-  as.factor(substring(level7_fi_fil_melt_abs$Var1,3,nchar(as.character(level7_fi_fil_melt_abs$Var1))))
level7_fi_fil_melt_abs$filtration_2 <-  as.factor(substring(level7_fi_fil_melt_abs$Var1,3,3))

# heatmap
p_filtration_abs_heatmap <- ggplot(level7_fi_fil_melt_abs, aes(x = filtration, y = Var2, fill = value)) +
  geom_tile() +
  labs(title = "", x = "Filtration", y = "Taxa") +
  scale_fill_gradient(name = "",
                      low = "white",
                      high = "#542788") +
  scale_y_discrete(limits = rev(levels(level7_fi_fil_melt_abs$Var2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, face = "plain", angle = 330, hjust = 0),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + 
  # guides(fill=guide_legend(nrow =3,byrow=TRUE)) +
  facet_wrap(
    nrow = 3,
    ncol = 2,
    ~category)

p_filtration_abs_heatmap
#// ggsave("out_fig/heatmap_abs_abundance_filtration_no_singleton.pdf", width =  9, height = 9, unit = "in", dpi = 300)



### Filtration + locally rare taxa (level7_id) ---------------------------------------------------------------------------------
# locally rare OTUs were defined as having an abundance of < 0.01% (0.0001) within a sample
# the number of rare taxa for each sample
#// write.xlsx(level7_id, file = "out_metaxa_level7_id.xlsx")
ls_rare_no <- as.data.frame(colnames(level7_id)[1:ncol(level7_id)])
ls_rare_no <- rename(ls_rare_no, "sample" = "colnames(level7_id)[1:ncol(level7_id)]")

for (ii in 1:ncol(level7_id)){
  ls_rare_no$abun_no[ii] <- sum(level7_id[ ,ii] >= 0.01 , na.rm=TRUE)
  ls_rare_no$moderate_no[ii] <- sum(level7_id[ ,ii] < 0.01 & level7_id[ ,ii] >= 0.001, na.rm=TRUE)
  ls_rare_no$rare_1_no[ii] <- sum(level7_id[ ,ii] < 0.001 & level7_id[ ,ii] > 0, na.rm=TRUE)
  ls_rare_no$rare_2_no[ii] <- sum(level7_id[ ,ii] < 0.0001 & level7_id[ ,ii] > 0, na.rm=TRUE)
  ls_rare_no$no_zero_no[ii] <- sum(level7_id[ ,ii] != 0, na.rm=TRUE)
  }
  

which(ls_rare_no$abun_no + ls_rare_no$moderate_no != ls_rare_no$no_zero_no)  # integer(0)
# all abundant or moderate, no rare taxa

ls_rare_no$sample_2 <- sapply(str_split(ls_rare_no$sample, "[.]", n = 2), `[`, 1)

## t-test for total number of taxa
ls_total_taxa_no_fil <- ls_rare_no %>%
  filter(sample_2 %in% c('NNF',
                         'NNN',
                         'NPF',
                         'NPN',
                         'INF',
                         'INN',
                         'IPF',
                         'IPN',
                         'SNF',
                         'SNN',
                         'SPF',
                         'SPN'))

ls_total_taxa_no_fil$category <- substring(ls_total_taxa_no_fil$sample_2,1,2)
ls_total_taxa_no_fil$filtration <- substring(ls_total_taxa_no_fil$sample_2,3,nchar(as.character(ls_total_taxa_no_fil$sample_2)))

ttest.ls_taxa_no_fil <- ls_total_taxa_no_fil %>%
  group_by(category) %>%
  t_test(no_zero_no ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


# make plots
ls_rare_no_fil <- ls_rare_no %>%
  filter(sample_2 %in% c('NNF',
                         'NNN',
                         'NPF',
                         'NPN',
                         'INF',
                         'INN',
                         'IPF',
                         'IPN',
                         'SNF',
                         'SNN',
                         'SPF',
                         'SPN'))


ls_rare_no_fil <- ls_rare_no_fil %>%
  select(sample, abun_no, moderate_no) %>%
  pivot_longer(!sample, names_to = "abundance_category", values_to = "taxa_no")

# add category and filtration columns for facet_wrap()
ls_rare_no_fil$category <- substring(ls_rare_no_fil$sample,1,1)
ls_rare_no_fil$filtration <- substring(ls_rare_no_fil$sample,3,nchar(as.character(ls_rare_no_fil$sample)))
ls_rare_no_fil$filtration_2 <- substring(ls_rare_no_fil$sample,3,3)
ls_rare_no_fil$pma_group <- substring(ls_rare_no_fil$sample,2,2)

ls_rare_no_fil <- ls_rare_no_fil %>% mutate_if(is.character, as.factor)


# stacked barplot
p_taxa_no_fil <- ggplot(ls_rare_no_fil, aes(x = filtration, y = taxa_no)) +
    geom_bar(stat = "identity", colour="black", aes(fill = abundance_category)) +
    labs(title = "", x = "Filtration", y = "Number of taxa") +
    theme_bw() + 
    scale_fill_jco(breaks = c("abun_no", "moderate_no"), labels = c("Locally abundant taxa", "Locally moderate taxa")) +
    theme(axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
          plot.title = element_text(size = 15, hjust = 0, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          strip.text = element_text(size=12)) +
  facet_grid(labeller = labeller(category = c(I = "Sample + internal standard", N = "Sample", S = "External standard"), pma_group = c(N = "PMA-", P = "PMA+")), rows = vars(category), cols = vars(pma_group)) 



p_taxa_no_fil
#// ggsave("out_fig/taxa_no_filtration_bar.pdf", width =  6, height = 6, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/taxa_no_filtration_bar.pdf", width =  6, height = 8, unit = "in", dpi = 300)


## plot: p_abun_no_fil
# t test
ttest.ls_abun_no_fil <- subset(ls_rare_no_fil, abundance_category=="abun_no") %>%
  group_by(category, pma_group) %>%
  t_test(taxa_no ~ filtration_2,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "filtration_2")

# boxplot
p_abun_no_fil <- ggboxplot(subset(ls_rare_no_fil, abundance_category=="abun_no"), x = "filtration_2", y = "taxa_no", facet.by = c("category","pma_group"), panel.labs = list(category = c("Sample + internal standard", "Sample", "External standard"), pma_group = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "#0073C2FF", add = "jitter", legend = "none", xlab = "Filtration", ylab = "Number of locally abundant taxa", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))

p_abun_no_fil + stat_pvalue_manual(ttest.ls_abun_no_fil, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")
#// ggsave("out_fig/abun_no_filtration_box.pdf", width =  6, height = 6, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/abun_no_filtration_box.pdf", width =  6, height = 8, unit = "in", dpi = 300)


## plot: p_moderate_no_fil
# t test
ttest.ls_moderate_no_fil <- subset(ls_rare_no_fil, abundance_category=="moderate_no") %>%
  group_by(category, pma_group) %>%
  t_test(taxa_no ~ filtration_2,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "filtration_2")

p_moderate_no_fil <- ggboxplot(subset(ls_rare_no_fil, abundance_category=="moderate_no"), x = "filtration_2", y = "taxa_no", facet.by = c("category","pma_group"), panel.labs = list(category = c("Sample + internal standard", "Sample", "External standard"), pma_group = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "#EFC000FF", add = "jitter", legend = "none", xlab = "Filtration", ylab = "Number of locally moderate taxa", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))

p_moderate_no_fil + stat_pvalue_manual(ttest.ls_moderate_no_fil, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")
#// ggsave("out_fig/moderate_no_filtration_box.pdf", width =  6, height = 6, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/moderate_no_filtration_box.pdf", width =  6, height = 8, unit = "in", dpi = 300)




### Profile change (with filters) ----------------------------------------------------------------------------------------------
# * NN ----
level7_nn <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_nn <- column_to_rownames(level7_nn, var = "lowest_level")

level7_nn<- as.data.frame(t(level7_nn))
level7_nn <- rownames_to_column(level7_nn, var = "sample")

level7_nn$sample_2 <- substring(level7_nn$sample, 1,2)

level7_nn <- level7_nn %>%
  filter(sample_2 == 'NN')

level7_nn <- select(level7_nn, -sample_2)

level7_nn <- column_to_rownames(level7_nn,var = "sample")

level7_nn <- as.data.frame(t(level7_nn))

# order the taxa according to average abundance
level7_nn$avg <- apply(level7_nn, 1, mean)

level7_nn <- level7_nn[order(level7_nn$avg, decreasing = T),]

apply(level7_nn, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_nn$avg > 0.001))  # 12 (selected)
length(which(level7_nn$avg > 0.01))  # 6

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_nn_mod <- level7_nn[which(level7_nn$avg > 0.001), ]

level7_nn_mod = rbind(level7_nn_mod,
                       1-apply(level7_nn_mod, 2, sum))
rownames(level7_nn_mod)[nrow(level7_nn_mod)] = c("Others")

# delete avg column
level7_nn_mod <- select(level7_nn_mod, -avg)

# melt
level7_nn_melt <- melt(t(level7_nn_mod))

# add columns for facet_wrap()
level7_nn_melt$sample_2 <- sapply(str_split(level7_nn_melt$Var1, "[.]", n = 2), `[`, 1)

# plot, with expanded color palette
colourCount = length(unique(level7_nn_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

(p_nn <- ggplot(subset(level7_nn_melt, sample_2=="NN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")))

#// ggsave("out_fig/filter-detail/nn.pdf", width =  4, height = 4, unit = "in", dpi = 300)

filter_levels <- c('NNF100.1',
                   'NNF100.2',
                   'NNF100.swab',
                   'NNF80.1',
                   'NNF80.2',
                   'NNF41.1',
                   'NNF41.2',
                   'NNF41.3',
                   'NNF5.1',
                   'NNF5.2',
                   'NNF5.3')


p_nn_filter <- ggplot(subset(level7_nn_melt, grepl("\\d", level7_nn_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

#// ggsave("out_fig/filter-detail/nn_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_nnf <- ggplot(subset(level7_nn_melt, sample_2=="NNF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))
#// ggsave("out_fig/filter-detail/nnf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_nnn <- ggplot(subset(level7_nn_melt, sample_2=="NNN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))
#// ggsave("out_fig/filter-detail/nnn.pdf", width =  4, height = 4, unit = "in", dpi = 300)

# * NP ----
level7_np <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_np <- column_to_rownames(level7_np, var = "lowest_level")

level7_np<- as.data.frame(t(level7_np))
level7_np <- rownames_to_column(level7_np, var = "sample")

level7_np$sample_2 <- substring(level7_np$sample, 1,2)

level7_np <- level7_np %>%
  filter(sample_2 == 'NP')


level7_np <- select(level7_np, -sample_2)

level7_np <- column_to_rownames(level7_np,var = "sample")

level7_np <- as.data.frame(t(level7_np))

# order the taxa according to average abundance
level7_np$avg <- apply(level7_np, 1, mean)

level7_np <- level7_np[order(level7_np$avg, decreasing = T),]

apply(level7_np, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_np$avg > 0.001))  # 9 (selected)
length(which(level7_np$avg > 0.01))  # 6

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_np_mod <- level7_np[which(level7_np$avg>0.001), ]

level7_np_mod = rbind(level7_np_mod,
                      1-apply(level7_np_mod, 2, sum))
rownames(level7_np_mod)[nrow(level7_np_mod)] = c("Others")

# delete avg column
level7_np_mod <- select(level7_np_mod, -avg)

# melt
level7_np_melt <- melt(t(level7_np_mod))

# add columns for facet_wrap()
level7_np_melt$sample_2 <- sapply(str_split(level7_np_melt$Var1, "[.]", n = 2), `[`, 1)

# plot, with expanded color palette
colourCount = length(unique(level7_np_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_np <- ggplot(subset(level7_np_melt, sample_2=="NP"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) 

ggsave("out_fig/filter-detail/np.pdf", width =  4, height = 4, unit = "in", dpi = 300)


unique(subset(level7_np_melt, grepl("\\d", level7_np_melt$sample_2))$Var1)

filter_levels <- c('NPF100.1',
                   'NPF100.2',
                   'NPF100.3',
                   'NPF100.swab',
                   'NPF80.1',
                   'NPF80.2',
                   'NPF80.3',
                   'NPF41.1',
                   'NPF41.2',
                   'NPF41.3',
                   'NPF5.1',
                   'NPF5.2',
                   'NPF5.3')


p_np_filter <- ggplot(subset(level7_np_melt, grepl("\\d", level7_np_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

ggsave("out_fig/filter-detail/np_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_npf <- ggplot(subset(level7_np_melt, sample_2=="NPF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

ggsave("out_fig/filter-detail/npf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_npn <- ggplot(subset(level7_np_melt, sample_2=="NPN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))

ggsave("out_fig/filter-detail/npn.pdf", width =  4, height = 4, unit = "in", dpi = 300)






# * IN ----
level7_in <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_in <- column_to_rownames(level7_in, var = "lowest_level")

level7_in<- as.data.frame(t(level7_in))
level7_in <- rownames_to_column(level7_in, var = "sample")

level7_in$sample_2 <- substring(level7_in$sample, 1,2)

level7_in <- level7_in %>%
  filter(sample_2 == 'IN')

level7_in <- select(level7_in, -sample_2)

level7_in <- column_to_rownames(level7_in,var = "sample")

level7_in <- level7_in[rowSums(is.na(level7_in)) != ncol(level7_in), ]  # delete rows with all data as NA

level7_in <- as.data.frame(t(level7_in))

# order the taxa according to average abundance
level7_in$avg <- apply(level7_in, 1, mean)

level7_in <- level7_in[order(level7_in$avg, decreasing = T),]

apply(level7_in, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_in$avg > 0.001))  # 23 (selected)
length(which(level7_in$avg > 0.01))  # 8

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_in_mod <- level7_in[which(level7_in$avg>0.001), ]

level7_in_mod = rbind(level7_in_mod,
                      1-apply(level7_in_mod, 2, sum))
rownames(level7_in_mod)[nrow(level7_in_mod)] = c("Others")

# delete avg column
level7_in_mod <- select(level7_in_mod, -avg)

# melt
level7_in_melt <- melt(t(level7_in_mod))

# add columns for facet_wrap()
level7_in_melt$sample_2 <- sapply(str_split(level7_in_melt$Var1, "[.]", n = 2), `[`, 1)

# plot, with expanded color palette
colourCount = length(unique(level7_in_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_in <- ggplot(subset(level7_in_melt, sample_2=="IN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(ncol =1))

ggsave("out_fig/filter-detail/in.pdf", width =  4, height = 4, unit = "in", dpi = 300)


unique(subset(level7_in_melt, grepl("\\d", level7_in_melt$sample_2))$Var1)

filter_levels <- c('INF100.1',
                   'INF100.2',
                   'INF100.3',
                   'INF100.swab',
                   'INF80.1',
                   'INF80.2',
                   'INF41.1',
                   'INF41.2',
                   'INF41.3',
                   'INF5.1',
                   'INF5.2',
                   'INF5.3')


p_in_filter <- ggplot(subset(level7_in_melt, grepl("\\d", level7_in_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

ggsave("out_fig/filter-detail/in_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_inf <- ggplot(subset(level7_in_melt, sample_2=="INF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))
ggsave("out_fig/filter-detail/inf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_inn <- ggplot(subset(level7_in_melt, sample_2=="INN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(ncol =1,bycol=TRUE))

p_inn 
ggsave("out_fig/filter-detail/inn.pdf", width =  4, height = 4, unit = "in", dpi = 300)

# * IP ----
level7_ip <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_ip <- column_to_rownames(level7_ip, var = "lowest_level")

level7_ip<- as.data.frame(t(level7_ip))
level7_ip <- rownames_to_column(level7_ip, var = "sample")

level7_ip$sample_2 <- substring(level7_ip$sample, 1,2)

level7_ip <- level7_ip %>%
  filter(sample_2 == 'IP')


level7_ip <- select(level7_ip, -sample_2)

level7_ip <- column_to_rownames(level7_ip,var = "sample")

level7_ip <- level7_ip[rowSums(is.na(level7_ip)) != ncol(level7_ip), ]  # delete rows with all data as NA

level7_ip <- as.data.frame(t(level7_ip))

# order the taxa according to average abundance
level7_ip$avg <- apply(level7_ip, 1, mean)

level7_ip <- level7_ip[order(level7_ip$avg, decreasing = T),]

apply(level7_ip, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_ip$avg > 0.001))  # 7 (selected)
length(which(level7_ip$avg > 0.01))  # 5

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_ip_mod <- level7_ip[which(level7_ip$avg>0.001), ]

level7_ip_mod = rbind(level7_ip_mod,
                      1-apply(level7_ip_mod, 2, sum))
rownames(level7_ip_mod)[nrow(level7_ip_mod)] = c("Others")

# delete avg column
level7_ip_mod <- select(level7_ip_mod, -avg)

# melt
level7_ip_melt <- melt(t(level7_ip_mod))

# add columns for facet_wrap()
level7_ip_melt$sample_2 <- sapply(str_split(level7_ip_melt$Var1, "[.]", n = 2), `[`, 1)


p_ip <- ggplot(subset(level7_ip_melt, sample_2=="IP"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_brewer(palette="Set1") +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) 

ggsave("out_fig/filter-detail/ip.pdf", width =  4, height = 4, unit = "in", dpi = 300)


unique(subset(level7_ip_melt, grepl("\\d", level7_ip_melt$sample_2))$Var1)

filter_levels <- c('IPF100.1',
                   'IPF100.2',
                   'IPF100.3',
                   'IPF100.swab',
                   'IPF80.1',
                   'IPF80.2',
                   'IPF80.3',
                   'IPF41.1',
                   'IPF41.2',
                   'IPF41.3',
                   'IPF5.1',
                   'IPF5.3')


p_ip_filter <- ggplot(subset(level7_ip_melt, grepl("\\d", level7_ip_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_brewer(palette="Set1") +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

ggsave("out_fig/filter-detail/ip_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_ipf <- ggplot(subset(level7_ip_melt, sample_2=="IPF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_brewer(palette="Set1") +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))
ggsave("out_fig/filter-detail/ipf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_ipn <- ggplot(subset(level7_ip_melt, sample_2=="IPN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_brewer(palette="Set1") +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))
ggsave("out_fig/filter-detail/ipn.pdf", width =  4, height = 4, unit = "in", dpi = 300)


# * SN ----
level7_sn <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_sn <- column_to_rownames(level7_sn, var = "lowest_level")

level7_sn<- as.data.frame(t(level7_sn))
level7_sn <- rownames_to_column(level7_sn, var = "sample")

level7_sn$sample_2 <- substring(level7_sn$sample, 1,2)

level7_sn <- level7_sn %>%
  filter(sample_2 == 'SN')


level7_sn <- select(level7_sn, -sample_2)

level7_sn <- column_to_rownames(level7_sn,var = "sample")

level7_sn <- level7_sn[rowSums(is.na(level7_sn)) != ncol(level7_sn), ]  # delete rows with all data as NA

level7_sn <- as.data.frame(t(level7_sn))

# order the taxa according to average abundance
level7_sn$avg <- apply(level7_sn, 1, mean)

level7_sn <- level7_sn[order(level7_sn$avg, decreasing = T),]

apply(level7_sn, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_sn$avg > 0.001))  # 24 (selected)
length(which(level7_sn$avg > 0.01))  # 17

# keep only the taxa with relative abundance larger than 0.1%; the rest grouped to "Others"
level7_sn_mod <- level7_sn[which(level7_sn$avg>0.001), ]

level7_sn_mod = rbind(level7_sn_mod,
                      1-apply(level7_sn_mod, 2, sum))
rownames(level7_sn_mod)[nrow(level7_sn_mod)] = c("Others")

# delete avg column
level7_sn_mod <- select(level7_sn_mod, -avg)

# melt
level7_sn_melt <- melt(t(level7_sn_mod))

# add columns for facet_wrap()
level7_sn_melt$sample_2 <- sapply(str_split(level7_sn_melt$Var1, "[.]", n = 2), `[`, 1)

# plot, with expanded color palette
colourCount = length(unique(level7_sn_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_sn <- ggplot(subset(level7_sn_melt, sample_2=="SN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(ncol =1))

ggsave("out_fig/filter-detail/sn.pdf", width =  4, height = 4, unit = "in", dpi = 300)


unique(subset(level7_sn_melt, grepl("\\d", level7_sn_melt$sample_2))$Var1)

filter_levels <- c('SNF100.2',
                   'SNF100.3',
                   'SNF80.2',
                   'SNF41.1',
                   'SNF41.2',
                   'SNF41.3',
                   'SNF5.1',
                   'SNF5.2')


p_sn_filter <- ggplot(subset(level7_sn_melt, grepl("\\d", level7_sn_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))

ggsave("out_fig/filter-detail/sn_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_snf <- ggplot(subset(level7_sn_melt, sample_2=="SNF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(nrow =4,byrow=TRUE))
ggsave("out_fig/filter-detail/snf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_snn <- ggplot(subset(level7_sn_melt, sample_2=="SNN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))
ggsave("out_fig/filter-detail/snn.pdf", width =  4, height = 4, unit = "in", dpi = 300)



# * SP ----
level7_sp <- level7_fi[, c(2, 10:ncol(level7_fi)) ]

level7_sp <- column_to_rownames(level7_sp, var = "lowest_level")

level7_sp<- as.data.frame(t(level7_sp))
level7_sp <- rownames_to_column(level7_sp, var = "sample")

level7_sp$sample_2 <- substring(level7_sp$sample, 1,2)

level7_sp <- level7_sp %>%
  filter(sample_2 == 'SP')


level7_sp <- select(level7_sp, -sample_2)

level7_sp <- column_to_rownames(level7_sp,var = "sample")

level7_sp <- level7_sp[rowSums(is.na(level7_sp)) != ncol(level7_sp), ]  # delete rows with all data as NA

level7_sp <- as.data.frame(t(level7_sp))

# order the taxa according to average abundance
level7_sp$avg <- apply(level7_sp, 1, mean)

level7_sp <- level7_sp[order(level7_sp$avg, decreasing = T),]

apply(level7_sp, 2, sum)  # check if the sum of abundance = 1

# check the number of taxa
length(which(level7_sp$avg > 0.001))  # 3 (selected)
length(which(level7_sp$avg > 0.01))  # 3

# delete avg column
level7_sp <- select(level7_sp, -avg)

# delete rows all equal to 0.0
level7_sp <- level7_sp[rowSums(level7_sp) != 0, ]  

# melt
level7_sp_melt <- melt(t(level7_sp))

# add columns for facet_wrap()
level7_sp_melt$sample_2 <- sapply(str_split(level7_sp_melt$Var1, "[.]", n = 2), `[`, 1)

# plot, with expanded color palette
colourCount = length(unique(level7_sp_melt$Var2))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # set1 is the best among Set1, accent, paired, set3, and spectral

p_sp <- ggplot(subset(level7_sp_melt, sample_2=="SP"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), units = "cm")) + guides(fill=guide_legend(ncol =1))

ggsave("out_fig/filter-detail/sp.pdf", width =  4, height = 4, unit = "in", dpi = 300)


unique(subset(level7_sp_melt, grepl("\\d", level7_sp_melt$sample_2))$Var1)



p_sp_filter <- ggplot(subset(level7_sp_melt, grepl("\\d", level7_sp_melt$sample_2)), aes(x = factor(Var1,levels = filter_levels), y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 330, hjust = 0),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "right",
        plot.margin = unit(c(1,1,1,1), units = "cm")) 

ggsave("out_fig/filter-detail/sp_filter.pdf", width =  8, height = 4, unit = "in", dpi = 300)

p_spf <- ggplot(subset(level7_sp_melt, sample_2=="SPF"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))
ggsave("out_fig/filter-detail/spf.pdf", width =  4, height = 4, unit = "in", dpi = 300)

p_spn <- ggplot(subset(level7_sp_melt, sample_2=="SPN"), aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Sample", y = "Relative Abundance") +
  theme_bw() +  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "cm"))
ggsave("out_fig/filter-detail/spn.pdf", width =  4, height = 4, unit = "in", dpi = 300)
