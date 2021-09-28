# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 
# purpose: 



### Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(openxlsx) # handle xlsx files
library(tidyr)
library(rstatix)
library(ggpubr)
library(stringr)
library(tibble)

### Functions ----




### Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/knead/filtration/")



### Import data ----
# check the data type of columns of interest by str()
reads_homo <- read.csv("human_contaminant_calculation.csv")[,3:5]
reads_homo_online <- read.xlsx("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/knead/filtration/knead_log_other_studies/human_contaminant_calculation_online.xlsx")[,c(1,2,4)]

silva132_level1 <- read.xlsx("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/compo_analysis/input/merged_SILVA132_level1.xlsx", sheet = "merged_SILVA132_level1")

# remove singleton
# for (ii in 2:ncol(silva132_level1)) {silva132_level1[which(silva132_level1[,ii]<2),ii] <- 0}

rownames(silva132_level1) <- silva132_level1$Taxa
silva132_level1 <- silva132_level1[,-1]
silva132_level1 <- as.data.frame(t(silva132_level1))
silva132_level1 <- rownames_to_column(silva132_level1, var = "sample")

###  my samples----
reads_homo <- pivot_wider(reads_homo, names_from = category, values_from = reads_no)

reads_homo <- mutate(reads_homo, reads_eukaryote = before - after, percent_eukaryote = reads_eukaryote/before*100)  # calculation

reads_homo_sample <- filter(reads_homo, substring(sample_id_short_dot, 1, 1) != "S")  # exclude all standard only samples

## differences between filtered and not filtered samples
reads_homo_2 <- reads_homo[,c(1,4,5)]

reads_homo_2$sample_id_short_dot <- as.character(reads_homo_2$sample_id_short_dot)
reads_homo_2$sample_2 <- sapply(str_split(reads_homo_2$sample_id_short_dot, "[.]", n = 2), `[`, 1)

reads_homo_3 <- reads_homo_2 %>%
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

reads_homo_3$category <- substring(reads_homo_3$sample_2,1,2)
reads_homo_3$filtration <- substring(reads_homo_3$sample_2,3,3)

reads_homo_3$category <- as.factor(reads_homo_3$category)
reads_homo_3$filtration <- as.factor(reads_homo_3$filtration)
reads_homo_3$filtration <- factor(reads_homo_3$filtration, levels = c("F","N"))

ttest.reads_homo_3 <- reads_homo_3 %>%
  group_by(is_cat, pma) %>%
  t_test(percent_eukaryote ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

reads_homo_3$is_cat <- substring(reads_homo_3$category, 1, 1) 
reads_homo_3$pma <- substring(reads_homo_3$category, 2, 2)
reads_homo_3 <- reads_homo_3 %>% mutate_if(is.character, as.factor)

p_reads_filtration <- ggboxplot(reads_homo_3, x = "filtration", y = "percent_eukaryote", facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "filtration", palette = "npg", legend = "none", xlab = "Filtration", ylab = "Percentage of human-associated reads (%)", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))


ttest.reads_homo_3 <- ttest.reads_homo_3 %>% add_xy_position(x = "filtration")

p_reads_filtration + 
  stat_pvalue_manual(ttest.reads_homo_3, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")

#// ggsave("fig/reads_filtration.pdf", width =  5, height = 8, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/reads_filtration.pdf", width =  5, height = 8, unit = "in", dpi = 300)

## all filtered samples VS filters VS non-filtered samples
#// write.xlsx(reads_homo, file = "out_reads_homo.xlsx")

### my samples: after kneaddata; taxa ----
silva132_level1$sample_2 <- sapply(str_split(silva132_level1$sample, "[.]", n = 2), `[`, 1)

silva132_level1_2 <- silva132_level1 %>%
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


silva132_level1_2$reads_sum = rowSums(silva132_level1_2[, 2:8])
silva132_level1_2$Eukaryota_percentage = silva132_level1_2$Eukaryota / silva132_level1_2$reads_sum *100

silva132_level1_2$category <- substring(silva132_level1_2$sample_2,1,2)
silva132_level1_2$filtration <- substring(silva132_level1_2$sample_2,3,3)

silva132_level1_2$category <- as.factor(silva132_level1_2$category)
silva132_level1_2$filtration <- as.factor(silva132_level1_2$filtration)
silva132_level1_2$filtration <- factor(silva132_level1_2$filtration, levels = c("F","N"))

silva132_level1_2$is_cat <- substring(silva132_level1_2$category, 1, 1) 
silva132_level1_2$pma <- substring(silva132_level1_2$category, 2, 2)
silva132_level1_2 <- silva132_level1_2 %>% mutate_if(is.character, as.factor)

ttest.silva132_level1_2 <- silva132_level1_2 %>%
  group_by(is_cat, pma) %>%
  t_test(Eukaryota_percentage ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


p_reads_filtration_taxa <- ggboxplot(silva132_level1_2, x = "filtration", y = "Eukaryota_percentage", facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "filtration", palette = "npg", legend = "none", xlab = "Filtration", ylab = "Percentage of eukaryotes unassociated to humans (%)", ggtheme = theme_pubr(border = TRUE)) + 
  yscale(.scale="log2", .format = FALSE) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))


p_reads_filtration_taxa + 
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")

#// ggsave("fig/reads_filtration_non_human.pdf", width =  5, height = 8, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/reads_filtration_non_human.pdf", width =  5, height = 8, unit = "in", dpi = 300)


### add up percentage ----
percentage_sum <- reads_homo_3[, c(1,3:6)]
percentage_sum <-dplyr::rename(percentage_sum, sample = sample_id_short_dot)
percentage_sum$human <- "Y"
silva132_level1_2$human <- "N"
silva132_level1_3 <- silva132_level1_2[ ,c(1,12,10,13,14,17)]
colnames(silva132_level1_3) <- colnames(percentage_sum)
percentage_sum <- rbind(percentage_sum,silva132_level1_3)

percentage_sum_2 <- percentage_sum %>%
  group_by(sample, category, filtration) %>%
  summarise(percent_eukaryote_sum = sum(percent_eukaryote))

percentage_sum_3$is_cat <- substring(percentage_sum_3$category, 1, 1) 
percentage_sum_3$pma <- substring(percentage_sum_3$category, 2, 2)
percentage_sum_3 <- percentage_sum_3 %>% mutate_if(is.character, as.factor)


ttest.percentage_sum <- percentage_sum_2 %>%
  group_by(category) %>%
  t_test(percent_eukaryote_sum ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

percentage_sum_3 <- percentage_sum %>%
  group_by(category, filtration, human) %>%
  summarise(percent_eukaryote_mean = mean(percent_eukaryote))

p_reads_filtration_sum <- ggbarplot(percentage_sum_3, x = "filtration", y = "percent_eukaryote_mean", facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "human", palette = "npg", xlab = "Filtration", ylab = "Mean percentage of eukaryotes (%)", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))

p_reads_filtration_sum
# ttest.percentage_sum <- ttest.percentage_sum  %>% add_xy_position(x = "filtration")
# 
# p_reads_filtration_sum + 
#   stat_pvalue_manual(ttest.percentage_sum, label = "p.adj.signif")

ggpar(p_reads_filtration_sum,
      legend = "bottom",
      legend.title = "Human-associated eukaryotes",
      font.legend = c(12),
      font.tickslab = c(12),
      font.x = c(14),
      font.y = c(14),
      ylim = c(0,6.5)
      )

#// ggsave("fig/reads_filtration_sum.pdf", width =  5, height = 8, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/reads_filtration_sum.pdf", width =  5, height = 8, unit = "in", dpi = 300)



### online samples ----
reads_homo_online <- reads_homo_online[!duplicated(reads_homo_online), ]

reads_homo_online <- pivot_wider(reads_homo_online, names_from = category, values_from = reads_no)

# two samples were not included in the text extraction step by bash, so they are manually added here
# SRR5420295 is not found in log files (data found in "reads_brooks.txt")
# ERR3928927 is missed becasue it is the last line in the list text (data found in log files)
reads_homo_online[nrow(reads_homo_online) + 1, ] = list("SRR5420295",3412582, 3381712)
reads_homo_online[nrow(reads_homo_online) + 1, ] = list("SERR3928927",35063937,35053206)

reads_homo_online <- mutate(reads_homo_online, reads_eukaryote = before - after, percent_eukaryote = reads_eukaryote/before*100)  # calculation

(p_online <- ggplot(reads_homo_online, aes(x = percent_eukaryote)) +
    geom_histogram(binwidth =  5, color = "black", fill = "gray") +
    # geom_vline(aes(xintercept = median(file$diversity)), 
              #  linetype = "dashed", size = 0.6) +
    scale_y_continuous(trans='log10') +
    labs(x="Percentage of eukaryote reads (%)", y="Count of samples", caption = "Binwidth = 5%")+
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "bottom",
          panel.grid = element_blank(),  # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"))) # Adding a 1cm margin around the plot

#// ggsave("fig/reads_homo_online_distribution.pdf", width =  7, height = 5, unit = "in", dpi = 300)
#// ggsave("fig/reads_homo_online_distribution_binwidth_5.pdf", width =  7, height = 5, unit = "in", dpi = 300)

sum(reads_homo_online$percent_eukaryote <= 1)
sum(reads_homo_online$percent_eukaryote == 1)
sum(reads_homo_online$percent_eukaryote > 1)
sum(reads_homo_online$percent_eukaryote > 5)
sum(reads_homo_online$percent_eukaryote > 10)
sum(reads_homo_online$percent_eukaryote > 25)
sum(reads_homo_online$percent_eukaryote > 50)

'''
> sum(reads_homo_online$percent_eukaryote > 1)
[1] 132
> sum(reads_homo_online$percent_eukaryote > 5)
[1] 69
> sum(reads_homo_online$percent_eukaryote > 10)
[1] 51
> sum(reads_homo_online$percent_eukaryote > 25)
[1] 36
> sum(reads_homo_online$percent_eukaryote > 50)
[1] 20
'''
