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
library(Cairo)



# Functions ----




# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/qpcr")



# Import data ----
# check the data type of columns of interest by str()
qpcr <- read.xlsx("merged_cdc_qpcr_result.xlsx",
          sheet = 1, startRow = 1, colNames = TRUE,
          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
          skipEmptyCols = TRUE, rows = NULL, cols = NULL,
          check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

### manipulate the entire imported dataset ----
# change CDC to AG
qpcr$sample[qpcr$sample == "CDC_1"] <- "AG_1"
qpcr$sample[qpcr$sample == "CDC_2"] <- "AG_2"
qpcr$sample[qpcr$sample == "CDC_3"] <- "AG_3"

qpcr <- separate(qpcr, "sample", c("sample","replicate"), sep = "_")  # separate the replicate info

qpcr$sample <- ifelse((qpcr$replicate == "swab"), paste(qpcr$sample,qpcr$replicate,sep="_"), qpcr$sample)  # make "swab" samples distinguishable

# change the class of columns
qpcr$sample <- as.factor(qpcr$sample)
qpcr$category <- as.factor(qpcr$category)
qpcr[,4:7] <- lapply(qpcr[,4:7], as.numeric) # BLQ was converted to NA

### make the "copy number per ng DNA" plot (w/ t test) ----
# # version 1 (abandoned)
# qpcr_cp_ng_summary <- qpcr %>%
#   group_by(sample) %>%
#   summarise(copy_number_per_ng_DNA_mean = mean(copy_number_per_ng_DNA),
#             copy_number_per_ng_DNA_se = sd(copy_number_per_ng_DNA) / sqrt(length(copy_number_per_ng_DNA))) %>%
#   filter(sample %in% c('NNF',
#                        'NNN',
#                        'NPF',
#                        'NPN',
#                        'INF',
#                        'INN',
#                        'IPF',
#                        'IPN',
#                        'SNF',
#                        'SNN',
#                        'SPF',
#                        'SPN'))
# (p1 <- ggplot(qpcr_cp_ng_summary, aes(x=sample, y=copy_number_per_ng_DNA_mean)) +
#     geom_bar(stat = "identity", width = 0.3, fill = c('#d73027',	'#d73027',
#                                                       '#fc8d59',	'#fc8d59',
#                                                       '#fee090',	'#fee090',
#                                                       '#e0f3f8',	'#e0f3f8',
#                                                       '#91bfdb',	'#91bfdb',
#                                                       '#4575b4',	'#4575b4') ) +
#     geom_errorbar(aes(ymin=copy_number_per_ng_DNA_mean-copy_number_per_ng_DNA_se, ymax=copy_number_per_ng_DNA_mean + copy_number_per_ng_DNA_se), width=.2, position=position_dodge(0.05)) +
#     labs(title=NULL, x="Sample", y="Copy number per ng DNA") +
#     theme_bw() +
#     theme(axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 330),  # vjust/hjust is between 0 and 1 
#           axis.text.y = element_text(size = 12),
#           axis.title = element_text(size = 14, face = "plain"),  
#           plot.title = element_text(lineheight=.8, face="bold"),
#           legend.position = "none",
#           panel.grid = element_blank()))  # Removing the background grid lines      
# 
# # qpcr_no_is <- filter(qpcr, category == "no_is")
# # qpcr_with_is <- filter(qpcr, sample %in% sample_list$sample_with_is)
# # qpcr_s_only <- filter(qpcr, sample %in% sample_list$sample_s_only)
# # qpcr_neg <- filter(qpcr, sample %in% sample_list$neg_ctrl)

## t test 
qpcr_cp_ng <- qpcr %>%
  select(sample, copy_number_per_ng_DNA) %>%
  filter(sample %in% c('NNF',
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

qpcr_cp_ng$category <- substring(qpcr_cp_ng$sample,1,2)
qpcr_cp_ng$filtration <- substring(qpcr_cp_ng$sample,3,3)
qpcr_cp_ng <- qpcr_cp_ng[,-1]

# create facet category
qpcr_cp_ng$is_cat <- substring(qpcr_cp_ng$category, 1, 1) 
qpcr_cp_ng$pma <- substring(qpcr_cp_ng$category, 2, 2)
qpcr_cp_ng <- qpcr_cp_ng %>% mutate_if(is.character, as.factor)

ttest.stat <- qpcr_cp_ng %>%
  group_by(is_cat, pma) %>%
  t_test(copy_number_per_ng_DNA ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


# Create the box plot (version 1)
p_qpcr_cp_ng <- ggboxplot(qpcr_cp_ng, x = "filtration", y = "copy_number_per_ng_DNA", fill = "filtration", palette = "npg", legend = "none", xlab = "Filtration", ylab = "Copy number per ng DNA", ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(nrow = 3,
             ncol = 2,
             ~category) 

p_qpcr_cp_ng + 
  stat_pvalue_manual(ttest.stat, label = "p.adj.signif") +
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  font("xy.text", size = 14, face = "plain")

#// ggsave("fig/copy_number_ng_dna.pdf", width =  5, height = 8, unit = "in", dpi = 300)



# bar plot (version 2; ggpubr)
p_qpcr_cp_ng_bar <- ggbarplot(qpcr_cp_ng, x = "filtration", y = "copy_number_per_ng_DNA",  add = c("mean_se", "jitter"), facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "filtration", palette = "jco", legend = "none", xlab = "", ylab = "16S rRNA gene copy number per ng DNA", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))

p_qpcr_cp_ng_bar
# Add statistical test p-values
ttest.stat <- ttest.stat %>% add_xy_position(x = "filtration")


p_qpcr_cp_ng_bar + stat_pvalue_manual(ttest.stat, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")

#// ggsave("fig/copy_number_ng_dna_bar.pdf", width =  5, height = 8, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/copy_number_ng_dna.pdf", width =  6, height = 8, unit = "in", dpi = 300)

# bar plot + log transformnation (version 3; ggpubr)
qpcr_cp_ng <- qpcr_cp_ng %>%
  mutate(copy_number_per_ng_DNA_log10 = log10(copy_number_per_ng_DNA))



p_qpcr_cp_ng_bar_log10 <- ggbarplot(qpcr_cp_ng, x = "filtration", y = "copy_number_per_ng_DNA_log10",  add = c("mean_se", "jitter"), facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "filtration", palette = "jco", legend = "none", xlab = "", ylab = "16S rRNA gene copy number per ng DNA (log10 transformed)", ggtheme = theme_pubr(border = TRUE)) + scale_x_discrete(labels=c("F" = "Filtration+", "N" = "Filtration-"))

p_qpcr_cp_ng_bar_log10
# Add statistical test p-values
ttest.stat <- qpcr_cp_ng %>%
  group_by(is_cat, pma) %>%
  t_test(copy_number_per_ng_DNA_log10 ~ filtration,  paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()  %>% 
  add_xy_position(x = "filtration")


p_qpcr_cp_ng_bar_log10 + stat_pvalue_manual(ttest.stat, label = "p.adj.signif") +
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")



### biomass loss after filtration ----
## import data
qpcr_2 <- read.xlsx("cdc_qpcr_sample.xlsx",
                  sheet = 1, startRow = 1, colNames = TRUE,
                  rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                  skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                  check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

## manipulate the entire imported dataset
# rename col name
qpcr_2 <-dplyr::rename(qpcr_2, "copy_number_per_250uL_sample" = "copy_number_per_250uL_sample.(account.for.concentration)")

# change CDC to AG
qpcr_2$sample[qpcr_2$sample == "CDC_1"] <- "AG_1"
qpcr_2$sample[qpcr_2$sample == "CDC_2"] <- "AG_2"
qpcr_2$sample[qpcr_2$sample == "CDC_3"] <- "AG_3"

qpcr_2$sample_2 <- sapply(strsplit(qpcr_2$sample, "_"), `[`, 1)
qpcr_2$sample_2 <- substring(qpcr_2$sample_2, 1, 3)

qpcr_2 <- qpcr_2 %>%
  select(sample, category, category_2, copy_number_per_250uL_sample, sample_2) %>%
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

# add up all biomass retended by filters
qpcr_2_filter <- filter(qpcr_2, category_2 == "filter")

qpcr_2_filter$replicate <- sapply(strsplit(qpcr_2_filter$sample, "_"), `[`, 2)

qpcr_2_filter_sum <- qpcr_2_filter %>%
  group_by(sample_2, replicate) %>%
  summarise(copy_number_per_250uL_sample = sum(copy_number_per_250uL_sample))

# add the swab value
for (i in c(1,5,9,13) ){
qpcr_2_filter_sum[i, 3] = qpcr_2_filter_sum[i, 3] + qpcr_2_filter_sum[i+3, 3]
qpcr_2_filter_sum[i+1, 3] = qpcr_2_filter_sum[i+1, 3] + qpcr_2_filter_sum[i+3, 3]
qpcr_2_filter_sum[i+2, 3] = qpcr_2_filter_sum[i+2, 3] + qpcr_2_filter_sum[i+3, 3]
}

qpcr_2_filter_sum <- filter(qpcr_2_filter_sum, replicate != "swab")  # delete the swab rows

# organize the dataframe
qpcr_2_filter_sum$category <- substring(qpcr_2_filter_sum$sample_2, 1, 2)
qpcr_2_filter_sum$filtration <- "retention"
qpcr_2_filter_sum <- qpcr_2_filter_sum[, c(4,5,3)]

# get the liquid samples ready
qpcr_2_liquid <- filter(qpcr_2, category_2 != "filter")
qpcr_2_liquid$category <- substring(qpcr_2_liquid$sample_2,1,2)
qpcr_2_liquid$filtration <- substring(qpcr_2_liquid$sample_2,3,3)
qpcr_2_liquid <- qpcr_2_liquid[, c(2,6,4)]

# combine liquid samples and filter retentions
qpcr_2_merge <- rbind(qpcr_2_filter_sum, qpcr_2_liquid)

# create facet category
qpcr_2_merge$is_cat <- substring(qpcr_2_merge$category, 1, 1) 
qpcr_2_merge$pma <- substring(qpcr_2_merge$category, 2, 2)
qpcr_2_merge <- qpcr_2_merge %>% mutate_if(is.character, as.factor)

# make plots
p_qpcr_2_merge <- ggbarplot(qpcr_2_merge, x = "filtration", y = "copy_number_per_250uL_sample",  add = c("mean_se", "jitter"), facet.by = c("is_cat","pma"), panel.labs = list(is_cat = c("Sample + internal standard", "Sample", "External standard"), pma = c("PMA-", "PMA+") ), panel.labs.font = list(size = 12), fill = "filtration", palette = "jco", legend = "none", xlab = "", ylab = "16S rRNA gene copy number per 250 \u03BCL sample", ggtheme = theme_pubr(border = TRUE))  + scale_x_discrete(labels=c("F" = "Filtration\n+", "N" = "Filtration\n-", "retention" = "Filter\nretention"))


# calculate percentage loss
qpcr_2_merge_summarise <- qpcr_2_merge %>%
  group_by(category, filtration) %>%
  summarise(mean = mean(copy_number_per_250uL_sample))

qpcr_2_merge_summarise <- pivot_wider(qpcr_2_merge_summarise, names_from = filtration, values_from = mean)  # convert to wide format

qpcr_2_merge_summarise <- qpcr_2_merge_summarise %>%
  mutate(loss = N - F, 
         loss_percent = loss/N*100,
         retention_percent = retention/N*100)
# This data frame contains percentage loss & retention

# output this table
#// write.xlsx(qpcr_2_merge_summarise, file="out_percentage_loss_retention.xlsx")


## Tukey post-hoc test
anova.stat <- qpcr_2_merge %>%
  group_by(category) %>%
  tukey_hsd(copy_number_per_250uL_sample ~ filtration) %>%
  add_significance()

get_anova_table(anova.stat, correction = "auto")  # this outputs the test results

p_qpcr_2_merge + 
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12, face = "plain")

#// ggsave("fig/sample_loss.pdf", width =  5, height = 7, unit = "in", dpi = 300)
#// ggsave(filename="/Users/jiaxianshen/Box/Publication/CDC_paper/figures/b_sample_loss.pdf", plot=p_qpcr_2_merge, device=cairo_pdf, width = 6, height = 8.2, units = "in", dpi = 300)

### reference -----



