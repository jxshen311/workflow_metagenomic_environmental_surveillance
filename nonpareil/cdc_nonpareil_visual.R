# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 2020-07-13
# purpose: 



# Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(nlme)  # For mixed effects models
library(gtools)
library(reshape2)
library(goeveg)
library(openxlsx) # handle xlsx files
library(multcomp)
library(car)
library(FSA)
library(tidyverse)
library(rstatix)


# Functions ----




# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/nonpareil")


# Import data ----
# check the data type of columns of interest by str()
file <- read.xlsx("out_parameter_all.xlsx",
                  sheet = 1, startRow = 1, colNames = TRUE,
                  rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                  skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                  check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)

file_pool <- read.xlsx("out_parameter_pool.xlsx",
                  sheet = 1, startRow = 1, colNames = TRUE,
                  rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                  skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                  check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
# change data type
# character to factor
for (ii in c(8:12,14:17)) {
  file[,ii] <-  as.factor(file[,ii])  
}


# density plot/overall (distribution) ----
(p1 <- ggplot(file, aes(x = diversity)) +
  geom_density(aes(y = ..count.., color = study)) +
  geom_density(aes(y = ..count..)) +
  geom_vline(aes(xintercept = median(file$diversity)), 
             linetype = "dashed", size = 0.6) +
  scale_colour_discrete(name="Study",
                      breaks=c("brooks", "chng", "constantinides","lax","mahnert","ohara","shen"),
                      labels=c("Brooks", "Chng", "Constantinides","Lax","Mahnert","O'Hara","Shen")) + 
  labs(x="Nonpareil diversity (Nd)", y="Count of samples")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),    
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"),  
        plot.title = element_text(lineheight=.8, face="bold", size = 15),
        legend.title = element_text(size=12, face="plain"),
        legend.text = element_text(size = 12, face = "plain"),
        legend.position = "right",
        panel.grid = element_blank(),  # Removing the background grid lines      
        )) # Adding a 1cm margin around the plot

# Modification of p1 (keep for reference)
# title="Distribution of nonpareil diversity of hospital-associated\nenvironmental samples"


#// ggsave("fig_selected/density_plot_study_group_eposter.pdf", width =  7, height = 7, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/density_plot_study_group.pdf", width =  7, height = 4, unit = "in", dpi = 300)

# * empirical cumulative distribution ----
# version 1: ggplot2
(p2 <- ggplot(file, aes(x = diversity)) + 
    stat_ecdf(geom = "step") + 
    theme_bw() +
    labs(x="Nonpareil diversity (Nd)", y="Cumulative probability of samples", title = "Empirical cumulative distribution", caption="") +
    geom_hline(yintercept=0.95, linetype="dashed", color = "blue") + 
    geom_hline(yintercept=0.05, linetype="dashed", color = "blue") + 
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          panel.grid = element_blank(),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none")) 

#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/ecf_diversity.pdf", width =  4, height = 4, unit = "in", dpi = 300)


# version 2: plot() in base R
diver_sort <- sort(file$diversity)
e_cdf <- 1:length(diver_sort) / length(diver_sort)

# plot
# #// pdf("fig_selected/ecf_diversity_eposter.pdf", width =  6, height = 6)
# pdf("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/ecf_diversity.pdf", width =  4, height = 4)
# 
# p2 <- plot(diver_sort, e_cdf, type = "s",
#            ylab="Cumulative frequency of samples",
#            xlab="Nonpareil diversity (Nd)",
#            # cex.main = 1.1,
#            cex.lab = 0.8,
#            cex.axis = 0.6) +
#   # title("Empirical cumulative distribution of hospital-associated \nenvironmental samples", adj = 0, line = 1) +
#   abline(h = 0.95, lty = 2, col="blue") +
#   abline(h = 0.05, lty = 2, col="blue")
# 
# dev.off()

# calculate diversity of 95% and 5%, 25% and 75%
diver_sort[which(e_cdf >= 0.05)[1]]
diver_sort[which(e_cdf >= 0.95)[1]]

diver_sort[which(e_cdf >= 0.25)[1]]  # output: 16.19031
diver_sort[which(e_cdf >= 0.75)[1]]  # output: 18.92985
diver_sort[which(e_cdf >= 0.5)[1]]   # output: 17.31595

quantile(file$diversity, 0.25)[[1]]  # 16.19048
quantile(file$diversity, 0.75)[[1]]  # 18.93076
IQR(file$diversity)  # 2.740287

e_cdf[[which(diver_sort >= 15.0)[1]]]
e_cdf[[which(diver_sort >= 15.5)[1]]]
e_cdf[[which(diver_sort >= 20.0)[1]]]



# * test if diversity follows normal distribution + best fit distribution ----
shapiro.test(file$diversity)  # W = 0.95665, p-value = 3.338e-16 -> not follow normal distribution

library(fitdistrplus)
library(logspline)

# extract left and right part of the distribution
# (p_test <- ggplot(file, aes(x = diversity)) +
#     geom_density(aes(y = ..count..))) # Adding a 1cm margin around the plot
# 
# test_1 <- ggplot_build(p_test)
# head(test_1$data[[1]], 3)
# write.csv(test_1$data[[1]], file = "/Users/jiaxianshen/Box/Publication/CDC_paper/figures/test.csv")
# 
# test_g1 <- file$diversity[which(file$diversity <= 16.1087920172407)]
# test_g2 <- file$diversity[which(file$diversity >= 19.3309453207092)]
# descdist(test_g1, discrete = FALSE)

descdist(file$diversity, discrete = FALSE)
descdist(file$diversity, discrete = FALSE, boot = 1000)


fit.norm <- fitdist(file$diversity, "norm")
plot(fit.norm)
summary(fit.norm)
fit.norm$aic

fit.logi <- fitdist(file$diversity, "logis")
plot(fit.logi)
summary(fit.logi)
fit.logi$aic



# significant effect of factor "study" on diversity? (testing) ----
file_study <- filter(file, study != "mahnert")  # exclude mahnert as it only has 1 point

## visualize by box plot
(p_study <- ggplot(file_study, aes(x=study, y=diversity, fill = study)) + geom_boxplot() +
    theme_bw() +
    labs(x="Study", y="Nonpareil diversity (Nd)") +
    scale_x_discrete(breaks=c("brooks", "chng", "constantinides","lax","ohara","shen"),
                     labels=c("Brooks", "Chng", "Constantinides","Lax","O'Hara","Shen")) +
    scale_fill_discrete(name="Study",
                          breaks=c("brooks", "chng", "constantinides","lax","ohara","shen"),
                          labels=c("Brooks", "Chng", "Constantinides","Lax","O'Hara","Shen")) + 
    theme(axis.title = element_text(size = 18, face = "plain"),
          axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 0, angle = 300),    
          axis.text.y = element_text(size = 14),
          plot.title = element_text(lineheight=.8, face="bold", size = 19),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=16, face="plain"),
          legend.text = element_text(size = 16, face = "plain"),
          legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~study_eposter.pdf", width =  6, height = 6, unit = "in", dpi = 300)

## violin plot
(p_study_vio <- ggplot(file_study, aes(x=study, y=diversity, fill = study)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
    labs(x="Study", y="Nonpareil diversity (Nd)") +
    scale_x_discrete(breaks=c("brooks", "chng", "constantinides","lax","ohara","shen"),
                     labels=c("Brooks", "Chng", "Constantinides","Lax","O'Hara","Shen")) +
    scale_fill_discrete(name="Study",
                        breaks=c("brooks", "chng", "constantinides","lax","ohara","shen"),
                        labels=c("Brooks", "Chng", "Constantinides","Lax","O'Hara","Shen")) + 
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 315),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~study_violin.pdf", width =  6, height = 6, unit = "in", dpi = 300)

## ANOVA test (one-way)
res.aov.study <- aov(diversity ~ study, data = file_study)
summary(res.aov.study) # re: significant

# Tukey multiple pairwise-comparisons
res.aov.study.pair <- TukeyHSD(res.aov.study)
res.aov.study.pair <- as.data.frame(res.aov.study.pair[["study"]])
write.csv(res.aov.study.pair, file = "test_study.csv")

# Check ANOVA assumptions
leveneTest(diversity ~ study, data = file_study)  # significant -> violate homogeneity of variance assumption
plot(res.aov.study, 1)

shapiro.test(x = residuals(object = res.aov.study ) )  # p= 1.33e-15 (significant) -> violate normality assumption

## so I use kruskal.test
kruskal.test(diversity ~ study, data = file_study)  # significant

# post-hoc test
DT.study = dunnTest(diversity ~ study,
              data=file_study,
              method="bh")  
write.csv(DT.study[["res"]], file = "fig_selected/test_study.csv")  # export p values

### check the density plot of the separated groups ----
file_study$study_group <- ifelse(file_study$study %in% c("brooks", "chng", "constantinides","lax"), "H", "L")


df_study_group <- file_study %>% 
  group_by(study_group) %>%
  summarise(grp.med = median(diversity))

## density plot after grouping the studies
(p_study_2 <- ggplot(file_study, aes(x = diversity)) +
  geom_density(aes(y = ..count.., color = study_group)) +
  geom_density(aes(y = ..count..)) +
  geom_vline(aes(xintercept = grp.med, color = study_group),
             data = df_study_group, linetype = "dashed") +
  scale_colour_discrete(name="Study group",
    breaks=c("H", "L"),
  labels=c("H: Brooks, Chng, Constantinides, Lax", "L: O'Hara, Shen")) +
  labs(x="Nonpareil diversity (Nd)", y="Count of samples")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),    
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"),  
        plot.title = element_text(lineheight=.8, face="bold"),
        panel.grid = element_blank(),  # Removing the background grid lines      
        plot.margin = unit(c(1,1,1,1), units = , "cm"),
        legend.text = element_text(size = 12),
        legend.position = "bottom")) # Adding a 1cm margin around the plot


## empirical cumulative distribution of separated study groups
# calculate diversity of 95% and 5%
file_study_high <- filter(file_study, study %in% c("brooks", "chng", "constantinides","lax"))
file_study_low <- filter(file_study, study %in% c("ohara", "shen"))

diver_sort_h <- sort(file_study_high$diversity)
e_cdf_h <- 1:length(diver_sort_h) / length(diver_sort_h)
diver_sort_h[which(e_cdf_h >= 0.05)[1]]  # re: 16.4
diver_sort_h[which(e_cdf_h >= 0.95)[1]]  # re: 20.1

diver_sort_l <- sort(file_study_low$diversity)
e_cdf_l <- 1:length(diver_sort_l) / length(diver_sort_l)
diver_sort_l[which(e_cdf_l >= 0.05)[1]]  # re: 15.2
diver_sort_l[which(e_cdf_l >= 0.95)[1]]  #re: 18.2

# plot ECD
(p_study_3 <- ggplot(file_study, aes(x = diversity)) +
  stat_ecdf(aes(color = study_group), 
              geom = "step", size = 1.5) +
  geom_hline(aes(yintercept = 0.05), color = "blue", linetype = "dashed") +
  geom_hline(aes(yintercept = 0.95), color = "blue", linetype = "dashed") +
  labs(x = "Nonpareil diversity (Nd)", y = "Cumulative frequency of samples") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,1,1), units = , "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom"))

### significant effect of factor "sample type" on diversity? (w/testing) ----
## box plot
(p_sample.type_box <- ggplot(file, aes(x=sample_type, y=diversity, fill = sample_type)) + geom_boxplot() +
   theme_bw() +
   labs(x="Sample type", y="Nonpareil diversity (Nd)") +
   theme(axis.title = element_text(size = 14, face = "plain"),
         axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
         axis.text.y = element_text(size = 12),
         plot.title = element_text(lineheight=.8, face="bold", size = 15),
         plot.margin = unit(c(1,1,1,1), units = , "cm"),
         # panel.grid = element_blank(), 
         legend.title = element_text(size=12, face="plain"),
         legend.text = element_text(size = 12, face = "plain"),
         legend.position = "none")) 

#// ggsave("fig_selected/diversity~sample type_2.pdf", width =  4, height = 4, unit = "in", dpi = 300)

## kruskal.test (the dataset violates normality upon quick tests)
kruskal.test(diversity ~ sample_type, data = file)  # not significant (Kruskal-Wallis chi-squared = 0.00012736, df = 1, p-value = 0.991)

t.test(diversity ~ sample_type, data = file)  # t = -1.3363, df = 218.99, p-value = 0.1828

## density plot
df_sample.type_median <- file %>% 
  group_by(sample_type) %>%
  summarise(grp.med = median(diversity))

(p_sample.type_density <- ggplot(file, aes(x = diversity)) +
    geom_density(aes(y = ..count.., color = sample_type)) +
    geom_density(aes(y = ..count..)) +
    geom_vline(aes(xintercept = grp.med, color = sample_type),
               data = df_sample.type_median, linetype = "dashed") +
    labs(x="Nonpareil diversity (Nd)", y="Count of samples")+
    scale_color_discrete(name="Sample type") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = c(0.88,0.85),
          legend.background=element_blank(),
          panel.grid = element_blank(),  # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"))) # Adding a 1cm margin around the plot

#// ggsave("fig_selected/diversity~sample type_density.pdf", width =  6, height = 4, unit = "in", dpi = 300)


# before ggplot: calculate test; run with ggplot
stat.test <- file %>%
  t_test(diversity ~ sample_type, paired = FALSE) %>%
  add_significance() %>%  
  add_xy_position(x = "sample_type")


## violin plot
(p_sample.type_vio <- ggplot(file, aes(x=sample_type, y=diversity)) +
    geom_violin(trim=FALSE, aes(fill = sample_type))+
    geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
    labs(x="Sample type", y="Nonpareil diversity (Nd)") +
    scale_x_discrete(breaks=c("sink", "surface"),
                     labels=c("Sink", "Surface")) +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none") +
    stat_pvalue_manual(stat.test, label = "p = {p}"))


#// ggsave("fig_selected/diversity~sample type_violin.pdf", width =  4, height = 4, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/diversity~sample type_violin.pdf", width =  4, height = 5, unit = "in", dpi = 300)



### significant effect of factor "touch frequency" on diversity? (w/testing) ----
## violin plot
(p_touch.frequency_vio <- ggplot(file, aes(x=touch_frequency, y=diversity, fill = touch_frequency)) +
   geom_violin(trim=FALSE)+
   geom_boxplot(width=0.1, fill = "white") +
   theme_bw() +
   scale_x_discrete(breaks=c("averaged", "high", "low","sink"),
                    labels=c("Averaged", "High", "Low","Sink")) +
   labs(x="Touch frequency", y="Nonpareil diversity (Nd)") +
   theme(axis.title = element_text(size = 14, face = "plain"),
         axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
         axis.text.y = element_text(size = 12),
         plot.title = element_text(lineheight=.8, face="bold", size = 15),
         plot.margin = unit(c(1,1,1,1), units = , "cm"),
         legend.title = element_text(size=12, face="plain"),
         legend.text = element_text(size = 12, face = "plain"),
         legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~touch frequency_violin.pdf", width =  4, height = 4, unit = "in", dpi = 300)

## kruskal.test (the dataset violates normality upon quick tests)
kruskal.test(diversity ~ touch_frequency, data = file)  # significant (Kruskal-Wallis chi-squared = 120.59, df = 3, p-value < 2.2e-16)

# post-hoc test
file$touch_frequency <- as.factor(file$touch_frequency)
DT.touch.frequency = dunnTest(diversity ~ touch_frequency,
                    data=file,
                    method="bh")  
write.csv(DT.touch.frequency[["res"]], file = "fig_selected/test_touch.csv")  # export p values


### significant effect of factor "sampling method" on diversity? (w/testing) ----
file$type_method <- paste(file$sample_type,file$sampling_method,sep=" ")
file$type_method <- as.factor(file$type_method)

# anova + t test
# t test
file %>%
  group_by(sample_type) %>%
  t_test(diversity ~ sampling_method, paired = FALSE, p.adjust.method = "BH") %>%
  add_significance() %>%  
  add_xy_position(x = "sampling_method")


# anova + post hoc
file_sink <- filter(file, sample_type == "sink")
summary(aov(diversity ~ sampling_method, data = file_sink))
aov(diversity ~ sampling_method, data = file_sink) %>% tukey_hsd()

# examine the number of samples for each category
table(file_sink$sampling_method)


file_surface <- filter(file, sample_type == "surface")
table(file_surface$sampling_method)
t.test(diversity ~ sampling_method, data = file_surface, paired = FALSE) %>% add_significance()

## violin plot
p_type.method_vio <- ggplot(file, aes(x=sampling_method, y=diversity)) +
    geom_violin(trim=FALSE, aes(fill = type_method)) +
    # geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
   scale_x_discrete(breaks=c("aspirate","sink","swab","swab","wipe"),
                    labels=c("Aspirate","Sink","Swab","Swab","Wipe")) +
    labs(x="Sampling method", y="Nonpareil diversity (Nd)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none",
          panel.spacing=unit(0, "lines")) + 
  facet_grid(~ sample_type, scales = "free_x", labeller = labeller(sample_type = c(sink = "Sink", surface = "Surface")))

p_type.method_vio


#// ggsave("fig_selected/diversity~type+method_violin.pdf", width =  5, height = 4, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/diversity~type+method_violin.pdf", width =  5, height = 4, unit = "in", dpi = 300)

## kruskal.test (the dataset violates normality upon quick tests)
kruskal.test(diversity ~ type_method, data = file)  # slightly significant (Kruskal-Wallis chi-squared = 19.077, df = 4, p-value = 0.000759)

# post-hoc test
DT.type_method = dunnTest(diversity ~ type_method,
                              data=file,
                              method="bh")  
write.csv(DT.type_method[["res"]], file = "fig_selected/type_method.csv")  # export p values


file %>%
  group_by(sample_type) %>%
  dunn_test(diversity ~ sampling_method,
            p.adjust.method = "BH") 




### significant effect of factor "geography" on diversity? (w/testing) ----
## violin plot
(p_geography_vio <- ggplot(file, aes(x=geography, y=diversity, fill = geography)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
    # scale_x_discrete(breaks=c("sink aspirate","sink sink","sink swab","surface swab","surface wipe"),
    #                  labels=c("sink\naspirate","sink","sink\nswab","surface\nswab","surface\nwipe")) +
    labs(x="Geography", y="Nonpareil diversity (Nd)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 315),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~geography_violin.pdf", width =  6, height = 4, unit = "in", dpi = 300)

## group geography into country
file$country <- file$geography
file$country <- ifelse(file$country %in% c("Chicago", "Pittsburgh", "w_coast", "s_w_w_coast", "w", "s_e","e_coast"), "US", file$country )

file_country <- file %>%
  filter(country != "Austria")

file_country$country <- as.factor(file_country$country)

(p_country_vio <- ggplot(file_country, aes(x=country, y=diversity, fill = country)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
    # scale_x_discrete(breaks=c("sink aspirate","sink sink","sink swab","surface swab","surface wipe"),
    #                  labels=c("sink\naspirate","sink","sink\nswab","surface\nswab","surface\nwipe")) +
    labs(x="Geography", y="Nonpareil diversity (Nd)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~country_violin.pdf", width =  5, height = 4, unit = "in", dpi = 300)

kruskal.test(diversity ~ country, data = file_country)
DT.country = dunnTest(diversity ~ country,
                          data=file_country,
                          method="bh")  
write.csv(DT.country[["res"]], file = "fig_selected/country.csv")

## within US
file_us <- file %>%
  filter(country == "US")

file_us$geography <- as.factor(file_us$geography)

(p_us_vio <- ggplot(file_us, aes(x=geography, y=diversity, fill = geography)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill = "white") +
    theme_bw() +
    # scale_x_discrete(breaks=c("sink aspirate","sink sink","sink swab","surface swab","surface wipe"),
    #                  labels=c("sink\naspirate","sink","sink\nswab","surface\nswab","surface\nwipe")) +
    labs(x="Geography", y="Nonpareil diversity (Nd)") +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 315),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          plot.margin = unit(c(1,1,1,1), units = , "cm"),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~us_geography_violin.pdf", width =  5, height = 4, unit = "in", dpi = 300)

kruskal.test(diversity ~ geography, data = file_us)
DT.us = dunnTest(diversity ~ geography,
                      data=file_us,
                      method="bh")  
write.csv(DT.us[["res"]], file = "fig_selected/us.csv")


### significant effect of factor "pooling" on diversity? (w/testing) ----
## violin plot
(p_pool_vio <- ggplot(file, aes(x=pooled, y=diversity, fill = pooled)) +
   geom_violin(trim=FALSE)+
   geom_boxplot(width=0.1, fill = "white") +
   theme_bw() +
   scale_x_discrete(breaks=c("averaged","N","pooled"),
                    labels=c("Averaged","Not pooled","Pooled")) +
   labs(x="Pooling", y="Nonpareil diversity (Nd)") +
   theme(axis.title = element_text(size = 14, face = "plain"),
         axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
         axis.text.y = element_text(size = 12),
         plot.title = element_text(lineheight=.8, face="bold", size = 15),
         plot.margin = unit(c(1,1,1,1), units = , "cm"),
         legend.title = element_text(size=12, face="plain"),
         legend.text = element_text(size = 12, face = "plain"),
         legend.position = "none")) # exported // text size was changed for e-poster

#// ggsave("fig_selected/diversity~pool_violin.pdf", width =  5, height = 5, unit = "in", dpi = 300)

kruskal.test(diversity ~ pooled, data = file)
DT.pool = dunnTest(diversity ~ pooled,
                 data=file,
                 method="bh")  
write.csv(DT.pool[["res"]], file = "fig_selected/pool.csv")

# *  control the location, compare pooling ----
file_pool$loc_adj <- file_pool$location
file_pool$loc_adj[which(file_pool$location == "Sink basin")] <- "Sink"
file_pool$loc_adj[which(file_pool$location == "Sink trap")] <- "Sink"
file_pool$loc_adj[which(file_pool$location == "Rear cabinets/counters")] <- "Counter"
file_pool$loc_adj[which(file_pool$location == "Rear handles/rails")] <- "Hand rail"

# t test
file_pool %>%
  filter(loc_adj %in% c("Sink"  ,    "Monitor"  , "Counter")) %>%
  group_by(loc_adj) %>%
  t_test(diversity ~ pooled, paired = FALSE) %>%
  add_significance()

# box plot
(p_pool_control <- ggplot(file_pool, aes(x=loc_adj, y=diversity, fill = pooled)) + geom_boxplot() +
    theme_bw() +
    labs(x="Location", y="Nonpareil diversity (Nd)") +
    scale_fill_discrete(name="Sample pooling",
                        breaks=c("N","pooled"),
                        labels=c("No","Yes")) +
    theme(axis.title = element_text(size = 14, face = "plain"),
          axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
          axis.text.y = element_text(size = 12),
          plot.title = element_text(lineheight=.8, face="bold", size = 15),
          legend.title = element_text(size=12, face="plain"),
          legend.text = element_text(size = 12, face = "plain"),
          legend.position = "bottom")) 


#// ggsave("fig_selected/diversity~pool_location.pdf", width =  7, height = 5, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/d_diversity~pool_location.pdf", width =  5, height = 4.55, unit = "in", dpi = 300)

# test (previous)
file_pool  %>% 
  group_by(loc_adj) %>%
  summarise(p.value=wilcox.test(diversity[pooled=="N"], diversity[pooled=="pooled"],exact = NULL, alternative = "two.sided")$p.value)


### plot the relationship between diversity and LRstar ----
(ggplot(file, aes(x = diversity, y = LRstar)) +
  geom_point() + scale_y_log10() +
  stat_smooth(method = lm)) + 
  theme_bw() +
  labs(x="Nonpareil diversity (Nd)", y="Estimated sequencing effort\nat 95% coverage (bp)") +
  theme(axis.title = element_text(size = 14, face = "plain"),
        axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
        axis.text.y = element_text(size = 12),
        plot.title = element_text(lineheight=.8, face="bold", size = 15),
        legend.title = element_text(size=12, face="plain"),
        legend.text = element_text(size = 12, face = "plain")
        ) 


#// ggsave("fig_selected/linear-diversity~LRstar.pdf", width =  7, height = 5, unit = "in", dpi = 300)
#// ggsave("/Users/jiaxianshen/Box/Publication/CDC_paper/figures/c_linear-diversity~LRstar.pdf", width =  5, height = 5, unit = "in", dpi = 300)


cor(log(file$LRstar), file$diversity)  # out: 0.7756707
linear_reg <- lm(log(LRstar) ~ diversity, data = file)  # log() is natural log
summary(linear_reg)

pridiction_error_rate = sigma(linear_reg)*100/mean(log(file$LRstar))

### R-squared for all factors (get the explanatory power ranking) ----
file$country <- as.factor(file$country)

r2 <- data.frame("factor" = c("location", 
                              "study", 
                              "building", 
                              "country", 
                              "touch_frequency", 
                              "pooled", 
                              "sampling_method", 
                              "sample_type"), "r.squared" = NA, "adj.r.squared" = NA)
# location
r2[1,2] = summary(lm(diversity ~ location, data=file))$r.squared
r2[1,3] = summary(lm(diversity ~ location, data=file))$adj.r.squared

# building
r2[2,2] = summary(lm(diversity ~ building, data=file))$r.squared
r2[2,3] = summary(lm(diversity ~ building, data=file))$adj.r.squared

# study
r2[3,2] = summary(lm(diversity ~ study, data=file))$r.squared
r2[3,3] = summary(lm(diversity ~ study, data=file))$adj.r.squared

# country
r2[4,2] = summary(lm(diversity ~ country, data=file))$r.squared
r2[4,3] = summary(lm(diversity ~ country, data=file))$adj.r.squared

# touch_frequency
r2[5,2] = summary(lm(diversity ~ touch_frequency, data=file))$r.squared
r2[5,3] = summary(lm(diversity ~ touch_frequency, data=file))$adj.r.squared

# sample_type
r2[6,2] = summary(lm(diversity ~ sample_type, data=file))$r.squared
r2[6,3] = summary(lm(diversity ~ sample_type, data=file))$adj.r.squared

# sampling_method
r2[7,2] = summary(lm(diversity ~ sampling_method, data=file))$r.squared
r2[7,3] = summary(lm(diversity ~ sampling_method, data=file))$adj.r.squared

# pooled
r2[8,2] = summary(lm(diversity ~ pooled, data=file))$r.squared
r2[8,3] = summary(lm(diversity ~ pooled, data=file))$adj.r.squared

## plot the varible importance (explanatory power) explained by R2
(ggplot(r2, aes(x = r.squared, y = reorder(factor,r.squared))) +
    geom_point() + 
    geom_segment(aes(x=0,xend=r.squared,y=factor,yend=factor)) +
    scale_x_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
    theme_bw() +
    labs(x="R squared", y="") +
    theme(axis.title = element_text(size = 14, face = "plain"),
        axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 0),    
        axis.text.y = element_text(size = 12),
        plot.title = element_text(lineheight=.8, face="bold", size = 15),
        plot.margin = unit(c(1,1,1,1), units = , "cm"),
        legend.title = element_text(size=12, face="plain"),
        legend.text = element_text(size = 12, face = "plain"))) 


#// ggsave("fig_selected/factor~R squared.pdf", width =  7, height = 5, unit = "in", dpi = 300)


                                                   