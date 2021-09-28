# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 2020-03-24
# purpose: 



# Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(nlme)  # For mixed effects models
library(gtools)
library(reshape2)
#library(goeveg)
library(openxlsx) # handle xlsx files
library(tidyr)



# Functions ----


# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/knead")

# import data ----
reads <- read.csv("reads_no.csv")[ ,c(1,2,3,5)]

# separate sample_id column
reads <- separate(reads, "sample_id", c("sampling_event","project","sample","dna_extraction","replicate","dilution"), sep = "-")


# plot negative controls (TODO: change mean calculation // search R TODO // nc section is not reproducible for running again)----
nc <- filter(reads, quality_sequenced_reads_no == "#N/A")


nc <- mutate(nc, sample_group = ifelse(sample == "NMCA"| sample == "NMCB", "NMC",
   ifelse(sample == "NFCA"| sample == "NFCB", "NFC", sample)))
                 

nc2 <- data_summary(nc, varname="R1_reads_no", 
                    groupnames=c("sample_group"))

# plot
target <- c("NFC100", "NFC80", "NFC41", "NFC5", "NKC", "NMC", "NFC"  )

p_nc <- ggplot(nc2, aes(x=factor(sample_group,level=target), y=R1_reads_no)) + 
  geom_line() +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=R1_reads_no-sd, ymax=R1_reads_no+sd), width=.2,
                position=position_dodge(0.05))+
  labs(title="Read number of negative controls", x="control sample", y = "number of reads")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),    
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "plain"),  
        plot.title = element_text(lineheight=.8, face="bold"),
        # panel.grid = element_blank(),                                # Removing the background grid lines      
        plot.margin = unit(c(1,1,1,1), units = , "cm")                # Adding a 1cm margin around the plot
        )
  
p_nc
  

# plot AG with 3 dilutions ----
ag <- filter(reads, sample == "AG")


# calculate quality_sequenced_ratio
ag$quality_sequenced_reads_no <- as.numeric(levels(ag$quality_sequenced_reads_no))[ag$quality_sequenced_reads_no]

ag$quality_sequenced_ratio = ag$quality_sequenced_reads_no/ag$R1_reads_no*100


ag$sample <- as.factor(ag$sample)
ag$dilution <- as.factor(ag$dilution)

ag <- ag %>%
  mutate(sample_dilution = paste(sample, dilution, sep = "-"))

ag$sample_dilution <- as.factor(ag$sample_dilution)

ag2 <- ag %>%
  group_by(sample_dilution) %>%
  summarise(R1_reads_no_mean = mean(R1_reads_no),
            R1_reads_no_sd = sd(R1_reads_no),
            quality_sequenced_reads_no_mean = mean(quality_sequenced_reads_no),
            quality_sequenced_reads_no_sd = sd(quality_sequenced_reads_no),
            quality_sequenced_ratio_mean = mean(quality_sequenced_ratio),
            quality_sequenced_ratio_sd = sd(quality_sequenced_ratio)
  )


target <- c("AG-1", "AG-5", "AG-10")

# plot
(p1 <- ggplot(ag2, aes(x=factor(sample_dilution, levels = target))) +
    geom_point(aes(y= R1_reads_no_mean, color = "#FF8C00"), size =2) +
    geom_line(aes(y= R1_reads_no_mean, color = "#FF8C00"), group = 1) +
    geom_errorbar(aes(ymin=R1_reads_no_mean-R1_reads_no_sd, ymax=R1_reads_no_mean+R1_reads_no_sd, color = "#FF8C00"), width=.1, position=position_dodge(0.05)) +
    
    geom_point(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), size =2) +
    geom_line(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), group = 1) +
    geom_errorbar(aes(ymin=quality_sequenced_reads_no_mean-quality_sequenced_reads_no_sd, ymax=quality_sequenced_reads_no_mean+quality_sequenced_reads_no_sd, color = "#21D19F"), width=.1) +
    
    geom_bar(aes(y = quality_sequenced_ratio_mean*3000), fill = "#DB99999E",  width = 0.3, stat = 'identity') +
    geom_errorbar(aes(ymin=(quality_sequenced_ratio_mean-quality_sequenced_ratio_sd)*3000, ymax=(quality_sequenced_ratio_mean+quality_sequenced_ratio_sd)*3000, width=.1)) + 
    scale_y_continuous(sec.axis = sec_axis(~./3000, name = "Quality-sequenced ratio (%)")) +
    
    scale_colour_discrete(name  ="",
                          breaks=c("#FF8C00", "#21D19F"),
                          labels=c("Raw reads", "Quality-sequenced reads")) +
    
    labs(title="Reads number of aggregate samples", x="Aggregate sample at three dilution levels", y="Number of reads") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 0.2, hjust = 0.5),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),                                # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"),  # Adding a 1cm margin around the plot
          
          legend.title = element_text(colour="black", size=9, face="plain"),
          legend.text = element_text(colour="black", size = 9, face = "plain"),
          legend.position=c(1,1.15),
          legend.background=element_blank()
    ))             

# plot reads of all samples except negative controls (dilution = 1) ----
sp <- filter(reads, dilution == "1" & quality_sequenced_reads_no != "#N/A")

# calculate quality_sequenced_ratio
sp$quality_sequenced_reads_no <- as.numeric(levels(sp$quality_sequenced_reads_no))[sp$quality_sequenced_reads_no]

sp$quality_sequenced_ratio = sp$quality_sequenced_reads_no/sp$R1_reads_no*100


sp$sample <- as.factor(sp$sample)
sp$dilution <- as.factor(sp$dilution)



sp_temp <- filter(sp, replicate == "swab")
sp <- filter(sp, replicate != "swab")

sp_temp <- sp_temp %>%
  mutate(sample2 = paste(sample, replicate, sep = "-"))

sp_temp <- select(sp_temp, 11, 8:10)
sp_temp <- rename(sp_temp, sample = sample2)

sp <- select(sp, 3, 8:10)

sp <- rbind(sp, sp_temp)

sp2 <- sp %>%
  group_by(sample) %>%
  summarise(R1_reads_no_mean = mean(R1_reads_no),
            R1_reads_no_sd = sd(R1_reads_no),
            quality_sequenced_reads_no_mean = mean(quality_sequenced_reads_no),
            quality_sequenced_reads_no_sd = sd(quality_sequenced_reads_no),
            quality_sequenced_ratio_mean = mean(quality_sequenced_ratio),
            quality_sequenced_ratio_sd = sd(quality_sequenced_ratio)
  )

sp2 <- sp2[order(-sp2$R1_reads_no_mean), ]


# plot
(p2 <- ggplot(sp2, aes(x=reorder(sample,-R1_reads_no_mean))) +
    geom_point(aes(y= R1_reads_no_mean, color = "#FF8C00"), size =2) +
    geom_line(aes(y= R1_reads_no_mean, color = "#FF8C00"), group = 1) +
    geom_errorbar(aes(ymin=R1_reads_no_mean-R1_reads_no_sd, ymax=R1_reads_no_mean+R1_reads_no_sd, color = "#FF8C00"), width=.1, position=position_dodge(0.05)) +
    
    geom_point(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), size =2) +
    geom_line(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), group = 1) +
    geom_errorbar(aes(ymin=quality_sequenced_reads_no_mean-quality_sequenced_reads_no_sd, ymax=quality_sequenced_reads_no_mean+quality_sequenced_reads_no_sd, color = "#21D19F"), width=.1) +
    
    scale_colour_discrete(name  ="",
                          breaks=c("#FF8C00", "#21D19F"),
                          labels=c("Raw reads", "Quality-sequenced reads")) +
    
    labs(title="Reads number of samples", x="Sample", y="Number of reads") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, vjust = 0.2, hjust = 1, angle = 90),    
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 12, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),                                # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = "cm"),  # Adding a 1cm margin around the plot
          
          legend.title = element_text(colour="black", size=9, face="plain"),
          legend.text = element_text(colour="black", size = 9, face = "plain"),
          legend.position="bottom",
          legend.background=element_blank()
    )) 




(p3 <- ggplot(sp2, aes(x=reorder(sample,-R1_reads_no_mean))) +
  # geom_point(aes(y= quality_sequenced_ratio_mean), size =2) +
  # geom_line(aes(y= quality_sequenced_ratio_mean), group = 1, size =0.5) +
  geom_bar(aes(y = quality_sequenced_ratio_mean), fill = "#DB99999E",  width = 0.3, stat = 'identity') +
  geom_errorbar(aes(ymin=(quality_sequenced_ratio_mean-quality_sequenced_ratio_sd), ymax=(quality_sequenced_ratio_mean+quality_sequenced_ratio_sd)), width=.3, size = 0.2) +
  labs(title="Quality-sequenced ratio of samples", x="Sample", y="Quality-sequenced ratio (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, vjust = 0.2, hjust = 1, angle = 90),    
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12, face = "plain"),  
        plot.title = element_text(lineheight=.8, face="bold"),
        # panel.grid = element_blank(),                                # Removing the background grid lines      
        plot.margin = unit(c(1,1,1,1), units = "cm")))  # Adding a 1cm margin around the plot
      

## plot reads number of standard only samples (text for e-poster)
sp2$sample_group <- substr(sp2$sample,1,1)

sp2_std.only <- filter(sp2, sp2$sample_group == "S")

(p4 <- ggplot(sp2_std.only, aes(x=reorder(sample,-R1_reads_no_mean))) +
    geom_point(aes(y= R1_reads_no_mean, color = "#FF8C00"), size =2) +
    # geom_line(aes(y= R1_reads_no_mean, color = "#FF8C00"), group = 1) +
    geom_errorbar(aes(ymin=R1_reads_no_mean-R1_reads_no_sd, ymax=R1_reads_no_mean+R1_reads_no_sd, color = "#FF8C00"), width=.1, position=position_dodge(0.05)) +
    
    geom_point(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), size =2) +
    # geom_line(aes(y = quality_sequenced_reads_no_mean, color = "#21D19F"), group = 1) +
    geom_errorbar(aes(ymin=quality_sequenced_reads_no_mean-quality_sequenced_reads_no_sd, ymax=quality_sequenced_reads_no_mean+quality_sequenced_reads_no_sd, color = "#21D19F"), width=.1) +
    
    scale_colour_discrete(name  ="",
                          breaks=c("#FF8C00", "#21D19F"),
                          labels=c("Raw reads", "Quality-sequenced reads")) +
    
    labs(title="Reads number of internal standard samples", x="Sample", y="Number of reads") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 0.2, hjust = 1, angle = 90),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),                                # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = "cm"),  # Adding a 1cm margin around the plot
          
          legend.title = element_text(colour="black", size=12, face="plain"),
          legend.text = element_text(colour="black", size = 12, face = "plain"),
          legend.position=c(0.75,0.95),
          legend.background=element_blank()
    )) 

ggsave("fig/reads number~standards_no line.pdf", width =  6, height = 5, unit = "in", dpi = 300)
