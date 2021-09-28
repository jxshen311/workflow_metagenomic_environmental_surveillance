# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 
# purpose: 



# Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(nlme)  # For mixed effects models
library(gtools)
library(reshape2)
# library(goeveg)
library(openxlsx) # handle xlsx files
library(tidyr)



# Functions ----


# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/knead")

# import data ----
seq <- read.csv("reads_seq_success.csv")


# plot 
(p1 <- ggplot(seq, aes(x=dna_con_sub)) +
    geom_point(aes(y= raw, color = success), size =2) +
    
    geom_vline(xintercept = 0.75, linetype="longdash", 
               color = "blue", size=0.6) + 

    
    geom_hline(yintercept = 1e5, linetype="longdash", 
               color = "blue", size=0.6) + 
    labs(title="Relationship between reads number and submitted DNA concentration", x= expression(paste("submitted DNA concentration (ng/", mu, "L)")), y="Number of raw reads") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 0.2, hjust = 0.5),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 10, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),                                # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"),  # Adding a 1cm margin around the plot
          
          legend.title = element_text(colour="black", size=12, face="plain"),
          legend.text = element_text(colour="black", size = 12, face = "plain"),
          legend.position=c(0.9,0.9),
          legend.background=element_blank()
    ))             

ggsave("fig/fig_reads number~DNA concentration_eposter.pdf", width = 7.5, height = 6, units = "in", dpi = 300)

# (text size for e-poster)
(p2 <- ggplot(seq, aes(x=dna_input_sub)) +
    geom_point(aes(y= raw, color = success), size =2) +
    
    geom_vline(xintercept = 11.2, linetype="longdash", 
               color = "blue", size=0.6) + 
    
    
    geom_hline(yintercept = 1e5, linetype="longdash", 
               color = "blue", size=0.6) + 
    
    labs(title="Relationship between reads number and submitted DNA input", x="submitted DNA input (ng)", y="Number of raw reads") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 0.2, hjust = 0.5),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),                                # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"),  # Adding a 1cm margin around the plot
          
          legend.title = element_text(colour="black", size=12, face="plain"),
          legend.text = element_text(colour="black", size = 12, face = "plain"),
          legend.position=c(0.9,0.9),
          legend.background=element_blank()
    ))     

ggsave("fig/fig_reads number~DNA input_eposter.pdf", width = 7, height = 6, units = "in", dpi = 300)
