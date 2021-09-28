# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 2020-03-29
# purpose: 



# Libraries ----
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(nlme)  # For mixed effects models
library(gtools)
library(reshape2)
library(Nonpareil)


# Functions ----




# Set the working directory ----
setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/nonpareil")

# output nonpareil parameters ----
# kmer mode
meta_kmer <- read.table('metadata_kmer.txt', sep='\t', h=T);
attach(meta_kmer);
np_kmer <- Nonpareil.curve.batch(file, libnames=name, plot.observed = FALSE);
kmer <- print(np_kmer)
df_kmer <- as.data.frame(kmer)
detach(meta_kmer);

# align mode // part 1
meta_align_1 <- read.table('metadata_align_2kb&1kb.txt', sep='\t', h=T);
np_align_1 <- Nonpareil.curve.batch(meta_align_1$file, libnames=meta_align_1$name, modelOnly=TRUE);
align_1 <- print(np_align_1)
df_align_1 <- as.data.frame(align_1)

# align mode // part 2
# C-CDC-SPF80-Lu-2-1.npo is excluded because it is bugging
meta_align_2 <- read.table('metadata_align_below.txt', sep='\t', h=T);
np_align_2 <- Nonpareil.curve.batch(meta_align_2$file, libnames=meta_align_2$name, modelOnly=TRUE);
align_2 <- print(np_align_2)
df_align_2 <- as.data.frame(align_2)

df_align <- rbind(df_align_1, df_align_2)

df_align$kappa <- as.numeric(df_align$kappa)
df_align$C <- as.numeric(df_align$C)
df_align$LR <- as.numeric(df_align$LR)
df_align$modelR <- as.numeric(df_align$modelR)
df_align$LRstar <- as.numeric(df_align$LRstar)
df_align$diversity <- as.numeric(df_align$diversity)

# output nonpareil outputs
write.csv(df_kmer, file = "out_parameter_kmer.csv")
write.csv(df_align, file = "out_parameter_align.csv")

# plot ----
# kmer // coverage
df_kmer <- tibble::rownames_to_column(df_kmer, var = "sample")

# group replicates
df_kmer$sample_s <- strsplit(df_kmer$sample,"-")

for (i in 1:nrow(df_kmer)){
  ifelse(df_kmer$sample_s[[i]][6] == 1,
         df_kmer$sample_s[i] <- df_kmer$sample_s[[i]][3],
         df_kmer$sample_s[i] <- paste(df_kmer$sample_s[[i]][3], df_kmer$sample_s[[i]][6], sep = "-")
         )
}

# manually change for swab sample
df_kmer$sample_s <- ifelse(grepl("swab", df_kmer$sample), paste(df_kmer$sample_s, "swab",sep = "-"), df_kmer$sample_s)

df_kmer$sample_s <- as.factor(as.character(df_kmer$sample_s))

kmer_sum <- df_kmer %>%
  group_by(sample_s) %>%
  summarise(C_mean = mean(C),
            C_sd = sd(C),
            LR_mean = mean(LR),
            LR_sd = sd(LR),
            LRstar_mean = mean(LRstar),
            LRstar_sd = sd(LRstar),
            diversity_mean = mean(diversity),
            diversity_sd = sd(diversity)
  )

# plot
(p1 <- ggplot(kmer_sum, aes(x=reorder(sample,-C_mean))) +
    geom_point(aes(y=C_mean),size = 2) +
    geom_errorbar(aes(ymin=C_mean-C_sd, ymax=C_mean+C_sd), width=.2, position=position_dodge(0.05))+
  
    labs(title="", x="sample", y="mean coverage")+
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 90),    
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14, face = "plain"),  
            plot.title = element_text(lineheight=.8, face="bold"),
            # panel.grid = element_blank(),  # Removing the background grid lines      
            plot.margin = unit(c(1,1,1,1), units = , "cm"))) # Adding a 1cm margin around the plot

   
(p2 <- ggplot(kmer_sum, aes(x=reorder(sample,-C_mean))) +
  geom_point(aes(y=diversity_mean),size = 2) +
  geom_errorbar(aes(ymin=diversity_mean-diversity_sd, ymax=diversity_mean+diversity_sd), width=.2) + 
  
  labs(title="", x="sample", y="mean diversity")+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 90),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),  # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"))) # Adding a 1cm margin around the plot

  
    
# concentration & diversity ----
con <- read.csv("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/dna_concentration.csv")[,c(2,7)]

# group replicates
con$sample_id_dash <- as.character(con$sample_id_dash)
con$sample <- strsplit(con$sample_id_dash,"-")

for (i in 1:nrow(con)){
  ifelse(con$sample[[i]][6] == 1,
         con$sample[i] <- con$sample[[i]][3],
         con$sample[i] <- paste(con$sample[[i]][3], con$sample[[i]][6], sep = "-")
  )}

# manually change for swab sample
con$sample <- ifelse(grepl("swab", con$sample_id_dash), paste(con$sample, "swab",sep = "-"), con$sample)

con$sample <- as.factor(as.character(con$sample))
con$dna_con <- as.numeric(levels(con$dna_con))[con$dna_con]

con_sum <- con %>%
  group_by(sample) %>%
  summarise(dna_con_mean = mean(dna_con),
            dna_con_sd = sd(dna_con))

kmer_sum <- rename(kmer_sum, sample = sample_s)

kmer_sum <- left_join(kmer_sum, con_sum, by = "sample")

(p3 <- ggplot(kmer_sum, aes(x=reorder(sample,-C_mean))) +
    geom_point(aes(y=dna_con_mean),size = 2, color = "orange") +
    geom_errorbar(aes(ymin=dna_con_mean-dna_con_sd, ymax=dna_con_mean+dna_con_sd), width=.2) +
    
    geom_point(aes(y=diversity_mean),size = 2, color = "green") +
    geom_errorbar(aes(ymin=diversity_mean-diversity_sd, ymax=diversity_mean+diversity_sd), width=.2) +
    
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 90),    
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"),  
          plot.title = element_text(lineheight=.8, face="bold"),
          # panel.grid = element_blank(),  # Removing the background grid lines      
          plot.margin = unit(c(1,1,1,1), units = , "cm"))) # Adding a 1cm margin around the plot
)
    #scale_y_continuous(sec.axis = sec_axis(~./, name = ""))
    
# //CONCLUSION: did not see a relaitonship between diversity and biomass (at least for this dataset, but i still think it is worth examing the relationship for the Rush samples)

# other studies ----
# brooks 2017
m_brooks <- read.table('metadata_brooks_kmer.txt', sep='\t', h=T);
np_brooks <- Nonpareil.curve.batch(m_brooks$file, libnames=m_brooks$name, plot.observed = FALSE);
df_brooks <- as.data.frame(print(np_brooks))
write.csv(df_brooks, file = "out_parameter_brooks.csv")

# mahnert
m_mahnert <- read.table('metadata_mahnert_kmer.txt', sep='\t', h=T);
np_mahnert <- Nonpareil.curve.batch(m_mahnert$file, libnames=m_mahnert$name, plot.observed = FALSE);
df_mahnert <- as.data.frame(print(np_mahnert))
write.csv(df_mahnert, file = "out_parameter_mahnert.csv")

# lax
m_lax <- read.table('metadata_lax_kmer.txt', sep='\t', h=T);
np_lax <- Nonpareil.curve.batch(m_lax$file, libnames=m_lax$name, plot.observed = FALSE);
df_lax <- as.data.frame(print(np_lax))
write.csv(df_lax, file = "out_parameter_lax.csv")

# lax_multi
m_lax_multi<- read.table('output_lax_multi/metadata_lax_multi.txt', sep='\t', h=T);
np_lax_multi  <- Nonpareil.curve.batch(m_lax_multi$file, libnames=m_lax_multi$name, plot.observed = FALSE);
df_lax_multi <- as.data.frame(print(np_lax_multi))
write.csv(df_lax_multi, file = "out_parameter_lax_multi.csv")

# ohara
m_ohara<- read.table('output_ohara_kmer/metadata_ohara.txt', sep='\t', h=T);
np_ohara  <- Nonpareil.curve.batch(m_ohara$file, libnames=m_ohara$name, plot.observed = FALSE);
df_ohara <- as.data.frame(print(np_ohara))
write.csv(df_ohara, file = "out_parameter_ohara.csv")

#chng
## old
m_chng<- read.table('output_chng_kmer/metadata_chng.txt', sep='\t', h=T);
np_chng  <- Nonpareil.curve.batch(m_chng$file, libnames=m_chng$name, plot.observed = FALSE);
df_chng <- as.data.frame(print(np_chng))
write.csv(df_chng, file = "out_parameter_chng.csv")

## new
m_chng_new<- read.table('output_chng_kmer/metadata_chng_new.txt', sep='\t', h=T);
np_chng_new  <- Nonpareil.curve.batch(m_chng_new$file, libnames=m_chng_new$name, plot.observed = FALSE);
df_chng_new <- as.data.frame(print(np_chng_new))
write.csv(df_chng_new, file = "out_parameter_chng_new.csv")

# constantinides
m_constantinides<- read.table('output_constantinides_kmer/metadata_constantinides.txt', sep='\t', h=T);
np_constantinides  <- Nonpareil.curve.batch(m_constantinides$file, libnames=m_constantinides$name, plot.observed = FALSE);
df_constantinides <- as.data.frame(print(np_constantinides))
write.csv(df_constantinides, file = "out_parameter_constantinides.csv")


### draft codes for predict sequencing effort from coverage and diversity + predict coverage from sequencing effort ----
# nonpareil 
# ?Nonpareil.curve
# ?Nonpareil.curve.batch
# ?Nonpareil.legend
# ?Nonpareil.f
# ?Nonpareil.antif
# ?Nonpareil.coverageFactor

np_test <- Nonpareil.curve('C-CDC-AG-Lu-1-1.npo')
predict(np_test,1.711330e+08) # Predict coverage from sequencing effort

## get sequencing effort from diversity and targeted coverage

coef(np_test$model)[1]
coef(np_test$model)[2]

(coef(np_test$model)[1]-1)/coef(np_test$model)[2] # diversity
Nonpareil.antif(0.9, coef(np_test$model)[1], coef(np_test$model)[2])  # predict sequencing effort given coefficient and target coverage
# the above operations only work for Formal class 'Nonpareil.Curve', not work for Formal class 'Nonpareil.Set'

# more organized version
pa <- coef(np_example$model)[1]
pb <- coef(np_example$model)[2]

diversity <- (pa-1)/pb 
LR_0.9 <- Nonpareil.antif(0.9, pa, pb)  # sequencing effort at 90% coverage

# reference code copied from "Nonpareil.R"
if(np$has.model){
  pa <- coef(np$model)[1]
  pb <- coef(np$model)[2]
  if(pa > 1) np$diversity <- (pa-1)/pb
  np$LRstar <- Nonpareil.antif(np$star/100, pa, pb)
  np$modelR <- cor(data$y, predict(np, lr=data$x))
}else{
  np$warnings <- c(np$warnings,
                   "Model didn't converge. Try modifying the values of weights.exp.")
}


###  see if similar diversity produce similar curve ----
# diversity = 15.4
m_test_similar_diversity_15.4 <- read.table('metadata_test_similar_diversity_15.4.txt', sep='\t', h=T);
m_test_similar_diversity_15.4$name <- as.character(m_test_similar_diversity_15.4$name)
pdf("fig_selected/np_curve_diversity_15.4.pdf")
np_test_similar_diversity_15.4 <- Nonpareil.set(m_test_similar_diversity_15.4$file, labels=m_test_similar_diversity_15.4$name, plot.opts=list(plot.observed=FALSE)); #plot only the model curve
dev.off()


# diversity = 19.3
m_test_similar_diversity_19.3 <- read.table('metadata_test_similar_diversity_19.3.txt', sep='\t', h=T);
m_test_similar_diversity_19.3$name <- as.character(m_test_similar_diversity_19.3$name)
pdf("fig_selected/np_curve_diversity_19.3.pdf")
np_test_similar_diversity_19.3 <- Nonpareil.set(m_test_similar_diversity_19.3$file, labels=m_test_similar_diversity_19.3$name, plot.opts=list(plot.observed=FALSE)); #plot only the model curve
dev.off()

# diversity = 15.9
m_test_similar_diversity_15.9 <- read.table('metadata_test_similar_diversity_15.9.txt', sep='\t', h=T);
m_test_similar_diversity_15.9$name <- as.character(m_test_similar_diversity_15.9$name)
pdf("fig_selected/np_curve_diversity_15.9.pdf")
np_test_similar_diversity_15.9 <- Nonpareil.set(m_test_similar_diversity_15.9$file, labels=m_test_similar_diversity_15.9$name, plot.opts=list(plot.observed=FALSE)); #plot only the model curve
dev.off()




# reserve alternative codes ----
# # keep the code for reference (one-by-one, easy to troubleshoot)
# # part 2 // one by one
# meta_align_3 <- read.table('metadata_align_below_test.txt', sep='\t', h=T);
# 
# meta_align_3$file <- as.character(meta_align_3$file)
# meta_align_3$name <- as.character(meta_align_3$name)
# 
# np_align_3 <- vector("list", length = nrow(meta_align_3))
# df_align_3 <- data.frame("parameter" = c("kappa","C", "LR", "modelR", "LRstar", "diversity"))
# align_3 <- vector("list", length = nrow(meta_align_3))
# 
# for (ii in 1:nrow(meta_align_3)){
#   np_align_3[[ii]] <- Nonpareil.curve(meta_align_3[ii,1])
#   align_3[[ii]] <- as.numeric(print(np_align_3[[ii]]))
#   align_3[[ii]] <- as.data.frame(align_3[[ii]])
#   colnames(align_3[[ii]]) <- meta_align_3$name[ii]
#   df_align_3 <- cbind(df_align_3, align_3[[ii]])
# }