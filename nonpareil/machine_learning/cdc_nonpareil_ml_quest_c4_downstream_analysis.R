library(reshape)
library(caret)
library(randomForest)  # rf
library(ggplot2)

setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/nonpareil/ml")

load(file = "cdc_nonpareil_ml_quest_c4.RData")

# plot variable importance ----
## rf
rfImp <- varImp(mod.rf)  # only non-fomula format works for generating variable importance
rfImp$importance

sink("data_importance/rf_Imp_div_c4.txt")
print(rfImp$importance)
sink()

# plot
pdf("rfIMP_div_c4.pdf", width =  6, height = 4)
plot(rfImp)
dev.off()
