library(reshape)
library(caret)
library(randomForest)  # rf

setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/nonpareil/ml")

load(file = "cdc_nonpareil_ml_quest_c3.RData")



### estimate skill of model on the validation dataset ----
cm_list <- list()

cm_list[[1]] <- confusionMatrix(predict(mod.rf, test.data), test.data$div_c3)  # statistics

cm_list[[1]]

cm_list[[2]] <- confusionMatrix(predict(mod.rf.no.study, test_1), test_1$div_c3)  # statistics

cm_list[[2]]






# alternatively, manually construct CM and calculate accuracy ----
if(FALSE){
  cm_test = as.matrix(table(Actual = test.data$div_c3, Predicted = predict(models[[1]], test.data)))
  
  sum(diag(cm_test))/length(test.data$div_c3)  # overall accuracy
  
}












# ### plot variable importance ----
# ## rf
# rfImp <- varImp(mod.rf.1)  # only non-fomula format works for generating variable importance
# rfImp$importance
# 
# # plot
# pdf("rfIMP_div_c2_no_study.pdf", width =  6, height = 4)
# plot(rfImp)
# dev.off()