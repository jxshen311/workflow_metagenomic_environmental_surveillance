# Introduction ----
# author: Jiaxian Shen (jiaxianshen2022@u.northwestern.edu)
# date: 
# purpose: 
# no preprocessing



# Libraries ----
library(easypackages)  # to load multiple packages at once
library(methods) # quest does not load this automatically
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualisation
library(openxlsx) # handle xlsx files
library(caret)  # machine learning
library(doParallel)
library(MLmetrics)  # required for random forest
library(randomForest)  # rf
library(e1071)
libraries("glmnet", "Matrix")  # glmnet
library(RRF)  # RRFglobal
library(gbm)  # gbm
library(C50)  # C5.0Tree
library(pls)  # pls
library(kernlab)  # svmLinear & svmRadial
library(kknn)  # kknn




# Set the working directory
setwd("/projects/p30892/cdc/nonpareil/ml")


# Import data 
# check the data type of columns of interest by str()
file <- read.xlsx("out_parameter_all.xlsx",
                  sheet = 1, startRow = 1, colNames = TRUE,
                  rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                  skipEmptyCols = TRUE, rows = NULL, cols = NULL,
                  check.names = FALSE, namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)[,c(7:12,14:16)]

## add country column
file$country <- file$geography
file$country <- ifelse(file$country %in% c("Chicago", "Pittsburgh", "w_coast", "s_w_w_coast", "w", "s_e","e_coast"), "US", file$country )


# change character to factor
for (i in 2:10){
  file[,i] <- as.factor(file[,i])
}

## check data
anyNA(file)  # no missing value

## categorize diversity (0.5 increment)
file$div_c1[file$diversity <= 15.5 ] = 1


for (ii in seq(0.5,4.5,by=0.5)) {
  file$div_c1[file$diversity > (15+ii) & file$diversity <= (ii+15.5)] = (ii*2+1)
}

file$div_c1[file$diversity > 20 ] = 11

file$div_c1 <- as.factor(file$div_c1) # change to factor

file <- file %>% 
  mutate(div_c1 = factor(div_c1, 
                        labels = make.names(levels(div_c1))))


# caret parallel
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

### make 80/20 training/testing split ----
set.seed(87) 
train.index <- createDataPartition(file$div_c1, p=0.80, list=FALSE)

## Run algorithms using repeated 5-fold cross validation (5 times)
tr_ctrl <- trainControl(method="repeatedcv", number=5, repeats=5,
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        #sampling='down',
                        allowParallel=TRUE)


### ml algorithms ----
## case 1: no study
train_1 <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("study",colnames(file)))] 
test_1 <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("study",colnames(file)))]

X1 <- train_1[, !(names(train_1) %in% "div_c1")]
Y1 <- train_1$div_c1

## random forest: rf: nno study
set.seed(87)
mod.rf.1 <- train(x = X1, y = Y1,               
                method="rf", 
                tuneLength = 5, 
                trControl=tr_ctrl) 

cat("=============================\n")
cat("case 1: random forest (rf) + no [study]\n")
print(mod.rf.1)

rf.1.Imp <- varImp(mod.rf.1)  # only non-fomula format works for generating variable importance
rf.1.Imp$importance

# plot
pdf("fig/rfIMP_div_c1_no_study.pdf", width =  6, height = 4)
plot(rf.1.Imp)
dev.off()


## case 2: eliminate "sampling_method + pooled", keep "study"
train_2 <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("sampling_method",colnames(file)), grep("pooled",colnames(file)))] 
test_2 <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("sampling_method",colnames(file)), grep("pooled",colnames(file)))]

X2 <- train_2[, !(names(train_2) %in% "div_c1")]
Y2 <- train_2$div_c1

## random forest: rf: nno study
set.seed(87)
mod.rf.2 <- train(x = X2, y = Y2,               
                  method="rf", 
                  tuneLength = 5, 
                  trControl=tr_ctrl) 

cat("=============================\n")
cat("case 2: random forest (rf) + no [sampling_method + pooled] + keep [study]\n")
print(mod.rf.2)

rf.2.Imp <- varImp(mod.rf.2)  # only non-fomula format works for generating variable importance
rf.2.Imp$importance

# plot
pdf("fig/rfIMP_div_c1_no sampling_method & pooled_keep study.pdf", width =  6, height = 4)
plot(rf.2.Imp)
dev.off()

## case 3: eliminate "sampling_method + pooled + study"
train_3 <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("sampling_method",colnames(file)), grep("pooled",colnames(file)), grep("study",colnames(file)))] 
test_3 <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("sampling_method",colnames(file)), grep("pooled",colnames(file)), grep("study",colnames(file)))]

X3 <- train_3[, !(names(train_3) %in% "div_c1")]
Y3 <- train_3$div_c1

## random forest: rf: nno study
set.seed(87)
mod.rf.3 <- train(x = X3, y = Y3,               
                  method="rf", 
                  tuneLength = 5, 
                  trControl=tr_ctrl) 

cat("=============================\n")
cat("case 3: random forest (rf) + no [sampling_method + pooled + study]\n")
print(mod.rf.3)

rf.3.Imp <- varImp(mod.rf.3)  # only non-fomula format works for generating variable importance
rf.3.Imp$importance

# plot
pdf("fig/rfIMP_div_c1_no sampling_method & pooled & study.pdf", width =  6, height = 4)
plot(rf.3.Imp)
dev.off()


stopCluster(cl)  # stop parellel

# assess model performance on training/validation set CV
models <- list(mod.rf.1 = mod.rf.1,
               mod.rf.2 = mod.rf.2,
               mod.rf.3 = mod.rf.3)

resample_models <- resamples(models)

cat("=============================\n")
cat("summary of models:\n")
summary(resample_models, metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))

# plot all models
pdf('fig/ml_div_c1_elimvar.pdf',height=6,width=13)
bwplot(resample_models,metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))
dev.off()

### estimate skill of model on the validation dataset ----
cm_list <- list()
test_mean_balanced_accuracy <- numeric(3)  


for (kk in 1:3){
  cm_list[[kk]] <- confusionMatrix(predict(models[[kk]], get(paste("test",kk,sep = "_"))), get(paste("test",kk,sep = "_"))$div_c1)
  test_mean_balanced_accuracy[kk] <- mean(cm_list[[kk]]$byClass[ ,"Balanced Accuracy"])  # mean balanced accuracy across diversity categories
}

cat("=============================\n")
cat("test_mean_balanced_accuracy:\n")
test_mean_balanced_accuracy

### save all variables ----
save.image(file = "rdata/cdc_nonpareil_ml_quest_c1_elimvar.RData")

