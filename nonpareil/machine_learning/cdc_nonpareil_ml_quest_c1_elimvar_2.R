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
## case 4: 4 predictors: location, building, country, touch_frequency
train_4 <- file[train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("country",colnames(file)), grep("touch_frequency",colnames(file)), grep("div_c1",colnames(file)))] 
test_4 <- file[-train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("country",colnames(file)), grep("touch_frequency",colnames(file)), grep("div_c1",colnames(file)))]

X4 <- train_4[, !(names(train_4) %in% "div_c1")]
Y4 <- train_4$div_c1

## random forest: 
set.seed(87)
mod.rf.4 <- train(x = X4, y = Y4,               
                method="rf", 
                tuneLength = 5, 
                trControl=tr_ctrl) 

cat("=============================\n")
cat("case 4: mod.rf.4\n")
print(mod.rf.4)

rf.4.Imp <- varImp(mod.rf.4)  # only non-fomula format works for generating variable importance
rf.4.Imp$importance

# plot
pdf("fig/div_c1_rf.4.Imp.pdf", width =  6, height = 4)
plot(rf.4.Imp)
dev.off()

## case 5: 3 predictors: location, building, country
train_5 <- file[train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("country",colnames(file)), grep("div_c1",colnames(file)))] 
test_5 <- file[-train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("country",colnames(file)), grep("div_c1",colnames(file)))]

X5 <- train_5[, !(names(train_5) %in% "div_c1")]
Y5 <- train_5$div_c1

## random forest: 
set.seed(87)
mod.rf.5 <- train(x = X5, y = Y5,               
                        method="rf", 
                        tuneLength = 5, 
                        trControl=tr_ctrl) 

cat("=============================\n")
cat("case 5: mod.rf.5\n")
print(mod.rf.5)

rf.5.Imp <- varImp(mod.rf.5)  # only non-fomula format works for generating variable importance
rf.5.Imp$importance

# plot
pdf("fig/div_c1_rf.5.Imp.pdf", width =  6, height = 4)
plot(rf.5.Imp)
dev.off()


## case 6: 2 predictors: location, building
train_6 <- file[train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("div_c1",colnames(file)))] 
test_6 <- file[-train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("div_c1",colnames(file)))]

X6 <- train_6[, !(names(train_6) %in% "div_c1")]
Y6 <- train_6$div_c1

## random forest: 
set.seed(87)
mod.rf.6 <- train(x = X6, y = Y6,               
                        method="rf", 
                        tuneLength = 5, 
                        trControl=tr_ctrl) 

cat("=============================\n")
cat("case 6: mod.rf.6\n")
print(mod.rf.6)

rf.6.Imp <- varImp(mod.rf.6)  # only non-fomula format works for generating variable importance
rf.6.Imp$importance

# plot
pdf("fig/div_c1_rf.6.Imp.pdf", width =  6, height = 4)
plot(rf.6.Imp)
dev.off()


## case 7: 1 predictor: location
train_7 <- file[train.index, c(grep("location",colnames(file)), grep("div_c1",colnames(file)))] 
test_7 <- file[-train.index, c(grep("location",colnames(file)), grep("div_c1",colnames(file)))]

X7 <- train_7[, !(names(train_7) %in% "div_c1"), drop=FALSE]  # to keep X as a dataframe, otherwise bugging in train()
Y7 <- train_7$div_c1

## random forest: 
set.seed(87)
mod.rf.7 <- train(x = X7, y = Y7,               
                  method="rf", 
                  tuneLength = 5, 
                  trControl=tr_ctrl) 

cat("=============================\n")
cat("case 7: mod.rf.7\n")
print(mod.rf.7)

rf.7.Imp <- varImp(mod.rf.7)  # only non-fomula format works for generating variable importance
rf.7.Imp$importance

# plot
pdf("fig/div_c1_rf.7.Imp.pdf", width =  6, height = 4)
plot(rf.7.Imp)
dev.off()

stopCluster(cl)  # stop parellel

### assess model performance on training/validation set CV ----
models <- list(mod.rf.4 = mod.rf.4,
               mod.rf.5 = mod.rf.5,
               mod.rf.6 = mod.rf.6,
               mod.rf.7 = mod.rf.7)

resample_models <- resamples(models)

cat("=============================\n")
cat("summary of models:\n")
summary(resample_models, metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))

# plot all models
pdf('fig/ml_div_c1_elimvar_2.pdf',height=6,width=13)
bwplot(resample_models,metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))
dev.off()

### estimate skill of model on the validation dataset ----
cm_list <- list()
test_mean_balanced_accuracy <- numeric(4)  


for (kk in 1:4){
  cm_list[[kk]] <- confusionMatrix(predict(models[[kk]], get(paste("test",kk+3,sep = "_"))), get(paste("test",kk+3,sep = "_"))$div_c1)
  test_mean_balanced_accuracy[kk] <- mean(cm_list[[kk]]$byClass[ ,"Balanced Accuracy"])  # mean balanced accuracy across diversity categories
}

cat("=============================\n")
cat("test_mean_balanced_accuracy:\n")
test_mean_balanced_accuracy

### save all variables ----
save.image(file = "rdata/cdc_nonpareil_ml_quest_c1_elimvar_2.RData")

