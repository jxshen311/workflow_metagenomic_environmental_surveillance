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

# make 80/20 training/testing split ----
set.seed(87) 

train.index <- createDataPartition(file$div_c1, p=0.80, list=FALSE)
train.data <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)))] 
test.data <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)))]


### ml algorithms ----
## Run algorithms using repeated 5-fold cross validation (5 times)
tr_ctrl <- trainControl(method="repeatedcv", number=5, repeats=5,
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        #sampling='down',
                        allowParallel=TRUE)

X <- train.data[, !(names(train.data) %in% "div_c1")]
Y <- train.data$div_c1

## random forest: rf
set.seed(87)
mod.rf <- train(x = X, y = Y,               
                method="rf", 
                tuneLength = 5, 
                trControl=tr_ctrl) 

cat("=============================")
cat("model 1: random forest: rf")
print(mod.rf)

# rfImp <- varImp(mod.rf)  # only non-fomula format works for generating variable importance
# rfImp$importance
# 
# # plot
# pdf("fig/rfIMP.pdf", width =  6, height = 4) 
# plot(rfImp)
# dev.off()

## glmnet
set.seed(87)
param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05)) # this tuning parameter is from ohara paper
mod.glmnet <- train(div_c1 ~ .,  # non-fomula version is bugging (Something is wrong; all the Accuracy metric values are missing)
                    data=train.data, 
                    method='glmnet',  # Lasso and Elastic-Net Regularized Generalized Linear Models
                    trControl=tr_ctrl,
                    tuneGrid=param_sweep,
                    maxit=1000000)

cat("=============================\n")
cat("model 2: Lasso and Elastic-Net Regularized Generalized Linear Models: glmnet\n")
print(mod.glmnet)

# glmnetImp <- varImp(mod.glmnet)
# glmnetImp$importance
# 
# # plot
# pdf("fig/glmnetImp.pdf", width =  6, height = 4) 
# plot(glmnetImp)
# dev.off()


## regularized random forest
set.seed(87)
param_sweep <- expand.grid(mtry=round(seq(2,ncol(X),length=4)),coefReg=c(.1,.5,.9))
mod.rrf <- train(x = X, y = Y,               
                 method='RRFglobal',  # regularized rf
                 trControl=tr_ctrl,
                 tuneGrid=param_sweep)

cat("=============================")
cat("model 3: regularized random forest: RRFglobal")
print(mod.rrf)

# rrfImp <- varImp(mod.rrf)
# rrfImp$importance

## Stochastic Gradient Boosting
set.seed(87)
mod.gbm <- train(x = X, y = Y,               
                 method='gbm',
                 trControl=tr_ctrl,
                 tuneLength = 5)

cat("=============================")
cat("model 4: Stochastic Gradient Boosting: gbm")
print(mod.gbm)

# gbmImp <- varImp(mod.gbm)
# gbmImp$importance


## Single C5.0 Tree
set.seed(87)
mod.c50 <- train(x = X, y = Y,               
                 method='C5.0Tree',
                 trControl=tr_ctrl)  # No tuning parameters for this model

cat("=============================")
cat("model 5: Single C5.0 Tree: C5.0Tree")
print(mod.c50)

# c50Imp <- varImp(mod.c50)
# c50Imp$importance

## partial least squares
set.seed(87)
mod.pls <- train(div_c1 ~ .,  # non-fomula version is bugging (Something is wrong; all the Accuracy metric values are missing)
                 data=train.data,               
                 method='pls',
                 trControl=tr_ctrl,
                 tuneLength = 5)

cat("=============================")
cat("model 6: partial least squares: pls")
print(mod.pls)

# plsImp <- varImp(mod.pls)  # A model-specific variable importance metric is available.
# plsImp$importance


## linear svm (Support Vector Machines with Linear Kernel)
set.seed(87)
mod.svmlinear <- train(div_c1 ~ .,  # non-fomula version is bugging (Something is wrong; all the Accuracy metric values are missing)
                       data=train.data,              
                       method='svmLinear',
                       trControl=tr_ctrl,
                       tuneLength = 5)

cat("=============================")
cat("model 7: Support Vector Machines with Linear Kernel: svmLinear")
print(mod.svmlinear)

# svmlinearImp <- varImp(mod.svmlinear)  # no "A model-specific variable importance metric is available." mentioned
# svmlinearImp$importance

## rbf svm (Support Vector Machines with Radial Basis Function Kernel)
set.seed(87)
mod.svmradial <- train(div_c1 ~ .,  # non-fomula version is bugging
                       data=train.data,
                       method='svmRadial',
                       trControl=tr_ctrl,
                       tuneLength = 5)

cat("=============================")
cat("model 8: Support Vector Machines with Radial Basis Function Kernel: svmRadial")
print(mod.svmradial)

# svmradialImp <- varImp(mod.svmradial)  # no "A model-specific variable importance metric is available." mentioned
# svmradialImp$importance


## k nearest neighbors
set.seed(87)
mod.kknn <- train(x = X, y = Y, 
                 method='kknn',
                 trControl=tr_ctrl,
                 tuneLength = 5)

cat("=============================")
cat("model 9: k nearest neighbors: kknn")
print(mod.kknn)

# kknnImp <- varImp(mod.kknn)  # no "A model-specific variable importance metric is available." mentioned
# kknnImp$importance

stopCluster(cl)  # stop parellel

# assess model performance on training/validation set CV
models <- list(mod.rf=mod.rf,
               mod.glmnet=mod.glmnet,
               mod.rrf=mod.rrf,
               mod.gbm=mod.gbm,
               mod.c50=mod.c50,
               mod.pls=mod.pls,
               mod.svmlinear=mod.svmlinear,
               mod.svmradial=mod.svmradial,
               mod.kknn=mod.kknn)

resample_models <- resamples(models)

cat("=============================")
cat("summary of models:")
summary(resample_models, metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))

# plot all models
pdf('fig/ml_div_c1_models.pdf',height=6,width=13)
bwplot(resample_models,metric=c("AUC","Kappa", "Mean_Balanced_Accuracy"))
dev.off()

### estimate skill of model on the validation dataset ----
cm_list <- list()
test_mean_balanced_accuracy <- numeric(9)  # 9 models

for (kk in 1:9){
  cm_list[[kk]] <- confusionMatrix(predict(models[[kk]], test.data), test.data$div_c1)
  test_mean_balanced_accuracy[kk] <- mean(cm_list[[kk]]$byClass[ ,"Balanced Accuracy"])  # mean balanced accuracy across diversity categories
}

cat("=============================")
cat("test_mean_balanced_accuracy:")
test_mean_balanced_accuracy

### save all variables ----
save.image(file = "rdata/cdc_nonpareil_ml_quest_c1.RData")

