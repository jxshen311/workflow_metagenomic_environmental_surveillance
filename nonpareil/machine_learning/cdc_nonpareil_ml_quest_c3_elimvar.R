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

## categorize diversity (2.5 increment)
file$div_c3[file$diversity <= 17.5 ] = 1
file$div_c3[file$diversity > 17.5 ] = 2


file$div_c3 <- as.factor(file$div_c3) # change to factor

file <- file %>% 
  mutate(div_c3 = factor(div_c3, 
                        labels = make.names(levels(div_c3))))


# caret parallel
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# make 80/20 training/testing split ----
set.seed(87) 

train.index <- createDataPartition(file$div_c3, p=0.80, list=FALSE)
train.data <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)))] 
test.data <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)))]


### ml algorithms ----
## Run algorithms using repeated 5-fold cross validation (5 times)
tr_ctrl <- trainControl(method="repeatedcv", number=5, repeats=5,
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        #sampling='down',
                        allowParallel=TRUE)

X <- train.data[, !(names(train.data) %in% "div_c3")]
Y <- train.data$div_c3

## case 6: 2 predictors: location, building ----
train_6 <- file[train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("div_c3",colnames(file)))] 
test_6 <- file[-train.index, c(grep("location",colnames(file)), grep("building",colnames(file)), grep("div_c3",colnames(file)))]

X6 <- train_6[, !(names(train_6) %in% "div_c3")]
Y6 <- train_6$div_c3

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
pdf("fig/div_c3_rf.6.Imp.pdf", width =  6, height = 4)
plot(rf.6.Imp)
dev.off()


## case 7: 1 predictor: location ----
train_7 <- file[train.index, c(grep("location",colnames(file)), grep("div_c3",colnames(file)))] 
test_7 <- file[-train.index, c(grep("location",colnames(file)), grep("div_c3",colnames(file)))]

X7 <- train_7[, !(names(train_7) %in% "div_c3"), drop=FALSE]  # to keep X as a dataframe, otherwise bugging in train()
Y7 <- train_7$div_c3

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
pdf("fig/div_c3_rf.7.Imp.pdf", width =  6, height = 4)
plot(rf.7.Imp)
dev.off()

stopCluster(cl)  # stop parellel


### save all variables ----
save.image(file = "rdata/cdc_nonpareil_ml_quest_c3_elimvar.RData")

