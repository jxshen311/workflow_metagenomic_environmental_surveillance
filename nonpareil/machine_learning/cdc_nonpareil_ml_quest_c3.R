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

## random forest: rf
set.seed(87)
mod.rf <- train(x = X, y = Y,               
                method="rf", 
                tuneLength = 5, 
                trControl=tr_ctrl) 

cat("=============================")
cat("model 1: random forest: rf")
print(mod.rf)

rfImp <- varImp(mod.rf)  # only non-fomula format works for generating variable importance
rfImp$importance

# plot
pdf("fig/rfIMP_div_c3.pdf", width =  6, height = 4)
plot(rfImp)
dev.off()



## case 1: no study ----
train_1 <- file[train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("study",colnames(file)))] 
test_1 <- file[-train.index, -c(grep("diversity",colnames(file)), grep("geography",colnames(file)), grep("study",colnames(file)))]

X1 <- train_1[, !(names(train_1) %in% "div_c3")]
Y1 <- train_1$div_c3

## random forest: rf: no study
set.seed(87)
mod.rf.no.study <- train(x = X1, y = Y1,               
                  method="rf", 
                  tuneLength = 5, 
                  trControl=tr_ctrl) 

cat("=============================\n")
cat("case 1: random forest (rf) + no [study]\n")
print(mod.rf.no.study)

rf.no.study.Imp <- varImp(mod.rf.no.study)  # only non-fomula format works for generating variable importance
rf.no.study.Imp$importance

# plot
pdf("fig/rfIMP_div_c3_no_study.pdf", width =  6, height = 4)
plot(rf.no.study.Imp)
dev.off()


stopCluster(cl)  # stop parellel


### save all variables ----
save.image(file = "rdata/cdc_nonpareil_ml_quest_c3.RData")

