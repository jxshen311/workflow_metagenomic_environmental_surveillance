library(caret)

setwd("/Users/jiaxianshen/Documents/HartmannLab/CDC_project/metaseq/output/nonpareil/ml")

load(file = "cdc_nonpareil_ml_quest_c1.RData")

cm_manual = as.matrix(table(Predicted = predict(models[[1]], test.data), Actual = test.data$div_c1))

class_no = ncol(cm_manual)
sens <- numeric(class_no)
spec <- numeric(class_no)
bal_accu <- numeric(class_no)
# class 1 (account for neibours)----
tp = cm_manual[1,1] + cm_manual[2,1]
tn = sum(cm_manual[-1,-1])
fp = sum(cm_manual[1,]) - cm_manual[1,1]
fn = sum(cm_manual[ ,1]) - tp

sum(cm_manual) == tp + tn + fp + fn

sens[1] = tp / (tp + fn)
spec[1] = tn / (fp + tn)


# class middle (account for neibours)----
for (ii in 2:(class_no-1)){
  tp = cm_manual[ii,ii] + cm_manual[ii-1,ii] + cm_manual[ii+1,ii]
  tn = sum(cm_manual[-ii,-ii])
  fp = sum(cm_manual[ii,]) - cm_manual[ii,ii]
  fn = sum(cm_manual[ ,ii]) - tp
  
  sum(cm_manual) == tp + tn + fp + fn
  
  sens[ii] = tp / (tp + fn)
  spec[ii] = tn / (fp + tn)
  
}

# class last one (account for neibours)----
tp = cm_manual[class_no,class_no] + cm_manual[class_no - 1,class_no]
tn = sum(cm_manual[-class_no,-class_no])
fp = sum(cm_manual[class_no,]) - cm_manual[class_no,class_no]
fn = sum(cm_manual[ ,class_no]) - tp

sum(cm_manual) == tp + tn + fp + fn

sens[class_no] = tp / (tp + fn)
spec[class_no] = tn / (fp + tn)

for (ii in 1:class_no){bal_accu[ii] = (sens[ii] + spec[ii])/2}

mean(bal_accu)  # 0.7717405

# Save for reference ----
# * test codes for calculating sensitivity and specificity ----
if(FALSE){
  sensitivity(cm_manual, rownames(cm_manual)[2])
  specificity(cm_manual, rownames(cm_manual)[c(1,3:7)])
  
  # class 1 ----
  tp = cm_manual[1,1]
  tn = sum(cm_manual[-1,-1])
  fp = sum(cm_manual[1,]) - cm_manual[1,1]
  fn = sum(cm_manual[ ,1]) - cm_manual[1,1]
  
  sens = tp / (tp + fn)
  spec = tn / (fp + tn)
  
  # class 2 ----
  tp = cm_manual[2,2]
  tn = sum(cm_manual[-2,-2])
  fp = sum(cm_manual[2,]) - cm_manual[2,2]
  fn = sum(cm_manual[ ,2]) - cm_manual[2,2]
  
  sum(cm_manual) == tp + tn + fp + fn
  
  sens = tp / (tp + fn)
  spec = tn / (fp + tn)
  
}