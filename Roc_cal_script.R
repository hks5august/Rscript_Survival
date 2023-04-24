install.packages("pROC")

library(pROC)

#Load required libraries
library(car)
library(caret)
library(dplyr)
library(ggfortify)
library(ggplot2)
library(kernlab)
library(e1071)
library(randomForest)
library(DataExplorer)
library(ROSE)
library(skimr)
library(RANN)
library(fastAdaboost)
library(xgboost)
library(caretEnsemble)
library(C50)
library(earth)
library(gbm)
library(rpart)
library(bnclassify)
library(RSNNS)

test1_pred <- read.table("test1_pred", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#test2_pred <- read.table("test2_pred", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


colnames(test1_pred)

roc_data <- test1_pred
E= length(roc_data )
for(i in seq(from=2, to=E ,by=1))
{
  #roc1 = roc(as.factor(test1_pred[,1]), as.factor(test1_pred[,i]))
  #ROC curve
   roc<- roc.curve(as.factor(roc_data[,1]), as.factor(roc_data[,i]), plotit = T)
   roc_val <- round(roc$auc,2)
   #print(roc_val )
   write.table(cbind(colnames(roc_data[i]),roc_val), file="roc_results.txt",row.names=F,col.names=F,sep = '\t',append = T)
}


?roc.curve
roc.curve(as.factor(test1_pred[,1]), as.factor(test1_pred[,2]), plotit = T, add=TRUE, main= paste0("ROC curve", names(test1_pred[,2])) )

roc.curve(as.factor(test1_pred[,1]), as.factor(test1_pred[,2]), plotit = F, add=F, main= "ROC curve" )


