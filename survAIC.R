

#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(MASS)
library(survivalROC)
library(glmnet)
args <- commandArgs(TRUE)

set.seed(7)
setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/ssGSEA/Reactome/")
#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")


#data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
data <- read.table("ssgsea_train_reactome_with_Clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

data1 <- subset(data, OS_month!="NA")
data1_te <- subset(data_test, OS_month!="NA")

head(data1)
dim(data1)

##############

######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
features<- read.table("16_lasso_features_list", header =TRUE,  sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)

head(features)
rownames(features) <- features$ID
head(features)

######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
sel_features_results<-read.table("16_Lasso_key_variables_recatome_from_116.txt",header =TRUE, sep = "\t", check.names = FALSE, row.names=1)

head(sel_features_results)
dim(sel_features_results)


#data preparation with selected features 
sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
head(sel_train,1)
dim(sel_train)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data1["OS_month"],data1["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


# create survival object
surv_object<- Surv(time = as.numeric(train_feature_mat1$OS_month), event = train_feature_mat1$Death_Status)

train_feature_mat2 <- train_feature_mat1[3:ncol(train_feature_mat1)]
train_feature_mat3 <- cbind(surv_object, train_feature_mat2)
head(train_feature_mat3,4)
#Cox multivariate model
#cox_multivariate <- coxph(surv_object ~ train_feature_mat1[3:ncol(train_feature_mat1)],   data=train_feature_mat1 )
cox_multivariate <- coxph(surv_object ~ .,   data=train_feature_mat3 )

summary(cox_multivariate)


stepAIC(cox_multivariate , direction="both")


fwd <- stepAIC(cox_multivariate , direction="forward")

fwd

stepAIC(cox_multivariate , direction="forward")

bckwd <- stepAIC(cox_multivariate , direction="backward")

bckwd

stepAIC(cox_multivariate , direction="backward")

both <- stepAIC(cox_multivariate , direction="both")

both 

stepAIC(cox_multivariate , direction="both")

