######## cox lasso model for feature selection ###
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

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Primary_grade4_survival/")

set.seed(7)

data <- read.table("185_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

head(data[1:25],2)
dim(data)

Clin_data <- data[1:25]
head(Clin_data)
dim(Clin_data)

head(data[25],2)

data1 <- subset(data, OS_time!="NA")
dim(data1)
data1 <- subset(data1, OS_month > 0)
dim(data1)
data1 <- subset(data1, Death_Status!="NA")
dim(data1)

#data1 <- subset(data1, OS < 0)

surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
set.seed(10)
cvfit1 <- cv.glmnet(as.matrix(data1[26:ncol(data1)]), 
                    surv_object2, # create survival object from the data
                    family = "cox", # specify Cox PH model
                    type.measure = "C", 
                    nfolds = 5, 
                    alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                    maxit = 1000)

lambda_min <- cvfit1$lambda.min
lambda_min

plot(cvfit1)
#plot(cvfit1)

#jpeg(file="Cox_Lasso_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=350)
#plot(cvfit1)
#dev.off()


est.coef = coef(cvfit1, s = cvfit1$lambda.min) # returns the p length coefficient vector
head(est.coef)
# of the solution corresponding to lambda 
active.k = which(est.coef != 0)
active.k

#Extract coefficinet values
active.k.vals = est.coef[active.k]
active.k.vals

key_variables <- as.data.frame(est.coef[est.coef[,1]!=0,])
colnames(key_variables) <- c("coeff")
key_variables <- round(key_variables,3)
key_variables 
dim(key_variables)

#write.table(key_variables,file="Lasso_key_variables.txt",row.names=T,col.names=T,sep = '\t', quote = F)
write.table(cbind("ID"=rownames(key_variables), key_variables),file="Lasso_11key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",sep="\t",quote=F, row.names=F)




############## PI Model #######


######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
sel_features_results<-read.table("Lasso_11key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)
head(sel_features_results)
dim(sel_features_results)


#data preparation with selected features 
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(sel_features_results$ID), ])
sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(sel_features_results)), ])
head(sel_train,2)
dim(sel_train)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data1["OS_month"],data1["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


#train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

head(tr,4)
tr[,2]
sel_features_results[2,1]
head(sel_features_results,2)

PI_tr=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,1]))
}



head(tr,3)
head(sel_features_results)

LUSC_C_tr$PI<-PI_tr
head(LUSC_C_tr)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_with_11_genes_PI.txt",sep="\t",quote=F, row.names=F)

tr_PI <- as.data.frame(LUSC_C_tr$PI)
rownames(tr_PI) <- rownames(LUSC_C_tr)
colnames(tr_PI) <- c("PI")
write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="tr_PI",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_tr <- Surv(time = LUSC_C_tr$OS_month, event = LUSC_C_tr$Death_Status)
surv_object_tr
dim(surv_object_tr )


dim(tr)
head(tr_PI)
mean(tr_PI$PI)

# ###survival analysis: fits cox ph model to find HR for PI
fit_tr <- survfit(surv_object_tr~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_tr)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_tr <- coxph(surv_object_tr ~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit.coxph_tr)

fit.coxph_tr$linear.predictors


tr_res <- summary(fit.coxph_tr)
coeff <- round(tr_res$coefficients[1],2)
HR <- round(tr_res$conf.int[1],2)
int1 <- round(tr_res$conf.int[3],2)
int2 <- round(tr_res$conf.int[4],2)
CI <- round (tr_res$concordance[1],2)
pval <- tr_res$sctest[3]
pval1 <- format(pval, scientific=T)
pval1

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

tr_res1 <- cbind(coeff, HR1, CI, pval1)
tr_res1
#save results as a file
write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res_PI.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI,
                  #pval=TRUE,
                  risk.table=TRUE,
                  tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  #title= paste0("KM plot based on PI/Risk score"),
                  surv.median.line = "hv", # Add medians survival
                  palette=c("dodgerblue2", "red"), #add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  #legend.title = paste0(pathway),
                  legend.labs = c("Less than mean " , " more than Mean"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp

# customised the plot
pp$plot <- pp$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    #label = "HR = 0.9 \n p < 0.001",
    #label =  paste0("HR = ", HR1, \n, "p-val <", pval),
    label =  paste0("HR = ", HR1, "\n",  "p-val = ", pval1,  "\n",  "C-Index = ", CI),
    size = 5)

# now plot
pp


jpeg("KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()


#ROC plots

LUSC_C_tr$PI


##Create ROC plot for 1-,3- and 5-years survival time prediction
tr_roc1 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
tr_roc1

tr_roc3 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 36,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc3

tr_roc4 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc4
tr_roc5 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5
tr_roc6 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 72,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")

tr_roc6
jpeg(file="ROC_train.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()
