
#Load required libraries
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

#library(pca3d)

setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/top_10pathways_survival/16_genes_models/final_Risk_score_model/")
#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")
#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/top_10pathways_survival/Apoptosis")
set.seed(7)



tr_clin_new <- read.table("tr_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te1_clin_new <- read.table("te1_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te2_clin_new <- read.table("te2_clin_PI", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

dim(tr_clin_new)
dim(te1_clin_new)
dim(te2_clin_new)







# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_event)
#Cox multivariate model
#model_tr <- coxph(surv_object_p_tr ~ Gender + Grade + IDH_mutation_status +  codel_1p19q_status +  PI ,   data=tr_clin_new  )
#model_tr <- coxph(surv_object_p_tr ~ Gender + Grade + Class +IDH_mutation_status +  codel_1p19q_status + RT + chemo + PI ,   data=tr_clin_new  )
model_tr <- coxph(surv_object_p_tr ~ Gender + Grade + Class +IDH_mutation_status +  codel_1p19q_status + PI ,   data=tr_clin_new  )

summary(model_tr )
#jpeg(file="train_multivariate.jpeg", units="in", width=10, height=10, res=300)
#jpeg(file="tr_new_multivariate.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="tr_new2_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(model_tr,data=tr_clin_new)
dev.off()


install.packages("survivalROC")
library(survivalROC)



tr_roc <- survivalROC(Stime        = tr_clin_new$OS_month,
            status       = tr_clin_new$Death_event,
            marker       = tr_clin_new$PI,
            predict.time = 12,
            method       = "NNE", lambda = 0.055, span = NULL, window ="symmetric")


plot(tr_roc$FP, tr_roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(tr_roc$AUC,2)),
     ylab="TP",main="AUC for 1 year prediction on Training Data with NNE")
abline(0,1)



#Calculated  Lmbda min value = 0.002102803

tr_roc1 <- survivalROC(Stime        = tr_clin_new$OS_month,
                      status       = tr_clin_new$Death_event,
                      marker       = tr_clin_new$PI,
                      predict.time = 12,
                      method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")


tr_roc3 <- survivalROC(Stime        = tr_clin_new$OS_month,
                       status       = tr_clin_new$Death_event,
                       marker       = tr_clin_new$PI,
                       predict.time = 36,
                       method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")

tr_roc5 <- survivalROC(Stime        = tr_clin_new$OS_month,
                       status       = tr_clin_new$Death_event,
                       marker       = tr_clin_new$PI,
                       predict.time = 60,
                       method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")
tr_roc3
tr_roc3$AUC
?plot
jpeg(file="Training_ROC.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", cex = 2, xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",  cex.lab=1.5, cex.axis=1.5, cex.main=2,
     ylab="TP",main="AUC Curve Training Data")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, cex = 1.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)

dev.off()



# create survival object
surv_object_p_te1 <- Surv(time = as.numeric(te1_clin_new$OS_month), event = te1_clin_new$Death_event)
#Cox multivariate model
#model_te1 <- coxph(surv_object_p_te1 ~ Gender + Grade + IDH_mutation_status +  codel_1p19q_status +  PI ,   data=te1_clin_new  )
#model_te1 <- coxph(surv_object_p_te1 ~ Gender + Grade + Class + IDH_mutation_status +  codel_1p19q_status + RT + chemo + PI ,   data=te1_clin_new  )
model_te1 <- coxph(surv_object_p_te1 ~ Gender + Grade + Class + IDH_mutation_status +  codel_1p19q_status+ PI ,   data=te1_clin_new  )

summary(model_te1 )
#jpeg(file="test1_multivariate.jpeg", units="in", width=10, height=10, res=300)
#jpeg(file="te1_new_multivariate.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="te1_new2_multivariate.jpeg", units="in", width=10, height=10, res=300)

#ggforest(cox_multivariate,data=LGG_clin )
ggforest(model_te1,data=te1_clin_new)
dev.off()




te1_roc1 <- survivalROC(Stime        = te1_clin_new$OS_month,
                       status       = te1_clin_new$Death_event,
                       marker       = te1_clin_new$PI,
                       predict.time = 12,
                       method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")


te1_roc3 <- survivalROC(Stime        = te1_clin_new$OS_month,
                       status       = te1_clin_new$Death_event,
                       marker       = te1_clin_new$PI,
                       predict.time = 36,
                       method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")

te1_roc5 <- survivalROC(Stime        = te1_clin_new$OS_month,
                       status       = te1_clin_new$Death_event,
                       marker       = te1_clin_new$PI,
                       predict.time = 60,
                       method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")

jpeg(file="Test1_ROC.jpeg", units="in", width=10, height=10, res=300)
plot(te1_roc1$FP, te1_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red", cex.lab=1.5, cex.axis=1.5, cex.main=2,
     ylab="TP",main="AUC Curve Test1 Data")
lines(te1_roc3$FP, te1_roc3$TP, type="l", lty=2, col="blue")
lines(te1_roc5$FP, te1_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, cex=1.5, legend=c( paste( "1 Year AUC = ",round(te1_roc1$AUC,2)),  paste("3 Years AUC = ",round(te1_roc3$AUC,2)),  paste("5 Years AUC = ",round(te1_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()


# create survival object
surv_object_p_te2 <- Surv(time = as.numeric(te2_clin_new$OS_month), event = te2_clin_new$Death_event)
#Cox multivariate model
#Cox multivariate model
#model_te2 <- coxph(surv_object_p_te2 ~ Gender + Grade + IDH_mutation_status +  codel_1p19q_status +  PI ,   data=te2_clin_new  )
#model_te2 <- coxph(surv_object_p_te2 ~ Gender + Grade + Class + IDH_mutation_status +  codel_1p19q_status + RT + chemo + PI  ,   data=te2_clin_new  )
model_te2 <- coxph(surv_object_p_te2 ~ Gender + Grade + Class + IDH_mutation_status +  codel_1p19q_status + PI  ,   data=te2_clin_new  )

summary(model_te2 )
#jpeg(file="test2_multivariate.jpeg", units="in", width=10, height=10, res=300)
#jpeg(file="te2_new_multivariate.jpeg", units="in", width=10, height=10, res=300)
jpeg(file="te2_new2_multivariate.jpeg", units="in", width=10, height=10, res=300)

ggforest(model_te2,data=te2_clin_new)
dev.off()





te2_roc1 <- survivalROC(Stime        = te2_clin_new$OS_month,
                        status       = te2_clin_new$Death_event,
                        marker       = te2_clin_new$PI,
                        predict.time = 12,
                        method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")


te2_roc3 <- survivalROC(Stime        = te2_clin_new$OS_month,
                        status       = te2_clin_new$Death_event,
                        marker       = te2_clin_new$PI,
                        predict.time = 36,
                        method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")

te2_roc5 <- survivalROC(Stime        = te2_clin_new$OS_month,
                        status       = te2_clin_new$Death_event,
                        marker       = te2_clin_new$PI,
                        predict.time = 60,
                        method       = "NNE", lambda = 0.0021, span = NULL, window ="symmetric")


jpeg(file="Test2_ROC.jpeg", units="in", width=10, height=10, res=300)
plot(te2_roc1$FP, te2_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red", cex.lab=1.5, cex.axis=1.5, cex.main=2,
     ylab="TP",main="AUC Curve Validation Data2")
lines(te2_roc3$FP, te2_roc3$TP, type="l", lty=2, col="blue")
lines(te2_roc5$FP, te2_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, cex=1.5, legend=c( paste( "1 Year AUC = ",round(te2_roc1$AUC,2)),  paste("3 Years AUC = ",round(te2_roc3$AUC,2)),  paste("5 Years AUC = ",round(te2_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()




######### COX regrwssion #####
