#https://rpubs.com/clayford/nomogram
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
library(rms)
library(SurvMetrics)
args <- commandArgs(TRUE)

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Radiation_responsive_GBM/")

set.seed(7)

#tr_clin_new <- read.table("tr_PI_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tr_clin_new <- read.table("Train_with_PI_and_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(tr_clin_new,2)
dim(tr_clin_new)
tr_clin_new$Age_group
# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_Status)


d <-  cbind(tr_clin_new["OS_month"], tr_clin_new["Death_Status"], tr_clin_new["Gender"],tr_clin_new["Age"], tr_clin_new["IDH_mutation_status"], tr_clin_new["Codel_1p19q_status"], tr_clin_new["MGMTp_meth_status"], tr_clin_new["PI"])
d_te <-  cbind(te_clin_new["OS_month"], te_clin_new["Death_Status"], te_clin_new["Grade"], te_clin_new["Gender"],te_clin_new["Age"], te_clin_new["IDH_mutation_status"], te_clin_new["Codel_1p19q_status"], te_clin_new["PI"])

dim(d)
ddist <- datadist(d)
options(datadist='ddist')

ddist_te <- datadist(d_te)
options(datadist='ddist_te')

### logistic regression models
f <- lrm(Death_Status ~ Age + Gender +  IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status + PI , data = d)
f1 <- lrm(Death_Status ~  PI , data = d,  x = T,y = T)
f2 <- lrm(Death_Status ~  PI + Grade , data = d,  x = T,y = T)
f3 <- lrm(Death_Status ~  PI + Grade +  IDH_mutation_status , data = d,  x = T,y = T)
f4 <- lrm(Death_Status ~  PI + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)
f5 <- lrm(Death_Status ~  PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)
f6 <- lrm(Death_Status ~  PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)

#nom <- nomogram(f, fun=plogis, funlabel="Risk of Death")
nom <- nomogram(f, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom1 <- nomogram(f1, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom2 <- nomogram(f2, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom3 <- nomogram(f3, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom4 <- nomogram(f4, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom5 <- nomogram(f5, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")

plot(nom)

jpeg("Nomogram_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom)
dev.off()

jpeg("Nomogram_PI_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom1)
dev.off()

jpeg("Nomogram_PI_Grade_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom2)
dev.off()

jpeg("Nomogram_PI_Grade_IDH_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom3)
dev.off()

jpeg("Nomogram_PI_Grade_IDH_codel_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom4)
dev.off()

jpeg("Nomogram_PI_Grade_IDH_codel_age_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom5)
dev.off()

f
Survival(f)
surv1<- Survival(f)

#coxph models
cox <-cph(Surv(OS_month,Death_Status==1) ~ Age + Gender +  IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  + PI,  x = T,y = T, data = d, surv = T)

cox1 <-cph(Surv(OS_month,Death_Status==1) ~ PI ,  x = T,y = T, data = d, surv = T)
cox2 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade ,  x = T,y = T, data = d, surv = T)
cox3 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status ,  x = T,y = T, data = d, surv = T)
cox4 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
cox5 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
cox6 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
cox7 <-cph(Surv(OS_month,Death_Status==1) ~ PI +  Grade +   Codel_1p19q_status,  x = T,y = T, data = d, surv = T)



#################
surv<- Survival(cox)
risk <- function(x)1/(1+exp(-x))
surv_1<- function(x)surv(1*12,lp=x) # defined time.inc,1 year OS
surv_2<- function(x)surv(1*36,lp=x) # defined time.inc,3 year OS
surv_3<- function(x)surv(1*60,lp=x) # defined time.inc,5 year OS

nom_cox<-nomogram(cox,fun = list(risk, surv_1,surv_2,surv_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom_cox),xfrac = .7)

jpeg("Nomogram_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom_cox)
dev.off()




###################
surv1<- Survival(cox1)
risk_1 <- function(x)1/(1+exp(-x))
surv1_1<- function(x)surv1(1*12,lp=x) # defined time.inc,1 year OS
surv1_2<- function(x)surv1(1*36,lp=x) # defined time.inc,3 year OS
surv1_3<- function(x)surv1(1*60,lp=x) # defined time.inc,5 year OS

nom1_cox<-nomogram(cox1,fun = list(risk_1, surv1_1,surv1_2,surv1_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom1_cox),xfrac = .7)

jpeg("Nomogram_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom1_cox)
dev.off()


surv2<- Survival(cox2)
risk_2 <- function(x)1/(1+exp(-x))
surv2_1<- function(x)surv2(1*12,lp=x) # defined time.inc,1 year OS
surv2_2<- function(x)surv2(1*48,lp=x) # defined time.inc,4 year OS
surv2_3<- function(x)surv2(1*60,lp=x) # defined time.inc,5 year OS

nom2_cox<-nomogram(cox2,fun = list(risk_2 , surv2_1,surv2_2,surv2_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom2_cox),xfrac = .7)


jpeg("Nomogram_PI_Grade_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom2_cox)
dev.off()


surv3<- Survival(cox3)
risk_3 <- function(x)1/(1+exp(-x))
surv3_1<- function(x)surv3(1*12,lp=x) # defined time.inc,1 year OS
surv3_2<- function(x)surv3(1*48,lp=x) # defined time.inc,4 year OS
surv3_3<- function(x)surv3(1*60,lp=x) # defined time.inc,5 year OS

nom3_cox<-nomogram(cox3,fun = list(risk_3 , surv3_1,surv3_2,surv3_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom3_cox),xfrac = .7)

jpeg("Nomogram_PI_grade_IDH_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom3_cox)
dev.off()


#nom4_C <- nomogram(cox4, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")


surv4<- Survival(cox4)
risk_4 <- function(x)1/(1+exp(-x))
surv4_1<- function(x)surv4(1*12,lp=x) # defined time.inc,1 year OS
surv4_2<- function(x)surv4(1*48,lp=x) # defined time.inc,4 year OS
surv4_3<- function(x)surv4(1*60,lp=x) # defined time.inc,5 year OS

nom4_cox<-nomogram(cox4,fun = list(risk_4, surv4_1,surv4_2,surv4_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom4_cox),xfrac = .7)

jpeg("Nomogram_PI_grade_IDH_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom4_cox)
dev.off()



surv5<- Survival(cox5)
risk_5 <- function(x)1/(1+exp(-x))
surv5_1<- function(x)surv5(1*12,lp=x) # defined time.inc,1 year OS
surv5_2<- function(x)surv5(1*48,lp=x) # defined time.inc,4 year OS
surv5_3<- function(x)surv5(1*60,lp=x) # defined time.inc,5 year OS

nom5_cox<-nomogram(cox5,fun = list(risk_5 , surv5_1,surv5_2,surv5_3),lp = F,
                   funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom5_cox),xfrac = .7)

jpeg("Nomogram_PI_grade_IDH_codel_age_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom5_cox)
dev.off()



surv6<- Survival(cox6)
risk6 <- function(x)1/(1+exp(-x))
surv6_1<- function(x)surv6(1*12,lp=x) # defined time.inc,1 year OS
surv6_2<- function(x)surv6(1*48,lp=x) # defined time.inc,4 year OS
surv6_3<- function(x)surv6(1*60,lp=x) # defined time.inc,5 year OS

nom6_cox<-nomogram(cox6,fun = list(risk6 , surv6_1,surv6_2,surv6_3),lp = F,
                   funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom6_cox),xfrac = .7)

jpeg("Nomogram_PI_grade_IDH_codel_age_gender_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom6_cox)
dev.off()


surv7<- Survival(cox7)
risk7 <- function(x)1/(1+exp(-x))
surv7_1<- function(x)surv7(1*12,lp=x) # defined time.inc,1 year OS
surv7_2<- function(x)surv7(1*48,lp=x) # defined time.inc,4 year OS
surv7_3<- function(x)surv7(1*70,lp=x) # defined time.inc,5 year OS

nom7_cox<-nomogram(cox7,fun = list(risk7 , surv7_1,surv7_2,surv7_3),lp = F,
                   funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom7_cox),xfrac = .7)

jpeg("Nomogram_PI_grade_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom7_cox)
dev.off()


########## C-index and p-value ######

#cox1$coefficients
f<-coxph(Surv(OS_month,Death_Status==1) ~  Age + Gender +  IDH_mutation_status + Codel_1p19q_status  + MGMTp_meth_status + PI , data = d)
f1_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI  , data = d)
f2_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade , data = d)
f3_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status  , data = d)
f4_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
f5_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
f6_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
f7_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  Codel_1p19q_status , data = d)
f8_c<-coxph(Surv(OS_month,Death_Status==1) ~  Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
f9_c<-coxph(Surv(OS_month,Death_Status==1) ~   Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)

sum.surv<-summary(f)
sum.surv1<-summary(f1_c)
sum.surv2<-summary(f2_c)
sum.surv3<-summary(f3_c)
sum.surv4<-summary(f4_c)
sum.surv5<-summary(f5_c)
sum.surv6<-summary(f6_c)
sum.surv7<-summary(f7_c)
sum.surv8<-summary(f8_c)
sum.surv9<-summary(f9_c)
sum.surv
sum.surv1
sum.surv2
sum.surv3
sum.surv4
sum.surv5
sum.surv6
sum.surv7
sum.surv8
sum.surv9
c_index<-sum.surv$concordance
c_index1<-sum.surv1$concordance
c_index2<-sum.surv2$concordance
c_index3<-sum.surv3$concordance
c_index4<-sum.surv4$concordance
c_index5<-sum.surv5$concordance
c_index6<-sum.surv6$concordance
c_index7<-sum.surv7$concordance
c_index8<-sum.surv8$concordance
c_index9<-sum.surv9$concordance
c_index
c_index1
c_index2
c_index3
c_index4
c_index5
c_index6
c_index7
c_index8
c_index9

################### Test data ########
########## C-index and p-value ######

#cox1$coefficients
f1_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI  , data = d_te)
f2_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade , data = d_te)
f3_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status  , data = d_te)
f4_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d_te)
f5_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d_te)
f6_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d_te)
f7_c_te<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  Codel_1p19q_status , data = d_te)
f8_c_te<-coxph(Surv(OS_month,Death_Status==1) ~  Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d_te)
f9_c_te<-coxph(Surv(OS_month,Death_Status==1) ~   Grade +  IDH_mutation_status + Codel_1p19q_status , data = d_te)

sum.surv1_te<-summary(f1_c_te)
sum.surv2_te<-summary(f2_c_te)
sum.surv3_te<-summary(f3_c_te)
sum.surv4_te<-summary(f4_c_te)
sum.surv5_te<-summary(f5_c_te)
sum.surv6_te<-summary(f6_c_te)
sum.surv7_te<-summary(f7_c_te)
sum.surv8_te<-summary(f8_c_te)
sum.surv9_te<-summary(f9_c_te)
sum.surv1_te
sum.surv2_te
sum.surv3_te
sum.surv4_te
sum.surv5_te
sum.surv6_te
sum.surv7_te
sum.surv8_te
sum.surv9_te
c_index1_te<-sum.surv1_te$concordance
c_index2_te<-sum.surv2_te$concordance
c_index3_te<-sum.surv3_te$concordance
c_index4_te<-sum.surv4_te$concordance
c_index5_te<-sum.surv5_te$concordance
c_index6_te<-sum.surv6_te$concordance
c_index7_te<-sum.surv7_te$concordance
c_index8_te<-sum.surv8_te$concordance
c_index9_te<-sum.surv9_te$concordance
c_index1_te
c_index2_te
c_index3_te
c_index4_te
c_index5_te
c_index6_te
c_index7_te
c_index8_te
c_index9_te

######################### ROC plots  #######

library(ROCR) 
cox7$sformula
cox2$linear.predictors
predvalue2 <- predict(cox2)
pred2 <- prediction(predvalue2 , d$Death_Status)
pred2
auc <- performance(pred2,"auc")
auc

perf<- performance(pred2,"tpr","fpr")
perf
plot(perf)
abline(0,1, col = 3, lty = 2)
##Create ROC plot for 1-,3- and 5-years survival time

tr_roc1 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
tr_roc1

tr_roc3 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc3
tr_roc5 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5


tr_roc6
jpeg(file="ROC_train_PI_Grade.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction based on PI+Grade Model")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("4 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()




##########################

?calibrate
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u = 12,  B = 100)
plot(cal1)
cal1 <- calibrate(f1, cmethod="KM", method="boot",  u = 365, m = 50, B = 100)
plot(cal1)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 1-Year",
     ylab="Actual Probability of 1-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))


surv1<- Survival(cox1)
surv2<- Survival(cox2)
#surv1_1<- function(x)surv(1*12,lp=x) # defined time.inc,1 year OS
#surv2_2<- function(x)surv(1*36,lp=x) # defined time.inc,3 year OS
#surv3_3<- function(x)surv(1*60,lp=x) # defined time.inc,5 year OS

nom2 <- nomogram(cox2, fun=list(function(x) surv2(12, x),
                                function(x) surv2(36, x),
                                function(x) surv2(60, x),
                                function(x) surv2(120, x)),
                 funlabel=c("1-Year Survival Probability", 
                            "3-Year Survival Probability",
                            "5-Year Survival Probability",
                            "10-Year Survival Probability"))
plot(nom2, xfrac=.7)

nom_cox1<-nomogram(cox1,fun = list(surv1,surv2,surv3),lp = F,funlabel = c("1-Year OS", "3-Year OS","5-Year OS)",maxscale = 100, fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')))
nom_cox1
plot(nom_cox1)

nom_cox2<-nomogram(cox2,fun=list(function(x) surv2(12, x),
                                 function(x) surv2(36, x),
                                 function(x) surv2(60, x),
                                 function(x) surv2(120, x)),
                   funlabel = c("1-Year OS", "3-Year OS","5-Year OS", "10-Year OS"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1') )



nom_cox2
plot(nom_cox2, xfrac=.7)


predvalue <- predict(cox1) 
predvalue2 <- predict(cox2) 
predvalue2


head(d)

?cindex
Cindex(Surv(OS_month, Death_status) ~ Age + Grade +  IDH_mutation_status + Codel_1p19q_status + PI, data = d)
pec::cindex(f_2, Surv(OS_month, Death_Status ) ~  + Grade + PI, data = d)
#concordance.index(Surv(OS_month, Death_status) ~ Age + Grade +  IDH_mutation_status + Codel_1p19q_status + PI, data = d)

plot(nom)

jpeg("Nomogram.jpg", units="in", width=10, height=10, res=300)
plot(nom)
dev.off()



####### Nomogram without Age #######
f1 <- lrm(Death_Status ~  Grade +  IDH_mutation_status + Codel_1p19q_status + PI , data = d)

?lrm
ddist1 <- datadist(d)
options(datadist='ddist1')


nom1 <- nomogram(f1, fun=plogis, funlabel="Risk of Death")
nom1
plot(nom1)

f1$coef[1]
f1$coefficients

jpeg("Nomogram_without_age.jpg", units="in", width=10, height=10, res=300)
plot(nom1)
dev.off()



# predict probability of Death based on vlue

Predict(f, IDH_mutation_status="Mutant", Grade='WHO IV', Codel_1p19q_status="Non-codel", PI=-1.3, fun=plogis)


lrm(Death_Status ~  Grade +  IDH_mutation_status + Codel_1p19q_status + PI , data = d)

head(d)
f2 <- psm(Surv(OS_month,Death_Status) ~  Grade +  IDH_mutation_status + Codel_1p19q_status +PI, data=d, dist='lognormal')
med  <- Quantile(f2)
surv1 <- Survival(f2)  # This would also work if f was from cph
surv1
plot(nomogram(f2, fun=function(x) med(lp=x), funlabel="Median Survival Time"))
nom2 <- nomogram(f2, fun=list(function(x) surv1(12, x),
                              function(x) surv1(36, x),
                              function(x) surv1(60, x),
                              function(x) surv1(120, x)),
                 funlabel=c("12-Month Survival Probability", 
                            "36-month Survival Probability",
                            "60-month Survival Probability",
                            "120-month Survival Probability"))
plot(nom2, xfrac=.7)

f2$score
f2$dist
f2$non.slopes
surv1(12, x)

f2$scale.pred



###########################
surv_object3 <- Surv(time = data3$OS_month, event = data3$Death_Status)
fit3 <- coxph(surv_object3 ~ Grade  +PI, data=data3)

fit3 <- lrm(Death_Status ~  Grade +  IDH_mutation_status + Codel_1p19q_status + PI , data = d)

f$
  # type of predicted value
  predict(fit3,type="lp") #linear predictor ("lp")
predict(fit3,type="expected") #expected number of events given the covariates and follow-up time 
predict(fit3,type="risk",se.fit=TRUE) # risk score exp(lp)
predict(fit3,type="terms",se.fit=TRUE) #terms of the linear predictor
predict(fit3, newdata=new_data, type="survival",se.fit=TRUE) #terms of the linear predictor



tr_CI<-concordance.index(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, method="noether")
#### x=predictions

tr_CI

D_index_tr<-D.index(x=train$pred1, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)
D_index_tr

HR_tr<-hazard.ratio(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)

HR_tr



plot(survfit(f, newdata=d,  xscale=365.25, xlab="Years", ylab="Survival", conf.int=F) )
# also plot the predicted survival for a 70 year old
lines(survfit(fit3, newdata=new_data), xscale=365.25, xlab="Years", ylab="Survival") 
