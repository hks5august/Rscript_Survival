library(rms)
library(caret)

set.seed(7)

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/")

train_data <- read.table("tr_PI_data_with_clin_without_peadetric_sample3.txt", header=T, sep="\t", row.names=1, check.names = F)
#train_data <- read.table("Train_Clin_data_with_11genes_PI.txt", header=T, sep="\t", row.names=1, check.names = F)

head(train_data,2)
test_data <- read.table("External_validation/Final_ext_test_with_PI_data2.txt", header=T, sep="\t", row.names=1, check.names = F)
head(test_data,2)

?cph
cox_PI <- cph(Surv(OS_month, Death_Status) ~ PI, data = train_data, x = T, y = T)
cox_PI


train_C <- as.numeric((cox_PI$stats[9] + 1)/2)
round(train_C, 2) # 

train_cv  <- validate(cox_PI, method = "crossvalidation", B = 5)
?validate
train_cv 


train_CV_C <- (train_cv[1, 3] + 1)/2
round(train_CV_C, 2) # 0.7167
train_CV_C


#'Test' Concordance:
  
actuals <- Surv(test_data$OS_month, test_data$Death_event)
estimates <- survest(cox_PI, newdata = test_data, times = 60)$surv # can pick an arbitrary time here, same results
test_C <- as.numeric(rcorr.cens(x = estimates, S = actuals)[1])
round(test_C, 2) 



########## 
head(train_data,2)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ Age, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ Gender, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ Grade, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ IDH_mutation_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ Codel_1p19q_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ MGMTp_meth_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Grade, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Grade + IDH_mutation_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Grade + IDH_mutation_status + Codel_1p19q_status, data = train_data, x = T, y = T)

cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Grade + IDH_mutation_status + Codel_1p19q_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Grade + IDH_mutation_status + Codel_1p19q_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Age + Grade + IDH_mutation_status + Codel_1p19q_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI + Gender + Age + Grade + IDH_mutation_status + Codel_1p19q_status, data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ PI +  Age  , data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~ Age  , data = train_data, x = T, y = T)
cox_2 <- cph(Surv(OS_month, Death_Status) ~   Age + Grade + Gender , data = train_data, x = T, y = T)

cox_2 <- cph(Surv(OS_month, Death_Status) ~   Gender + Age + Grade + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status , data = train_data, x = T, y = T)

cox_2 <- cph(Surv(OS_month, Death_Status) ~   Age + Grade + IDH_mutation_status + Codel_1p19q_status , data = train_data, x = T, y = T)

cox_2 <- cph(Surv(OS_month, Death_Status) ~   Grade + IDH_mutation_status + Codel_1p19q_status , data = train_data, x = T, y = T)

#cox_2 <- cph(Surv(OS_month, Death_Status) ~  ATOH8  + CP  + CRLF1  + DRAXIN + GLIS1 + H19  + HIST1H2AG +` RP11-189B4.6` + `RP11-300M24.1` + TNFSF14 , data = train_data, x = T, y = T)

head(train_data,2)



cox_2
an <- anova(cox_2 )

plot(Predict(cox_2), anova=an, pval=TRUE)
# Plot effects of all 4 predictors with test statistics from anova, and P
jpeg(file="Age_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="Gender_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="Grade_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="IDH_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="Codel_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="MGMT_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Grade_Annova.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Grade_IDH.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Grade_IDH_codel.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Grade_IDH_codel_age.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="PI_Grade_IDH_codel_age_gender.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="Grade_IDH_codel_age_gender.jpeg", units="in", width=10, height=10, res=350)

plot(Predict(cox_2), anova=an, pval=TRUE)
dev.off()


train_C2 <- as.numeric((cox_2$stats[9] + 1)/2)
round(train_C2, 2) # 

train_cv2  <- validate(cox_2, method = "crossvalidation", B = 5)

#surv_obj <- Surv(train_data$OS_month, train_data$Death_Status)
#fit <- glmnet(Surv(OS_month, Death_Status), Gender + Age + Grade + IDH_mutation_status + Codel_1p19q_status +MGMTp_meth_status, family = "cox", n=5, data=train_data)
#fit <- glmnet(surv_obj,  train_data$Sex1 + train_data$Age + train_data$Tumor_Grade + train_data$IDH_mut + train_data$codel_1p19q + train_data$MGMTp_meth, family = "cox", n=5, data=train_data, na.rm=T)

train_cv2

train_CV_C2 <- (train_cv2[1, 3] + 1)/2
round(train_CV_C2, 2) # 0.7167
train_CV_C2




#'Test' Concordance:

actuals2 <- Surv(test_data$OS_month, test_data$Death_event)
estimates2 <- survest(cox_2, newdata = test_data, times = 60)$surv # can pick an arbitrary time here, same results
estimates2
rcorr.cens(x = estimates2, S = actuals2)
test_C2 <- as.numeric(rcorr.cens(x = estimates2, S = actuals2)[1])
round(test_C2, 2) 






#Calibration curve
#fit1<-lrm( Death_Status ~  Gender + Age + Grade + IDH_mutation_status + Codel_1p19q_status +MGMTp_meth_status , data = train_data, x = T, y = T) 

fit1<-lrm( Survival_class ~ PI , data = train_data, x = T, y = T) 
fit1<-lrm( Death_Status ~ PI , data = train_data, x = T, y = T) 


cal2 <- calibrate(fit1, method ="boot", B = 100)
plot(cal2,xlim = c(0,1.0),ylim = c(0,1.0))









################ Validation on 11 sets ########
setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/New_11_sets/set10_survival/")

train_data1 <- read.table("Train_PI_Clin", header=T, sep="\t", row.names=1, check.names = F)
head(train_data1,2)
test_data1 <- read.table("Test_PI_Clin", header=T, sep="\t", row.names=1, check.names = F)
head(test_data1,2)

?cph
cox_PI_1 <- cph(Surv(OS_month, Death_Status) ~ PI, data = train_data1, x = T, y = T)
cox_PI_1


train_C_1 <- as.numeric((cox_PI_1$stats[9] + 1)/2)
round(train_C_1, 2) # 

train_cv_1  <- validate(cox_PI_1, method = "crossvalidation", B = 5)

train_cv

train_CV_C1 <- (train_cv_1[1, 3] + 1)/2
round(train_CV_C1, 2) # 0.7167
train_CV_C1


#'Test' Concordance:

actuals1 <- Surv(test_data1$OS_month, test_data1$Death_event)
estimates1 <- survest(cox_PI_1, newdata = test_data1, times = 60)$surv # can pick an arbitrary time here, same results
test_C1 <- as.numeric(rcorr.cens(x = estimates1, S = actuals1)[1])
round(test_C1, 2) 


