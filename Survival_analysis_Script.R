
#install.packages(c("dplyr", "ggplot2",  "ggfortify", "survival", "survminer", "pca3d")
#install.packages("rgl")

#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

#library(pca3d)


setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/GBM/Primary_GBM/survival")

set.seed(7)

Primary_tr_GBM_data <- read.table("Final_tr_clin_23171exp_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

dim(Primary_tr_GBM_data)

#remove samples where OS.time is NA

head(Primary_tr_GBM_data[1:24],3)
head(Primary_tr_GBM_data[18:24],3)

gbm_clin <- Primary_tr_GBM_data[1:24]
gbm_exp <- Primary_tr_GBM_data[25:ncol(Primary_tr_GBM_data)]


#Data Exploration using PCA 

#PCA object
pca_res1 <- prcomp(gbm_exp )

#extract variance explained by PCA componnents
var_explained1 <- pca_res1$sdev^2/sum(pca_res1$sdev^2)


#barplot for first 10 components
jpeg(file="PCA_bar.jpeg", units="in", width=10, height=10, res=300)
barplot(var_explained1[1:10], xlab="PCA components", ylab="Variability")
dev.off()

jpeg(file="subtype_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="PRS_type")
dev.off()


jpeg(file="Gender_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Gender")
dev.off()

jpeg(file="Age_grpup_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Age_group")
dev.off()

jpeg(file="IDH1_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="IDH_mutation_status")
dev.off()

jpeg(file="Codel_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="codel_1p19q_status")
dev.off()


jpeg(file="mgmt_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="MGMTp_meth_status")
dev.off()

jpeg(file="RT_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="RT")
dev.off()

jpeg(file="Chemo_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Chemo_TMZ")
dev.off()

jpeg(file="Both_trt_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Trt_both")
dev.off()

jpeg(file="Survival_class_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Survival_class")
dev.off()

jpeg(file="Grade_PCA.jpeg", units="in", width=10, height=10, res=300)
autoplot(pca_res1, data=Primary_tr_GBM_data, colour="Grade")
dev.off()



#Survival analysis

GBM_p_surv_data<-subset(Primary_tr_GBM_data ,OS_month!="NA")
dim(Primary_tr_GBM_data)
dim(GBM_p_surv_data)

# create survival object
surv_object_p <- Surv(time = as.numeric(gbm_clin$OS_month), event = gbm_clin$Death_event)


#Cox multivariate model
cox_multivariate <- coxph(surv_object_p ~  Gender + as.numeric(Age) + IDH_mutation_status + codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ ,  data=gbm_clin )
summary(cox_multivariate )
jpeg(file="multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=gbm_clin )
ggforest(cox_multivariate,data=gbm_clin)
dev.off()

covariates <- c( "Gender" ,  "Age", "IDH_mutation_status", "codel_1p19q_status", "MGMTp_meth_status", "RT" , "Chemo_TMZ")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('surv_object_p ~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = gbm_clin )})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #p.value<-signif(x$wald["pvalue"], digits=2) 
                         #wald<-signif(x$wald["test"], digits=2)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         CI1 <- round(as.numeric(signif(x$concordance[1])),2);#concordance
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, CI1, logrank, p.value)
                         #res<-c(beta, HR,  wald, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "CI", "logrank",  "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res_uni <- t(as.data.frame(univ_results))

write.table(res_uni , file = "Univariate_results_traing_data.txt", sep="\t", quote=F, row.names = T)


# Gender based
fit_gender_p = survfit(surv_object_p ~ Gender, data=gbm_clin)

?ggsurvplot
# we produce a Kaplan Meier plot
jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_gender_p, data=gbm_clin, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()
fit.coxph_gender_p <- coxph(surv_object_p ~ Gender, data=gbm_clin)
coef(summary(fit.coxph_gender_p))



#remove samples with no IDH mutation information
GBM_data_idh_p<-subset(GBM_p_surv_data, IDH_mutation_status!="NA")

#Create survival object
surv_object_idh_p <- Surv(time = as.numeric(GBM_data_idh_p$OS_month), event = GBM_data_idh_p$Death_event)

#fit model
fit_idh_p = survfit(surv_object_idh_p ~ IDH_mutation_status, data=GBM_data_idh_p)

# Kaplan Meier plot
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="IDH1mut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_idh_p <- coxph(surv_object_p ~ IDH_mutation_status, data=GBM_p_surv_data)
coef(summary(fit.coxph_idh_p))

#remove samples with no MGMT information
GBM_data_mgmt_p<-subset(GBM_p_surv_data, MGMTp_meth_status!="NA")

#Create survival object
surv_object_mgmt_p <- Surv(time = as.numeric(GBM_data_mgmt_p$OS_month), event = GBM_data_mgmt_p$Death_event)

#fit model
fit_mgmt_p = survfit(surv_object_mgmt_p ~ MGMTp_meth_status, data=GBM_data_mgmt_p)

# Kaplan Meier plot
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="mgmt.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_mgmt_p <- coxph(surv_object_p ~ MGMTp_meth_status, data=GBM_p_surv_data)
coef(summary(fit.coxph_mgmt_p))



#remove samples with no codeletion mutation information
GBM_data_codel_p<-subset(GBM_p_surv_data, codel_1p19q_status!="NA")

#Create survival object
surv_object_codel_p <- Surv(time = as.numeric(GBM_data_codel_p$OS_month), event = GBM_data_codel_p$Death_event)

#fit model
fit_codel_p = survfit(surv_object_codel_p ~ codel_1p19q_status, data=GBM_data_codel_p)

# Kaplan Meier plot
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="codelmut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_codel_p <- coxph(surv_object_codel_p ~ codel_1p19q_status, data=GBM_data_codel_p)
coef(summary(fit.coxph_codel_p))

#remove samples with no RT information
GBM_data_RT_p<-subset(GBM_p_surv_data, RT!="NA")

#Create survival object
surv_object_RT_p <- Surv(time = as.numeric(GBM_data_RT_p$OS_month), event = GBM_data_RT_p$Death_event)

#fit model
fit_RT_p = survfit(surv_object_RT_p ~ RT, data=GBM_data_RT_p)

# Kaplan Meier plot
ggsurvplot(fit_RT_p, data=GBM_data_RT_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")

jpeg(file="RT.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_RT_p, data=GBM_data_RT_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_RT_p <- coxph(surv_object_RT_p ~ RT, data=GBM_data_RT_p)
coef(summary(fit.coxph_RT_p))

#remove samples with no Chemo information
GBM_data_chemo_p<-subset(GBM_p_surv_data, Chemo_TMZ!="NA")

#Create survival object
surv_object_chemo_p <- Surv(time = as.numeric(GBM_data_chemo_p$OS_month), event = GBM_data_chemo_p$Death_event)

#fit model
fit_chemo_p = survfit(surv_object_chemo_p ~ Chemo_TMZ, data=GBM_data_chemo_p)

# Kaplan Meier plot
ggsurvplot(fit_chemo_p, data=GBM_data_chemo_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")

jpeg(file="Chemo.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_chemo_p, data=GBM_data_chemo_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_chemo_p <- coxph(surv_object_chemo_p ~ Chemo_TMZ, data=GBM_data_chemo_p)
coef(summary(fit.coxph_chemo_p))


#remove samples with no Chemo information
GBM_data_both_p<-subset(GBM_p_surv_data, Trt_both!="NA")

#Create survival object
surv_object_both_p <- Surv(time = as.numeric(GBM_data_both_p$OS_month), event = GBM_data_both_p$Death_event)

#fit model
fit_both_p = survfit(surv_object_both_p ~ Trt_both, data=GBM_data_both_p)

# Kaplan Meier plot
ggsurvplot(fit_both_p, data=GBM_data_both_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")

jpeg(file="Both_trt.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_both_p, data=GBM_data_both_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_both_p <- coxph(surv_object_both_p ~ Trt_both, data=GBM_data_both_p)
coef(summary(fit.coxph_both_p))


#remove samples with no Survival class information
GBM_data_survC_p<-subset(GBM_p_surv_data, Survival_class!="NA")

#Create survival object
surv_object_survC_p <- Surv(time = as.numeric(GBM_data_survC_p$OS_month), event = GBM_data_survC_p$Death_event)

#fit model
fit_survC_p = survfit(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)

#summary(fit_survC_p)

# Kaplan Meier plot
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="IDH1mut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_survC_p <- coxph(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)
coef(summary(fit.coxph_survC_p))


#GBM_p_surv_data_age <- subset(GBM_p_surv_data, Age_group!="pediatric")
GBM_p_surv_data_age = GBM_p_surv_data
#Create survival object
surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_event)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)
coef(summary(fit.coxph_surv_age_p))


surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_event)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 6,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age, data=GBM_p_surv_data)
coef(summary(fit.coxph_surv_age_p))

clin_factors_Coeff <- rbind(coef(summary(fit.coxph_gender_p)), coef(summary(fit.coxph_surv_age_p)), coef(summary(fit.coxph_idh_p)), coef(summary(fit.coxph_codel_p)),  coef(summary(fit.coxph_mgmt_p)), coef(summary(fit.coxph_RT_p)), coef(summary(fit.coxph_chemo_p)), coef(summary(fit.coxph_both_p)))
clin_factors_Coeff 


#Write into files
write.table(clin_factors_Coeff, file = "clin_factors_Coeff_train_COXPH.txt", sep="\t", quote=F, row.names = T)



dim(GBM_p_surv_data)
head(GBM_p_surv_data[17],2)


# create a file to store results
write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="Primary_GBM_tr_Survival_results_for_genes.csv",row.names=F,col.names=F,sep = ',');

#Here features to compute survival start from 17th column onwards, which is transcriptomics features

for(i in seq(from=18, to=length(GBM_p_surv_data), by=1))
{
  surv_object2 <- Surv(time = GBM_p_surv_data$OS_month, event = GBM_p_surv_data$Death_event)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object2 ~ (GBM_p_surv_data[,i])>(median(GBM_p_surv_data[1,i])), data=GBM_p_surv_data);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (GBM_p_surv_data[,i])>(median(GBM_p_surv_data[1,i])), data = GBM_p_surv_data)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  #check whether the pvalue is significant (< 0.05) or not
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(GBM_p_surv_data[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="Primary_GBM_tr_Survival_results_for_genes.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}


Surv_res_p <- read.table("Primary_GBM_tr_Survival_results_for_genes.csv", sep=",", header=TRUE, row.names=1)
dim(Surv_res_p)
head(Surv_res_p)

# Sort features based on Corcordance
CI_features_p <- round(as.data.frame(Surv_res_p[order(Surv_res_p[,"Concordance"], decreasing =  TRUE),]),3)

write.table(CI_features_p , file = "CI_features_p.txt", sep="\t", quote=F, row.names = T)

dim(CI_features_p)
head(CI_features_p, 10)


#Top significant prognostic features
sig_surv_features_p <- CI_features_p[CI_features_p$Concordance >= 0.545, ]
dim(sig_surv_features_p)
head(sig_surv_features_p)

survival_unfavorable_p <- sig_surv_features_p[sig_surv_features_p$HR >= 1.5, ]

dim(survival_unfavorable_p)
head(survival_unfavorable_p, 10)

#Good Prognostic features based on inverse of HR >= 1.5
survival_favorable_p <- sig_surv_features_p[sig_surv_features_p$Hr.Inv.lst >= 1.5, ]
dim(survival_favorable_p)
head(survival_favorable_p,10)


Survival_168fav_377_unfav_features <- rbind(survival_favorable_p,survival_unfavorable_p)
dim(Survival_168fav_377_unfav_features)



#Write into files
write.table(Survival_168fav_377_unfav_features, file = "Primary_GBM_Survival_168fav_377_unfav_features_pval_0_05_CI_55.txt", sep="\t", quote=F, row.names = T)


#Top significant prognostic features
sig_surv_features_p2 <- CI_features_p[CI_features_p$Concordance >= 0.595, ]
dim(sig_surv_features_p2)
head(sig_surv_features_p2)

write.table(sig_surv_features_p2 , file = "sig_surv_features_CI_60.txt", sep="\t", quote=F, row.names = T)


survival_unfavorable_p2 <- sig_surv_features_p2[sig_surv_features_p2$HR >= 1.5, ]

dim(survival_unfavorable_p2)
head(survival_unfavorable_p2, 16)

#Good Prognostic features based on inverse of HR >= 1.5
survival_favorable_p2 <- sig_surv_features_p2[sig_surv_features_p2$Hr.Inv.lst >= 1.5, ]
dim(survival_favorable_p2)
head(survival_favorable_p2,16)


Survival_6fav_16unfav_features <- rbind(survival_favorable_p2,survival_unfavorable_p2)
dim(Survival_6fav_16unfav_features)


#Write into files
write.table(Survival_6fav_16unfav_features, file = "Primary_GBM_Survival_6fav_16unfav_features.txt", sep="\t", quote=F, row.names = T)

#Good Prognostic features based on inverse of HR >= 1.5
survival_favorable_p3 <- sig_surv_features_p2[sig_surv_features_p2$Hr.Inv.lst >= 1.95, ]
dim(survival_favorable_p3)
survival_unfavorable_p3 <- sig_surv_features_p2[sig_surv_features_p2$HR >= 1.95, ]
dim(survival_unfavorable_p3 )


Survival_5fav_10unfav_features <- rbind(survival_favorable_p3,survival_unfavorable_p3)
dim(Survival_5fav_10unfav_features )
Survival_5fav_10unfav_features 


#Write into files
write.table(Survival_5fav_10unfav_features, file = "Survival_5fav_10unfav_feature_HR_2_CI_6.txt", sep="\t", quote=F, row.names = T)


# create survival object
surv_object_p <- Surv(time = as.numeric(GBM_p_surv_data$OS_month), event = GBM_p_surv_data$Death_event)


#Cox multivariate model
cox_multivariate1 <- coxph(surv_object_p ~  Gender + as.numeric(Age) + IDH_mutation_status + 
                             codel_1p19q_status + MGMTp_meth_status + Chemo_TMZ + RT +
                             as.numeric(AC073133.1) + as.numeric(AC140481.1 ) + as.numeric(`CTC-329D1.3`) + as.numeric(TSSC2 )
                           + as.numeric(DTX4)  + as.numeric(TSSC2 )  + as.numeric(AC138430.4)  + as.numeric(LARS2)
                           + as.numeric(FAM127A )  + as.numeric(GTPBP2)  + as.numeric(MIR3654) + as.numeric(ARSK)
                           + as.numeric(RPS4XP16)  + as.numeric(PISD)  + as.numeric(MPPED2) + as.numeric(AP001625.6 )  ,  
                           data=GBM_p_surv_data )
summary(cox_multivariate1 )
ggforest(cox_multivariate1, data=GBM_p_surv_data)

cox_multivariate2 <- coxph(surv_object_p ~  Gender + as.numeric(Age) + IDH_mutation_status + 
                             RT + Chemo_TMZ + 
                             as.numeric(AC073133.1) + as.numeric(AC140481.1 ) + as.numeric(`CTC-329D1.3`) + as.numeric(TSSC2 )
                           + as.numeric(DTX4)  + as.numeric(TSSC2 )  + as.numeric(AC138430.4)  + as.numeric(LARS2)
                           + as.numeric(FAM127A )  + as.numeric(GTPBP2)  + as.numeric(MIR3654)  + as.numeric(ARSK)
                           + as.numeric(RPS4XP16)  + as.numeric(PISD)  + as.numeric(MPPED2) + as.numeric(AP001625.6 ) ,  
                           data=GBM_p_surv_data )
  
  summary(cox_multivariate2 )
  ggforest(cox_multivariate2, data=GBM_p_surv_data)
  
  
  
  
  cox_multivariate3 <- coxph(surv_object_p ~  Gender + as.numeric(Age) + IDH_mutation_status + 
                               codel_1p19q_status + MGMTp_meth_status +
                               RT + Chemo_TMZ +  as.numeric(AC140481.1 ) + as.numeric(TSSC2 )+ as.numeric(LARS2)
                            + as.numeric(RPS4XP16)  ,   data=GBM_p_surv_data )
  
  summary(cox_multivariate3 )
  ggforest(cox_multivariate3, data=GBM_p_surv_data)
  
  
  

jpeg(file="Train_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=gbm_clin )
ggforest(cox_multivariate2,data=gbm_clin)
dev.off()

Primary_te_GBM_data <- read.table("Final_te_clin_23171exp_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
Primary_ext_GBM_data <- read.table("Final_ext_clin_23171exp_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)




GBM_p_surv_data_te<-subset(Primary_te_GBM_data ,OS_month!="NA")
dim(GBM_p_surv_data_te)
head(GBM_p_surv_data_te[1:10],2)


#Cox multivariate model

# create survival object
surv_object_p_te <- Surv(time = as.numeric(GBM_p_surv_data_te$OS_month), event = GBM_p_surv_data_te$Death_event)#Create survival object

cox_multivariate_te1 <- coxph(surv_object_p_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  MGMTp_meth_status +
                                 as.numeric(AC073133.1) + as.numeric(AC140481.1 ) + as.numeric(`CTC-329D1.3`) + as.numeric(TSSC2 )
                              + as.numeric(DTX4)  + as.numeric(TSSC2 )  + as.numeric(AC138430.4)  + as.numeric(LARS2)
                              + as.numeric(FAM127A )  + as.numeric(GTPBP2)  + as.numeric(MIR3654)  + as.numeric(ARSK)
                              + as.numeric(RPS4XP16)  +
                                as.numeric(PISD)  + as.numeric(MPPED2) + as.numeric(AP001625.6 ),  data=GBM_p_surv_data_te  )

summary(cox_multivariate_te1 )
ggforest(cox_multivariate_te1, data=GBM_p_surv_data_te )



cox_multivariate_te <- coxph(surv_object_p_te ~   as.numeric(Age) + IDH_mutation_status + 
                             RT + Chemo_TMZ + 
                             as.numeric(AC073133.1) + as.numeric(AC140481.1 ) + as.numeric(`CTC-329D1.3`) + as.numeric(TSSC2 )
                           + as.numeric(DTX4)  + as.numeric(TSSC2 )  + as.numeric(AC138430.4)  + as.numeric(LARS2)
                           + as.numeric(FAM127A )  + as.numeric(GTPBP2)  + as.numeric(MIR3654)  + as.numeric(ARSK)
                           + as.numeric(RPS4XP16)  + as.numeric(PISD)  + as.numeric(MPPED2) + as.numeric(AP001625.6 ) ,  
                           data=GBM_p_surv_data_te )

summary(cox_multivariate_te )
ggforest(cox_multivariate_te, data=GBM_p_surv_data_te)


cox_multivariate_te2 <- coxph(surv_object_p_te ~ Gender + as.numeric(Age) + IDH_mutation_status + 
  codel_1p19q_status + MGMTp_meth_status +
  RT + Chemo_TMZ +  as.numeric(AC140481.1 ) + as.numeric(TSSC2 )+ as.numeric(LARS2)
+ as.numeric(RPS4XP16) ,  
data=GBM_p_surv_data_te )


summary(cox_multivariate_te2 )
ggforest(cox_multivariate_te2, data=GBM_p_surv_data_te)


cox_multivariate_te3 <- coxph(surv_object_p_te ~  Gender + as.numeric(Age) + IDH_mutation_status + 
                                 RT + Chemo_TMZ +  as.numeric(AC140481.1 ) + as.numeric(TSSC2 )+ as.numeric(LARS2)
                               + as.numeric(RPS4XP16) , data=GBM_p_surv_data_te )


summary(cox_multivariate_te3 )
ggforest(cox_multivariate_te3, data=GBM_p_surv_data_te)




jpeg(file="Test_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=gbm_clin )
ggforest(cox_multivariate_te, data=GBM_p_surv_data_te)
dev.off()



GBM_p_surv_data_ext<-subset(Primary_ext_GBM_data ,OS_month!="NA")
dim(GBM_p_surv_data_ext)
head(GBM_p_surv_data_ext[1:10],2)

# create survival object
surv_object_p_ext <- Surv(time = as.numeric(GBM_p_surv_data_ext$OS_month), event = GBM_p_surv_data_ext$Death_event)#Create survival object

cox_multivariate_ext <- coxph(surv_object_p_ext ~   as.numeric(Age) + IDH_mutation_status + 
                               RT + Chemo_TMZ + 
                               as.numeric(AC073133.1) + as.numeric(AC140481.1 ) + as.numeric(`CTC-329D1.3`) + as.numeric(TSSC2 )
                             + as.numeric(DTX4)  + as.numeric(TSSC2 )  + as.numeric(AC138430.4)  + as.numeric(LARS2)
                             + as.numeric(FAM127A )  + as.numeric(GTPBP2)  + as.numeric(MIR3654)  + as.numeric(ARSK)
                             + as.numeric(RPS4XP16)  + as.numeric(PISD)  + as.numeric(MPPED2) + as.numeric(AP001625.6 ) ,  
                             data=GBM_p_surv_data_ext )

summary(cox_multivariate_ext )
ggforest(cox_multivariate_ext, data=GBM_p_surv_data_ext)


cox_multivariate_ext2 <- coxph(surv_object_p_ext ~   Gender + as.numeric(Age) + IDH_mutation_status + 
  RT + Chemo_TMZ +  as.numeric(AC140481.1 ) + as.numeric(TSSC2 )+ as.numeric(LARS2)
+ as.numeric(RPS4XP16) ,  data=GBM_p_surv_data_ext )

summary(cox_multivariate_ext2 )
ggforest(cox_multivariate_ext2, data=GBM_p_surv_data_ext)



jpeg(file="External_Test_multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=gbm_clin )
ggforest(cox_multivariate_ext, data=GBM_p_surv_data_ext)
dev.off()

