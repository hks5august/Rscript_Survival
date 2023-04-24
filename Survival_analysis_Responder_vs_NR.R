
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

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/")

#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")

set.seed(7)

#data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
data <- read.table("Primary_grades3_4_2_LTS_STS.clin2_expression.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
dim(data)
head(Clin_data)
Clin_data <- data[1:25]
head(data[1:30])

data1 <- subset(data, OS_month!="NA")
dim(data1)

# create a file to store results
write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');

#Here features to compute survival start from 16th column onwards, which is transcriptomics features

for(i in seq(from=25, to=length(data1), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
 
}


### Significant results 

write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');


for(i in seq(from=25, to=length(data), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(median(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  #check whether the pvalue is significant (< 0.05) or not
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
  }
}




##############
# create survival object
surv_object_p <- Surv(time = data1$OS_month, event = data1$Death_Status)
head(data1[1:13])

#Cox multivariate model
cox_multivariate <- coxph(surv_object_p ~  Gender + Grade +Group +  as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  ,  data=data1 )
summary(cox_multivariate )
jpeg(file="multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=data1 )
ggforest(cox_multivariate,data=data1)
dev.off()

covariates <- c( "Gender" , "Grade", "Group",  "Age", "IDH_mutation_status", "Codel_1p19q_status", "MGMTp_meth_status")
univ_formulas <- sapply(covariates,  function(x) as.formula(paste('surv_object_p ~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data1 )})
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

univ_results
univ_results
res_uni <- as.matrix(univ_results)
res_uni

write.table(res_uni , file = "Univariate_results_traing_data.txt", sep="\t", quote=F, row.names = T)


# Gender based
fit_gender_p = survfit(surv_object_p ~ Gender, data=data1)

?ggsurvplot
# we produce a Kaplan Meier plot
jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_gender_p, data=data1, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           legend.labs=c("Male", "Female"), legend.title="Gender",  
           palette=c("dodgerblue2", "red"), 
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()
fit.coxph_gender_p <- coxph(surv_object_p ~ Gender, data=data1)
coef(summary(fit.coxph_gender_p))





# Class - R/P based
fit_subtype_p = survfit(surv_object_p ~ Group, data=data1)

# we produce a Kaplan Meier plot
jpeg(file="Subtype.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_subtype_p, data=data1, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           legend.labs=c( "Non_Responder", "Responder"), legend.title="Group",  
           palette=c("dodgerblue2", "red"), 
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()
fit.coxph_subtype_p <- coxph(surv_object_p ~ Group, data=data1)
coef(summary(fit.coxph_subtype_p))



#remove samples with no IDH mutation information
GBM_data_idh_p<-subset(data1, IDH_mutation_status!="NA")

#Create survival object
surv_object_idh_p <- Surv(time = GBM_data_idh_p$OS_month, event = GBM_data_idh_p$Death_Status)

#fit model
fit_idh_p = survfit(surv_object_idh_p ~ IDH_mutation_status, data=GBM_data_idh_p)

# Kaplan Meier plot
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="IDH1mut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           legend.labs=c( "Mutant", "Wildtype"), legend.title="IDH1_mutation",  
           palette=c("dodgerblue2", "red"), 
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_idh_p <- coxph(surv_object_idh_p ~ IDH_mutation_status, data=GBM_data_idh_p)
coef(summary(fit.coxph_idh_p))

#remove samples with no MGMT information
GBM_data_mgmt_p<-subset(data1, MGMTp_meth_status!="NA")

#Create survival object
surv_object_mgmt_p <- Surv(time = as.numeric(GBM_data_mgmt_p$OS_month), event = GBM_data_mgmt_p$Death_Status)

#fit model
fit_mgmt_p = survfit(surv_object_mgmt_p ~ MGMTp_meth_status, data=GBM_data_mgmt_p)

# Kaplan Meier plot
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="mgmt.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           legend.labs=c( "Methylated", "Un-methylyalted"), legend.title="MGMTp methyl status",  
           palette=c("dodgerblue2", "red"), 
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_mgmt_p <- coxph(surv_object_mgmt_p ~ MGMTp_meth_status, data=GBM_data_mgmt_p)
coef(summary(fit.coxph_mgmt_p))



#remove samples with no codeletion mutation information
GBM_data_codel_p<-subset(data1, Codel_1p19q_status!="NA")

#Create survival object
surv_object_codel_p <- Surv(time = as.numeric(GBM_data_codel_p$OS_month), event = GBM_data_codel_p$Death_Status)

#fit model
fit_codel_p = survfit(surv_object_codel_p ~ Codel_1p19q_status, data=GBM_data_codel_p)

# Kaplan Meier plot
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="codelmut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           legend.labs=c( "Codel", "Non-Codel"), legend.title="1p19q Codeletion",  
           palette=c("dodgerblue2", "red"), 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_codel_p <- coxph(surv_object_codel_p ~ Codel_1p19q_status, data=GBM_data_codel_p)
coef(summary(fit.coxph_codel_p))


#remove samples with no Survival class information
GBM_data_survC_p<-subset(data1, Survival_class!="NA")

#Create survival object
surv_object_survC_p <- Surv(time = as.numeric(GBM_data_survC_p$OS_month), event = GBM_data_survC_p$Death_Status)

#fit model
fit_survC_p = survfit(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)

#summary(fit_survC_p)

# Kaplan Meier plot
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="Survival_class.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           legend.labs=c( "LTS", "STS"), legend.title="Survival Class",  
           palette=c("dodgerblue2", "red"), 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_survC_p <- coxph(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)
coef(summary(fit.coxph_survC_p))


#GBM_p_surv_data_age <- subset(GBM_p_surv_data, Age_group!="pediatric")
GBM_p_surv_data_age = subset(data1, Age!="NA")
#Create survival object
surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_Status)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           legend.labs=c( "Mid","old", "Pediatric", "Young"), 
           legend.title="Age Group",  
           palette=c("dodgerblue2", "red", "green", "Cyan"), 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)
coef(summary(fit.coxph_surv_age_p))


surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_Status)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age, data=GBM_p_surv_data_age)
coef(summary(fit.coxph_surv_age_p))

clin_factors_Coeff <- rbind(coef(summary(fit.coxph_gender_p)), coef(summary(fit.coxph_subtype_p)),  coef(summary(fit.coxph_surv_age_p)), coef(summary(fit.coxph_idh_p)), coef(summary(fit.coxph_codel_p)),  coef(summary(fit.coxph_mgmt_p)))
clin_factors_Coeff 


#Write into files
write.table(clin_factors_Coeff, file = "clin_factors_Coeff_train_COXPH.txt", sep="\t", quote=F, row.names = T)




######## cox modellasso model ###


#data1 <- read.table("282_Sig_genes_Quantile_mat_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data1 <- read.table("171_Sig_genes_Quantile_mat_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
cvfit1 <- cv.glmnet(as.matrix(data1[25:ncol(data1)]), 
                    surv_object2, # create survival object from the data
                    family = "cox", # specify Cox PH model
                    type.measure = "C", 
                    nfolds = 5, 
                    alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                    maxit = 1000)

lambda_min <- cvfit1$lambda.min
lambda_min

plot(cvfit1)


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
key_variables


write.table(key_variables,file="Lasso_key_variables.txt",row.names=T,col.names=T,sep = '\t', quote = F)



PI_model <- (AKR1C3 * (-0.148) + ATOH8 *(-0.127) + CP * (0.129) + CRLF1 * (-0.157) + DRAXIN * (0.325) + GLIS1 * (-0.48) + H19 *(0.0186) + HIST1H2AG * (0.2258) +` RP11-189B4.6`* (0.3927) + `RP11-300M24.1` * (0.1277) + TNFSF14 * (0.1389), data=data1)
PI_model <- data1$AKR1C3 * (-0.148) + data1$ATOH8 * (-0.127) + data1$CP * (0.129) +data1$CRLF1 * (-0.157) + data1$DRAXIN * (0.325) + data1$GLIS1 * (-0.48) +  data1$H19 *(0.0186) + data1$HIST1H2AG * (0.2258) +data1$`RP11-189B4.6`* (0.3927) + data1$`RP11-300M24.1` * (0.1277) + data1$TNFSF14 * (0.1389)


PI_model <- data1$AKR1C3 * (-0.148) + data1$ATOH8 * (-0.127) + data1$CP * (0.129) +data1$CRLF1 * (-0.157) + data1$DRAXIN * (0.325) + data1$GLIS1 * (-0.48) +  data1$H19 *(0.0186) + data1$HIST1H2AG * (0.2258) +data1$`RP11-189B4.6`* (0.3927) + data1$`RP11-300M24.1` * (0.1277) + data1$TNFSF14 * (0.1389)


data2 <-read.table("11_lasso_Sig_genes_Quantile_mat_data_with_clin_copy.txt",header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
dim(data2)
head(data2)
sel_features_results<-read.table("11_key_features.txt",header =TRUE, sep = "\t", check.names = FALSE)
head(sel_features_results)
dim(sel_features_results)
sel_ftr_surv <- as.data.frame(sel_features_results[,1])
names(sel_ftr_surv ) <- c("ID")
head(sel_ftr_surv )
dim(sel_ftr_surv )

#### prepare training, test and external validation data with selected features having significant value in univariate analysis ##########
sel_train <- as.data.frame(data2[,colnames(data2) %in% c(sel_ftr_surv$ID), ])
head(sel_train,2)
dim(sel_train)

data2$OS_month
######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data2["OS_month"],data2["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

PI_tr=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,2]))
}


LUSC_C_tr$PI<-PI_tr
head(LUSC_C_tr)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_with_PI.txt",sep="\t",quote=F, row.names=F)

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

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

tr_res1 <- cbind(coeff, HR1, CI, pval)
tr_res1
#save results as a file
write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res1.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file= paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg"), units="in", width=10, height=10, res=300)
#paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg")
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI, pval=TRUE,
                  risk.table=TRUE, tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  legend.labs=c( "Less than Mean of PI", "Greater than Mean of PI"), legend.title="PI Score",  
                  palette=c("dodgerblue2", "red"), 
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp
jpeg("PI_KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()


###############

# multivariate
#tr_clin_new <- read.table("tr_PI_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tr_clin_new <- read.table("tr_PI_data_with_clin_without_peadetric_sample3.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(tr_clin_new,2)
dim(tr_clin_new)
tr_clin_new$Age_group
# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_Status)
#Cox multivariate model
cox_multivariate_tr <- coxph(surv_object_p_tr ~    Gender + Group + Grade + as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  ,  data=tr_clin_new  )
cox_multivariate_tr <- coxph(surv_object_p_tr ~  Gender + Group + Grade + Age_group + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  + PI,  data=tr_clin_new  )

summary(cox_multivariate_tr )
ggforest(cox_multivariate_tr,data=tr_clin_new)
cox_multivariate_tr2 <- coxph(surv_object_p_tr ~  Grade + Age_group + IDH_mutation_status + Codel_1p19q_status + + PI,  data=tr_clin_new  )
summary(cox_multivariate_tr2 )
ggforest(cox_multivariate_tr2,data=tr_clin_new)

cox_multivariate_tr3 <- coxph(surv_object_p_tr ~  Sex + Grade + Age_group + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  ,  data=tr_clin_new  )
cox_multivariate_tr3 <- coxph(surv_object_p_tr ~  Tumor_Grade + Sex +   as.numeric(Age) + IDH_mut + codel_1p19q + MGMTp_meth ,  data=tr_clin_new  )

summary(cox_multivariate_tr3 )
jpeg(file="train_multivariate_only_clin.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_tr3,data=tr_clin_new, fontsize = 1.6, main = "Training Data")
dev.off()


data1$MGMTp_meth
cox_multivariate_tr4 <- coxph(surv_object_p_tr ~  Histology + Tumor_Grade + Sex +   as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status +PI  ,  data=tr_clin_new  )
cox_multivariate_tr4 <- coxph(surv_object_p_tr ~  Tumor_Grade + Sex +   as.numeric(Age) + IDH_mut + codel_1p19q + MGMTp_meth +PI  ,  data=tr_clin_new  )

summary(cox_multivariate_tr4 )
jpeg(file="train_multivariate_with_PI.jpeg", units="in", width=15, height=12, res=350)
ggforest(cox_multivariate_tr4,data=tr_clin_new, fontsize = 1.6, main = "Training Data")
dev.off()

jpeg(file="train_multivariate.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_tr,data=tr_clin_new, fontsize = 1.6, main = "Training Data")
dev.off()




#######################################################################################################

############# Responder v non-responders ###########

train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat2), train_feature_mat2),file="sel_train2.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_tr2=train_feature_mat2
tr2 <- LUSC_C_tr2[3:ncol(LUSC_C_tr2)]

E= length(tr2)
E

PI_tr2=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr2= PI_tr2+((tr2[,i])*(sel_features_results[i,2]))
}


LUSC_C_tr2$PI<-PI_tr2
head(LUSC_C_tr2)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_tr2), LUSC_C_tr2),file="train2_with_PI.txt",sep="\t",quote=F, row.names=F)

tr_PI2 <- as.data.frame(LUSC_C_tr2$PI)
rownames(tr_PI2) <- rownames(LUSC_C_tr2)
colnames(tr_PI2) <- c("PI")
write.table(cbind("ID"=rownames(tr_PI2 ), tr_PI2 ),file="tr_PI2",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_tr2 <- Surv(time = LUSC_C_tr2$OS_month, event = LUSC_C_tr2$Group)
surv_object_tr2
dim(surv_object_tr2 )


dim(tr2)
head(tr_PI2)
mean(tr_PI2$PI)

# ###survival analysis: fits cox ph model to find HR for PI
fit_tr2 <- survfit(surv_object_tr2~(tr_PI2$PI>mean(tr_PI2$PI)), data=tr_PI2)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_tr2)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_tr2 <- coxph(surv_object_tr2 ~(tr_PI2$PI>mean(tr_PI2$PI)), data=tr_PI2)
summary(fit.coxph_tr2)


tr_res2 <- summary(fit.coxph_tr2)
coeff2 <- round(tr_res2$coefficients[1],2)
HR2<- round(tr_res2$conf.int[1],2)
int12 <- round(tr_res2$conf.int[3],2)
int22 <- round(tr_res2$conf.int[4],2)
CI2 <- round (tr_res2$concordance[1],2)
pval2 <- tr_res2$sctest[3]

HR1_2 <- paste0(HR2 , " (",  int12,  " - ", int22, ") " )
HR1_2

tr_res1_2 <- cbind(coeff2, HR1_2, CI2, pval2)
tr_res1_2
#save results as a file
write.table(cbind("ID"=rownames(tr_res1_2), tr_res1_2),file="tr_res1_2.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file= paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg"), units="in", width=10, height=10, res=300)
#paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg")
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp2 <-  ggsurvplot(fit_tr2, data=tr_PI2, pval=TRUE,
                   risk.table=TRUE, tables.height = 0.3, #add risk table & height
                   xlab="Time in Months",
                   risk.table.col="strata", break.time.by = 12,
                   conf.int = F, censor = TRUE,
                   surv.median.line = "hv", # Add medians survival
                   #palette = c("red","blue"),#add desired color
                   size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                   font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp2
jpeg("KM_plot2.jpg", units="in", width=10, height=10, res=300)
print(pp2, newpage = FALSE)
dev.off()


#AKR1C3	ATOH8	CP	CRLF1	DRAXIN	GLIS1	H19	HIST1H2AG	RP11-189B4.6	RP11-300M24.1	TNFSF14

# AKR1C3 based
fit_AKR1C3 = survfit(surv_object_p ~ (data1$AKR1C3)>(mean(data1$AKR1C3)), data=data1)

plot_AKR1C3 <-  ggsurvplot(fit_AKR1C3, data=data1, pval=TRUE,
                   risk.table=TRUE, tables.height = 0.3, #add risk table & height
                   xlab="Time in Months",
                   risk.table.col="strata", break.time.by = 12,
                   conf.int = F, censor = TRUE,
                   surv.median.line = "hv", # Add medians survival
                   #palette = c("red","blue"),#add desired color
                   size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                   font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_AKR1C3
jpeg("KM_plot_AKR1C3.jpg", units="in", width=10, height=10, res=300)
print(plot_AKR1C3, newpage = FALSE)
dev.off()


# ATOH8 based
fit_ATOH8 = survfit(surv_object_p ~ (data1$ATOH8)>(median(data1$ATOH8)), data=data1)

plot_ATOH8 <-  ggsurvplot(fit_ATOH8, data=data1, pval=TRUE,
                           risk.table=TRUE, tables.height = 0.3, #add risk table & height
                           xlab="Time in Months",
                           risk.table.col="strata", break.time.by = 12,
                           conf.int = F, censor = TRUE,
                           surv.median.line = "hv", # Add medians survival
                           #palette = c("red","blue"),#add desired color
                           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_ATOH8
jpeg("KM_plot_ATOH8.jpg", units="in", width=10, height=10, res=300)
print(plot_ATOH8, newpage = FALSE)
dev.off()


# CP based
fit_CP = survfit(surv_object_p ~ (data1$CP)>(mean(data1$CP)), data=data1)

plot_CP <-  ggsurvplot(fit_CP, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_CP
jpeg("KM_plot_CP.jpg", units="in", width=10, height=10, res=300)
print(plot_CP, newpage = FALSE)
dev.off()




# ATOH8 based
fit_CRLF1 = survfit(surv_object_p ~ (data1$CRLF1)>(mean(data1$CRLF1)), data=data1)

plot_CRLF1 <-  ggsurvplot(fit_CRLF1, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_CRLF1
jpeg("KM_plot_CRLF1.jpg", units="in", width=10, height=10, res=300)
print(plot_CRLF1, newpage = FALSE)
dev.off()




# DRAXIN based
fit_DRAXIN = survfit(surv_object_p ~ (data1$DRAXIN)>(mean(data1$DRAXIN)), data=data1)

plot_DRAXIN <-  ggsurvplot(fit_DRAXIN, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_DRAXIN
jpeg("KM_plot_DRAXIN.jpg", units="in", width=10, height=10, res=300)
print(plot_DRAXIN, newpage = FALSE)
dev.off()




# GLIS1 based
fit_GLIS1 = survfit(surv_object_p ~ (data1$GLIS1)>(mean(data1$GLIS1)), data=data1)

plot_GLIS1 <-  ggsurvplot(fit_GLIS1, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_GLIS1
jpeg("KM_plot_GLIS1.jpg", units="in", width=10, height=10, res=300)
print(plot_GLIS1, newpage = FALSE)
dev.off()



# H19 based
fit_H19 = survfit(surv_object_p ~ (data1$H19)>(mean(data1$H19)), data=data1)

plot_H19 <-  ggsurvplot(fit_H19, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_H19
jpeg("KM_plot_H19.jpg", units="in", width=10, height=10, res=300)
print(plot_H19, newpage = FALSE)
dev.off()




# HIST1H2AG based
fit_HIST1H2AG = survfit(surv_object_p ~ (data1$HIST1H2AG)>(mean(data1$HIST1H2AG)), data=data1)

plot_HIST1H2AG <-  ggsurvplot(fit_HIST1H2AG, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_HIST1H2AG
jpeg("KM_plot_HIST1H2AG.jpg", units="in", width=10, height=10, res=300)
print(plot_HIST1H2AG, newpage = FALSE)
dev.off()




# RP11-189B4.6 based
fit_RP11_189B4_6 = survfit(surv_object_p ~ (data1$`RP11-189B4.6`)>(mean(data1$`RP11-189B4.6`)), data=data1)
fit_RP11_189B4_6 
plot_RP11_189B4_6 <- ggsurvplot(fit_RP11_189B4_6, pval=TRUE, break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

plot_RP11_189B4_6 <-  ggsurvplot(fit_RP11_189B4_6, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_RP11_189B4_6
jpeg("KM_plot_RP11_189B4_6.jpg", units="in", width=10, height=10, res=300)
print(plot_RP11_189B4_6, newpage = FALSE)
dev.off()




# RP11-300M24.1 based
fit_RP11_300M24_1 = survfit(surv_object_p ~ (data1$`RP11-300M24.1` )>(mean(data1$`RP11-300M24.1`)), data=data1)

plot_RP11_300M24_1 <-  ggsurvplot(fit_RP11_300M24_1, data=data1, pval=TRUE,
                        #  risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_RP11_300M24_1
jpeg("KM_plot_RP11_300M24_1.jpg", units="in", width=10, height=10, res=300)
print(plot_RP11_300M24_1, newpage = FALSE)
dev.off()




# TNFSF14 based
fit_TNFSF14 = survfit(surv_object_p ~ (data1$TNFSF14)>(mean(data1$TNFSF14)), data=data1)

plot_TNFSF14 <-  ggsurvplot(fit_TNFSF14, data=data1, pval=TRUE,
                          risk.table=TRUE, tables.height = 0.3, #add risk table & height
                          xlab="Time in Months",
                          risk.table.col="strata", break.time.by = 12,
                          conf.int = F, censor = TRUE,
                          surv.median.line = "hv", # Add medians survival
                          #palette = c("red","blue"),#add desired color
                          size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                          font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
plot_TNFSF14
jpeg("KM_plot_TNFSF14.jpg", units="in", width=10, height=10, res=300)
print(plot_TNFSF14, newpage = FALSE)
dev.off()




#risk prediction
data3 <- read.table("tr_PI_data_with_clin_without_peadetric_sample3.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

head(data3)
surv_object3 <- Surv(time = data3$OS_month, event = data3$Death_Status)
fit3 <- coxph(surv_object3 ~ Grade  +PI, data=data3)

# type of predicted value
predict(fit3,type="lp") #linear predictor ("lp")
predict(fit3,type="expected") #expected number of events given the covariates and follow-up time 
predict(fit3,type="risk",se.fit=TRUE) # risk score exp(lp)
predict(fit3,type="terms",se.fit=TRUE) #terms of the linear predictor
predict(fit3, newdata=new_data, type="survival",se.fit=TRUE) #terms of the linear predictor

ggforest(fit3,data=data3, fontsize = 1.6, main = "Training Data")

new_data <- head(data3,6)

plot(survfit(fit3, newdata=new_data),  xscale=365.25, xlab="Years", ylab="Survival", conf.int=F) 
# also plot the predicted survival for a 70 year old
lines(survfit(fit3, newdata=new_data), xscale=365.25, xlab="Years", ylab="Survival") 



# Derive model in the training data (after feature selection - I believe that in the paper you mentioned they use LASSO: R has a good package for this: glmnet)
#cox_model = coxph(Surv(training_data$Survival,training_data$Status) ~ ., data=training_data)

#BiocManager::install("survcomp")
library(survcomp)

head(new_data)
#risk prediction
data3_tr <- read.table("train_with_PI.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


surv_object3 <- Surv(time = data3_tr$OS_month, event = data3_tr$Death_Status)
fit3 <- coxph(surv_object3 ~ PI, data=data3_tr)
fit3 <- coxph(surv_object3 ~ IDH_mutation_status + PI, data=data3)
# Create survival estimates on validation data
new_data <- as.data.frame(data3_tr$PI)
colnames(new_data) <- c("PI")
head(new_data )
pred_train = predict (fit3, newdata = new_data)


# Determine concordance
cindex_train_val = concordance.index (pred_train, surv.time = data3_tr$OS_month, surv.event=data3_tr$Death_Status, method = "noether")
cindex_train_val

cindex_train_val$c.index
cindex_train_val$upper
cindex_train_val$lower
cindex_train_val$p.value
cindex_train_val$

head(data3)
#cross validation model
cvfit1 <- cv.glmnet(as.matrix(data3_tr[3:ncol(data3_tr)]), surv_object3, family = "cox", type.measure = "C", nfolds = 5)
lambda_min <- cvfit1$lambda.min
lambda_min
### COX model plot
jpeg(file="Cox_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=300)
plot(cvfit1)
dev.off()


##Create ROC plot for 1-,3- and 5-years survival time prediction
tr_roc1 <- survivalROC(Stime        = data3_tr$OS_month,
                       status       = data3_tr$Death_Status,
                       marker       = data3_tr$PI,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
tr_roc1

tr_roc3 <- survivalROC(Stime        = data3_tr$OS_month,
                       status       = data3_tr$Death_Status,
                       marker       = data3_tr$PI,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc3
tr_roc5 <- survivalROC(Stime        = data3_tr$OS_month,
                       status       = data3_tr$Death_Status,
                       marker       = data3_tr$PI,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5
tr_roc6 <- survivalROC(Stime        = data3_tr$OS_month,
                       status       = data3_tr$Death_Status,
                       marker       = data3_tr$PI,
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
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("4 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()







############################## External Validation ##############################

data3 <-read.table("External_validation/11gene_ext_test_QN_with_clin_mat.txt",header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
dim(data3)
head(data3)
sel_features_results<-read.table("11_key_features.txt",header =TRUE, sep = "\t", check.names = FALSE)
head(sel_features_results)
dim(sel_features_results)
sel_ftr_surv <- as.data.frame(sel_features_results[,1])
names(sel_ftr_surv ) <- c("ID")
head(sel_ftr_surv )
dim(sel_ftr_surv )



##### prepare training, test and external validation data with selected features having significant value in univariate analysis ##########
sel_test <- as.data.frame(data3[,colnames(data3) %in% c(sel_ftr_surv$ID), ])
head(sel_test,2)
dim(sel_test)

data3$OS_month
######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
test_feature_mat<- cbind(data3["OS_month"],data3["Death_event"],sel_test)
head(test_feature_mat,2)
dim(test_feature_mat)

############ remove where OS.time=NA ############
test_feature_mat1<-subset(test_feature_mat,OS_month!="NA")
test_feature_mat1<-subset(test_feature_mat1,OS_month!=0)


#test_feature_mat2<- cbind(data3["OS_month"],data3["Group"],sel_test)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(test_feature_mat1), test_feature_mat1),file="sel_test_ext.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_te=test_feature_mat1
te <- LUSC_C_te[3:ncol(LUSC_C_te)]

E1= length(te)
E1

PI_te=0
for(i in seq(from=1, to=E1 ,by=1))
{
  PI_te= PI_te+((te[,i])*(sel_features_results[i,2]))
}


LUSC_C_te$PI<-PI_te
head(LUSC_C_te)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_te), LUSC_C_te),file="test_with_PI.txt",sep="\t",quote=F, row.names=F)

te_PI <- as.data.frame(LUSC_C_te$PI)
rownames(te_PI) <- rownames(LUSC_C_te)
colnames(te_PI) <- c("PI")
write.table(cbind("ID"=rownames(te_PI ), te_PI ),file="te_PI.txt",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_te <- Surv(time = LUSC_C_te$OS_month, event = LUSC_C_te$Death_event)
surv_object_te
dim(surv_object_te )


dim(te)
head(te_PI)
mean(te_PI$PI)

# ###survival analysis: fits cox ph model to find HR for PI
fit_te <- survfit(surv_object_te~(te_PI$PI>mean(te_PI$PI)), data=te_PI)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_te)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_te <- coxph(surv_object_te ~(te_PI$PI>mean(te_PI$PI)), data=te_PI)
summary(fit.coxph_te)

fit.coxph_te$linear.predictors


te_res <- summary(fit.coxph_te)
coeff_te <- round(te_res$coefficients[1],2)
HR_te <- round(te_res$conf.int[1],2)
int1_te <- round(te_res$conf.int[3],2)
int2_te <- round(te_res$conf.int[4],2)
CI_te <- round (te_res$concordance[1],2)
pval_te <- te_res$sctest[3]

HR1_te <- paste0(HR_te , " (",  int1_te,  " - ", int2_te, ") " )
HR1_te

te_res1 <- cbind(coeff_te, HR1_te, CI_te, pval_te)
te_res1
#save results as a file
write.table(cbind("ID"=rownames(te_res1), te_res1),file="te_res1.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file= paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg"), units="in", width=10, height=10, res=300)
#paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg")
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp_te <-  ggsurvplot(fit_te, data=te_PI, pval=TRUE,
                  risk.table=TRUE, tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  legend.labs=c( "Less than Mean of PI", "Greater than Mean of PI"), legend.title="PI Score",  
                  palette=c("dodgerblue2", "red"), 
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp_te
jpeg("PI_KM_plot_te.jpg", units="in", width=10, height=10, res=300)
print(pp_te, newpage = FALSE)
dev.off()

######## test validation #########
data4 <- read.table("Final_ext_test_with_PI_data2.txt", sep="\t", header=T, row.names = 1, check.names = F)

fit4 <- coxph(surv_object_te ~ PI, data=data4)
fit4 <- coxph(surv_object_te ~ IDH_mutation_status + PI, data=data4)
# Create survival estimates on validation data
pred_validation1 = predict (fit4, newdata = data4)

pred_validation1


# Determine concordance
cindex_validation1 = concordance.index (pred_validation1, surv.time = data4$OS_month, surv.event=data4$Death_event, method = "noether")
cindex_validation1$c.index

cindex_validation1$lower
cindex_validation1$upper
cindex_validation1$p.value
uu# create survival object
surv_object_p_te <- Surv(time = as.numeric(data4$OS_month), event = data4$Death_event)
data4$Survival_class
#Cox multivariate model
cox_multivariate_te <- coxph(surv_object_p_te ~    Gender +  Grade + as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  ,  data=data4  )
summary(cox_multivariate_te)

jpeg(file="test_multivariate_only_clin.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te,data=data4, fontsize = 1.6, main = "Test Data")
dev.off()


cox_multivariate_te1 <- coxph(surv_object_p_te ~  Gender + Grade + as.numeric(Age)+ IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  + PI,  data=data4 )
summary(cox_multivariate_te1)
ggforest(cox_multivariate_te1,data=data4)
jpeg(file="test_multivariate.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te1,data=data4, fontsize = 1.6, main = "Test Data")
dev.off()



cox_multivariate_te3 <- coxph(surv_object_p_te ~  Tumor_Grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth ,  data=data4 )

summary(cox_multivariate_te3 )
jpeg(file="test_multivariate_only_clin1.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te3,data=data4, fontsize = 1.6, main = "Test Data")
dev.off()

cox_multivariate_te4 <- coxph(surv_object_p_te ~  Tumor_Grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth +PI ,  data=data4 )

summary(cox_multivariate_te4 )

jpeg(file="test_multivariate_with PI_1.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te4,data=data4, fontsize = 1.6, main = "Test Data")
dev.off()




#######
covariates_te <- c( "Gender" , "Grade", "Survival_class",  "Age", "IDH_mutation_status", "Codel_1p19q_status", "MGMTp_meth_status")
univ_formulas_te <- sapply(covariates_te,  function(x) as.formula(paste('surv_object_te ~', x)))

univ_models_te <- lapply(univ_formulas_te, function(x){coxph(x, data = data4 )})
# Extract data 
univ_results_te <- lapply(univ_models_te,
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

univ_results_te
univ_results_te
res_uni_te<- as.matrix(univ_results_te)
res_uni_te

write.table(res_uni , file = "Univariate_results_traing_data.txt", sep="\t", quote=F, row.names = T)

data4$OS_month

###### 
##Create ROC plot for 1-,3- and 5-years survival time prediction
te_roc1 <- survivalROC(Stime        = data4$OS_month,
                       status       = data4$Death_event,
                       marker       = data4$PI,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
te_roc1

te_roc3 <- survivalROC(Stime        = data4$OS_month,
                       status       = data4$Death_event,
                       marker       = data4$PI,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
te_roc3
te_roc5 <- survivalROC(Stime        = data4$OS_month,
                       status       = data4$Death_event,
                       marker       = data4$PI,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5
jpeg(file="ROC_test.jpeg", units="in", width=10, height=10, res=300)
plot(te_roc1$FP, te_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival time Prediction")
lines(te_roc3$FP, te_roc3$TP, type="l", lty=2, col="blue")
lines(te_roc5$FP, te_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(te_roc1$AUC,2)),  paste("4 Years AUC = ",round(te_roc3$AUC,2)),  paste("5 Years AUC = ",round(te_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()


te_roc6 <- survivalROC(Stime        = data4$OS_month,
                       status       = data4$Death_event,
                       marker       = data4$Codel_1p19q_status,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
te_roc6





################ Multivariate with genes ###############
data_new_tr <- read.table("Train_Clin_data_with_11genes_PI.txt", sep="\t", header=T, row.names = 1, check.names = F)
data_new_te <- read.table("Ext_test_Clin_data_with_11genes_PI.txt", sep="\t", header=T, row.names = 1, check.names = F)

head(data_new_tr ,2)

# create survival object
surv_object_tr <- Surv(time = as.numeric(data_new_tr$OS_month), event = data_new_tr$Death_Status)

surv_object_te <- Surv(time = as.numeric(data_new_te$OS_month), event = data_new_te$Death_Status)

#Cox multivariate model

cox_multivariate_tr <- coxph(surv_object_tr ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth + AKR1C3 + ATOH8 + CP + CRLF1 + DRAXIN + GLIS1 +  H19 + HIST1H2AG + `RP11-189B4.6` + `RP11-300M24.1` + TNFSF14,  data=data_new_tr  )
cox_multivariate_tr_c <- coxph(surv_object_tr ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth ,  data=data_new_tr  )
#cox_multivariate_tr <- coxph(surv_object_tr ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth + as.numeric(AKR1C3) + as.numeric(ATOH8) + as.numeric(CP) + as.numeric(CRLF1) + as.numeric(DRAXIN) + as.numeric(GLIS1) +  as.numeric(H19) + as.numeric(HIST1H2AG) + as.numeric(`RP11-189B4.6`) + as.numeric(`RP11-300M24.1`) + as.numeric(TNFSF14),  data=data_new_tr  )


#cox_multivariate_tr <- coxph(surv_object_tr ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q +  AKR1C3 + ATOH8 + CP + CRLF1 + DRAXIN + GLIS1 +  H19 + HIST1H2AG + `RP11-189B4.6` + `RP11-300M24.1` + TNFSF14,  data=data_new_tr  )


summary(cox_multivariate_tr)

ggforest(cox_multivariate_tr,data=data_new_tr)
ggforest(cox_multivariate_tr_c,data=data_new_tr)

jpeg(file="Train_multivariate_with_genes.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_tr,data=data_new_tr, fontsize = 1.6, main = "Training Data")
dev.off()

jpeg(file="Train_multivariate_with_clin.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_tr_c,data=data_new_tr, fontsize = 1.6, main = "Training Data")
dev.off()

#cox_multivariate_te <- coxph(surv_object_te ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth + AKR1C3 + ATOH8 + CP + CRLF1 + DRAXIN + GLIS1 +  H19 + HIST1H2AG + `RP11-189B4.6` + `RP11-300M24.1` + TNFSF14,  data=data_new_te  )
#cox_multivariate_te <- coxph(surv_object_te ~    Grade + Gender +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth,  data=data_new_te  )
cox_multivariate_te <- coxph(surv_object_te ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth + AKR1C3 + ATOH8 + CP + CRLF1 + DRAXIN + GLIS1 +  H19 + HIST1H2AG + `RP11-189B4.6` + `RP11-300M24.1` + TNFSF14,  data=data_new_te  )

cox_multivariate_te_c <- coxph(surv_object_te ~    Tumor_grade + Sex +   as.numeric(Age) + IDH_mut + Codel_1p19q + MGMTp_meth ,  data=data_new_te  )
#cox_multivariate_te <- coxph(surv_object_te ~    Grade + Gender +   as.numeric(Age) + IDH_mut + Codel_1p19q +  MGMTp_meth + as.numeric(AKR1C3) + as.numeric(ATOH8) + as.numeric(CP) + as.numeric(CRLF1) + as.numeric(DRAXIN) + as.numeric(GLIS1) +  as.numeric(H19) + as.numeric(HIST1H2AG) + as.numeric(`RP11-189B4.6`) + as.numeric(`RP11-300M24.1`) + as.numeric(TNFSF14),  data=data_new_te  )

#cox_multivariate_te <- coxph(surv_object_te ~    Grade + Gender +   as.numeric(Age) + IDH_mut + Codel_1p19q +  MGMTp_meth + as.numeric(AKR1C3) + as.numeric(ATOH8) + as.numeric(CP) + as.numeric(CRLF1) + as.numeric(DRAXIN) + as.numeric(GLIS1) +  as.numeric(H19) + as.numeric(HIST1H2AG) + as.numeric(`RP11-189B4.6`) + as.numeric(`RP11-300M24.1`) + as.numeric(TNFSF14),  data=data_new_te  )

summary(cox_multivariate_te)
ggforest(cox_multivariate_te,data=data_new_te)

jpeg(file="Test_multivariate_with_genes.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te,data=data_new_te, fontsize = 1.6, main = "Test Data")
dev.off()

ggforest(cox_multivariate_te_c,data=data_new_te)
jpeg(file="Test_multivariate_with_clin.jpeg", units="in", width=15, height=12, res=350)
#ggforest(cox_multivariate,data=LGG_clin )
ggforest(cox_multivariate_te_c,data=data_new_te, fontsize = 1.6, main = "Test Data")
dev.off()