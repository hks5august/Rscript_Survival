#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

#library(pca3d)


#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")

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


tr_clin_new1<-subset(tr_clin_new,OS_month!="NA")
# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new1$OS_month), event = tr_clin_new1$Death_event)


#Survival analysis
  
  surv_fit_tr <- survfit(surv_object_p_tr  ~ (tr_clin_new1$PI>mean(tr_clin_new$PI)), data=tr_clin_new1)
  summary(surv_fit_tr);
  #ggsurvplot(surv_fit)
  fit1.coxph_tr <- coxph(surv_object_p_tr  ~ (tr_clin_new1$PI>mean(tr_clin_new$PI)), data=tr_clin_new1)
  # summary(fit1.coxph);
  first_tr <- coef(summary(fit1.coxph_tr))
  HR_tr <- round(first_tr[2],2)
  CI_tr<- round(fit1.coxph_tr$concordance[6],2)
  #CI_tr
  options(scipen = 1)
  options(digits = 3)
  
  pval_tr <- format(first[5], scientific = T)
  #?ggsurvplot
  # we produce a Kaplan Meier plot
  #jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
  
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(plot.title=element_text(hjust=0.5,face="bold", size=25))}
  
  surv_plot1_tr <- ggsurvplot(  surv_fit_tr, data=tr_clin_new1, pval=F, 
                                risk.table=TRUE, 
                                tables.height = 0.25, #add risk table & height
                                xlab="Time in months", 
                                legend.labs = c("(>)mean PI","(<) mean PI"), 
                                #risk.table.col="strata", 
                                break.time.by = 12,
                                title=paste0("KM Plot for Training Data", "\n"),
                                ggtheme = custom_theme(),
                                #conf.int = F, 
                                censor = TRUE,
                                surv.median.line = "hv", # Add medians survival
                                palette = c("red","blue"),#add desired color
                                size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  surv_plot1_tr
  #Add HR and P-values
  surv_plot_train  <- surv_plot1_tr$plot +
    ggplot2::annotate(
      "text",
      x = Inf, y = Inf,
      vjust = 2, hjust = 1,
      label =  paste0("HR = ", HR_tr, "\n", "p-val = ", pval_tr, "\n", "C-Index = ", CI_tr),
      size = 8 ) 
  
  #save plot
  
  jpeg(file= "Training_PI_KM_plot.jpeg", units="in", width=12, height=12, res =350)
  surv_plot_train
  dev.off()
  
  surv_plot1_tr2 <- ggsurvplot(  surv_fit_tr, data=tr_clin_new1, pval=TRUE, 
                                risk.table=TRUE, 
                                tables.height = 0.25, #add risk table & height
                                xlab="Time in months", 
                                legend.labs = c("(>)mean PI","(<) mean PI"), 
                                #risk.table.col="strata", 
                                break.time.by = 12,
                                title=paste0("KM Plot for Training Data", "\n"),
                                ggtheme = custom_theme(),
                                #conf.int = F, 
                                censor = TRUE,
                                surv.median.line = "hv", # Add medians survival
                                palette = c("red","blue"),#add desired color
                                size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  jpeg(file= "Training_PI_KM_plot2.jpeg", units="in", width=12, height=12, res =350)
  surv_plot1_tr2
  dev.off()
  

  
  
  #### Test 1 data 
  
  te1_clin_new1<-subset(te1_clin_new,OS_month!="NA")
  # create survival object
  surv_object_p_te1 <- Surv(time = as.numeric(te1_clin_new1$OS_month), event = te1_clin_new1$Death_event)
  
  
  #Survival analysis
  
  surv_fit_te1 <- survfit(surv_object_p_te1  ~ (te1_clin_new1$PI>mean(te1_clin_new$PI)), data=te1_clin_new1)
  summary(surv_fit_te1);
  #ggsurvplot(surv_fit)
  fit1.coxph_te1 <- coxph(surv_object_p_te1  ~ (te1_clin_new1$PI>mean(te1_clin_new$PI)), data=te1_clin_new1)
  # summary(fit1.coxph);
  first_te1 <- coef(summary(fit1.coxph_te1))
  HR_te1 <- round(first_te1[2],2)
  CI_te1<- round(fit1.coxph_te1$concordance[6],2)
  #CI_te1
  options(scipen = 1)
  options(digits = 3)
  
  pval_te1 <- format(first[5], scientific = T)
  #?ggsurvplot
  # we produce a Kaplan Meier plot
  #jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
  
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(plot.title=element_text(hjust=0.5,face="bold", size=25))}
  
  surv_plot1_te1 <- ggsurvplot(  surv_fit_te1, data=te1_clin_new1, pval=F, 
                                 risk.table=TRUE, 
                                 tables.height = 0.25, #add risk table & height
                                 xlab="Time in months", 
                                 legend.labs = c("(>)mean PI","(<) mean PI"), 
                                 #risk.table.col="strata", 
                                 break.time.by = 12,
                                 title=paste0("KM Plot for Test1 Data", "\n"),
                                 ggtheme = custom_theme(),
                                 #conf.int = F, 
                                 censor = TRUE,
                                 surv.median.line = "hv", # Add medians survival
                                 palette = c("red","blue"),#add desired color
                                 size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                 font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  
  #Add HR and P-values
  surv_plot_te1  <- surv_plot1_te1$plot + 
    ggplot2::annotate(
      "text",
      x = Inf, y = Inf,
      vjust = 2, hjust = 1,
      label =  paste0("HR = ", HR_te1, "\n", "p-val = ", pval_te1, "\n", "C-Index = ", CI_te1),
      size = 8 ) 
  
  #save plot
  
  jpeg(file= "test1_PI_KM_plot.jpeg", units="in", width=12, height=12, res =350)
  surv_plot_te1
  dev.off()

  
  surv_plot1_te1_2 <- ggsurvplot(  surv_fit_te1, data=te1_clin_new1, pval=TRUE, 
                                 risk.table=TRUE, 
                                 tables.height = 0.25, #add risk table & height
                                 xlab="Time in months", 
                                 legend.labs = c("(>)mean PI","(<) mean PI"), 
                                 #risk.table.col="strata", 
                                 break.time.by = 12,
                                 title=paste0("KM Plot for Test1 Data", "\n"),
                                 ggtheme = custom_theme(),
                                 #conf.int = F, 
                                 censor = TRUE,
                                 surv.median.line = "hv", # Add medians survival
                                 palette = c("red","blue"),#add desired color
                                 size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                 font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  jpeg(file= "Test1_PI_KM_plot2.jpeg", units="in", width=12, height=12, res =350)
  surv_plot1_te1_2
  dev.off()
  
  
  
  ####### Test2 data
  
  te2_clin_new1<-subset(te2_clin_new,OS_month!="NA")
  # create survival object
  surv_object_p_te2 <- Surv(time = as.numeric(te2_clin_new1$OS_month), event = te2_clin_new1$Death_event)
  
  
  #Survival analysis
  
  surv_fit_te2 <- survfit(surv_object_p_te2  ~ (te2_clin_new1$PI>mean(te2_clin_new$PI)), data=te2_clin_new1)
  summary(surv_fit_te2);
  #ggsurvplot(surv_fit)
  fit1.coxph_te2 <- coxph(surv_object_p_te2  ~ (te2_clin_new1$PI>mean(te2_clin_new$PI)), data=te2_clin_new1)
  # summary(fit1.coxph);
  first_te2 <- coef(summary(fit1.coxph_te2))
  HR_te2 <- round(first_te2[2],2)
  CI_te2<- round(fit1.coxph_te2$concordance[6],2)
  #CI_te2
  options(scipen = 1)
  options(digits = 3)
  
  pval_te2 <- format(first[5], scientific = T)
  #?ggsurvplot
  # we produce a Kaplan Meier plot
  #jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
  
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(plot.title=element_text(hjust=0.5,face="bold", size=25))}
  
  surv_plot1_te2 <- ggsurvplot(  surv_fit_te2, data=te2_clin_new1, pval=F, 
                                 risk.table=TRUE, 
                                 tables.height = 0.25, #add risk table & height
                                 xlab="Time in months", 
                                 legend.labs = c("(>)mean PI","(<) mean PI"), 
                                 #risk.table.col="strata", 
                                 break.time.by = 12,
                                 title=paste0("KM Plot for Validation Data2", "\n"),
                                 ggtheme = custom_theme(),
                                 #conf.int = F, 
                                 censor = TRUE,
                                 surv.median.line = "hv", # Add medians survival
                                 palette = c("red","blue"),#add desired color
                                 size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                 font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  
  #Add HR and P-values
  surv_plot_te2  <- surv_plot1_te2$plot +
    ggplot2::annotate(
      "text",
      x = Inf, y = Inf,
      vjust = 2, hjust = 1,
      label =  paste0("HR = ", HR_te2, "\n", "p-val = ", pval_te2, "\n", "C-Index = ", CI_te2),
      size = 8 ) 
  
  #save plot
  
  jpeg(file= "test2_PI_KM_plot.jpeg", units="in", width=12, height=12, res =350)
  surv_plot_te2
  dev.off()

  
  surv_plot1_te2_2 <- ggsurvplot(  surv_fit_te2, data=te2_clin_new1, pval=TRUE, 
                                   risk.table=TRUE, 
                                   tables.height = 0.25, #add risk table & height
                                   xlab="Time in months", 
                                   legend.labs = c("(>)mean PI","(<) mean PI"), 
                                   #risk.table.col="strata", 
                                   break.time.by = 12,
                                   title=paste0("KM Plot for Validation Data2", "\n"),
                                   ggtheme = custom_theme(),
                                   #conf.int = F, 
                                   censor = TRUE,
                                   surv.median.line = "hv", # Add medians survival
                                   palette = c("red","blue"),#add desired color
                                   size=2, font.tickslab = c(20,  "black"), font.y = c(20,  "black"), 
                                   font.x = c(20, "black"), font.legend = c(20, "black"), font.label = c(20, "bold", "black"))
  
  jpeg(file= "Test2_PI_KM_plot2.jpeg", units="in", width=12, height=12, res =350)
  surv_plot1_te2_2
  dev.off()
  