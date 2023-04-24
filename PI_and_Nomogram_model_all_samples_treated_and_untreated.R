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

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/Both_treated_untreated/")

#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")

set.seed(7)

    #data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
    data <- read.table("11_genes_LTS_STS_mat_with_clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
   
    
    data1 <- subset(data, OS_month!="NA")
    dim(data1)
    head(data1)
    
    
    ##############
    
    ######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
    ######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
    
    features<- read.table("genes_list", header =TRUE,  sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
    head(features)
    rownames(features) <- features$ID
    
    
    ######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
    sel_features_results<-read.table("11_key_features.txt",header =TRUE, sep = "\t", check.names = FALSE, row.names=1)
    head(sel_features_results)
    dim(sel_features_results)
    
    
    #data preparation with selected features 
    sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
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
    write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="LTS_STS_mat_with_PI.txt",sep="\t",quote=F, row.names=F)
    
    tr_PI <- as.data.frame(LUSC_C_tr$PI)
    rownames(tr_PI) <- rownames(LUSC_C_tr)
    colnames(tr_PI) <- c("PI")
    write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="tr_STS_LTS_PI",sep="\t",quote=F, row.names=F)
    
    
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
    #write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res1.txt",sep="\t",quote=F, row.names=F)
    
    
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
    #jpeg("Train_PI_KM_plot.jpg", units="in", width=10, height=10, res=300)
    jpeg("LTS_STS_PI_KM_plot.jpg", units="in", width=10, height=10, res=350)
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
                           predict.time = 48,
                           method       = "KM", 
                           #lambda = lambda_min, 
                           span = NULL, 
                           window ="symmetric")
    tr_roc3
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
    #jpeg(file="ROC_train.jpeg", units="in", width=10, height=10, res=300)
    jpeg(file="ROC_train_LTS_STS.jpeg", units="in", width=10, height=10, res=300)
    plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="FP", col="red",
         ylab="TP",main= "AUC Curve for Survival Prediction")
    lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
    lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
    legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("4 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")
    
    abline(0,1)
    dev.off()
    


####################################  STS/MTS/LTS with PI           ##############################################
    #data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
    data <- read.table("11_genes_LTS_MTS_STS_mat_with_clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    
    data1 <- subset(data, OS_month!="NA")
    head(data1)
    dim(data1)
    
    
    ##############
    
    ######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
    ######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
    
    features<- read.table("genes_list", header =TRUE,  sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)
    head(features)
    rownames(features) <- features$ID
    
    
    ######## Extract selected genes results from univariate analysis results, for the purpose of obtaining beta coefficient  ####################
    sel_features_results<-read.table("11_key_features.txt",header =TRUE, sep = "\t", check.names = FALSE, row.names=1)
    head(sel_features_results)
    dim(sel_features_results)
    
    
    #data preparation with selected features 
    sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
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
    write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train_STS_MTS_LTS.txt",sep="\t",quote=F, row.names=F)
    
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
    write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="STS_MTS_LTS_mat_with_PI.txt",sep="\t",quote=F, row.names=F)
    
    tr_PI <- as.data.frame(LUSC_C_tr$PI)
    rownames(tr_PI) <- rownames(LUSC_C_tr)
    colnames(tr_PI) <- c("PI")
    write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="tr_STS_MTS_LTS_PI",sep="\t",quote=F, row.names=F)
    
    
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
    #write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res1.txt",sep="\t",quote=F, row.names=F)
    
    
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
    #jpeg("Train_PI_KM_plot.jpg", units="in", width=10, height=10, res=300)
    jpeg("STS_MTS_LTS_PI_KM_plot.jpg", units="in", width=10, height=10, res=350)
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
    tr_roc5 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                           status       = LUSC_C_tr$Death_Status,
                           marker       = LUSC_C_tr$PI,
                           predict.time = 60,
                           method       = "KM", 
                           #lambda = lambda_min, 
                           span = NULL, 
                           window ="symmetric")
    tr_roc5
    tr_roc4 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                           status       = LUSC_C_tr$Death_Status,
                           marker       = LUSC_C_tr$PI,
                           predict.time = 48,
                           method       = "KM", 
                           #lambda = lambda_min, 
                           span = NULL, 
                           window ="symmetric")
    
    tr_roc4
    #jpeg(file="ROC_train.jpeg", units="in", width=10, height=10, res=300)
    jpeg(file="ROC_train_STS_MTS_LTS.jpeg", units="in", width=10, height=10, res=300)
    plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
         xlab="FP", col="red",
         ylab="TP",main= "AUC Curve for Survival Prediction")
    lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
    lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
    legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")
    
    abline(0,1)
    dev.off()
    
    
    
    
    ######### multivariate
    #Cox multivariate model
    data2 <- read.table("tr_STS_MTS_LTS_PI_with_clin_data", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    surv_object_tr2 <- Surv(time = as.numeric(data2$OS_month), event = data2$Death_Status)
    
    cox_multivariate <- coxph( surv_object_tr2 ~  Gender + Grade +  as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ,  data=data2 )
    cox_multivariate_2 <- coxph( surv_object_tr2 ~  Sex + Tumor_Grade+  as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ,  data=data2 )
    
    summary(cox_multivariate )
    jpeg(file="multivariate_STS_MTS_LTS_clin_without_PI.jpeg", units="in", width=10, height=10, res=300)
    #ggforest(cox_multivariate,data=data1 )
    ggforest(cox_multivariate,data=data2)
    dev.off()
    
    
    
    cox_multivariate2 <- coxph( surv_object_tr2 ~  Gender + Grade +  as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status + RT + Chemo_TMZ + PI,  data=data2 )
    summary(cox_multivariate2 )
    jpeg(file="multivariate_STS_MTS_LTS_clin_with_PI.jpeg", units="in", width=10, height=10, res=300)
    #ggforest(cox_multivariate,data=data1 )
    ggforest(cox_multivariate2,data=data2)
    dev.off()
    
    
    ######################### Nomogram STS and LTS #########
    
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
    library(foreign)
    args <- commandArgs(TRUE)
    
    
    setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/Both_treated_untreated/")

    
    #setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")
    
    set.seed(7)
    
    tr_clin_new <- read.table("tr_STS_MTS_LTS_PI_with_clin_data", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    tr_clin_new1 <- read.table("tr_STS_LTS_PI_with_clin_data", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    
    # te_clin_new <- read.table("te_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    head(tr_clin_new,2)
    dim(tr_clin_new)
    tr_clin_new$Age_group
    # create survival object
    surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_Status)
    
    
    head(tr_clin_new)
    
    
    d <-  cbind(tr_clin_new["OS_month"], tr_clin_new["Death_Status"], tr_clin_new["Grade"], tr_clin_new["Gender"],tr_clin_new["Age"], tr_clin_new["IDH_mutation_status"], tr_clin_new["Codel_1p19q_status"], tr_clin_new["MGMTp_meth_status"], tr_clin_new["RT"],  tr_clin_new["Chemo_TMZ"], tr_clin_new["Histology"], tr_clin_new["PI"])
   # d_te <-  cbind(te_clin_new["OS_month"], te_clin_new["Death_Status"], te_clin_new["Grade"], te_clin_new["Gender"],te_clin_new["Age"], te_clin_new["IDH_mutation_status"], te_clin_new["Codel_1p19q_status"], te_clin_new["PI"])
    
    dim(d)
    surv_object_tr <- Surv(time = as.numeric(d$OS_month), event = d$Death_Status)
    
    
    ddist <- datadist(d)
    options(datadist='ddist')
    
    #ddist_te <- datadist(d_te)
    #options(datadist='ddist_te')
    
    ### logistic regression models
    f1 <- lrm(Death_Status ~  PI , data = d,  x = T,y = T)
    f2 <- lrm(Death_Status ~  PI + Grade , data = d,  x = T,y = T)
    f3 <- lrm(Death_Status ~  PI + Grade +  IDH_mutation_status , data = d,  x = T,y = T)
    f4 <- lrm(Death_Status ~  PI + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)
    f5 <- lrm(Death_Status ~  PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)
    #f6 <- lrm(Death_Status ~  PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status , data = d,  x = T,y = T)
    f6 <- lrm(Death_Status ~  PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d,  x = T,y = T)
    f7 <- lrm(Death_Status ~  RT , data = d,  x = T,y = T)
    f8<- lrm(Death_Status ~  Chemo_TMZ , data = d,  x = T,y = T)
    #f9<- lrm(Death_Status ~  PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + RT + Chemo_TMZ + MGMTp_meth_status , data = d,  x = T,y = T)
    f9<- lrm(Death_Status ~  PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + RT + Chemo_TMZ  , data = d,  x = T,y = T)
    
    #nom <- nom3 <- nomogram(f, fun=plogis, funlabel="Risk of Death")
    nom1 <- nomogram(f1, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom2 <- nomogram(f2, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom3 <- nomogram(f3, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom4 <- nomogram(f4, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom5 <- nomogram(f5, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom6 <- nomogram(f6, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom7 <- nomogram(f7, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom8 <- nomogram(f8, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom9 <- nomogram(f9, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
    nom1
    plot(nom1)
    plot(nom7)
    plot(nom8)
    plot(nom9)
    
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
    
    f1
    Survival(f1)
    surv1<- Survival(f1)
    
    #coxph models
    cox1 <-cph(Surv(OS_month,Death_Status==1) ~ PI ,  x = T,y = T, data = d, surv = T)
    cox2 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade ,  x = T,y = T, data = d, surv = T)
    cox3 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status ,  x = T,y = T, data = d, surv = T)
    cox4 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox5 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox6 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox7 <-cph(Surv(OS_month,Death_Status==1) ~ PI +  Grade +   Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox8 <-cph(Surv(OS_month,Death_Status==1) ~ RT,  x = T,y = T, data = d, surv = T)
    cox9 <-cph(Surv(OS_month,Death_Status==1) ~ Chemo_TMZ,  x = T,y = T, data = d, surv = T)
    cox10 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + RT + Chemo_TMZ  ,  x = T,y = T, data = d, surv = T)
    
    
    surv1<- Survival(cox1)
    risk_1 <- function(x)1/(1+exp(-x))
    surv1_1<- function(x)surv1(1*12,lp=x) # defined time.inc,1 year OS
    surv1_2<- function(x)surv1(1*48,lp=x) # defined time.inc,4 year OS
    surv1_3<- function(x)surv1(1*60,lp=x) # defined time.inc,5 year OS
    
    nom1_cox<-nomogram(cox1,fun = list(risk_1, surv1_1,surv1_2,surv1_3),lp = F,
                       funlabel = c("Risk", "1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
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
    surv7_3<- function(x)surv7(1*60,lp=x) # defined time.inc,5 year OS
    
    nom7_cox<-nomogram(cox7,fun = list(risk7 , surv7_1,surv7_2,surv7_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom7_cox),xfrac = .7)
    
    jpeg("Nomogram_PI_grade_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom7_cox)
    dev.off()
    
    
    
    surv8<- Survival(cox8)
    risk8 <- function(x)1/(1+exp(-x))
    surv8_1<- function(x)surv8(1*12,lp=x) # defined time.inc,1 year OS
    surv8_2<- function(x)surv8(1*48,lp=x) # defined time.inc,4 year OS
    surv8_3<- function(x)surv8(1*60,lp=x) # defined time.inc,5 year OS
    
    nom8_cox<-nomogram(cox8,fun = list(risk8 , surv8_1,surv8_2,surv8_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom8_cox))
    
    jpeg("Nomogram_PI_grade_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom8_cox)
    dev.off()
    
    
    
    surv9<- Survival(cox9)
    risk9 <- function(x)1/(1+exp(-x))
    surv9_1<- function(x)surv9(1*12,lp=x) # defined time.inc,1 year OS
    surv9_2<- function(x)surv9(1*48,lp=x) # defined time.inc,4 year OS
    surv9_3<- function(x)surv9(1*60,lp=x) # defined time.inc,5 year OS
    
    nom9_cox<-nomogram(cox9,fun = list(risk9 , surv9_1,surv9_2,surv9_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot(nom9_cox)
    
    
    
    surv10<- Survival(cox10)
    risk10 <- function(x)1/(1+exp(-x))
    surv10_1<- function(x)surv10(1*12,lp=x) # defined time.inc,1 year OS
    surv10_2<- function(x)surv10(1*24,lp=x) # defined time.inc,4 year OS
    surv10_3<- function(x)surv10(1*36,lp=x) # defined time.inc,4 year OS
    surv10_4<- function(x)surv10(1*48,lp=x) # defined time.inc,5 year OS
    surv10_5<- function(x)surv10(1*60,lp=x) # defined time.inc,5 year OS
    
    nom10_cox<-nomogram(cox10,fun = list(risk10 , surv10_1,surv10_2,surv10_3, surv10_4, surv10_5),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "2-Year Survival Probability", "3-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('1.0', '0.95', '0.9','0.85','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom10_cox),xfrac = .9)
    
    jpeg("Nomogram_PI_Clin_features_STS_MTS_LTS__OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom10_cox)
    dev.off()
    ########## C-index and p-value ######
    
    #cox1$coefficients
    f1_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI  , data = d)
    f2_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade , data = d)
    f3_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status  , data = d)
    f4_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
    f5_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
    f6_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
    f7_c<-coxph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  Codel_1p19q_status , data = d)
    f8_c<-coxph(Surv(OS_month,Death_Status==1) ~  Age + Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
    f9_c<-coxph(Surv(OS_month,Death_Status==1) ~   Grade +  IDH_mutation_status + Codel_1p19q_status , data = d)
    f10_c<-coxph(Surv(OS_month,Death_Status==1) ~   PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + RT + Chemo_TMZ , data = d)
    
    sum.surv1<-summary(f1_c)
    sum.surv2<-summary(f2_c)
    sum.surv3<-summary(f3_c)
    sum.surv4<-summary(f4_c)
    sum.surv5<-summary(f5_c)
    sum.surv6<-summary(f6_c)
    sum.surv7<-summary(f7_c)
    sum.surv8<-summary(f8_c)
    sum.surv9<-summary(f9_c)
    sum.surv10<-summary(f10_c)
    sum.surv1
    sum.surv2
    sum.surv3
    sum.surv4
    sum.surv5
    sum.surv6
    sum.surv7
    sum.surv8
    sum.surv9
    sum.surv10
    c_index1<-sum.surv1$concordance
    c_index2<-sum.surv2$concordance
    c_index3<-sum.surv3$concordance
    c_index4<-sum.surv4$concordance
    c_index5<-sum.surv5$concordance
    c_index6<-sum.surv6$concordance
    c_index7<-sum.surv7$concordance
    c_index8<-sum.surv8$concordance
    c_index9<-sum.surv9$concordance
    c_index10<-sum.surv10$concordance
    c_index1
    c_index2
    c_index3
    c_index4
    c_index5
    c_index6
    c_index7
    c_index8
    c_index9
    c_index10
    
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
    
   ################  STS and LTS #####
    
    
    tr_clin_new <- read.table("tr_STS_LTS_PI_with_clin_data", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    
    # te_clin_new <- read.table("te_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    head(tr_clin_new,2)
    dim(tr_clin_new)
    tr_clin_new$Age_group
    # create survival object
    surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_Status)
    
    
    head(tr_clin_new)
    
    
    d <-  cbind(tr_clin_new["OS_month"], tr_clin_new["Death_Status"], tr_clin_new["Grade"], tr_clin_new["Gender"],tr_clin_new["Age"], tr_clin_new["IDH_mutation_status"], tr_clin_new["Codel_1p19q_status"], tr_clin_new["MGMTp_meth_status"], tr_clin_new["RT"],  tr_clin_new["Chemo_TMZ"], tr_clin_new["Histology"], tr_clin_new["PI"])
    # d_te <-  cbind(te_clin_new["OS_month"], te_clin_new["Death_Status"], te_clin_new["Grade"], te_clin_new["Gender"],te_clin_new["Age"], te_clin_new["IDH_mutation_status"], te_clin_new["Codel_1p19q_status"], te_clin_new["PI"])
    
    dim(d)
    surv_object_tr <- Surv(time = as.numeric(d$OS_month), event = d$Death_Status)
    
    
    ddist <- datadist(d)
    options(datadist='ddist')
    
    
    #coxph models
    cox1 <-cph(Surv(OS_month,Death_Status==1) ~ PI ,  x = T,y = T, data = d, surv = T)
    cox2 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade ,  x = T,y = T, data = d, surv = T)
    cox3 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status ,  x = T,y = T, data = d, surv = T)
    cox4 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox5 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox6 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox7 <-cph(Surv(OS_month,Death_Status==1) ~ PI +  Grade +   Codel_1p19q_status,  x = T,y = T, data = d, surv = T)
    cox8 <-cph(Surv(OS_month,Death_Status==1) ~ RT,  x = T,y = T, data = d, surv = T)
    cox9 <-cph(Surv(OS_month,Death_Status==1) ~ Chemo_TMZ,  x = T,y = T, data = d, surv = T)
    cox10 <-cph(Surv(OS_month,Death_Status==1) ~ PI + Age + Gender + Grade +  IDH_mutation_status + Codel_1p19q_status + RT + Chemo_TMZ  ,  x = T,y = T, data = d, surv = T)
    
    
    surv1<- Survival(cox1)
    risk_1 <- function(x)1/(1+exp(-x))
    surv1_1<- function(x)surv1(1*12,lp=x) # defined time.inc,1 year OS
    surv1_2<- function(x)surv1(1*48,lp=x) # defined time.inc,4 year OS
    surv1_3<- function(x)surv1(1*60,lp=x) # defined time.inc,5 year OS
    
    nom1_cox<-nomogram(cox1,fun = list(risk_1, surv1_1,surv1_2,surv1_3),lp = F,
                       funlabel = c("Risk", "1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
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
    surv7_3<- function(x)surv7(1*60,lp=x) # defined time.inc,5 year OS
    
    nom7_cox<-nomogram(cox7,fun = list(risk7 , surv7_1,surv7_2,surv7_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom7_cox),xfrac = .7)
    
    jpeg("Nomogram_PI_grade_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom7_cox)
    dev.off()
    
    
    
    surv8<- Survival(cox8)
    risk8 <- function(x)1/(1+exp(-x))
    surv8_1<- function(x)surv8(1*12,lp=x) # defined time.inc,1 year OS
    surv8_2<- function(x)surv8(1*48,lp=x) # defined time.inc,4 year OS
    surv8_3<- function(x)surv8(1*60,lp=x) # defined time.inc,5 year OS
    
    nom8_cox<-nomogram(cox8,fun = list(risk8 , surv8_1,surv8_2,surv8_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom8_cox))
    
    jpeg("Nomogram_PI_grade_codel_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom8_cox)
    dev.off()
    
    
    
    surv9<- Survival(cox9)
    risk9 <- function(x)1/(1+exp(-x))
    surv9_1<- function(x)surv9(1*12,lp=x) # defined time.inc,1 year OS
    surv9_2<- function(x)surv9(1*48,lp=x) # defined time.inc,4 year OS
    surv9_3<- function(x)surv9(1*60,lp=x) # defined time.inc,5 year OS
    
    nom9_cox<-nomogram(cox9,fun = list(risk9 , surv9_1,surv9_2,surv9_3),lp = F,
                       funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                       maxscale = 100,
                       fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot(nom9_cox)
    
    
    
    surv10<- Survival(cox10)
    risk10 <- function(x)1/(1+exp(-x))
    surv10_1<- function(x)surv10(1*12,lp=x) # defined time.inc,1 year OS
    surv10_2<- function(x)surv10(1*24,lp=x) # defined time.inc,4 year OS
    surv10_3<- function(x)surv10(1*36,lp=x) # defined time.inc,4 year OS
    surv10_4<- function(x)surv10(1*48,lp=x) # defined time.inc,5 year OS
    surv10_5<- function(x)surv10(1*60,lp=x) # defined time.inc,5 year OS
    
    nom10_cox<-nomogram(
      #cox10,fun = list(risk10 , surv10_1,surv10_2,surv10_3, surv10_4, surv10_5),
      cox10,fun = list(risk10 , surv10_1, surv10_4, surv10_5),
      lp = F,
                        funlabel = c("Risk","1-Year Survival Probability", "4-Year Survival Probability","5-Year Survival Probability"),
                        maxscale = 100,
                        fun.at = c('1.0', '0.95', '0.9','0.85','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1')  )
    plot((nom10_cox),xfrac = .9)
    
    jpeg("Nomogram_PI_Clin_features_STS_LTS__OS_time_COX.jpg", units="in", width=15, height=10, res=300)
    plot(nom10_cox)
    dev.off()
    
