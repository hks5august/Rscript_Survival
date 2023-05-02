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
setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/updated_folder2/ssGSEA_oncogenic/Hallmark_Pathways/")


path <- paste0("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/updated_folder2/ssGSEA_oncogenic/Hallmark_Pathways/")



setwd(path)
#select cancer type
cancer <- as.character("GBM")
cancer <- as.character("LGG")
cancer <- as.character(args[1])
cancer

setwd(paste0(path,cancer, "/"))

#getwd()

#data<-read.table(paste0(cancer,".csv"),header =TRUE, sep = ",",  check.names = FALSE)
data<-read.table(paste0(cancer,".tsv"),header =TRUE, sep = "\t",  check.names = FALSE, row.names = 1)
print("data")
dim(data)
head(data)


head(data[19:21],2)
#extract clinical data
Clin_data <- data[1:19]
print("dim of clin data:")
dim(Clin_data)
Clin_data1 <- Clin_data %>% mutate(OS_month = round(OS.time/30.417, digit=0))

#extract expression data
exp_data <- data[20:ncol(data)]
print("dim of exp data:")
dim(exp_data)

#combine new clinical data with exp data
data2 <- cbind(Clin_data1, exp_data)
data3 <- subset(data2, OS!="NA")
train_mat <- subset(data3, OS_month!="NA")
train_mat <- subset(train_mat, OS_month >0)
dim(train_mat)




#select/provide pathway
pathway <- as.character("APOPTOSIS")
pathway <- as.character(args[2])
print("pathway name:")
pathway

#features<- read.csv(paste0(path2,pathway, ".txt"), header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)

#create directory for each pathway then create prognostic index analysis for each pathway
folder1<-dir.create(paste0(path, cancer,"_temp_", pathway)) ### create new directory to store results
print("folder1:")
folder1

setwd(paste0(path, cancer,"_temp_", pathway))


# Create Survival object
surv_object_tr <- Surv(time = train_mat$OS_month, event = train_mat$OS)


# ###survival analysis: fits cox ph model to find HR for PI
#fit_tr <- survfit(surv_object_tr~(train_mat$Pathway>mean(train_mat$Pathway)), data=train_mat)

fit_tr <- survfit(surv_object_tr~(train_mat[, pathway]>mean(train_mat[, pathway])), data=train_mat)

summary(fit_tr)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_tr <- coxph(surv_object_tr ~(train_mat[, pathway]>mean(train_mat[, pathway])), data=train_mat)
summary(fit.coxph_tr)


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
write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="train_result.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=train_mat, 
                  #pval=TRUE,
                  risk.table=TRUE, 
                  tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  title= paste0("KM plot based on ssGSEA score of ", pathway , " in ", cancer),
                  surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  legend.title = paste0(pathway),
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
    #label =  paste0("HR = ", HR1, "\n",  "p-val = ", pval1),
    label =  paste0("HR = ", HR1, "\n",  "p-val = ", pval1,  "\n",  "C-Index = ", CI),
size = 5)

# now plot
pp


jpeg("KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()





#Create ROC plot for 1-,3- and 5-years survival time prediction

#Calculated  Lmbda min value = 0.002102803

tr_roc1 <- survivalROC(Stime        = train_mat$OS_month,
                       status       = train_mat$OS,
                       marker       = train_mat[, pathway],
                       predict.time = 12,
                       method       = "KM",
                       #method       = "NNE",
                       #lambda = lambda_min ,
                       span = NULL, window ="symmetric")
tr_roc1 

tr_roc3 <- survivalROC(Stime        = train_mat$OS_month,
                       status       = train_mat$OS,
                       marker       = train_mat[, pathway],
                       predict.time = 36,
                       method       = "KM",
                       #method       = "NNE",
                       #lambda = lambda_min,
                       span = NULL, window ="symmetric")

tr_roc5 <- survivalROC(Stime        = train_mat$OS_month,
                       status       = train_mat$OS,
                       marker       = train_mat[, pathway],
                       predict.time = 60,
                       method       = "KM",
                       #method       = "NNE",
                       #lambda = lambda_min,
                       span = NULL, window ="symmetric")

jpeg(file="ROC.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",
     #main= "AUC Curve for Survival Prediction")
     main= paste0("AUC Curve for Survival Prediction in ", cancer, " based on ssGSEA of ",pathway))
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)

dev.off()
