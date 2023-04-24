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
setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/survival/16_genes_Survival_Plots")
set.seed(7)


tr_data <- read.table("16_genes_train.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te_data <- read.table("16_genes_test.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
ext_data <- read.table("16_genes_ext.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

tr_data  <- tr_data 
#Survival analysis

tr_data1<-subset(tr_data,OS_month!="NA")
dim(tr_data1)
head(tr_data1)

# create survival object
surv_object<- Surv(time = as.numeric(tr_data1$OS_month), event = tr_data1$Death_event)


dim(tr_data1)
tr_data1[2]
#for(i in seq(from=3, to=length(tr_data1), by=1))

for(i in seq(from=3, to=18, by=1))

  {
gene <-  colnames(tr_data1[i])
  
surv_fit <- survfit(surv_object ~ (tr_data1[,i])>(median(tr_data1[1,i])), data=tr_data1);
summary(surv_fit);
#ggsurvplot(surv_fit)
fit1.coxph <- coxph(surv_object ~ (tr_data1[,i])>(median(tr_data1[1,i])), data = tr_data1)
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
HR <- round(first[2],2)
options(scipen = 1)
options(digits = 3)

pval <- format(first[5], scientific = T)
#?ggsurvplot
# we produce a Kaplan Meier plot
#jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(plot.title=element_text(hjust=0.5,face="bold", size=30))}

surv_plot1 <- ggsurvplot(surv_fit, data=tr_data1, pval=F, 
           #risk.table=TRUE, 
           risk.table=F, 
           tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           #risk.table.col="strata", break.time.by = 12,
           title=paste0("KM Plot for ", gene, "\n"),
           ggtheme=custom_theme(),
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           legend.labs=c(paste0(">Median Exp of ", gene), paste0("<Median Exp of ", gene)),
           #palette = c("red","blue"),#add desired color
           size=2.5, font.tickslab = c(30,  "black"), font.y = c(30,  "black"), 
           font.x = c(30, "black"), font.legend = c(25, "black"), font.label = c(30, "bold", "black"))



#Add HR and P-values
surv_plot2  <- surv_plot1$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 2, hjust = 1,
    label =  paste0("HR = ", HR, "\n", "p-val = ", pval),
    size = 10
  ) 

print(gene)
surv_plot2
#save plot

#jpeg(file= paste0(gene,"_KM_plot.jpeg"), units="in", width=10, height=10, res=300)
#surv_plot2
#dev.off()

ggsave(surv_plot2, file=paste0(gene,"_KM_plot.jpeg"), units="in", width=12, height=12, dpi =350)


}
gene



