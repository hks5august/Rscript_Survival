library(survival)
library(glmnet)

setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/primary_vs_Recurrent/DGE/top_10pathways_survival/16_genes_models/")

set.seed(7)


library(survival)
data(CoxExample)
x <- CoxExample$x
y <- CoxExample$y
#y1 <- CoxExample$y

head(x,3)

head(y,3)



surv_object_tr


train <- read.table("sel_train.txt", sep="\t", header = TRUE, row.names=1, check.names = FALSE)
test1<- read.table("sel_test1.txt", sep="\t", header = TRUE, row.names=1, check.names = FALSE)
head(train[1:6],2)


y_tr <- train[1:2]
colnames(y_tr) <- c("time", "status")

# create survival object
surv_object_tr <- Surv(time = as.numeric(train$OS_month), event = train$Death_event)

data_tr <- train[3:ncol(train)]
data_te1 <- test1[3:ncol(test1)]
#data(CoxExample)
#x <- CoxExample$x ## x is expression matrix where samples in rows and co-variate in columns
#y <- CoxExample$y ## y is y is an n × 2 matrix, with a column "time" of failure/censoring times, and "status" a 0/1 indicator, with 1 meaning the time is a failure (death) time, and 0 a censoring time

#apply the glmnet function to compute the solution path under default settings:
  #fit <- glmnet(x, y, family = "cox")
  
fit <- glmnet(data_tr, y_tr, family = "cox")
  
  
# plot the coefficients with the plot method:
plot(fit)
head(fit)

# extract the coefficients at certain values of λ:
coef(fit, s = 0.05)

#cross validation model
set.seed(1)
cvfit1 <- cv.glmnet(as.matrix(data_tr), as.matrix(y_tr), family = "cox", type.measure = "C", nfolds = 5)

#view the optimal λ value and a cross validated error plot to help evaluate our model.
plot(cvfit1)

#extract such optimal λ
cvfit$lambda.min

glmnet_fit <- glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", lambda = 0.03055984)


cvfit$lambda.1se



set.seed(1)
cvfit1 <- cv.glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", type.measure = "C", nfolds = 5)
plot(cvfit1)


coxph_fit1 <- coxph(surv_object_tr ~ as.matrix(data_tr))
coxph_fit1
plot(coef(glmnet_fit), coef(coxph_fit1))
abline(0, 1)




set.seed(1)
cvfit1 <- cv.glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", type.measure = "C", nfolds = 5)

cvfit1

#view the optimal λ value and a cross validated error plot to help evaluate our model.

plot(cvfit1)

jpeg(file="Cox_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=300)
plot(cvfit1)
dev.off()
#extract such optimal λ
cvfit1$lambda.min

coef(cvfit1, s = 0.055)


cv.fit <- cv.glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", nfolds = 5)
cv.fit

plot(cv.fit)
#extract such optimal λ
cvfit$lambda.min


glmnet_fit <- glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", lambda = 0)
glmnet_fit <- glmnet(as.matrix(data_tr), surv_object_tr, family = "cox", lambda = 0.055)
coxph_fit <- coxph(surv_object_tr ~ as.matrix(data_tr), ties = "breslow")
plot(coef(glmnet_fit), coef(coxph_fit))
abline(0, 1)

survival::survfit(glmnet_fit, s = 0.05, x = as.matrix(data_tr), y = surv_object_tr)


plot(survival::survfit(glmnet_fit, s = 0.055, x = as.matrix(data_tr), y = surv_object_tr))

plot(survival::survfit(glmnet_fit, s = 0.055, x = as.matrix(data_tr), y = surv_object_tr, newx = as.matrix(data_te1)))
