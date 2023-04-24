
set.seed(7)

#path where we have input data all files and folders
path <- paste0("/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/")

setwd(path)
#select cancer type
cancer <- as.character("LIHC")
cancer <- as.character(args[1])
cancer

setwd(paste0(path,cancer, "/"))

getwd()

data<-read.table(paste0(path,cancer,".csv"),header =TRUE, sep = ",", row.names=1,  check.names = FALSE)
dim(data)

######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
all_features_results<-read.csv("Significant_Survival_results_for_genes.txt",header =TRUE, sep = "\t", check.names = FALSE)
dim(all_features_results)
head(all_features_results)

#provide path of folder containing  files having list of genes for each pathways
path2 <- "/Users/kaurh8/Documents/Pathways_Survival_analysis_all_cancer/All_pathways/"

#select/provide pathway
#pathway <- as.character("Apoptosis")
pathway <- as.character("mTOR_signaling_pathway")

features<- read.csv(paste0(path2,pathway, ".txt"), header =TRUE, sep = "\t", dec = ".",stringsAsFactors=FALSE, check.names = FALSE)

dim(features)


sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]
sel_features_results

dim(sel_features_results)


#create data frame for selected features
sel_ftr_surv <- as.data.frame(sel_features_results[,1])
names(sel_ftr_surv ) <- c("ID")
head(sel_ftr_surv )

# calculate number of genes for enrichment
total_g <- 20501 #total number of genes in the data
path_g <- as.numeric(nrow(features)) # number of genes in a specific pathway
sig_g <- as.numeric(nrow(all_features_results)) #number of total significant prognostic genes in a specific cancer 
sel_g <- as.numeric(nrow(sel_features_results)[1]) ##number of total significant genes found common in pathway

sig_g1 <- sig_g - sel_g #no .pf significant genes which are not in specific pathway
path_g1 <- path_g - sel_g #pathway genes which are not in significant
rest_g <-  path_g1 + sel_g + sig_g1 #significant and genes of pathway
total2_g <- total_g - rest_g  ## genes which are not significant and not in pathway
perc_path <- (sel_g/path_g) * 100
perc_path


#prepare contigenecy table
contingencyTable <- data.frame(
  Sig_Yes=c(sel_g, sig_g1),
  Sig_No=c(path_g1, total2_g))

#row names
row.names(contingencyTable) <- c("PathwayYes", "PathwayNO")
#print contigenecy table
contingencyTable

#apply fisher exact test
Fisher_res <- fisher.test(contingencyTable, alternative = "greater") #one sided
pval <- Fisher_res$p.value
odd_ratio <- Fisher_res$estimate

#create dataframe considering important things in fisher test
path_Fisher_test <- cbind(path_g, sel_g, perc_path, pval, odd_ratio[1])
#colnames
colnames(path_Fisher_test) <- c("Total_path_genes" , "Path_sig_Yes", "Percentage", "p.value", "odd_ratio")
#row names
row.names(path_Fisher_test) <- NULL
path_Fisher_test


#######################################







sel_g = 100
sig_g1= 782
path_g1 = 54
total2_g = 19064

contingencyTable <- data.frame(
  Sig_Yes=c(sel_g, sig_g1),
  Sig_No=c(path_g1, total2_g))


row.names(contingencyTable) <- c("PathwayYes", "PathwayNO")
contingencyTable


Fisher_res <- fisher.test(contingencyTable, alternative = "greater") #one sided
Fisher_res
pval <- Fisher_res$p.value
odd_ratio <- Fisher_res$estimate
odd_ratio
#Fisher_res2 <- fisher.test(contingencyTable, alternative = "two.sided") #two sided test
#Fisher_res2

?fisher.test

pval <- Fisher_res$p.value
odd_ratio <- Fisher_res$estimate
odd_ratio

Fisher_res$null.value
?fisher.test
m <-  path_g      # Genes IN GO term
m
n <- total_g  - path_g      # Genes NOT IN GO term
n
k <-   sig_g     # Gene hits, that is, differentially expressed
k
x <- sel_g  # Genes both IN GO term and differentially expressed 'hits'
x
# Use the dhyper built-in function for hypergeometric density
probabilities <- dhyper(x, m, n, k, log = FALSE)

probabilities
pvalue <- sum(probabilities)
pvalue
probabilities
