setwd("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/03_PhylogeneticLeastSquareRegression/Molecular_data")

require(caper)
library("Hmisc")
library("corrplot")

#read in ultrametric tree and dataset and assign both for pgls
tree<-read.tree("Masked_18S_rerooted_ultrametric.treefile")
dt<-read.csv("CompiledProperties.csv")
dt_go_wt <- read.delim("GeneOrder_aligned_with_tRNA.tsv", row.names=NULL)
dt_go_wt$Gene_order <- as.numeric(as.factor(dt_go_wt$Gene_order))
dt_go_wo <- read.delim("GeneOrder_aligned_reducedMissingness.tsv", row.names=NULL)
dt_go_wo$Gene_order <- as.numeric(as.factor(dt_go_wo$Gene_order))

#generate list of column names
varList <- colnames(dt)

#determine phylogenetic signal for reverse distance with and without tRNA
dt_mitp_lambda1<-comparative.data(tree, dt, names.col="Taxon",warn.dropped=TRUE, na.omit=FALSE, vcv=TRUE, vcv.dim=3)
mod_lambda1<-pgls(Reverse_Dist_withRNA~1, data=dt_mitp_lambda1, lambda = 'ML')
summary(mod_lambda1)
lambda.mod1 <-pgls.profile(mod_lambda1, "lambda")
mod_lambda2<-pgls(Reverse_Dist_notRNA~1, data=dt_mitp_lambda1, lambda = 'ML')
summary(mod_lambda2)
lambda.mod2 <-pgls.profile(mod_lambda2, "lambda")

#determine phylogenetic signal directly for numbered gene orders with and without tRNA
dt_mitp_lambda2<-comparative.data(tree, dt_go_wt, names.col="Taxon",warn.dropped=TRUE, na.omit=FALSE, vcv=TRUE, vcv.dim=3)
mod_lambda3<-pgls(Gene_order~1, data=dt_mitp_lambda2, lambda = 'ML')
summary(mod_lambda3)
lambda.mod3 <-pgls.profile(mod_lambda3, "lambda")
dt_mitp_lambda3<-comparative.data(tree, dt_go_wo, names.col="Taxon",warn.dropped=TRUE, na.omit=FALSE, vcv=TRUE, vcv.dim=3)
mod_lambda4<-pgls(Gene_order~1, data=dt_mitp_lambda3, lambda = 'ML')
summary(mod_lambda4)
lambda.mod4 <-pgls.profile(mod_lambda4, "lambda")
par(mar=c(1,1,1,1))
par(mfrow =c(2, 2))
plot(lambda.mod1)
plot(lambda.mod2)
plot(lambda.mod3)
plot(lambda.mod4)


#gene order without tRNAs as response
#loop through each possible predictors, save only the AIC and p values and extract significant ones < 0.001
sets <- data.frame(response=character(0), predictor=character(0), AIC=numeric(0), P_value=numeric(0), adjusted_R2=numeric(0), stringsAsFactors=FALSE)
sig_sets <- data.frame(response=character(0), predictor=character(0), AIC=numeric(0), P_value=numeric(0), adjusted_R2=numeric(0), stringsAsFactors=FALSE)
cat("Predictor\nReverse_Dist_notRNA Response",file="lambda_fixed.txt",sep="\n")
for(i in 2:1203) {
 print(varList[i])
 dt_red <- subset(dt, select = c(1, i, 1205))
 dt_mitp<-comparative.data(tree, dt_red, names.col="Taxon",warn.dropped=FALSE, na.omit=TRUE, vcv=TRUE, vcv.dim=3)
 mod <- tryCatch(pgls(Reverse_Dist_notRNA ~ dt_mitp$data[,1], data = dt_mitp, lambda = 'ML'), 
                 error = function(e) {
                   print("varList[i]: fixed lambda used")
                   cat(varList[i], file="lambda_fixed.txt", sep="\n", append=TRUE)
                   pgls(Reverse_Dist_notRNA ~ dt_mitp$data[,1], data = dt_mitp, lambda = 1)
                 }
 )
 sum_mod <- summary(mod)
 res_mod <- cbind("Reverse_Dist_notRNA", varList[i], mod$aic, sum_mod$coefficients[2,4], sum_mod$adj.r.squared)
 sets <- rbind(sets, res_mod)
 if(sum_mod$coefficients[2,4] < 0.001) {sig_sets <- rbind(sig_sets, res_mod)}
}
colnames(sets) <- c("response", "predictor", "AIC", "P_Values", "adjusted_R2")
colnames(sig_sets) <- c("response", "predictor", "AIC", "P_Values", "adjusted_R2")
write.table(sets, "All_predictors_without_tRNAs.txt", sep="\t")
write.table(sig_sets, "Significant_predictors_without_tRNAs.txt", sep="\t")

#subset the set of predictor with significant p values
dt_sig <- subset(dt, select = c(sig_sets$predictor))

#Generate a hierarchical cluster and corellogarm of the correlation coefficients of the significant ones
corr_dt_sig <- rcorr(as.matrix(na.omit(dt_sig)))
corr_dt_sig_r <- corr_dt_sig$r
hc <- hclust(dist(corr_dt_sig_r), "ave")
dev.off()
plot(hc,  hang = -1, cex = 0.3)
corrplot.mixed(as.matrix(corr_dt_sig_r), number.cex = .7, tl.col = "black", tl.pos = "lt", tl.cex = .7, diag = "u", order="hclust", hclust.method = "average", addgrid.col = "NA")

#gene order with tRNAs as response
#loop through each possible predictors, save only the AIC and p values and extract significant ones < 0.001
sets <- data.frame(response=character(0), predictor=character(0), AIC=numeric(0), P_value=numeric(0), adjusted_R2=numeric(0), stringsAsFactors=FALSE)
sig_sets <- data.frame(response=character(0), predictor=character(0), AIC=numeric(0), P_value=numeric(0), adjusted_R2=numeric(0), stringsAsFactors=FALSE)
cat("Reverse_Dist_withRNA Response",file="lambda_fixed.txt",sep="\n", append=TRUE)
for(i in 2:1203) {
  print(varList[i])
  dt_red <- subset(dt, select = c(1, i, 1204))
  dt_mitp<-comparative.data(tree, dt_red, names.col="Taxon",warn.dropped=FALSE, na.omit=TRUE, vcv=TRUE, vcv.dim=3)
  mod <- tryCatch(pgls(Reverse_Dist_withRNA ~ dt_mitp$data[,1], data = dt_mitp, lambda = 'ML'), 
                  error = function(e) {
                    print("varList[i]: fixed lambda used")
                    cat(varList[i], file="lambda_fixed.txt", sep="\n", append=TRUE)
                    pgls(Reverse_Dist_withRNA ~ dt_mitp$data[,1], data = dt_mitp, lambda = 1)
                  }
  )
  sum_mod <- summary(mod)
  res_mod <- cbind("Reverse_Dist_withRNA", varList[i], mod$aic, sum_mod$coefficients[2,4], sum_mod$adj.r.squared)
  sets <- rbind(sets, res_mod)
  if(sum_mod$coefficients[2,4] < 0.001) {sig_sets <- rbind(sig_sets, res_mod)}
}
colnames(sets) <- c("response", "predictor", "AIC", "P_Values", "adjusted_R2")
colnames(sig_sets) <- c("response", "predictor", "AIC", "P_Values", "adjusted_R2")
write.table(sets, "All_predictors_with_tRNAs.txt", sep="\t")
write.table(sig_sets, "Significant_predictors_with_tRNAs.txt", sep="\t")

#subset the set of predictor with significant p values
dt_sig <- subset(dt, select = c(sig_sets$predictor))

#Generate a hierarchical cluster and corellogarm of the correlation coefficients of the significant ones
corr_dt_sig <- rcorr(as.matrix(na.omit(dt_sig)))
corr_dt_sig_r <- corr_dt_sig$r
hc <- hclust(dist(corr_dt_sig_r), "ave")
plot(hc,  hang = -1, cex = 0.3)
corrplot.mixed(as.matrix(corr_dt_sig_r), number.cex = .7, tl.col = "black", tl.pos = "lt", tl.cex = .7, diag = "u", order="hclust", hclust.method = "average", addgrid.col = "NA")

#generate correlation networks to determine connected predictors
library(tidyverse)  
library(corrr)

#gene order without tRNAs
Significant_predictors_without_tRNAS <- read.delim("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/03_PhylogeneticLeastSquareRegression/Molecular_data/Significant_predictors_without_tRNAs.txt")
dt_sig <- subset(dt, select = c(Significant_predictors_without_tRNAS$predictor))
dt_sig %>% correlate() %>% network_plot(min_cor = 0.7, repel = FALSE, curved = FALSE)

#gene order with tRNAs
Significant_predictors_with_tRNAS <- read.delim("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/03_PhylogeneticLeastSquareRegression/Molecular_data/Significant_predictors_with_tRNAs.txt")
dt_sig <- subset(dt, select = c(Significant_predictors_with_tRNAS$predictor))
dt_sig %>% correlate() %>% network_plot(min_cor = 0.7, repel = FALSE, curved = FALSE)

#testing the combined information of the best-fitting variable in the three categories (Dist, RCFV/Freq, Num) present in both with and without tRNA
dt_comb<-comparative.data(tree, dt, names.col="Taxon",warn.dropped=FALSE, na.omit=FALSE, vcv=TRUE, vcv.dim=3)
#best variable from the pool of shared variables, chosen on AIC
#without tRNAs
#including just the three variables
mod_without_comb <- pgls(Reverse_Dist_notRNA ~ (COX3_Nuc_Dist_TR + Whole_Nuc_Num_IG + NAD3_AA_Freq_H), data = dt_comb, lambda = 'ML')
summary(mod_without_comb)
#considering also the interaction between them
mod_without_comb_int <- pgls(Reverse_Dist_notRNA ~ (COX3_Nuc_Dist_TR * Whole_Nuc_Num_IG * NAD3_AA_Freq_H), data = dt_comb, lambda = 'ML')
summary(mod_without_comb_int)
#with tRNAs
#including just the three variables
mod_with_comb <- pgls(Reverse_Dist_withRNA ~ (ATP6_Nuc_Dist_TR + Whole_Nuc_Num_IG + NAD3_AA_Freq_H), data = dt_comb, lambda = 'ML')
summary(mod_with_comb)
#considering also the interaction between them
mod_with_comb_int <- pgls(Reverse_Dist_withRNA ~ (ATP6_Nuc_Dist_TR * Whole_Nuc_Num_IG * NAD3_AA_Freq_H), data = dt_comb, lambda = 'ML')
summary(mod_with_comb_int)
#AIC values of each
mod_without_comb$aic
mod_without_comb_int$aic
mod_with_comb$aic
mod_with_comb_int$aic

#best variable from each pool independently, chosen on AIC
#without tRNAs
#including just the three variables
#same as the shared ones
#considering also the interaction between them
#same as the shared ones
#with tRNAs
#including just the three variables
mod_with_comb <- pgls(Reverse_Dist_withRNA ~ (ATP6_Nuc_Dist_TR + Whole_Nuc_Num_IG + ATP8_Nuc_Skew_GC), data = dt_comb, lambda = 'ML')
summary(mod_with_comb)
#considering also the interaction between them
mod_with_comb_int <- pgls(Reverse_Dist_withRNA ~ (ATP6_Nuc_Dist_TR * Whole_Nuc_Num_IG * ATP8_Nuc_Skew_GC), data = dt_comb, lambda = 'ML')
summary(mod_with_comb_int)
#AIC values of each
#mod_without_comb$aic
#mod_without_comb_int$aic
mod_with_comb$aic
mod_with_comb_int$aic

