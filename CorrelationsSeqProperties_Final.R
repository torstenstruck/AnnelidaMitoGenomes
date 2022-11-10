#set working directory to the folder with the input files
#files needed:
#CompiledProperties.csv

library("ggpubr")
library("corrplot")
library("Hmisc")
require("fastcluster")
require("graphics")
library("lares")
library("corrr")
library("data.table")
library("tidyverse")  
library("dplyr")

#Load data
CompiledSeqProperties <- read.csv("CompiledProperties.csv", header = TRUE, row.names=1)
View(CompiledSeqProperties)

#calculate correlation coefficents and plot the distribution of correlation coefficents
res2 <- rcorr(as.matrix(na.omit(CompiledSeqProperties)))
resNA <- res2$r
diag(resNA) <- NA
View(resNA)
plot(density(na.omit(as.numeric(resNA))))

#extract all correlated pairs, which are highly correlated (>0.5 or <-0.5)
res0 <- res2$r
diag(res0) <- 0

cor.pair <- data.frame(v1=character(0), v2=character(0), cor=numeric(0), stringsAsFactors=FALSE)
for(i in 1:nrow(res0)) { #rownumber
  for(j in 1:ncol(res0)) {  #columnnumber 
    testvalue <- abs(res0[i,j])
    if(testvalue > 0.5) {
      cor.pair <- rbind(cor.pair, data.frame(v1=rownames(res0)[i], v2=colnames(res0)[j], cor=res0[i,j]))
      res0[,j] <- 0
      res0[j,i] <- 0
    } 
  }
}
View(res0)

#remove all columns with only 0 values
resReduced <- res0[,colSums(res0)!=0]
View(resReduced)

#write out table with highly correlated pairs found and determine different section parameters for pairs
write.table(cor.pair, "Correlation_HighPairs.txt", sep="\t")
sets <- data.frame(condition=character(0), unique_v1=numeric(0), unique_v2=numeric(0), intersect_v1_v2=numeric(0), diff_v1_v2=numeric(0), diff_v2_v1=numeric(0), abs_min_cor=numeric(0), stringsAsFactors=FALSE)
sets <- rbind(sets, data.frame(
  condition="HighCorrPairs",
  unique_v1=length(unique(cor.pair$v1)),
  unique_v2=length(unique(cor.pair$v2)),
  intersect_v1_v2=length(intersect(cor.pair$v1,cor.pair$v2)),
  diff_v1_v2=length(setdiff(cor.pair$v1,cor.pair$v2)),
  diff_v2_v1=length(setdiff(cor.pair$v2,cor.pair$v1)),
  abs_min_cor=min(abs(cor.pair$cor))
))
write.table(sets, "Correlation_Pair_Sets.txt", sep="\t")

#reduce original data matrix to the remaining columns
CompiledSeqProperties_Corr_Clean <- CompiledSeqProperties %>% select(colnames(resReduced))
View(CompiledSeqProperties_Corr_Clean)
write.table(CompiledSeqProperties_Corr_Clean, "CompiledSeqProperties_Corr_Clean.txt", sep="\t")

#calculate correlation coefficents plus significant values from first reduced data matrix
res3 <- rcorr(as.matrix(na.omit(CompiledSeqProperties_Corr_Clean)))
resClean <- res3$r

#calculate hierarchical cluster of correlation coefficients of the remaining ones after extracting highly correlated ones
hc <- hclust(dist(resClean), "ave")
#plot the original hierarchical cluster
plot(hc,  hang = -1, cex = 0.3)
#Generate a corellogarm for the same
corrplot(as.matrix(resClean), type="upper", order="hclust", hclust.method = "average", addgrid.col = "NA", tl.cex = 0.3)

#determine groups of clustered parameters in the hierarchical cluster based on a cutoff values for height(80% of maximal height), sort them and save all except one in a variable for removing
threshold <- (max(hc$height) * 0.8)
hclust_groups <- sort(cutree(hc, h = threshold))
RemoveHclust <- c()
for(i in 1:max(hclust_groups)) {
  RemoveHclustTemp <- hclust_groups[hclust_groups == i] 
  RemoveHclustTemp <- RemoveHclustTemp[-1]
  RemoveHclust <- c(RemoveHclust,RemoveHclustTemp)
}
RemoveHclust <- stack(RemoveHclust)

#reduce clustered parameters to one
resCleaning <- as.data.frame(resClean)
resCleaning[ ,c(as.character(RemoveHclust$ind))] <- list(NULL)
setDT(resCleaning, keep.rownames = TRUE)
setkey(resCleaning, rn)
Remove <- c(as.character(RemoveHclust$ind))
resCleaning <- resCleaning[!Remove]

#prepare data for further analyses
resReduced <- as.data.frame(resCleaning)
row.names(resReduced) <- resReduced$rn
resReduced <- resReduced[, -1]
resReduced <- resReduced[, order(names(resReduced))]
resReduced <- resReduced[order(rownames(resReduced)), ]

#calculate hierarchical cluster of correlation coefficients of the reduced dataset
hc <- hclust(dist(resReduced), "ave")
#plot the original hierarchical cluster
plot(hc,  hang = -1)
#Generate a corellogarm for the same
corrplot(as.matrix(resReduced), type="upper", order="hclust", hclust.method = "average", addgrid.col = "NA")

#reduce reduced data matrix further to the remaining columns
CompiledSeqProperties_Reduced <- CompiledSeqProperties %>% select(colnames(resReduced))
View(CompiledSeqProperties_Reduced)
write.table(CompiledSeqProperties_Reduced, "CompiledSeqProperties_Reduced.txt", sep="\t")
