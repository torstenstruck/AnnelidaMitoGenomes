#set working directory to the folder with the input files
#the following files are needed:
#Breakpoints_wtRNA.csv
#CommonInterval_wtRNA.csv
#ReverseDistance_wtRNA.csv

#load distance matrices
Breakpoints_wotRNA <- read.csv("Breakpoints_wtRNA.csv")
CommonInterval_wotRNA <- read.csv("CommonInterval_wtRNA.csv")
ReverseDistance_wotRNA <- read.csv("ReverseDistance_wtRNA.csv")

#convert to columns and merge in one matrix
Breakpoints_wotRNA_column <- unlist(Breakpoints_wotRNA)
CommonInterval_wotRNA_column <- unlist(CommonInterval_wotRNA)
ReverseDistance_wotRNA_column <- unlist(ReverseDistance_wotRNA)
CompiledData <- as.data.frame(cbind(Breakpoints_wotRNA_column, CommonInterval_wotRNA_column, ReverseDistance_wotRNA_column))
View(CompiledData)

#check correlation by plotting
library("ggpubr")
ggscatter(CompiledData, x = "Breakpoints_wotRNA_column", y = "CommonInterval_wotRNA_column", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "breakpoints", ylab = "common interval",
          ylim = c(0, 1500))
ggscatter(CompiledData, x = "ReverseDistance_wotRNA_column", y = "CommonInterval_wotRNA_column", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "reverse distance", ylab = "common interval",
          ylim = c(0, 1500))
ggscatter(CompiledData, x = "ReverseDistance_wotRNA_column", y = "Breakpoints_wotRNA_column", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "reverse distance", ylab = "breakpoints")

#Visual inspection of the data normality using Q-Q plots (quantile-quantile plots) as there are too many data points for the Shapiro-Wilk normality test
library("ggpubr")
ggqqplot(CompiledData$Breakpoints_wotRNA_column, ylab = "breakpoints")
ggqqplot(CompiledData$CommonInterval_wotRNA_column, ylab = "common interval")
ggqqplot(CompiledData$ReverseDistance_wotRNA_column, ylab = "reverse distance")

#Spearman's rho  rank correlation coefficient
cor.test(CompiledData$Breakpoints_wotRNA_column, CompiledData$CommonInterval_wotRNA_column,  method = "spearman")
cor.test(CompiledData$ReverseDistance_wotRNA_column, CompiledData$CommonInterval_wotRNA_column,  method = "spearman")
cor.test(CompiledData$ReverseDistance_wotRNA_column, CompiledData$Breakpoints_wotRNA_column,  method = "spearman")

