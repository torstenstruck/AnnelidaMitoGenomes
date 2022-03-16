#load distance matrices
Breakpoints_wotRNA <- read.csv("C:/Users/torsths/Documents/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/01_CorrelationAnalyses/GeneOrder/without_tRNAs/Breakpoints_wotRNA.csv")
CommonInterval_wotRNA <- read.csv("C:/Users/torsths/Documents/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/01_CorrelationAnalyses/GeneOrder/without_tRNAs/CommonInterval_wotRNA.csv")
ReverseDistance_wotRNA <- read.csv("C:/Users/torsths/Documents/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/04_Macroevolution/01_CorrelationAnalyses/GeneOrder/without_tRNAs/ReverseDistance_wotRNA.csv")

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
          ylim = c(0, 250))
ggscatter(CompiledData, x = "ReverseDistance_wotRNA_column", y = "CommonInterval_wotRNA_column", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "reverse distance", ylab = "common interval",
          ylim = c(0, 250))
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

