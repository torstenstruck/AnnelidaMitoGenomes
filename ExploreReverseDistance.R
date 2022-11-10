#set working directory to the folder with input files
#the following files are needed
#without_tRNAs/ReverseDistance_wotRNA.csv
#with_tRNAs/ReverseDistance_wtRNA.csv
#Family_Species_list.csv
#Nuc_supermatrix_partition_rooted_cladogram.txt.treefile
#PB_Nuc_chain2_rooted_cladogram.treefile
#Masked_18S_rooted_cladogram.treefile
#AA_supermatrix_partition_rooted_cladogram.txt.treefile
#PB_AA_chain1_rooted_cladogram.treefile
#Nuc_supermatrix_partition_rooted_cladogram.txt
#PB_Nuc_chain2_rooted_cladogram.txt
#Masked_18S_rooted_cladogram.txt
#AA_supermatrix_partition_rooted_cladogram.txt
#PB_AA_chain1_rooted_cladogram.txt

library("tidyr")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("ape")
library("TreeDist")
library("TreeTools")


#load distance matrices
ReverseDistance_wotRNA <- read.csv("without_tRNAs/ReverseDistance_wotRNA.csv")
ReverseDistance_wtRNA <- read.csv("with_tRNAs/ReverseDistance_wtRNA.csv")
Family_species <- read.csv("Family_Species_list.csv")

#compile data into single columns
RD_wo_long <- ReverseDistance_wotRNA %>%  pivot_longer(., cols = c(1:ncol(ReverseDistance_wotRNA)), names_to = "Species_wo", values_to = "RD_wo_tRNA")
RD_w_long <- ReverseDistance_wtRNA %>%  pivot_longer(., cols = c(1:ncol(ReverseDistance_wtRNA)), names_to = "Species_w", values_to = "RD_w_tRNA")
#compile with and without together
RD_both_long <- na.omit(cbind(RD_wo_long,RD_w_long))
#smaller subset columns based on only RD values
RD_both_values <- select(RD_both_long, matches("RD_"))
#make to single column
RD_both_values_long <- RD_both_values %>%  pivot_longer(., cols = c(1:ncol(RD_both_values)), names_to = "tRNA", values_to = "RD_value")

#generate boxplots of all Reverse Distances with and without tRNA
ggplot(RD_both_values_long, aes(x = tRNA, y = RD_value)) +
  geom_boxplot()

#generate data frame with mean values for reverse distances
Mean_RD_wo <- colMeans(ReverseDistance_wotRNA)
Mean_RD_w <- colMeans(ReverseDistance_wtRNA)
Mean_RDs <- as.data.frame(cbind(Mean_RD_wo,Mean_RD_w))

#make to single column of only mean values
Mean_RDs_long <- Mean_RDs %>%  pivot_longer(., cols = c(1:2), names_to = "tRNA", values_to = "RD_value")

#generate boxplots of mean Reverse Distances with and without tRNA
ggplot(Mean_RDs_long, aes(x = tRNA, y = RD_value)) +
  geom_boxplot()

#plot mean reverse distances with and without tRNA against each other
View(Mean_RDs)
ggscatter(Mean_RDs, x = "Mean_RD_w", y = "Mean_RD_wo", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_w", ylab = "Mean_RD_wo")

#add species as the first column to mean values
Mean_RDs_species <- cbind(rownames(Mean_RDs), data.frame(Mean_RDs, row.names=NULL))
colnames(Mean_RDs_species)[1] <- "Species"

#add families to the species
Mean_RDs_species_family <- merge(Family_species,Mean_RDs_species,by="Species")

#plot mean reverse distances with and without tRNA against each other and color by family
ggplot(data = Mean_RDs_species_family, mapping = aes(x = Mean_RD_w, y = Mean_RD_wo,color = Family)) +
  geom_line() +
  geom_point(aes(alpha = 0.5))

#compare mean reverse distances to difference in trees of nucleotide alignment 
tree_nuc_uncon <- read.tree(file = "Nuc_supermatrix_partition_rooted_cladogram.txt.treefile");
tree_nuc_PB <- read.tree(file = "PB_Nuc_chain2_rooted_cladogram.treefile");
tree_ref<- read.tree(file = "Masked_18S_rooted_cladogram.treefile");

#generate distance matrix of nuclear trees
Matrix_nuc_uncon <- cophenetic(tree_nuc_uncon)
Matrix_nuc_uncon <- Matrix_nuc_uncon[sort(rownames(Matrix_nuc_uncon)),sort(colnames(Matrix_nuc_uncon))] 
Matrix_nuc_PB <- cophenetic(tree_nuc_PB)
Matrix_nuc_PB <- Matrix_nuc_PB[sort(rownames(Matrix_nuc_PB)),sort(colnames(Matrix_nuc_PB))] 
Matrix_ref <- cophenetic(tree_ref)
Matrix_ref <- Matrix_ref[sort(rownames(Matrix_ref)),sort(colnames(Matrix_ref))] 

Matrix_nuc_diff_uncon <- abs(Matrix_nuc_uncon - Matrix_ref)
Matrix_nuc_diff_PB <- abs(Matrix_nuc_PB - Matrix_ref)

#generate data frame with mean values for tree difference and combine with mean of reverse distances
Mean_nuc_Tree_uncon <- colMeans(Matrix_nuc_diff_uncon)
Mean_RDs_nuc_Tree_uncon <- as.data.frame(cbind(Mean_RD_wo,Mean_RD_w,Mean_nuc_Tree_uncon))
View(Mean_RDs_nuc_Tree_uncon)
Mean_nuc_Tree_PB <- colMeans(Matrix_nuc_diff_PB)
Mean_RDs_nuc_Tree_PB <- as.data.frame(cbind(Mean_RD_wo,Mean_RD_w,Mean_nuc_Tree_PB))
View(Mean_RDs_nuc_Tree_PB)

#plot mean reverse distances with and without tRNA against mean tree difference
ggscatter(Mean_RDs_nuc_Tree_uncon, x = "Mean_RD_w", y = "Mean_nuc_Tree_uncon", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_w", ylab = "Mean_nuc_Tree_uncon")
ggscatter(Mean_RDs_nuc_Tree_uncon, x = "Mean_RD_wo", y = "Mean_nuc_Tree_uncon", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_wo", ylab = "Mean_nuc_Tree_uncon")
ggscatter(Mean_RDs_nuc_Tree_PB, x = "Mean_RD_w", y = "Mean_nuc_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_w", ylab = "Mean_nuc_Tree_PB")
ggscatter(Mean_RDs_nuc_Tree_PB, x = "Mean_RD_wo", y = "Mean_nuc_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_wo", ylab = "Mean_nuc_Tree_PB")

#compare mean reverse distances to difference in trees of amino acid alignment 
tree_aa_uncon <- read.tree(file = "AA_supermatrix_partition_rooted_cladogram.txt.treefile");
tree_aa_PB <- read.tree(file = "PB_AA_chain1_rooted_cladogram.treefile");

#generate distance matrix of nuclear trees
Matrix_aa_uncon <- cophenetic(tree_aa_uncon)
Matrix_aa_uncon <- Matrix_aa_uncon[sort(rownames(Matrix_aa_uncon)),sort(colnames(Matrix_aa_uncon))] 
Matrix_aa_PB <- cophenetic(tree_aa_PB)
Matrix_aa_PB <- Matrix_aa_PB[sort(rownames(Matrix_aa_PB)),sort(colnames(Matrix_aa_PB))] 

Matrix_aa_diff_uncon <- abs(Matrix_aa_uncon - Matrix_ref)
Matrix_aa_diff_PB <- abs(Matrix_aa_PB - Matrix_ref)

#generate data frame with mean values for tree difference and combine with mean of reverse distances
Mean_aa_Tree_uncon <- colMeans(Matrix_aa_diff_uncon)
Mean_RDs_aa_Tree_uncon <- as.data.frame(cbind(Mean_RD_wo,Mean_RD_w,Mean_aa_Tree_uncon))
View(Mean_RDs_aa_Tree_uncon)
Mean_aa_Tree_PB <- colMeans(Matrix_aa_diff_PB)
Mean_RDs_aa_Tree_PB <- as.data.frame(cbind(Mean_RD_wo,Mean_RD_w,Mean_aa_Tree_PB))
View(Mean_RDs_aa_Tree_PB)

#plot mean reverse distances with and without tRNA against mean tree difference
ggscatter(Mean_RDs_aa_Tree_uncon, x = "Mean_RD_w", y = "Mean_aa_Tree_uncon", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_w", ylab = "Mean_aa_Tree_uncon")
ggscatter(Mean_RDs_aa_Tree_uncon, x = "Mean_RD_wo", y = "Mean_aa_Tree_uncon", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_wo", ylab = "Mean_aa_Tree_uncon")
ggscatter(Mean_RDs_aa_Tree_PB, x = "Mean_RD_w", y = "Mean_aa_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_w", ylab = "Mean_aa_Tree_PB")
ggscatter(Mean_RDs_aa_Tree_PB, x = "Mean_RD_wo", y = "Mean_aa_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_RD_wo", ylab = "Mean_aa_Tree_PB")

#generate data frame with mean values for tree difference
Mean_aa_nuc_Tree <- as.data.frame(cbind(Mean_aa_Tree_uncon,Mean_nuc_Tree_uncon,Mean_aa_Tree_PB,Mean_nuc_Tree_PB))
View(Mean_aa_nuc_Tree)

#make to single column of only mean values
Mean_aa_nuc_Tree_long <- Mean_aa_nuc_Tree %>%  pivot_longer(., cols = c(1:4), names_to = "category", values_to = "value")

#generate boxplots of mean tree differences for different analyses
ggplot(Mean_aa_nuc_Tree_long, aes(x = category, y = value)) +
  geom_boxplot()

#plot different mean tree differences against each other
ggscatter(Mean_aa_nuc_Tree, x = "Mean_aa_Tree_uncon", y = "Mean_nuc_Tree_uncon", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_aa_Tree_uncon", ylab = "Mean_nuc_Tree_uncon")
ggscatter(Mean_aa_nuc_Tree, x = "Mean_aa_Tree_PB", y = "Mean_nuc_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_aa_Tree_PB", ylab = "Mean_nuc_Tree_PB")
ggscatter(Mean_aa_nuc_Tree, x = "Mean_aa_Tree_uncon", y = "Mean_aa_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_aa_Tree_uncon", ylab = "Mean_aa_Tree_PB")
ggscatter(Mean_aa_nuc_Tree, x = "Mean_nuc_Tree_uncon", y = "Mean_nuc_Tree_PB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean_nuc_Tree_uncon", ylab = "Mean_nuc_Tree_PB")

#compare mean reverse distances to difference in trees of nucleotide alignment 
tree_nuc_uncon <- read.tree(file = "Nuc_supermatrix_partition_rooted_cladogram.txt");
tree_nuc_PB <- read.tree(file = "PB_Nuc_chain2_rooted_cladogram.txt");
tree_ref<- read.tree(file = "Masked_18S_rooted_cladogram.txt");

#visualize similarity between trees using shared phylogenetic information showing matching clades
VisualizeMatching(RobinsonFouldsMatching, tree_nuc_uncon, tree_ref)
RobinsonFoulds(tree_nuc_uncon, tree_ref) # RF = 196
VisualizeMatching(RobinsonFouldsMatching, tree_nuc_PB, tree_ref)
RobinsonFoulds(tree_nuc_PB, tree_ref) # RF = 162

#compare mean reverse distances to difference in trees of amino acid alignment 
tree_aa_uncon <- read.tree(file = "AA_supermatrix_partition_rooted_cladogram.txt");
tree_aa_PB <- read.tree(file = "PB_AA_chain1_rooted_cladogram.txt");

#visualize similarity between trees using shared phylogenetic information showing matching clades
VisualizeMatching(RobinsonFouldsMatching, tree_aa_uncon, tree_ref)
RobinsonFoulds(tree_aa_uncon, tree_ref) # RF = 174
VisualizeMatching(RobinsonFouldsMatching, tree_aa_PB, tree_ref)
RobinsonFoulds(tree_aa_PB, tree_ref) # RF = 154

