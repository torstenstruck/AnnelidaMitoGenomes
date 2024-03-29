#set working directory to the folder with the input file
#the following file is needed:
#CompiledProperties.csv

library("dplyr")
library("ggplot2")
library("tidyr")
library("ggpubr")
library("ggpmisc")
library("data.table") 

#Load data
CompiledSeqProperties <- read.csv("CompiledProperties.csv", header = TRUE, row.names=1)
View(CompiledSeqProperties)

#Generate box plots of genome structure such as size, duplication, introns, intergenic regions
#subset columns based on name structure
structure <- select(CompiledSeqProperties, matches(c("_Size_","_Num_")))
View(structure)

#subset the complete genomes not lacking data (NA)
GeneOrder_aligned_NA <- na.omit (read.delim("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/03_MitochondrialProperties/StructuralInformation/GeneOrder_aligned_with_tRNA.txt", row.names=1))
View(GeneOrder_aligned_NA)
structure_complete <- subset(structure, rownames(structure) %in% rownames(GeneOrder_aligned_NA))
View(structure_complete)

#generate boxplots for the different datasets of similar size
structure_1 <- structure_complete %>%  pivot_longer(., cols = c(1), names_to = "Var", values_to = "Val") %>%
  ggplot(aes(x = Var, y = Val)) +
  geom_boxplot() +
  ylim(12000,22000)
structure_2 <- structure_complete %>%  pivot_longer(., cols = c(2:4), names_to = "Var", values_to = "Val") %>%
  ggplot(aes(x = Var, y = Val)) +
  geom_boxplot()
structure_3 <- structure_complete %>%  pivot_longer(., cols = c(5), names_to = "Var", values_to = "Val") %>%
  ggplot(aes(x = Var, y = Val)) +
  geom_boxplot()
#combine the three graphics into one
figure_boxplot_structure <- ggarrange(structure_1, structure_2, structure_3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure_boxplot_structure

#Generate plot of genome size to intergenic size, provide numbers of duplications as color, of introns as shape, of intergenic regions as size, transparency of 50%,
ggplot(data = structure_complete, mapping = aes(x = Whole_Nuc_Size_bp, y = Whole_Nuc_Num_IG_bp, size = Whole_Nuc_Num_IG)) +
  geom_point(aes(alpha = 0.5, shape = as.character(Whole_Nuc_Num_Intron), color = as.character(Whole_Nuc_Num_Dup))) + 
  scale_size(range = c(1, 3)) +
  geom_smooth(method = "lm") +
# include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)

#explore distance measurements in different aspects
#subset columns based on name distance
distance <- select(CompiledSeqProperties, matches("_Dist_"))
View(distance)
#smaller subset columns based on name specific distance
distance_LB <- select(distance, matches("_LB"))
distance_PD <- select(distance, matches("_PD"))
distance_TR <- select(distance, matches("_TR"))
#group smaller subsets into one column
distance_LB_long <- distance_LB %>%  pivot_longer(., cols = c(1:ncol(distance_LB)), names_to = "Var_LB", values_to = "Val_LB")
distance_PD_long <- distance_PD %>%  pivot_longer(., cols = c(1:ncol(distance_PD)), names_to = "Var_PD", values_to = "Val_PD")
distance_TR_long <- distance_TR %>%  pivot_longer(., cols = c(1:ncol(distance_TR)), names_to = "Var_TR", values_to = "Val_TR")
#combine the long formats
distance_long <- na.omit(cbind(distance_LB_long,distance_PD_long,distance_TR_long))
#generate plots comparing the three measurements to each, add R2 for each and arrange them in one panel
distance_1 <- ggplot(data = distance_long, mapping = aes(x = Val_LB, y = Val_PD)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
distance_2 <- ggplot(data = distance_long, mapping = aes(x = Val_LB, y = Val_TR)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
distance_3 <- ggplot(data = distance_long, mapping = aes(x = Val_PD, y = Val_TR)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
figure_plot_distance <- ggarrange(distance_1, distance_2, distance_3,
                                      labels = c("A", "B", "C"),
                                      ncol = 3, nrow = 1)
figure_plot_distance

#generate boxplots for PD distance and the individual genes, group based on genes
distance_PD_Single <- select(distance_PD, matches(c("ATP","CO","NAD","^rRNA")))
distance_PD_Single_Nuc <- select(distance_PD_Single, matches(c("_Nuc_")))
distance_PD_Single_AA <- select(distance_PD_Single, matches(c("_AA_","^rRNA"))) #rRNA included once more to achieve similar spacing between plots and some color range
distance_PD_Single_Nuc_long <- distance_PD_Single_Nuc %>%  pivot_longer(., cols = c(1:ncol(distance_PD_Single_Nuc)), names_to = "Var_PD", values_to = "Val_PD")
distance_PD_Single_AA_long <- distance_PD_Single_AA %>%  pivot_longer(., cols = c(1:ncol(distance_PD_Single_AA)), names_to = "Var_PD", values_to = "Val_PD")
distance_PD_Single_Nuc_long <- na.omit(distance_PD_Single_Nuc_long)
distance_PD_Single_AA_long <- na.omit(distance_PD_Single_AA_long)

distance_4 <- ggplot(distance_PD_Single_Nuc_long, aes(x = Var_PD, y = Val_PD, fill = Var_PD)) +
  geom_boxplot()  +
  ylim(0,7) #28 outlier values above 7 are not displayed due to visiblity
distance_5 <- ggplot(distance_PD_Single_AA_long, aes(x = Var_PD, y = Val_PD, fill = Var_PD)) +
  geom_boxplot()  +
  ylim(0,7) #50 outlier values above 7 are not displayed due to visiblity
figure_boxplot_distance <- ggarrange(distance_4, distance_5,
                                  labels = c("A", "B"),
                                  ncol = 1, nrow = 2)
figure_boxplot_distance

#explore skew measurements in different aspects
#subset columns based on name skew
skew <- select(CompiledSeqProperties, matches("_Skew_"))
View(skew)
#compare AT/GC skews for single genes
skew_AT <- select(skew, matches("_AT"))
skew_AT_Single <- select(skew_AT, matches(c("ATP","CO","NAD","^rRNA")))
skew_GC <- select(skew, matches("_GC"))
skew_GC_Single <- select(skew_GC, matches(c("ATP","CO","NAD","^rRNA")))
#group smaller subsets into one column
skew_AT_Single_long <- skew_AT_Single %>%  pivot_longer(., cols = c(1:ncol(skew_AT_Single)), names_to = "Var_AT", values_to = "Val_AT")
skew_GC_Single_long <- skew_GC_Single %>%  pivot_longer(., cols = c(1:ncol(skew_GC_Single)), names_to = "Var_GC", values_to = "Val_GC")
skew_AT_Single_long <- na.omit(skew_AT_Single_long)
skew_GC_Single_long <- na.omit(skew_GC_Single_long)
#combine the long formats
skew_AT_GC_long <- na.omit(cbind(skew_AT_Single_long,skew_GC_Single_long))
#generate plots comparing AT and GC to each other, add color for genes
ggplot(data = skew_AT_GC_long, mapping = aes(x = Val_AT, y = Val_GC, color = Var_AT)) +
  geom_point(aes(alpha = 0.5))

#generate boxplots for AT/GC skews and the individual genes, group based on genes
skew_AT_GC_1 <- ggplot(skew_AT_Single_long, aes(x = Var_AT, y = Val_AT, fill = Var_AT)) +
  geom_boxplot()
skew_AT_GC_2 <- ggplot(skew_GC_Single_long, aes(x = Var_GC, y = Val_GC, fill = Var_GC)) +
  geom_boxplot()
figure_boxplot_skew_AT_GC <- ggarrange(skew_AT_GC_1, skew_AT_GC_2,
                                     labels = c("A", "B"),
                                     ncol = 1, nrow = 2)
figure_boxplot_skew_AT_GC

#compare CT/AG skews for single genes
skew_CT <- select(skew, matches("_CT"))
skew_CT_Single <- select(skew_CT, matches(c("ATP","CO","NAD","^rRNA")))
skew_AG <- select(skew, matches("_AG"))
skew_AG_Single <- select(skew_AG, matches(c("ATP","CO","NAD","^rRNA")))
#group smaller subsets into one column
skew_CT_Single_long <- skew_CT_Single %>%  pivot_longer(., cols = c(1:ncol(skew_CT_Single)), names_to = "Var_CT", values_to = "Val_CT")
skew_AG_Single_long <- skew_AG_Single %>%  pivot_longer(., cols = c(1:ncol(skew_AG_Single)), names_to = "Var_AG", values_to = "Val_AG")
skew_CT_Single_long <- na.omit(skew_CT_Single_long)
skew_AG_Single_long <- na.omit(skew_AG_Single_long)
#combine the long formats
skew_CT_AG_long <- na.omit(cbind(skew_CT_Single_long,skew_AG_Single_long))
#generate plots comparing CT and AG to each other, add color for genes
ggplot(data = skew_CT_AG_long, mapping = aes(x = Val_CT, y = Val_AG, color = Var_CT)) +
  geom_point(aes(alpha = 0.5))

#generate boxplots for CT/AG skews and the individual genes, group based on genes
skew_CT_AG_1 <- ggplot(skew_CT_Single_long, aes(x = Var_CT, y = Val_CT, fill = Var_CT)) +
  geom_boxplot()
skew_CT_AG_2 <- ggplot(skew_AG_Single_long, aes(x = Var_AG, y = Val_AG, fill = Var_AG)) +
  geom_boxplot()
figure_boxplot_skew_CT_AG <- ggarrange(skew_CT_AG_1, skew_CT_AG_2,
                                       labels = c("A", "B"),
                                       ncol = 1, nrow = 2)
figure_boxplot_skew_CT_AG

#generate plot comparing AT and GC for the whole dataset and all protein-coding genes
skew_AT_GC_3 <- ggplot(data = skew, mapping = aes(x = Whole_Nuc_Skew_AT, y = Whole_Nuc_Skew_GC)) +
  geom_point(aes(alpha = 0.5)) +
  xlim(-0.4,0.1) +
  ylim(-0.5,0.5) 
skew_AT_GC_4 <- ggplot(data = skew, mapping = aes(x = AllPCG_Nuc_Skew_AT, y = AllPCG_Nuc_Skew_GC)) +
  geom_point(aes(alpha = 0.5)) +
  xlim(-0.4,0.1) +
  ylim(-0.5,0.5)
figure_plot_skew_AT_GC <- ggarrange(skew_AT_GC_3, skew_AT_GC_4,
                                       labels = c("A", "B"),
                                       ncol = 2, nrow = 1)
figure_plot_skew_AT_GC

#generate plot comparing CT and AG for the whole dataset and all protein-coding genes
skew_CT_AG_3 <- ggplot(data = skew, mapping = aes(x = Whole_Nuc_Skew_CT, y = Whole_Nuc_Skew_AG)) +
  geom_point(aes(alpha = 0.5)) +
  xlim(-0.7,0.1)  +
  ylim(-0.3,0.6) 
skew_CT_AG_4 <- ggplot(data = skew, mapping = aes(x = AllPCG_Nuc_Skew_CT, y = AllPCG_Nuc_Skew_AG)) +
  geom_point(aes(alpha = 0.5)) +
  xlim(-0.7,0.1)  +
  ylim(-0.3,0.6)
figure_plot_skew_CT_AG <- ggarrange(skew_CT_AG_3, skew_CT_AG_4,
                                    labels = c("A", "B"),
                                    ncol = 2, nrow = 1)
figure_plot_skew_CT_AG


#explore ntRCFV measurements in different aspects
#subset columns based on name RCFV
RCFV <- select(CompiledSeqProperties, matches("_RCFV_"))
View(RCFV)
#subset columns for nucleotides and amino acids
RCFV_Nuc <- select(RCFV, matches("_Nuc_"))
RCFV_AA <- select(RCFV, matches("_AA_"))

#compare ntRCFV values for single genes as nucleotides
RCFV_Nuc_all <- select(RCFV_Nuc, matches("_all"))
RCFV_Nuc_all_Single <- select(RCFV_Nuc_all, matches(c("ATP","CO","NAD","^rRNA")))
#group smaller subsets into one column
RCFV_Nuc_all_Single_long <- RCFV_Nuc_all_Single %>%  pivot_longer(., cols = c(1:ncol(RCFV_Nuc_all_Single)), names_to = "Var_Nuc_all", values_to = "Val_Nuc_all")
#remove missing values
RCFV_Nuc_long <- na.omit(RCFV_Nuc_all_Single_long)

#generate boxplots for comparing ntRCFV
RCFV_Nuc_box <- ggplot(RCFV_Nuc_long, aes(x = Var_Nuc_all, y = Val_Nuc_all, fill = Var_Nuc_all)) +
  geom_boxplot() + theme(legend.position="none")

#compare ntRCFV values for single genes as amino acids
RCFV_AA_all <- select(RCFV_AA, matches("_all"))
RCFV_AA_all_Single <- select(RCFV_AA_all, matches(c("ATP","CO","NAD","^rRNA")))
#add rRNA to achieve similar spacing between box plots and some color range
RCFV_AA_all_Single_RNA <- cbind(RCFV_AA_all_Single, RCFV_AA_all_Single$NAD5_AA_RCFV_all, RCFV_AA_all_Single$NAD6_AA_RCFV_all) #add two more columns to make up for the additional RNAs in the nucleotides
#group smaller subsets into one column
RCFV_AA_all_Single_RNA_long <- RCFV_AA_all_Single_RNA %>%  pivot_longer(., cols = c(1:ncol(RCFV_AA_all_Single_RNA)), names_to = "Var_AA_all", values_to = "Val_AA_all")
#remove missing values
RCFV_AA_long_RNA <- na.omit(RCFV_AA_all_Single_RNA_long)

#generate boxplots for comparing ntRCFV 
RCFV_AA_box <- ggplot(RCFV_AA_long_RNA, aes(x = Var_AA_all, y = Val_AA_all, fill = Var_AA_all)) +
  geom_boxplot() + theme(legend.position="none")
figure_boxplot_RCFV <- ggarrange(RCFV_Nuc_box, RCFV_AA_box,
                                     labels = c("A", "B"),
                                     ncol = 1, nrow = 2)
figure_boxplot_RCFV

#generate plot comparing RCFV of all, and purines for the whole dataset and all protein-coding genes
RCFV_Nuc_plot <- ggplot(data = RCFV, mapping = aes(x = Whole_Nuc_RCFV_all, y = AllPCG_Nuc_RCFV_all)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
RCFV_AA_Nuc_plot <- ggplot(data = RCFV, mapping = aes(x = AllPCG_Nuc_RCFV_all, y = AllPCG_AA_RCFV_all)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
figure_plot_RCFV <- ggarrange(RCFV_Nuc_plot, RCFV_AA_Nuc_plot,
                                    labels = c("A", "B"),
                                    ncol = 2, nrow = 1)
figure_plot_RCFV

#explore selected frequencies for single genes
#subset columns based on name Freq
Freq <- select(CompiledSeqProperties, matches("_Freq_"))
View(Freq)
#subset columns for nucleotides and amino acids
Freq_Nuc <- select(Freq, matches("_Nuc_"))
Freq_AA <- select(Freq, matches("_AA_"))
Freq_AT <- select(Freq_Nuc, matches("_AT"))
Freq_AT_Single <- select(Freq_AT, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Purine <- select(Freq_Nuc, matches("_Purine"))
Freq_Purine_Single <- select(Freq_Purine, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Gap <- select(Freq_Nuc, matches("_Gap"))
Freq_Gap_Single <- select(Freq_Gap, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Hydrophobic <- select(Freq_AA, matches("_Hydrophobic"))
Freq_Hydrophobic_Single <- select(Freq_Hydrophobic, matches(c("ATP","CO","NAD")))
Freq_Hydrophobic_Single_RNA <- cbind(Freq_Hydrophobic_Single, Freq_AT_Single$rRNAL_Nuc_Freq_AT, Freq_AT_Single$rRNAS_Nuc_Freq_AT) #rRNA included once more to achieve similar spacing between plots and some color range
#group smaller subsets into one column
Freq_AT_Single_long <- Freq_AT_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_AT_Single)), names_to = "Var_AT", values_to = "Val_AT")
Freq_Purine_Single_long <- Freq_Purine_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Purine_Single)), names_to = "Var_Purine", values_to = "Val_Purine")
Freq_Gap_Single_long <- Freq_Gap_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Gap_Single)), names_to = "Var_Gap", values_to = "Val_Gap")
Freq_Hydrophobic_Single_long <- Freq_Hydrophobic_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Hydrophobic_Single)), names_to = "Var_Hydrophobic", values_to = "Val_Hydrophobic")
Freq_Hydrophobic_Single_long_RNA <- Freq_Hydrophobic_Single_RNA %>%  pivot_longer(., cols = c(1:ncol(Freq_Hydrophobic_Single_RNA)), names_to = "Var_Hydrophobic", values_to = "Val_Hydrophobic")
#remove missing values
Freq_AT_Single_long <- na.omit(Freq_AT_Single_long)
Freq_Purine_Single_long <- na.omit(Freq_Purine_Single_long)
Freq_Gap_Single_long <- na.omit(Freq_Gap_Single_long)
Freq_Hydrophobic_Single_long <- na.omit(Freq_Hydrophobic_Single_long)
Freq_Hydrophobic_Single_long_RNA <- na.omit(Freq_Hydrophobic_Single_long_RNA)

#generate boxplots for subsets, group based on genes
freq_1 <- ggplot(Freq_Gap_Single_long, aes(x = Var_Gap, y = Val_Gap, fill = Var_Gap)) +
  geom_boxplot() + theme(legend.position="none") +
  ylim(0.075,0.6) #8 outlier values above 0.6 are not displayed due to visibility
freq_2 <- ggplot(Freq_AT_Single_long, aes(x = Var_AT, y = Val_AT, fill = Var_AT)) +
  geom_boxplot() + theme(legend.position="none")
freq_3 <- ggplot(Freq_Hydrophobic_Single_long_RNA, aes(x = Var_Hydrophobic, y = Val_Hydrophobic, fill = Var_Hydrophobic)) +
  geom_boxplot() + theme(legend.position="none")
figure_boxplot_freq <- ggarrange(freq_1, freq_2, freq_3,
                                     labels = c("A", "B", "c"),
                                     ncol = 1, nrow = 3)
figure_boxplot_freq

#combine the long formats
Freq_Nuc_long <- na.omit(cbind(Freq_Purine_Single_long,Freq_AT_Single_long))

#generate plots comparing purines and AT to each other
Freq_Nuc_plot <- ggplot(data = Freq_Nuc_long, mapping = aes(x = Val_Purine, y = Val_AT)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
Freq_Nuc_plot

#compare frequencies for positive, neutral and negative for single genes
Freq_Pos <- select(Freq_AA, matches("_Pos"))
Freq_Pos_Single <- select(Freq_Pos, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Neu <- select(Freq_AA, matches("_Neu"))
Freq_Neu_Single <- select(Freq_Neu, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Neg <- select(Freq_AA, matches("_Neg"))
Freq_Neg_Single <- select(Freq_Neg, matches(c("ATP","CO","NAD","^rRNA")))
#group smaller subsets into one column
Freq_Pos_Single_long <- Freq_Pos_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Pos_Single)), names_to = "Var_Pos", values_to = "Val_Pos")
Freq_Neu_Single_long <- Freq_Neu_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Neu_Single)), names_to = "Var_Neu", values_to = "Val_Neu")
Freq_Neg_Single_long <- Freq_Neg_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Neg_Single)), names_to = "Var_Neg", values_to = "Val_Neg")
#combine the long formats
Freq_Charge_long <- na.omit(cbind(Freq_Pos_Single_long,Freq_Neu_Single_long,Freq_Neg_Single_long))

#generate plots comparing positive, neutral and negative to each other
Freq_Charge_1 <- ggplot(data = Freq_Charge_long, mapping = aes(x = Val_Pos, y = Val_Neu)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
Freq_Charge_2 <- ggplot(data = Freq_Charge_long, mapping = aes(x = Val_Pos, y = Val_Neg)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
Freq_Charge_3 <- ggplot(data = Freq_Charge_long, mapping = aes(x = Val_Neu, y = Val_Neg)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
figure_plot_Freq_charge <- ggarrange(Freq_Charge_1, Freq_Charge_2, Freq_Charge_3,
                                  labels = c("A", "B", "C"),
                                  ncol = 3, nrow = 1)
figure_plot_Freq_charge

#compare RCFV skews for hydrophobic, polar and neutral for single genes
Freq_Polar <- select(Freq_AA, matches("_Polar"))
Freq_Polar_Single <- select(Freq_Polar, matches(c("ATP","CO","NAD","^rRNA")))
Freq_Hydrophobic <- select(Freq_AA, matches("_Hydrophobic"))
Freq_Hydrophobic_Single <- select(Freq_Hydrophobic, matches(c("ATP","CO","NAD")))
#group smaller subsets into one column
Freq_Polar_Single_long <- Freq_Polar_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Polar_Single)), names_to = "Var_Polar", values_to = "Val_Polar")
Freq_Hydrophobic_Single_long <- Freq_Hydrophobic_Single %>%  pivot_longer(., cols = c(1:ncol(Freq_Hydrophobic_Single)), names_to = "Var_Hydrophobic", values_to = "Val_Hydrophobic")
#combine the long formats
Freq_AA_long <- na.omit(cbind(Freq_Hydrophobic_Single_long,Freq_Polar_Single_long,Freq_Neu_Single_long))

#generate plots comparing all, hydrophobic, polar and neutral to each other
Freq_AA_1 <- ggplot(data = Freq_AA_long, mapping = aes(x = Val_Hydrophobic, y = Val_Polar)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
Freq_AA_2 <- ggplot(data = Freq_AA_long, mapping = aes(x = Val_Hydrophobic, y = Val_Neu)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
Freq_AA_3 <- ggplot(data = Freq_AA_long, mapping = aes(x = Val_Polar, y = Val_Neu)) +
  geom_point(aes(alpha = 0.5)) + 
  geom_smooth(method = "lm") +
  # include R2
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = y ~ x, parse = TRUE, size = 5)
figure_plot_AA <- ggarrange(Freq_AA_1, Freq_AA_2, Freq_AA_3,
                                     labels = c("A", "B", "C"),
                                     ncol = 3, nrow = 1)
figure_plot_AA


  