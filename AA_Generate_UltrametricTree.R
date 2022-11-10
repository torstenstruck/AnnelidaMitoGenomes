#set working directory to the folder with the inout files
#the following files are needed:
#AA_supermatrix_partition_rooted.txt.treefile

library(ape)

#load rooted tree in newick format
mytree <- read.tree("AA_supermatrix_partition_rooted.txt.treefile")

#generate time-calibrated ultrametric tree using the correlated method, lambda of 1 and relative dates from 1 to 0
mytimetree <- chronos(mytree, lambda = 1, model = "correlated", control = chronos.control())
plot(mytimetree)
write.tree(mytimetree, file = "AA_supermatrix_partition_Ultrametric.txt.treefile", append = FALSE,
           digits = 10, tree.names = FALSE)
