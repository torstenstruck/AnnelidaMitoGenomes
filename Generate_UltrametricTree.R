# set working directory to the appropriate folder
# the following files are needed Masked_18S.treefile and Unmasked_18S.treefile

library(ape)

#load rooted tree in newick format
mytree <- read.tree("Masked_18S.treefile")

#generate time-calibrated ultrametric tree using the correlated method and lambda of 1
mytimetree <- chronos(mytree, lambda = 1, model = "correlated")
plot(mytimetree)
write.tree(mytimetree, file = "Masked_18S_ultrametric.treefile", append = FALSE,
           digits = 10, tree.names = FALSE)

#load rooted tree in newick format
mytree <- read.tree("Unmasked_18S.treefile")

#generate time-calibrated ultrametric tree using the correlated method and lambda of 1
mytimetree <- chronos(mytree, lambda = 1, model = "correlated")
plot(mytimetree)
write.tree(mytimetree, file = "Unmasked_18S_ultrametric.treefile", append = FALSE,
           digits = 10, tree.names = FALSE)
