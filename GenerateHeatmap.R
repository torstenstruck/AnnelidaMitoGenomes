#set working directory to the folder with input file

#import data
PropertiesCorrelatedFactors <- read.delim("PropertiesCorrelatedFactors.txt", row.names=1)
View(PropertiesCorrelatedFactors)

library(RColorBrewer)
library(gplots)

#creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "orange", "yellow", "green", "cyan", "blue"))(n = 299) 
#defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,0.01,length=100),  # for white
#               seq(0.02,0.3,length=100),  # for orange
#               seq(0.31,0.5,length=100), # for green
#               seq(0.51,0.8,length=100), # for cyan
#               seq(0.81,1.0,length=100)) # for blue

#generate heatmap using the color scheme and add a key
heatmap.2(as.matrix(PropertiesCorrelatedFactors),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=c("white", "chocolate4", "chocolate", "sandybrown", 
                "orange", "orangered", "red", "coral1", "gold", "yellow", 
                "yellowgreen", "green", "cyan", "blue", "darkslateblue"),       # use on color palette defined earlier
          breaks=c(0,0.001,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),    # enable color transition at specified limits
          dendrogram="none",     # no dendrograms
          Colv="NA",            # turn off column clustering
          Rowv="NA")            # turn off row clustering
