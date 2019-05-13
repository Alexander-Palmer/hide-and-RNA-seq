#########################
#Installing Packages
#########################

source("http://bioconductor.org/biocLite.R")
biocLite("superheat", dependencies = TRUE)
biocLite("RColorBrewer", dependencies = TRUE)

library(superheat)
library(RColorBrewer)

#########################
#Input
#########################

n <- read.delim("Combined ALG1 abs.txt")

#########################
#Body
#########################

#Data insertion
o <- n
o$X <- NULL
colnames(o) <- c("Intestine", "Neuron", "Muscle")
rownames(o) <- (n[,1])
o <- as.matrix(o)

#Heatmap
superheat(o, left.label.size = 0.2, bottom.label.size = 0.05, scale = FALSE, 
          row.dendrogram = TRUE, title = "ALG1 miRNA logFC Tissue loading comparison", 
          row.title = "miRNA", column.title = "Tissue", left.label.text.size = 2, 
          grid.hline.col = "white", grid.vline.col = "white", grid.hline.size = 500, 
          grid.vline.size = 1, legend.height = 0.05, bottom.label.text.size = 3)
          
#########################
#Output
#########################

#Export manually - metafile is not working

#https://rlbarter.github.io/superheat/heatmap-colormap.html#heatmap-palette
