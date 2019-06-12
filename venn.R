library(VennDiagram)

z1 <- read.csv('IntALG1 vs N2 Top - Pos.csv')
z1 <- z1[,2]
z1 <- as.vector(z1)
z2 <- read.csv('IntALG2 vs N2 Top - Pos.csv')
z2 <- z2[,2]
z2 <- as.vector(z2)
z3 <- read.csv('NeuALG1 vs N2 Top - Pos.csv')
z3 <- z3[,2]
z3 <- as.vector(z3)
z4 <- read.csv('NeuALG2 vs N2 Top - Pos.csv')
z4 <- z4[,2]
z4 <- as.vector(z4)
z5 <- read.csv('MusALG1 vs N2 Top - Pos.csv')
z5 <- z5[,2]
z5 <- as.vector(z5)
z6 <- read.csv('MusALG2 vs N2 Top - Pos.csv')
z6 <- z6[,2]
z6 <- as.vector(z6)

zz1 <- union(z1, z2)
zz2 <- union(z3, z4)
zz3 <- union(z5, z6)

venn.plot <- venn.diagram(
  list("Intestine" = zz1, "Neuron" = zz2, "Muscle" = zz3), 
  fill = c("cornflowerblue", "green", "red"), "Venn.tiff", scaled=TRUE, 
  main = "Differential isomiR loading between tissues", cex = 1.3, main.cex = 1.3) 
       
venn.plot <- venn.diagram(
  list("ALG1" = z5, "ALG2" = z6), 
  fill = c("yellow", "violet"), "Venn2.tiff", scaled=TRUE, 
  main = "Muscle isomiRs / N2", cex = 1.3, main.cex = 1.3)
