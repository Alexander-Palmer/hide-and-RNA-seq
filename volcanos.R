biocLite("calibrate")
library(calibrate)

r <- read.delim("Neuron ALG2.txt")
colnames(r) <- c("miRNA", "logFC", "pval")



plot(r$logFC, -log10(r$pval),
     main="Neuron ALG2 Volcano Plot", xlim=c(0,11), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(r, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(r, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(r, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(r, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

#############################################################################################

s <- read.delim("Muscle Combined.txt")
colnames(s) <- c("miRNA", "logFC.ALG1", "pval.ALG1", "logFC.ALG2", "pval.ALG2")

plot(s$logFC.ALG1, -log10(s$pval.ALG1),
     main="Muscle Combined ALG Volcano Plot", xlim=c(0,11), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)") #, col = "red")
points(s$logFC.ALG2, -log10(s$pval.ALG2)) #, col = "blue")

with(subset(s), points(logFC.ALG1, -log10(pval.ALG1), pch=20, col="red"))
with(subset(s), points(logFC.ALG2, -log10(pval.ALG2), pch=20, col="green"))
