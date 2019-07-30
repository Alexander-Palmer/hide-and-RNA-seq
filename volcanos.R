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
#############################################################################################
##########################################Complex_form#######################################
#############################################################################################
#############################################################################################

biocLite("calibrate")
library(calibrate)

##############
#Data_input_1#
##############

comb <- read.csv("all_miRNAs_known.csv")
comb1 <- comb[, c(1:8, 12:20, 24:26)]
comb1 <- comb1[grep("cel-", comb1[,1]), ]

int_alg1 = comb1[, c(1:2, 3:5)]
colnames(int_alg1) <- c("miRNA", "seq", "logFC", "FDR", "pval")
int_alg2 = comb1[, c(1:2, 12:14)]
colnames(int_alg2) <- c("miRNA", "seq", "logFC", "FDR", "pval")
neu_alg1 = comb1[, c(1:2, 6:8)]
colnames(neu_alg1) <- c("miRNA", "seq", "logFC", "FDR", "pval")
neu_alg2 = comb1[, c(1:2, 15:17)]
colnames(neu_alg2) <- c("miRNA", "seq", "logFC", "FDR", "pval")
mus_alg1 = comb1[, c(1:2, 9:11)]
colnames(mus_alg1) <- c("miRNA", "seq", "logFC", "FDR", "pval")
mus_alg2 = comb1[, c(1:2, 18:20)]
colnames(mus_alg2) <- c("miRNA", "seq", "logFC", "FDR", "pval")

comb2 <- rbind(int_alg1, int_alg2, neu_alg1, neu_alg2, mus_alg1, mus_alg2)

########
#Plot_1#
########

plot(comb2$logFC, -log10(comb2$pval),
     main="Total miRNAs", xlim=c(0,11), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")
with(subset(comb2, pval<0.05 ), points(logFC, -log10(pval), pch=20, col="green"))

rm(list=ls())

##############
#Data_input_2#
##############

IntALG1 <- read.csv("Intestine ALG1.csv")
colnames(IntALG1) <- c("miRNA", "sequence", "logFC", "FDR", "pval")
IntALG2 <- read.csv("Intestine ALG2.csv")
colnames(IntALG2) <- c("miRNA", "sequence", "logFC", "FDR", "pval")
NeuALG1 <- read.csv("Neuron ALG1.csv")
colnames(NeuALG1) <- c("miRNA", "sequence", "logFC", "FDR", "pval")
NeuALG2 <- read.csv("Neuron ALG2.csv")
colnames(NeuALG2) <- c("miRNA", "sequence", "logFC", "FDR", "pval")
MusALG1 <- read.csv("Muscle ALG1.csv")
colnames(MusALG1) <- c("miRNA", "sequence", "logFC", "FDR", "pval")
MusALG2 <- read.csv("Muscle ALG2.csv")
colnames(MusALG2) <- c("miRNA", "sequence", "logFC", "FDR", "pval")

########
#Plot_2#
########
{
plot(IntALG1$logFC, -log10(IntALG1$pval),
     main="Intestine ALG1 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(IntALG1, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(IntALG1, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(IntALG1, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(IntALG1, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

plot(IntALG2$logFC, -log10(IntALG2$pval),
     main="Intestine ALG2 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(IntALG2, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(IntALG2, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(IntALG2, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(IntALG2, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

###
###
###

plot(NeuALG1$logFC, -log10(NeuALG1$pval),
     main="Neuron ALG1 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(NeuALG1, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(NeuALG1, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(NeuALG1, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(NeuALG1, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

plot(NeuALG2$logFC, -log10(NeuALG2$pval),
     main="Neuron ALG2 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(NeuALG2, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(NeuALG2, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(NeuALG2, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(NeuALG2, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

###
###
###

plot(MusALG1$logFC, -log10(MusALG1$pval),
     main="Muscle ALG1 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(MusALG1, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(MusALG1, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(MusALG1, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(MusALG1, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

plot(MusALG2$logFC, -log10(MusALG2$pval),
     main="Muscle ALG2 Volcano Plot", xlim=c(0,11.1), ylim=c(0,11),
     xlab="log(FC)", ylab="-log10 (p-value)")

with(subset(MusALG2, pval<9e-05 ), points(logFC, -log10(pval), pch=20, col="red"))
with(subset(MusALG2, abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="orange"))
with(subset(MusALG2, pval<9e-05 & abs(logFC)>5), points(logFC, -log10(pval), pch=20, col="green"))
with(subset(MusALG2, pval<9e-05 & abs(logFC)>5), textxy(logFC, -log10(pval), labs=miRNA, cex=0.55))

}
#############################################################################################
