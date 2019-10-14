biocLite("calibrate")
library(calibrate)
library(RColorBrewer)
library(wesanderson)
############################################################################################

r <- read.csv("IntALG1 vs MusALG1.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Int1_Mus1.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[1], c[2], c[5], c[6], c[7], c[8], 
                                                        c[9], c[10], c[11], c[12], c[14], c[15]))

############################################################################################

r <- read.csv("IntALG1 vs NeuALG1.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Int1_Neu1.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[2], c[3], c[4], c[5], c[8], c[9], 
                                                        c[10], c[12], c[13], c[14]))
############################################################################################

r <- read.csv("NeuALG1 vs MusALG1.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Neu1_Mus1.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[1], c[3], c[4], c[5], c[6], c[7], c[8], c[9], 
                                                        c[10], c[11], c[12], c[13], c[14], c[15]))


############################################################################################

r <- read.csv("IntALG2 vs MusALG2.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Int2_Mus2.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[1], c[5], c[7], c[8], c[9], c[10], 
                                                        c[11], c[12], c[14], c[15]))

############################################################################################

r <- read.csv("IntALG2 vs NeuALG2.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Int2_Neu2.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[3], c[4], c[5], c[6], c[8], c[9], 
                                                        c[12], c[13]))
############################################################################################

r <- read.csv("NeuALG2 vs MusALG2.csv")
colnames(r) <- c("ID", "sequence", "logFC", "logCPM", "F", "PValue", "FDR", "miRNA", "ID_2")
n <- read.delim("Neu2_Mus2.txt")

plot(r$logFC, jitter(-log10(r$FDR), 5),
     main="isomiRs diff. exp. between IntALG1 & MusALG1", xlim=c(-10,10), ylim=c(0,5),
     xlab="logFC", ylab="-log10 (adj. p-value)")

c <- rainbow(15)

with(subset(r, miRNA == "3p1"), points(logFC, -log10(FDR), pch=20, col=c[1]))
with(subset(r, miRNA == "3p243"), points(logFC, -log10(FDR), pch=20, col=c[2]))
with(subset(r, miRNA == "3p245"), points(logFC, -log10(FDR), pch=20, col=c[3]))
with(subset(r, miRNA == "3p249"), points(logFC, -log10(FDR), pch=20, col=c[4]))
with(subset(r, miRNA == "3p253"), points(logFC, -log10(FDR), pch=20, col=c[5]))
with(subset(r, miRNA == "3p54"), points(logFC, -log10(FDR), pch=20, col=c[6]))
with(subset(r, miRNA == "3p80"), points(logFC, -log10(FDR), pch=20, col=c[7]))
with(subset(r, miRNA == "3p90"), points(logFC, -log10(FDR), pch=20, col=c[8]))
with(subset(r, miRNA == "5p34"), points(logFC, -log10(FDR), pch=20, col=c[9]))
with(subset(r, miRNA == "5p48"), points(logFC, -log10(FDR), pch=20, col=c[10]))
with(subset(r, miRNA == "5p50"), points(logFC, -log10(FDR), pch=20, col=c[11]))
with(subset(r, miRNA == "5p71"), points(logFC, -log10(FDR), pch=20, col=c[12]))
with(subset(r, miRNA == "5p80"), points(logFC, -log10(FDR), pch=20, col=c[13]))
with(subset(r, miRNA == "5p86"), points(logFC, -log10(FDR), pch=20, col=c[14]))
with(subset(r, miRNA == "5p795"), points(logFC, -log10(FDR), pch=20, col=c[15]))
with(subset(r, FDR>0.05 ), points(logFC, -log10(FDR), pch=20, col="grey"))

abline(v=0, col="black", lty=1, lwd = 2)
abline(h=1.301, col="black", lty=2, lwd = 1)

legend("bottomright", legend=n$miRNA, cex=0.7, fill = c(c[1], c[3], c[4], c[5], c[6], c[7], c[8], c[9], 
                                                        c[10], c[11], c[12], c[13], c[14], c[15]))

