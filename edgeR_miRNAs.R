###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")
library(edgeR)

a <- read.delim('CPM Table and Raw Table wo ident 2 integer.txt')
a <- read.delim('ALL VARIANTS with at least 1cpmillion 2.txt')

a <- as.data.frame(a)
a[is.na(a)] <- 0

levels <- letters[1:18]
levels[c(1,2)] <- "A"
levels[c(3,4)] <- "B"
levels[c(5,6)] <- "C"
levels[c(7,8)] <- "D"
levels[c(9,10)] <- "E"
levels[c(11,12)] <- "F"
levels[c(13,14)] <- "G"
levels[c(15,16)] <- "H"
levels[c(17,18)] <- "I"


#levels <- letters[1:24]
#levels[c(1,2,3)] <- "A"
#levels[c(4,5,6)] <- "B"
#levels[c(7,8,9)] <- "C"
#levels[c(10,11,12)] <- "D"
#levels[c(13,14,15)] <- "E"
#levels[c(16,17,18)] <- "F"
#levels[c(19,20,21)] <- "G"
#levels[c(22,23,24)] <- "H"

b <- read.delim('CPM Table and Raw Table.txt')
b <- read.delim('ALL VARIANTS with at least 1cpmillion.txt')
sequences <- b$seq

a <- DGEList(counts=a, group=levels, genes=sequences)

#Filter out genes with less than 1CPM in less than 2 groups
a$samples
keep <- rowSums(cpm(a)>1) >= 2
a <- a[keep, , keep.lib.sizes=FALSE]
a$samples

#Normalisation not required, due to no variation in total read abundances

#Adjust sample for RNA composition bias (CPM affected by higher expression in dif. samples)
a <- calcNormFactors(a)

#Estimate dispersion between groups
a <- estimateDisp(a)

###Classic approach
et <- exactTest(a, pair=c("A","B"))
topTags(et)
write.table(et, file="AZ vs EB.csv", sep = ",", col.names = NA, qmethod = "double")

###GLM between groups
design <- model.matrix(~0+levels, data=a$samples)
colnames(design) <- levels(a$samples$levels)
#a <- estimateGLMCommonDisp(a,design)
#a <- estimateGLMTrendedDisp(a,design)
#a <- estimateGLMTagwiseDisp(a,design)
fit <- glmQLFit(a, design)

capture.output(fit$genes, append = TRUE, file = "fitgenes2.txt")
capture.output(fit$fitted.values, append = TRUE, file = "fitfittedvalues2.txt")

#Mus_N2
Mus_Mus.inac <- glmQLFTest(fit, contrast=c(0,0,0,0,0,1,0,-1))
a1 <- topTags(Mus_Mus.inac, n=Inf)
b1 <- topTags(Mus_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a1, file="Mus vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b1, file="Mus vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#Neu_N2
Neu_Mus.inac <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,1,0,0))
a2 <- topTags(Neu_Mus.inac, n=Inf)
b2 <- topTags(Neu_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a2, file="Neu vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b2, file="Neu vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#Int_N2
Int_Mus.inac <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,1,0,0))
a3 <- topTags(Int_Mus.inac, n=Inf)
b3 <- topTags(Int_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a3, file="Int vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b3, file="Int vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#N2_N2
N2_Mus.inac <- glmQLFTest(fit, contrast=c(0,0,0,0,0,-1,1,0))
a4 <- topTags(N2_Mus.inac, n=Inf)
b4 <- topTags(N2_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a4, file="N2 vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b4, file="N2 vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#Isp.1_N2
Isp.1_Mus.inac <- glmQLFTest(fit, contrast=c(0,0,-1,0,0,1,0,0))
a5 <- topTags(Isp.1_Mus.inac, n=Inf)
b5 <- topTags(Isp.1_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a5, file="Isp.1 vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b5, file="Isp.1 vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#AZ_N2
AZ_Mus.inac <- glmQLFTest(fit, contrast=c(-1,0,0,0,0,1,0,0))
a6 <- topTags(AZ_Mus.inac, n=Inf)
b6 <- topTags(AZ_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a6, file="AZ vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b6, file="AZ vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

#EB_N2
EB_Mus.inac <- glmQLFTest(fit, contrast=c(0,-1,0,0,0,1,0,0))
a7 <- topTags(EB_Mus.inac, n=Inf)
b7 <- topTags(EB_Mus.inac, sort.by = "logFC", n=Inf, p.value=0.05)
write.table(a7, file="EB vs N2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(b7, file="EB vs N2 Top.csv", sep = ",", col.names = NA, qmethod = "double")
