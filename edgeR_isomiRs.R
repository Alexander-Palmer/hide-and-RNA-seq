BiocManager::install("edgeR")
library(edgeR)

#######
#Input#
#######

a <- read.delim('All_variants_above_1CPM_clean.txt')
b <- read.delim('All_variants_above_1CPM_seq.txt')

a <- as.data.frame(a)
a[is.na(a)] <- 0
sequences <- b$seq
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

##############
#Processing_1#
##############

a <- DGEList(counts=a, group=levels, genes=sequences)

keep <- rowSums(cpm(a)>50) >= 2
a <- a[keep, , keep.lib.sizes=FALSE]

a <- calcNormFactors(a)
a <- estimateDisp(a)

design <- model.matrix(~0+levels, data=a$samples)
colnames(design) <- levels(a$samples$levels)
fit <- glmQLFit(a, design)

########
#Body_1#
########

#IntALG1_IntALG2
Plcehldr <- glmQLFTest(fit, contrast=c(1,-1,0,0,0,0,0,0,0))
a1 <- topTags(Plcehldr, n=Inf)
a1 <- as.data.frame(a1)
b1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
b1 <- as.data.frame(b1)

#NeuALG1_NeuALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,-1,0,0,0,0,0))
c1 <- topTags(Plcehldr, n=Inf)
c1 <- as.data.frame(c1)
d1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
d1 <- as.data.frame(d1)

#MusALG1_MusALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,0,0,0,1,-1,0))
e1 <- topTags(Plcehldr, n=Inf)
e1 <- as.data.frame(e1)
f1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
f1 <- as.data.frame(f1)

########
#Body_2#
########

#IntALG1_NeuALG1
Plcehldr <- glmQLFTest(fit, contrast=c(1,0,-1,0,0,0,0,0,0))
g1 <- topTags(Plcehldr, n=Inf)
g1 <- as.data.frame(g1)
h1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
h1 <- as.data.frame(h1)

#IntALG1_MusALG1
Plcehldr <- glmQLFTest(fit, contrast=c(1,0,0,0,0,0,-1,0,0))
i1 <- topTags(Plcehldr, n=Inf)
i1 <- as.data.frame(i1)
j1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
j1 <- as.data.frame(j1)

#NeuALG1_MusALG1
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,1,0,0,0,-1,0,0))
k1 <- topTags(Plcehldr, n=Inf)
k1 <- as.data.frame(k1)
l1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
l1 <- as.data.frame(l1)


########
#Body_3#
########

#IntALG2_NeuALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,-1,0,0,0,0,0))
m1 <- topTags(Plcehldr, n=Inf)
m1 <- as.data.frame(m1)
n1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
n1 <- as.data.frame(n1)


#IntALG2_MusALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,1,0,0,0,0,0,-1,0))
o1 <- topTags(Plcehldr, n=Inf)
o1 <- as.data.frame(o1)
p1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
p1 <- as.data.frame(p1)

#NeuALG2_MusALG2
Plcehldr <- glmQLFTest(fit, contrast=c(0,0,0,1,0,0,0,-1,0))
q1 <- topTags(Plcehldr, n=Inf)
q1 <- as.data.frame(q1)
r1 <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
r1 <- as.data.frame(r1)

##############
#Processing_2#
##############
aa1 <- read.csv("IntALG1 vs N2 Top - Delim.csv", header=F)$V2
aa1 <- as.character(aa1)
bb1 <- read.csv("IntALG2 vs N2 Top - Delim.csv", header=F)$V2
bb1 <- as.character(bb1)
cc1 <- read.csv("NeuALG1 vs N2 Top - Delim.csv", header=F)$V2
cc1 <- as.character(cc1)
dd1 <- read.csv("NeuALG2 vs N2 Top - Delim.csv", header=F)$V2
dd1 <- as.character(dd1)
ee1 <- read.csv("MusALG1 vs N2 Top - Delim.csv", header=F)$V2
ee1 <- as.character(ee1)
ff1 <- read.csv("MusALG2 vs N2 Top - Delim.csv", header=F)$V2
ff1 <- as.character(ff1)

aaa1 <- union(aa1, bb1)
bbb1 <- union(cc1, dd1)
ccc1 <- union(ee1, ff1)
ddd1 <- union(aa1, cc1)
eee1 <- union(aa1, ee1)
fff1 <- union(cc1, ee1)
ggg1 <- union(bb1, dd1)
hhh1 <- union(bb1, ff1)
iii1 <- union(dd1, ff1)

aInt1_Int2 <- a1[grep(paste(aaa1, collapse = "|"), a1[,1]), ]
bInt1_Int2 <- b1[grep(paste(aaa1, collapse = "|"), b1[,1]), ]

cNeu1_Neu2 <- c1[grep(paste(bbb1, collapse = "|"), c1[,1]), ]
dNeu1_Neu2 <- d1[grep(paste(bbb1, collapse = "|"), d1[,1]), ]

eMus1_Mus2 <- e1[grep(paste(ccc1, collapse = "|"), e1[,1]), ]
fMus1_Mus2 <- f1[grep(paste(ccc1, collapse = "|"), f1[,1]), ]

gInt1_Neu1 <- g1[grep(paste(ddd1, collapse = "|"), g1[,1]), ]
hInt1_Neu1 <- h1[grep(paste(ddd1, collapse = "|"), h1[,1]), ]

iInt1_Mus1 <- i1[grep(paste(eee1, collapse = "|"), i1[,1]), ]
jInt1_Mus1 <- j1[grep(paste(eee1, collapse = "|"), j1[,1]), ]

kNeu1_Mus1 <- k1[grep(paste(fff1, collapse = "|"), k1[,1]), ]
lNeu1_Mus1 <- l1[grep(paste(fff1, collapse = "|"), l1[,1]), ]

mInt2_Neu2 <- m1[grep(paste(ggg1, collapse = "|"), m1[,1]), ]
nInt2_Neu2 <- n1[grep(paste(ggg1, collapse = "|"), n1[,1]), ]

oInt2_Mus2 <- o1[grep(paste(hhh1, collapse = "|"), o1[,1]), ]
pInt2_Mus2 <- p1[grep(paste(hhh1, collapse = "|"), p1[,1]), ]

qNeu2_Mus2 <- q1[grep(paste(iii1, collapse = "|"), q1[,1]), ]
rNeu2_Mus2 <- r1[grep(paste(iii1, collapse = "|"), r1[,1]), ]


########
#Output#
########

write.table(aInt1_Int2, file="IntALG1 vs IntALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(bInt1_Int2, file="IntALG1 vs IntALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(cNeu1_Neu2, file="NeuALG1 vs NeuALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(dNeu1_Neu2, file="NeuALG1 vs NeuALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(eMus1_Mus2, file="MusALG1 vs MusALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(fMus1_Mus2, file="MusALG1 vs MusALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

write.table(gInt1_Neu1, file="IntALG1 vs NeuALG1.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(hInt1_Neu1, file="IntALG1 vs NeuALG1 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(iInt1_Mus1, file="IntALG1 vs MusALG1.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(jInt1_Mus1, file="IntALG1 vs MusALG1 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(kNeu1_Mus1, file="NeuALG1 vs MusALG1.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(lNeu1_Mus1, file="NeuALG1 vs MusALG1 Top.csv", sep = ",", col.names = NA, qmethod = "double")

write.table(mInt2_Neu2, file="IntALG2 vs NeuALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(nInt2_Neu2, file="IntALG2 vs NeuALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(oInt2_Mus2, file="IntALG2 vs MusALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(pInt2_Mus2, file="IntALG2 vs MusALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(qNeu2_Mus2, file="NeuALG2 vs MusALG2.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(rNeu2_Mus2, file="NeuALG2 vs MusALG2 Top.csv", sep = ",", col.names = NA, qmethod = "double")

###
