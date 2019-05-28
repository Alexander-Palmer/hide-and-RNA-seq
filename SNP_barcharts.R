a <- read.delim('total_meta2.txt')
#a <- a[!a$SNP == 0, ]
#a$SNP[a$SNP == TRUE] <- T
#a$SNP[a$SNP == '1'] <- T
a <- a[!a$SNP == FALSE, ]

a$Count <- NA
a$A284_1 <- as.numeric(a$A284_1)
a$A285_1 <- as.numeric(a$A285_1)
a$A286_1 <- as.numeric(a$A286_1)
a$A287_1 <- as.numeric(a$A287_1)
a$A290_1 <- as.numeric(a$A290_1)
a$A291_1 <- as.numeric(a$A291_1)
a$Count <- rowSums(a[,c(25,26,27,28,29,30)])
a$CountALG1 <- rowSums(a[,c(25,27,29)])
a$CountALG2 <- rowSums(a[,c(26,28,30)])
a$CountInt <- rowSums(a[,c(25,26)])
a$CountNeu <- rowSums(a[,c(27,28)])
a$CountMus <- rowSums(a[,c(29,30)])

#####
IntALG1_A <- read.csv("IntALG1 vs N2 Top - Pos.csv", header=F)$V2
IntALG1_A <- as.character(IntALG1_A)
IntALG1_A <- intersect(IntALG1_A, a[,2])
IntALG1_T <- read.csv("IntALG1 vs N2 Top - Pos.csv", header=F)$V2
IntALG1_T <- as.character(IntALG1_T)
IntALG1_T <- intersect(IntALG1_T, a[,2])
IntALG1_G <- read.csv("IntALG1 vs N2 Top - Pos.csv", header=F)$V2
IntALG1_G <- as.character(IntALG1_G)
IntALG1_G <- intersect(IntALG1_G, a[,2])
IntALG1_C <- read.csv("IntALG1 vs N2 Top - Pos.csv", header=F)$V2
IntALG1_C <- as.character(IntALG1_C)
IntALG1_C <- intersect(IntALG1_C, a[,2])

NeuALG1_A <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header=F)$V2
NeuALG1_A <- as.character(NeuALG1_A)
NeuALG1_A <- intersect(NeuALG1_A, a[,2])
NeuALG1_T <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header=F)$V2
NeuALG1_T <- as.character(NeuALG1_T)
NeuALG1_T <- intersect(NeuALG1_T, a[,2])
NeuALG1_G <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header=F)$V2
NeuALG1_G <- as.character(NeuALG1_G)
NeuALG1_G <- intersect(NeuALG1_G, a[,2])
NeuALG1_C <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header=F)$V2
NeuALG1_C <- as.character(NeuALG1_C)
NeuALG1_C <- intersect(NeuALG1_C, a[,2])

MusALG1_A <- read.csv("MusALG1 vs N2 Top - Pos.csv", header=F)$V2
MusALG1_A <- as.character(MusALG1_A)
MusALG1_A <- intersect(MusALG1_A, a[,2])
MusALG1_T <- read.csv("MusALG1 vs N2 Top - Pos.csv", header=F)$V2
MusALG1_T <- as.character(MusALG1_T)
MusALG1_T <- intersect(MusALG1_T, a[,2])
MusALG1_G <- read.csv("MusALG1 vs N2 Top - Pos.csv", header=F)$V2
MusALG1_G <- as.character(MusALG1_G)
MusALG1_G <- intersect(MusALG1_G, a[,2])
MusALG1_C <- read.csv("MusALG1 vs N2 Top - Pos.csv", header=F)$V2
MusALG1_C <- as.character(MusALG1_C)
MusALG1_C <- intersect(MusALG1_C, a[,2])


IntALG2_A <- read.csv("IntALG2 vs N2 Top - Pos.csv", header=F)$V2
IntALG2_A <- as.character(IntALG2_A)
IntALG2_A <- intersect(IntALG2_A, a[,2])
IntALG2_T <- read.csv("IntALG2 vs N2 Top - Pos.csv", header=F)$V2
IntALG2_T <- as.character(IntALG2_T)
IntALG2_T <- intersect(IntALG2_T, a[,2])
IntALG2_G <- read.csv("IntALG2 vs N2 Top - Pos.csv", header=F)$V2
IntALG2_G <- as.character(IntALG2_G)
IntALG2_G <- intersect(IntALG2_G, a[,2])
IntALG2_C <- read.csv("IntALG2 vs N2 Top - Pos.csv", header=F)$V2
IntALG2_C <- as.character(IntALG2_C)
IntALG2_C <- intersect(IntALG2_C, a[,2])

NeuALG2_A <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header=F)$V2
NeuALG2_A <- as.character(NeuALG2_A)
NeuALG2_A <- intersect(NeuALG2_A, a[,2])
NeuALG2_T <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header=F)$V2
NeuALG2_T <- as.character(NeuALG2_T)
NeuALG2_T <- intersect(NeuALG2_T, a[,2])
NeuALG2_G <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header=F)$V2
NeuALG2_G <- as.character(NeuALG2_G)
NeuALG2_G <- intersect(NeuALG2_G, a[,2])
NeuALG2_C <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header=F)$V2
NeuALG2_C <- as.character(NeuALG2_C)
NeuALG2_C <- intersect(NeuALG2_C, a[,2])

MusALG2_A <- read.csv("MusALG2 vs N2 Top - Pos.csv", header=F)$V2
MusALG2_A <- as.character(MusALG2_A)
MusALG2_A <- intersect(MusALG2_A, a[,2])
MusALG2_T <- read.csv("MusALG2 vs N2 Top - Pos.csv", header=F)$V2
MusALG2_T <- as.character(MusALG2_T)
MusALG2_T <- intersect(MusALG2_T, a[,2])
MusALG2_G <- read.csv("MusALG2 vs N2 Top - Pos.csv", header=F)$V2
MusALG2_G <- as.character(MusALG2_G)
MusALG2_G <- intersect(MusALG2_G, a[,2])
MusALG2_C <- read.csv("MusALG2 vs N2 Top - Pos.csv", header=F)$V2
MusALG2_C <- as.character(MusALG2_C)
MusALG2_C <- intersect(MusALG2_C, a[,2])


IntALG1 <- Reduce(union, list(IntALG1_A, IntALG1_T, IntALG1_G, IntALG1_C))
IntALG1 <- intersect(IntALG1, a[,2])
write.table(IntALG1, file="IntALG1 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
IntALG2 <- Reduce(union, list(IntALG2_A, IntALG2_T, IntALG2_G, IntALG2_C))
IntALG2 <- intersect(IntALG2, a[,2])
write.table(IntALG2, file="IntALG2 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
NeuALG1 <- Reduce(union, list(NeuALG1_A, NeuALG1_T, NeuALG1_G, NeuALG1_C))
NeuALG1 <- intersect(NeuALG1, a[,2])
write.table(NeuALG1, file="NeuALG1 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
NeuALG2 <- Reduce(union, list(NeuALG2_A, NeuALG2_T, NeuALG2_G, NeuALG2_C))
NeuALG2 <- intersect(NeuALG2, a[,2])
write.table(NeuALG2, file="NeuALG2 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
MusALG1 <- Reduce(union, list(MusALG1_A, MusALG1_T, MusALG1_G, MusALG1_C))
MusALG1 <- intersect(MusALG1, a[,2])
write.table(MusALG1, file="MusALG1 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
MusALG2 <- Reduce(union, list(MusALG2_A, MusALG2_T, MusALG2_G, MusALG2_C))
MusALG2 <- intersect(MusALG2, a[,2])
write.table(MusALG2, file="MusALG2 vs N2 Top - Delim.csv", sep = ",", col.names = NA, qmethod = "double")
#Tissue, ALG and nucleotide as above

ALG1_A <- Reduce(union, list(IntALG2_A, MusALG1_A, NeuALG1_A))
ALG1_A <- intersect(ALG1_A, a[,2])
ALG1_T <- Reduce(union, list(IntALG2_T, MusALG1_T, NeuALG1_T))
ALG1_T <- intersect(ALG1_T, a[,2])
ALG1_G <- Reduce(union, list(IntALG2_G, MusALG1_G, NeuALG1_G))
ALG1_G <- intersect(ALG1_G, a[,2])
ALG1_C <- Reduce(union, list(IntALG2_C, MusALG1_C, NeuALG1_C))
ALG1_C <- intersect(ALG1_C, a[,2])
ALG2_A <- Reduce(union, list(IntALG2_A, MusALG2_A, NeuALG2_A))
ALG2_A <- intersect(ALG2_A, a[,2])
ALG2_T <- Reduce(union, list(IntALG2_T, MusALG2_T, NeuALG2_T))
ALG2_T <- intersect(ALG2_T, a[,2])
ALG2_G <- Reduce(union, list(IntALG2_G, MusALG2_G, NeuALG2_G))
ALG2_G <- intersect(ALG2_G, a[,2])
ALG2_C <- Reduce(union, list(IntALG2_C, MusALG2_C, NeuALG2_C))
ALG2_C <- intersect(ALG2_C, a[,2])

Int_A <- Reduce(union, list(IntALG1_A, IntALG2_A))
Int_A <- intersect(Int_A, a[,2])
Int_T <- Reduce(union, list(IntALG1_T, IntALG2_T))
Int_T <- intersect(Int_T, a[,2])
Int_G <- Reduce(union, list(IntALG1_G, IntALG2_G))
Int_G <- intersect(Int_G, a[,2])
Int_C <- Reduce(union, list(IntALG1_C, IntALG2_C))
Int_C <- intersect(Int_C, a[,2])
Neu_A <- Reduce(union, list(NeuALG1_A, NeuALG2_A))
Neu_A <- intersect(Neu_A, a[,2])
Neu_T <- Reduce(union, list(NeuALG1_T, NeuALG2_T))
Neu_T <- intersect(Neu_T, a[,2])
Neu_G <- Reduce(union, list(NeuALG1_G, NeuALG2_G))
Neu_G <- intersect(Neu_G, a[,2])
Neu_C <- Reduce(union, list(NeuALG1_C, NeuALG2_C))
Neu_C <- intersect(Neu_C, a[,2])
Mus_A <- Reduce(union, list(MusALG1_A, MusALG2_A))
Mus_A <- intersect(Mus_A, a[,2])
Mus_T <- Reduce(union, list(MusALG1_T, MusALG2_T))
Mus_T <- intersect(Mus_T, a[,2])
Mus_G <- Reduce(union, list(MusALG1_G, MusALG2_G))
Mus_G <- intersect(Mus_G, a[,2])
Mus_C <- Reduce(union, list(MusALG1_C, MusALG2_C))
Mus_C <- intersect(Mus_C, a[,2])

ALG1 <- Reduce(union, list(ALG1_A, ALG1_T, ALG1_G, ALG1_C))
ALG2 <- Reduce(union, list(ALG2_A, ALG2_T, ALG2_G, ALG2_C))
Int <- (Reduce(union, list(Int_A, Int_T, Int_G, Int_C)))
Neu <- (Reduce(union, list(Neu_A, Neu_T, Neu_G, Neu_C)))
Mus <- Reduce(union, list(Mus_A, Mus_T, Mus_G, Mus_C))
Meta <- Reduce(union, list(ALG1, ALG2, Int, Neu, Mus))
Meta_A <- Reduce(union, list(ALG1_A, ALG2_A))
Meta_T <- Reduce(union, list(ALG1_T, ALG2_T))
Meta_G <- Reduce(union, list(ALG1_G, ALG2_G))
Meta_C <- Reduce(union, list(ALG1_C, ALG2_C))

a <- a[grep(paste(Meta, collapse = "|"), a[,2]), ]
aALG1 <- a[grep(paste(ALG1, collapse = "|"), a[,2]), ]
aALG2 <- a[grep(paste(ALG2, collapse = "|"), a[,2]), ]
aInt <- a[grep(paste(Int, collapse = "|"), a[,2]), ]
aNeu <- a[grep(paste(Neu, collapse = "|"), a[,2]), ]
aMus <- a[grep(paste(Mus, collapse = "|"), a[,2]), ]

aInt_A <- a[grep(paste(Int_A, collapse = "|"), a[,2]), ]
aInt_T <- a[grep(paste(Int_T, collapse = "|"), a[,2]), ]
aInt_G <- a[grep(paste(Int_G, collapse = "|"), a[,2]), ]
aInt_C <- a[grep(paste(Int_C, collapse = "|"), a[,2]), ]

aNeu_A <- a[grep(paste(Neu_A, collapse = "|"), a[,2]), ]
aNeu_T <- a[grep(paste(Neu_T, collapse = "|"), a[,2]), ]
aNeu_G <- a[grep(paste(Neu_G, collapse = "|"), a[,2]), ]
aNeu_C <- a[grep(paste(Neu_C, collapse = "|"), a[,2]), ]

aMus_A <- a[grep(paste(Mus_A, collapse = "|"), a[,2]), ]
aMus_T <- a[grep(paste(Mus_T, collapse = "|"), a[,2]), ]
aMus_G <- a[grep(paste(Mus_G, collapse = "|"), a[,2]), ]
aMus_C <- a[grep(paste(Mus_C, collapse = "|"), a[,2]), ]

aIntALG1 <- a[grep(paste(IntALG2, collapse = "|"), a[,2]), ]
aIntALG2 <- a[grep(paste(IntALG2, collapse = "|"), a[,2]), ]
aNeuALG1 <- a[grep(paste(NeuALG1, collapse = "|"), a[,2]), ]
aNeuALG2 <- a[grep(paste(NeuALG2, collapse = "|"), a[,2]), ]
aMusALG1 <- a[grep(paste(MusALG1, collapse = "|"), a[,2]), ]
aMusALG2 <- a[grep(paste(MusALG2, collapse = "|"), a[,2]), ]

aIntALG1_A <- a[grep(paste(IntALG1_A, collapse = "|"), a[,2]), ]
aIntALG1_T <- a[grep(paste(IntALG1_T, collapse = "|"), a[,2]), ]
aIntALG1_G <- a[grep(paste(IntALG1_G, collapse = "|"), a[,2]), ]
aIntALG1_C <- a[grep(paste(IntALG1_C, collapse = "|"), a[,2]), ]

aIntALG2_A <- a[grep(paste(IntALG2_A, collapse = "|"), a[,2]), ]
aIntALG2_T <- a[grep(paste(IntALG2_T, collapse = "|"), a[,2]), ]
aIntALG2_G <- a[grep(paste(IntALG2_G, collapse = "|"), a[,2]), ]
aIntALG2_C <- a[grep(paste(IntALG2_C, collapse = "|"), a[,2]), ]

aNeuALG1_A <- a[grep(paste(NeuALG1_A, collapse = "|"), a[,2]), ]
aNeuALG1_T <- a[grep(paste(NeuALG1_T, collapse = "|"), a[,2]), ]
aNeuALG1_G <- a[grep(paste(NeuALG1_G, collapse = "|"), a[,2]), ]
aNeuALG1_C <- a[grep(paste(NeuALG1_C, collapse = "|"), a[,2]), ]

aNeuALG2_A <- a[grep(paste(NeuALG2_A, collapse = "|"), a[,2]), ]
aNeuALG2_T <- a[grep(paste(NeuALG2_T, collapse = "|"), a[,2]), ]
aNeuALG2_G <- a[grep(paste(NeuALG2_G, collapse = "|"), a[,2]), ]
aNeuALG2_C <- a[grep(paste(NeuALG2_C, collapse = "|"), a[,2]), ]

aMusALG1_A <- a[grep(paste(MusALG1_A, collapse = "|"), a[,2]), ]
aMusALG1_T <- a[grep(paste(MusALG1_T, collapse = "|"), a[,2]), ]
aMusALG1_G <- a[grep(paste(MusALG1_G, collapse = "|"), a[,2]), ]
aMusALG1_C <- a[grep(paste(MusALG1_C, collapse = "|"), a[,2]), ]

aMusALG2_A <- a[grep(paste(MusALG2_A, collapse = "|"), a[,2]), ]
aMusALG2_T <- a[grep(paste(MusALG2_T, collapse = "|"), a[,2]), ]
aMusALG2_G <- a[grep(paste(MusALG2_G, collapse = "|"), a[,2]), ]
aMusALG2_C <- a[grep(paste(MusALG2_C, collapse = "|"), a[,2]), ]

######################################################################################################
#aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1, (1/NROW(aPlcehldr))/(NROW(Int)), 0)#####
######################################################################################################
###########################################################################################################################
#Meta nucleotide
#####
{
  datalist_Count = list()

  for (Filename in Meta_A) {
    aPlcehldr <- a[grep(paste(Filename, collapse = "|"), a[,2]), ]
    aPlcehldr$Count_1 <- NA
    aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_2 <- NA
    aPlcehldr$Count_2 <- ifelse(as.numeric(aPlcehldr[,24]) == 2 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_3 <- NA
    aPlcehldr$Count_3 <- ifelse(as.numeric(aPlcehldr[,24]) == 3 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_4 <- NA
    aPlcehldr$Count_4 <- ifelse(as.numeric(aPlcehldr[,24]) == 4 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_5 <- NA
    aPlcehldr$Count_5 <- ifelse(as.numeric(aPlcehldr[,24]) == 5 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_6 <- NA
    aPlcehldr$Count_6 <- ifelse(as.numeric(aPlcehldr[,24]) == 6 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_7 <- NA
    aPlcehldr$Count_7 <- ifelse(as.numeric(aPlcehldr[,24]) == 7 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_8 <- NA
    aPlcehldr$Count_8 <- ifelse(as.numeric(aPlcehldr[,24]) == 8 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_9 <- NA
    aPlcehldr$Count_9 <- ifelse(as.numeric(aPlcehldr[,24]) == 9 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_10 <- NA
    aPlcehldr$Count_10 <- ifelse(as.numeric(aPlcehldr[,24]) == 10 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_11 <- NA
    aPlcehldr$Count_11 <- ifelse(as.numeric(aPlcehldr[,24]) == 11 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_12 <- NA
    aPlcehldr$Count_12 <- ifelse(as.numeric(aPlcehldr[,24]) == 12 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_13 <- NA
    aPlcehldr$Count_13 <- ifelse(as.numeric(aPlcehldr[,24]) == 13 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_14 <- NA
    aPlcehldr$Count_14 <- ifelse(as.numeric(aPlcehldr[,24]) == 14 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_15 <- NA
    aPlcehldr$Count_15 <- ifelse(as.numeric(aPlcehldr[,24]) == 15 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_16 <- NA
    aPlcehldr$Count_16 <- ifelse(as.numeric(aPlcehldr[,24]) == 16 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_17 <- NA
    aPlcehldr$Count_17 <- ifelse(as.numeric(aPlcehldr[,24]) == 17 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_18 <- NA
    aPlcehldr$Count_18 <- ifelse(as.numeric(aPlcehldr[,24]) == 18 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_19 <- NA
    aPlcehldr$Count_19 <- ifelse(as.numeric(aPlcehldr[,24]) == 19 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_20 <- NA
    aPlcehldr$Count_20 <- ifelse(as.numeric(aPlcehldr[,24]) == 20 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_21 <- NA
    aPlcehldr$Count_21 <- ifelse(as.numeric(aPlcehldr[,24]) == 21 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_22 <- NA
    aPlcehldr$Count_22 <- ifelse(as.numeric(aPlcehldr[,24]) == 22 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_23 <- NA
    aPlcehldr$Count_23 <- ifelse(as.numeric(aPlcehldr[,24]) == 23 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_24 <- NA
    aPlcehldr$Count_24 <- ifelse(as.numeric(aPlcehldr[,24]) == 24 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_25 <- NA
    aPlcehldr$Count_25 <- ifelse(as.numeric(aPlcehldr[,24]) == 25 & aPlcehldr[,23] == "A", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)

    datalist_Count[[Filename]] <- aPlcehldr
  }

  nt_a1 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_1')), na.rm = TRUE)
  nt_a2 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_2')), na.rm = TRUE)
  nt_a3 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_3')), na.rm = TRUE)
  nt_a4 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_4')), na.rm = TRUE)
  nt_a5 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_5')), na.rm = TRUE)
  nt_a6 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_6')), na.rm = TRUE)
  nt_a7 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_7')), na.rm = TRUE)
  nt_a8 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_8')), na.rm = TRUE)
  nt_a9 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_9')), na.rm = TRUE)
  nt_a10 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_10')), na.rm = TRUE)
  nt_a11 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_11')), na.rm = TRUE)
  nt_a12 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_12')), na.rm = TRUE)
  nt_a13 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_13')), na.rm = TRUE)
  nt_a14 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_14')), na.rm = TRUE)
  nt_a15 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_15')), na.rm = TRUE)
  nt_a16 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_16')), na.rm = TRUE)
  nt_a17 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_17')), na.rm = TRUE)
  nt_a18 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_18')), na.rm = TRUE)
  nt_a19 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_19')), na.rm = TRUE)
  nt_a20 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_20')), na.rm = TRUE)
  nt_a21 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_21')), na.rm = TRUE)
  nt_a22 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_22')), na.rm = TRUE)
  nt_a23 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_23')), na.rm = TRUE)
  nt_a24 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_24')), na.rm = TRUE)
  nt_a25 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_25')), na.rm = TRUE)

  datalist_Count = list()

  for (Filename in Meta_T) {
    aPlcehldr <- a[grep(paste(Filename, collapse = "|"), a[,2]), ]
    aPlcehldr$Count_1 <- NA
    aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_2 <- NA
    aPlcehldr$Count_2 <- ifelse(as.numeric(aPlcehldr[,24]) == 2 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_3 <- NA
    aPlcehldr$Count_3 <- ifelse(as.numeric(aPlcehldr[,24]) == 3 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_4 <- NA
    aPlcehldr$Count_4 <- ifelse(as.numeric(aPlcehldr[,24]) == 4 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_5 <- NA
    aPlcehldr$Count_5 <- ifelse(as.numeric(aPlcehldr[,24]) == 5 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_6 <- NA
    aPlcehldr$Count_6 <- ifelse(as.numeric(aPlcehldr[,24]) == 6 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_7 <- NA
    aPlcehldr$Count_7 <- ifelse(as.numeric(aPlcehldr[,24]) == 7 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_8 <- NA
    aPlcehldr$Count_8 <- ifelse(as.numeric(aPlcehldr[,24]) == 8 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_9 <- NA
    aPlcehldr$Count_9 <- ifelse(as.numeric(aPlcehldr[,24]) == 9 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_10 <- NA
    aPlcehldr$Count_10 <- ifelse(as.numeric(aPlcehldr[,24]) == 10 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_11 <- NA
    aPlcehldr$Count_11 <- ifelse(as.numeric(aPlcehldr[,24]) == 11 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_12 <- NA
    aPlcehldr$Count_12 <- ifelse(as.numeric(aPlcehldr[,24]) == 12 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_13 <- NA
    aPlcehldr$Count_13 <- ifelse(as.numeric(aPlcehldr[,24]) == 13 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_14 <- NA
    aPlcehldr$Count_14 <- ifelse(as.numeric(aPlcehldr[,24]) == 14 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_15 <- NA
    aPlcehldr$Count_15 <- ifelse(as.numeric(aPlcehldr[,24]) == 15 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_16 <- NA
    aPlcehldr$Count_16 <- ifelse(as.numeric(aPlcehldr[,24]) == 16 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_17 <- NA
    aPlcehldr$Count_17 <- ifelse(as.numeric(aPlcehldr[,24]) == 17 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_18 <- NA
    aPlcehldr$Count_18 <- ifelse(as.numeric(aPlcehldr[,24]) == 18 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_19 <- NA
    aPlcehldr$Count_19 <- ifelse(as.numeric(aPlcehldr[,24]) == 19 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_20 <- NA
    aPlcehldr$Count_20 <- ifelse(as.numeric(aPlcehldr[,24]) == 20 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_21 <- NA
    aPlcehldr$Count_21 <- ifelse(as.numeric(aPlcehldr[,24]) == 21 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_22 <- NA
    aPlcehldr$Count_22 <- ifelse(as.numeric(aPlcehldr[,24]) == 22 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_23 <- NA
    aPlcehldr$Count_23 <- ifelse(as.numeric(aPlcehldr[,24]) == 23 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_24 <- NA
    aPlcehldr$Count_24 <- ifelse(as.numeric(aPlcehldr[,24]) == 24 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_25 <- NA
    aPlcehldr$Count_25 <- ifelse(as.numeric(aPlcehldr[,24]) == 25 & aPlcehldr[,23] == "T", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)

    datalist_Count[[Filename]] <- aPlcehldr
  }

  nt_t1 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_1')), na.rm = TRUE)
  nt_t2 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_2')), na.rm = TRUE)
  nt_t3 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_3')), na.rm = TRUE)
  nt_t4 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_4')), na.rm = TRUE)
  nt_t5 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_5')), na.rm = TRUE)
  nt_t6 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_6')), na.rm = TRUE)
  nt_t7 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_7')), na.rm = TRUE)
  nt_t8 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_8')), na.rm = TRUE)
  nt_t9 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_9')), na.rm = TRUE)
  nt_t10 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_10')), na.rm = TRUE)
  nt_t11 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_11')), na.rm = TRUE)
  nt_t12 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_12')), na.rm = TRUE)
  nt_t13 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_13')), na.rm = TRUE)
  nt_t14 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_14')), na.rm = TRUE)
  nt_t15 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_15')), na.rm = TRUE)
  nt_t16 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_16')), na.rm = TRUE)
  nt_t17 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_17')), na.rm = TRUE)
  nt_t18 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_18')), na.rm = TRUE)
  nt_t19 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_19')), na.rm = TRUE)
  nt_t20 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_20')), na.rm = TRUE)
  nt_t21 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_21')), na.rm = TRUE)
  nt_t22 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_22')), na.rm = TRUE)
  nt_t23 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_23')), na.rm = TRUE)
  nt_t24 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_24')), na.rm = TRUE)
  nt_t25 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_25')), na.rm = TRUE)

  datalist_Count = list()

  for (Filename in Meta_G) {
    aPlcehldr <- a[grep(paste(Filename, collapse = "|"), a[,2]), ]
    aPlcehldr$Count_1 <- NA
    aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_2 <- NA
    aPlcehldr$Count_2 <- ifelse(as.numeric(aPlcehldr[,24]) == 2 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_3 <- NA
    aPlcehldr$Count_3 <- ifelse(as.numeric(aPlcehldr[,24]) == 3 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_4 <- NA
    aPlcehldr$Count_4 <- ifelse(as.numeric(aPlcehldr[,24]) == 4 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_5 <- NA
    aPlcehldr$Count_5 <- ifelse(as.numeric(aPlcehldr[,24]) == 5 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_6 <- NA
    aPlcehldr$Count_6 <- ifelse(as.numeric(aPlcehldr[,24]) == 6 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_7 <- NA
    aPlcehldr$Count_7 <- ifelse(as.numeric(aPlcehldr[,24]) == 7 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_8 <- NA
    aPlcehldr$Count_8 <- ifelse(as.numeric(aPlcehldr[,24]) == 8 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_9 <- NA
    aPlcehldr$Count_9 <- ifelse(as.numeric(aPlcehldr[,24]) == 9 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_10 <- NA
    aPlcehldr$Count_10 <- ifelse(as.numeric(aPlcehldr[,24]) == 10 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_11 <- NA
    aPlcehldr$Count_11 <- ifelse(as.numeric(aPlcehldr[,24]) == 11 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_12 <- NA
    aPlcehldr$Count_12 <- ifelse(as.numeric(aPlcehldr[,24]) == 12 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_13 <- NA
    aPlcehldr$Count_13 <- ifelse(as.numeric(aPlcehldr[,24]) == 13 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_14 <- NA
    aPlcehldr$Count_14 <- ifelse(as.numeric(aPlcehldr[,24]) == 14 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_15 <- NA
    aPlcehldr$Count_15 <- ifelse(as.numeric(aPlcehldr[,24]) == 15 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_16 <- NA
    aPlcehldr$Count_16 <- ifelse(as.numeric(aPlcehldr[,24]) == 16 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_17 <- NA
    aPlcehldr$Count_17 <- ifelse(as.numeric(aPlcehldr[,24]) == 17 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_18 <- NA
    aPlcehldr$Count_18 <- ifelse(as.numeric(aPlcehldr[,24]) == 18 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_19 <- NA
    aPlcehldr$Count_19 <- ifelse(as.numeric(aPlcehldr[,24]) == 19 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_20 <- NA
    aPlcehldr$Count_20 <- ifelse(as.numeric(aPlcehldr[,24]) == 20 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_21 <- NA
    aPlcehldr$Count_21 <- ifelse(as.numeric(aPlcehldr[,24]) == 21 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_22 <- NA
    aPlcehldr$Count_22 <- ifelse(as.numeric(aPlcehldr[,24]) == 22 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_23 <- NA
    aPlcehldr$Count_23 <- ifelse(as.numeric(aPlcehldr[,24]) == 23 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_24 <- NA
    aPlcehldr$Count_24 <- ifelse(as.numeric(aPlcehldr[,24]) == 24 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_25 <- NA
    aPlcehldr$Count_25 <- ifelse(as.numeric(aPlcehldr[,24]) == 25 & aPlcehldr[,23] == "G", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)

    datalist_Count[[Filename]] <- aPlcehldr
  }

  nt_g1 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_1')), na.rm = TRUE)
  nt_g2 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_2')), na.rm = TRUE)
  nt_g3 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_3')), na.rm = TRUE)
  nt_g4 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_4')), na.rm = TRUE)
  nt_g5 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_5')), na.rm = TRUE)
  nt_g6 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_6')), na.rm = TRUE)
  nt_g7 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_7')), na.rm = TRUE)
  nt_g8 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_8')), na.rm = TRUE)
  nt_g9 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_9')), na.rm = TRUE)
  nt_g10 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_10')), na.rm = TRUE)
  nt_g11 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_11')), na.rm = TRUE)
  nt_g12 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_12')), na.rm = TRUE)
  nt_g13 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_13')), na.rm = TRUE)
  nt_g14 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_14')), na.rm = TRUE)
  nt_g15 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_15')), na.rm = TRUE)
  nt_g16 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_16')), na.rm = TRUE)
  nt_g17 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_17')), na.rm = TRUE)
  nt_g18 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_18')), na.rm = TRUE)
  nt_g19 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_19')), na.rm = TRUE)
  nt_g20 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_20')), na.rm = TRUE)
  nt_g21 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_21')), na.rm = TRUE)
  nt_g22 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_22')), na.rm = TRUE)
  nt_g23 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_23')), na.rm = TRUE)
  nt_g24 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_24')), na.rm = TRUE)
  nt_g25 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_25')), na.rm = TRUE)

  datalist_Count = list()

  for (Filename in Meta_C) {
    aPlcehldr <- a[grep(paste(Filename, collapse = "|"), a[,2]), ]
    aPlcehldr$Count_1 <- NA
    aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_2 <- NA
    aPlcehldr$Count_2 <- ifelse(as.numeric(aPlcehldr[,24]) == 2 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_3 <- NA
    aPlcehldr$Count_3 <- ifelse(as.numeric(aPlcehldr[,24]) == 3 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_4 <- NA
    aPlcehldr$Count_4 <- ifelse(as.numeric(aPlcehldr[,24]) == 4 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_5 <- NA
    aPlcehldr$Count_5 <- ifelse(as.numeric(aPlcehldr[,24]) == 5 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_6 <- NA
    aPlcehldr$Count_6 <- ifelse(as.numeric(aPlcehldr[,24]) == 6 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_7 <- NA
    aPlcehldr$Count_7 <- ifelse(as.numeric(aPlcehldr[,24]) == 7 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_8 <- NA
    aPlcehldr$Count_8 <- ifelse(as.numeric(aPlcehldr[,24]) == 8 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_9 <- NA
    aPlcehldr$Count_9 <- ifelse(as.numeric(aPlcehldr[,24]) == 9 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_10 <- NA
    aPlcehldr$Count_10 <- ifelse(as.numeric(aPlcehldr[,24]) == 10 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_11 <- NA
    aPlcehldr$Count_11 <- ifelse(as.numeric(aPlcehldr[,24]) == 11 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_12 <- NA
    aPlcehldr$Count_12 <- ifelse(as.numeric(aPlcehldr[,24]) == 12 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_13 <- NA
    aPlcehldr$Count_13 <- ifelse(as.numeric(aPlcehldr[,24]) == 13 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_14 <- NA
    aPlcehldr$Count_14 <- ifelse(as.numeric(aPlcehldr[,24]) == 14 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_15 <- NA
    aPlcehldr$Count_15 <- ifelse(as.numeric(aPlcehldr[,24]) == 15 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_16 <- NA
    aPlcehldr$Count_16 <- ifelse(as.numeric(aPlcehldr[,24]) == 16 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_17 <- NA
    aPlcehldr$Count_17 <- ifelse(as.numeric(aPlcehldr[,24]) == 17 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_18 <- NA
    aPlcehldr$Count_18 <- ifelse(as.numeric(aPlcehldr[,24]) == 18 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_19 <- NA
    aPlcehldr$Count_19 <- ifelse(as.numeric(aPlcehldr[,24]) == 19 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_20 <- NA
    aPlcehldr$Count_20 <- ifelse(as.numeric(aPlcehldr[,24]) == 20 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_21 <- NA
    aPlcehldr$Count_21 <- ifelse(as.numeric(aPlcehldr[,24]) == 21 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_22 <- NA
    aPlcehldr$Count_22 <- ifelse(as.numeric(aPlcehldr[,24]) == 22 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_23 <- NA
    aPlcehldr$Count_23 <- ifelse(as.numeric(aPlcehldr[,24]) == 23 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_24 <- NA
    aPlcehldr$Count_24 <- ifelse(as.numeric(aPlcehldr[,24]) == 24 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)
    aPlcehldr$Count_25 <- NA
    aPlcehldr$Count_25 <- ifelse(as.numeric(aPlcehldr[,24]) == 25 & aPlcehldr[,23] == "C", (1/NROW(aPlcehldr))/(NROW(Meta)), 0)

    datalist_Count[[Filename]] <- aPlcehldr
  }

  nt_c1 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_1')), na.rm = TRUE)
  nt_c2 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_2')), na.rm = TRUE)
  nt_c3 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_3')), na.rm = TRUE)
  nt_c4 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_4')), na.rm = TRUE)
  nt_c5 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_5')), na.rm = TRUE)
  nt_c6 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_6')), na.rm = TRUE)
  nt_c7 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_7')), na.rm = TRUE)
  nt_c8 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_8')), na.rm = TRUE)
  nt_c9 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_9')), na.rm = TRUE)
  nt_c10 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_10')), na.rm = TRUE)
  nt_c11 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_11')), na.rm = TRUE)
  nt_c12 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_12')), na.rm = TRUE)
  nt_c13 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_13')), na.rm = TRUE)
  nt_c14 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_14')), na.rm = TRUE)
  nt_c15 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_15')), na.rm = TRUE)
  nt_c16 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_16')), na.rm = TRUE)
  nt_c17 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_17')), na.rm = TRUE)
  nt_c18 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_18')), na.rm = TRUE)
  nt_c19 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_19')), na.rm = TRUE)
  nt_c20 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_20')), na.rm = TRUE)
  nt_c21 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_21')), na.rm = TRUE)
  nt_c22 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_22')), na.rm = TRUE)
  nt_c23 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_23')), na.rm = TRUE)
  nt_c24 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_24')), na.rm = TRUE)
  nt_c25 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_25')), na.rm = TRUE)


  total_SNP_a <- c(nt_a1, nt_a2, nt_a3, nt_a4, nt_a5, nt_a6, nt_a7, nt_a8, nt_a9, nt_a10, nt_a11, nt_a12,
                    nt_a13, nt_a14, nt_a15, nt_a16, nt_a17, nt_a18, nt_a19, nt_a20, nt_a21, nt_a22, nt_a23,
                    nt_a24, nt_a25)
  total_SNP_t <- c(nt_t1, nt_t2, nt_t3, nt_t4, nt_t5, nt_t6, nt_t7, nt_t8, nt_t9, nt_t10, nt_t11, nt_t12,
                    nt_t13, nt_t14, nt_t15, nt_t16, nt_t17, nt_t18, nt_t19, nt_t20, nt_t21, nt_t22, nt_t23,
                    nt_t24, nt_t25)
  total_SNP_g <- c(nt_g1, nt_g2, nt_g3, nt_g4, nt_g5, nt_g6, nt_g7, nt_g8, nt_g9, nt_g10, nt_g11, nt_g12,
                    nt_g13, nt_g14, nt_g15, nt_g16, nt_g17, nt_g18, nt_g19, nt_g20, nt_g21, nt_g22, nt_g23,
                    nt_g24, nt_g25)
  total_SNP_c <- c(nt_c1, nt_c2, nt_c3, nt_c4, nt_c5, nt_c6, nt_c7, nt_c8, nt_c9, nt_c10, nt_c11, nt_c12,
                    nt_c13, nt_c14, nt_c15, nt_c16, nt_c17, nt_c18, nt_c19, nt_c20, nt_c21, nt_c22, nt_c23,
                    nt_c24, nt_c25)

  total_pos <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
                 "17", "18", "19", "20", "21", "22", "23", "24", "25")

  xlab <- list(title = "Nucleotide position",
               categoryorder = "array",
               categoryarray = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                                 "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25"))
  ylab <- list(title = paste0(c(rep("&nbsp;", 20),
                                "# isomiRs/miRNA",
                                rep("&nbsp;", 20),
                                rep("\n&nbsp;", 3)),
                              collapse = ""))

  p <- plot_ly(x = total_pos, y = total_SNP_c, type = "bar", name = "C", text = total_SNP_c) %>%
    add_trace(y = total_SNP_g, name = "G", text = total_SNP_g) %>%
    add_trace(y = total_SNP_t, name = "T", text = total_SNP_t) %>%
    add_trace(y = total_SNP_a, name = "A", text = total_SNP_a) %>%
    layout(xaxis=xlab, yaxis=ylab, barmode = 'stack', title = "Proportion of isomiR single nucleotide variants (Meta)")
  p
}
