a <- read.delim("Reference_no_SNP_meta.txt", header = F)

Int1 <- read.csv("IntALG1 vs N2 Top - Pos.csv", header=F)$V2
Int2 <- read.csv("IntALG2 vs N2 Top - Pos.csv", header=F)$V2
Neu1 <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header=F)$V2
Neu2 <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header=F)$V2
Mus1 <- read.csv("MusALG1 vs N2 Top - Pos.csv", header=F)$V2
Mus2 <- read.csv("MusALG2 vs N2 Top - Pos.csv", header=F)$V2

Int1Mus1 <- read.csv("IntALG1 vs MusALG1 Top.csv", header=F)$V2
Int2Mus2 <- read.csv("IntALG2 vs MusALG2 Top.csv", header=F)$V2
Int1Neu1 <- read.csv("IntALG1 vs NeuALG1 Top.csv", header=F)$V2
Int2Neu2 <- read.csv("IntALG2 vs NeuALG2 Top.csv", header=F)$V2
Neu1Mus1 <- read.csv("NeuALG1 vs MusALG1 Top.csv", header=F)$V2
Neu2Mus2 <- read.csv("NeuALG2 vs MusALG2 Top.csv", header=F)$V2

###############################################################

Int <- Reduce(union, list(Int1, Int2))
Int <-paste0("^",Int,"$")
Neu <- Reduce(union, list(Neu1, Neu2))
Neu <-paste0("^",Neu,"$")
Mus <- Reduce(union, list(Mus1, Mus2))
Mus <-paste0("^",Mus,"$")

IntMus <- Reduce(union, list(Int1Mus1, Int2Mus2))
IntMus <-paste0("^",IntMus,"$")
IntNeu <- Reduce(union, list(Int1Neu1, Int2Neu2))
IntNeu <-paste0("^",IntNeu,"$")
NeuMus <- Reduce(union, list(Neu1Mus1, Neu2Mus2))
NeuMus <-paste0("^",NeuMus,"$")

#####
test <- intersect(Int1, a[,1])
test <- intersect(Int2, a[,1])
test <- intersect(Neu1, a[,1])
test <- intersect(Neu2, a[,1])
test <- intersect(Mus1, a[,1])
test <- intersect(Mus2, a[,1])
#####

aInt <- a[grep(paste(Int, collapse = "|"), a[,1]), ]
aNeu <- a[grep(paste(Neu, collapse = "|"), a[,1]), ]
aMus <- a[grep(paste(Mus, collapse = "|"), a[,1]), ]

aIntMus <- a[grep(paste(IntMus, collapse = "|"), a[,1]), ]
aIntNeu <- a[grep(paste(IntNeu, collapse = "|"), a[,1]), ]
aNeuMus <- a[grep(paste(NeuMus, collapse = "|"), a[,1]), ]

################################################################

Hrpn_Int <- aInt[,21]
Hrpn_Neu <- aNeu[,21]
Hrpn_Mus <- aMus[,21]

Hrpn_IntMus <- aIntMus[,21]
Hrpn_IntNeu <- aIntNeu[,21]
Hrpn_NeuMus <- aNeuMus[,21]

################################################################

write.table(Hrpn_Int, file="Hairpin seq intestine.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(Hrpn_Neu, file="Hairpin seq neuron.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(Hrpn_Mus, file="Hairpin seq muscle.csv", sep = ",", col.names = NA, qmethod = "double")

write.table(Hrpn_IntMus, file="Hairpin seq intestine vs muscle.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(Hrpn_IntNeu, file="Hairpin seq intestine vs neuron.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(Hrpn_NeuMus, file="Hairpin seq neuron vs muscle.csv", sep = ",", col.names = NA, qmethod = "double")

################################################################

write.table(aInt, file="Hairpin table intestine.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(aNeu, file="Hairpin table neuron.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(aMus, file="Hairpin table muscle.csv", sep = ",", col.names = NA, qmethod = "double")

write.table(aIntMus, file="Hairpin table intestine vs muscle.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(aIntNeu, file="Hairpin table intestine vs neuron.csv", sep = ",", col.names = NA, qmethod = "double")
write.table(aNeuMus, file="Hairpin table neuron vs muscle.csv", sep = ",", col.names = NA, qmethod = "double")

#######################################
#######################################

miR90_Int <- read.delim("miR-90 Intestine.txt", header=FALSE)$V1
miR90_Int <- paste0("^",miR90_Int, "$")
miR90_Neu <- read.delim("miR-90 Neuron.txt", header=FALSE)$V1
miR90_Neu <- paste0("^",miR90_Neu, "$")
miR90_Mus <- read.delim("miR-90 Muscle.txt", header=FALSE)$V1
miR90_Mus <- paste0("^",miR90_Mus, "$")

amiR90_Int <- a[grep(paste(miR90_Int, collapse = "|"), a[,1]), ]
amiR90_Neu <- a[grep(paste(miR90_Neu, collapse = "|"), a[,1]), ]
amiR90_Mus <- a[grep(paste(miR90_Mus, collapse = "|"), a[,1]), ]
