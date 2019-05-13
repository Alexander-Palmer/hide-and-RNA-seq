###Find sum of columns containing original SNP type

fileNames <- Sys.glob("*.txt")

#Create average columns
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample[paste('av_284')] <- (sample[["A284_1"]] + sample[["A284_2"]])/2
  sample[paste('av_285')] <- (sample[["A285_1"]] + sample[["A285_2"]])/2
  sample[paste('av_286')] <- (sample[["A286_1"]] + sample[["A286_2"]])/2
  sample[paste('av_287')] <- (sample[["A287_1"]] + sample[["A287_2"]])/2
  sample[paste('av_290')] <- (sample[["A290_1"]] + sample[["A290_2"]])/2
  sample[paste('av_291')] <- (sample[["A291_1"]] + sample[["A291_2"]])/2
  sample[paste('Intcomb')] <- (sample[["av_284"]] + sample[["av_285"]])/2
  sample[paste('Neucomb')] <- (sample[["av_286"]] + sample[["av_287"]])/2
  sample[paste('Muscomb')] <- (sample[["av_290"]] + sample[["av_291"]])/2
  sample[paste('av_N2')] <- (sample[["N2_1"]] + sample[["N2_2"]])/2
  sample[paste('av_ALG1')] <- (sample[["av_284"]] + sample[["av_286"]] + sample[["av_290"]])/3
  sample[paste('av_ALG2')] <- (sample[["av_285"]] + sample[["av_287"]] + sample[["av_291"]])/3
  sample["SNP_A"] <- 0
  sample["SNP_T"] <- 0
  sample["SNP_G"] <- 0
  sample["SNP_C"] <- 0
  sample["N2_A"] <- 0
  sample["N2_T"] <- 0
  sample["N2_G"] <- 0
  sample["N2_C"] <- 0


  write.table(sample,
              fileName,
              append = FALSE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
}

IntAlg1_A <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG1_A.txt", header=F)$V1
IntAlg1_A <- as.character(IntAlg1_A)
IntAlg1_T <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG1_U.txt", header=F)$V1
IntAlg1_T <- as.character(IntAlg1_T)
IntAlg1_G <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG1_G.txt", header=F)$V1
IntAlg1_G <- as.character(IntAlg1_G)
IntAlg1_C <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG1_C.txt", header=F)$V1
IntAlg1_C <- as.character(IntAlg1_C)

NeuAlg1_A <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG1_A.txt", header=F)$V1
NeuAlg1_A <- as.character(NeuAlg1_A)
NeuAlg1_T <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG1_U.txt", header=F)$V1
NeuAlg1_T <- as.character(NeuAlg1_T)
NeuAlg1_G <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG1_G.txt", header=F)$V1
NeuAlg1_G <- as.character(NeuAlg1_G)
NeuAlg1_C <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG1_C.txt", header=F)$V1
NeuAlg1_C <- as.character(NeuAlg1_C)

MusAlg1_A <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG1_A.txt", header=F)$V1
MusAlg1_A <- as.character(MusAlg1_A)
MusAlg1_T <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG1_U.txt", header=F)$V1
MusAlg1_T <- as.character(MusAlg1_T)
MusAlg1_G <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG1_G.txt", header=F)$V1
MusAlg1_G <- as.character(MusAlg1_G)
MusAlg1_C <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG1_C.txt", header=F)$V1
MusAlg1_C <- as.character(MusAlg1_C)


IntAlg2_A <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG2_A.txt", header=F)$V1
IntAlg2_A <- as.character(IntAlg2_A)
IntAlg2_T <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG2_U.txt", header=F)$V1
IntAlg2_T <- as.character(IntAlg2_T)
IntAlg2_G <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG2_G.txt", header=F)$V1
IntAlg2_G <- as.character(IntAlg2_G)
IntAlg2_C <- read.delim("X_List_of_sig_SNP_miRNA_int_ALG2_C.txt", header=F)$V1
IntAlg2_C <- as.character(IntAlg2_C)


NeuAlg2_A <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG2_A.txt", header=F)$V1
NeuAlg2_A <- as.character(NeuAlg2_A)
NeuAlg2_T <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG2_U.txt", header=F)$V1
NeuAlg2_T <- as.character(NeuAlg2_T)
NeuAlg2_G <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG2_G.txt", header=F)$V1
NeuAlg2_G <- as.character(NeuAlg2_G)
NeuAlg2_C <- read.delim("X_List_of_sig_SNP_miRNA_neu_ALG2_C.txt", header=F)$V1
NeuAlg2_C <- as.character(NeuAlg2_C)


MusAlg2_A <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG2_A.txt", header=F)$V1
MusAlg2_A <- as.character(MusAlg2_A)
MusAlg2_T <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG2_U.txt", header=F)$V1
MusAlg2_T <- as.character(MusAlg2_T)
MusAlg2_G <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG2_G.txt", header=F)$V1
MusAlg2_G <- as.character(MusAlg2_G)
MusAlg2_C <- read.delim("X_List_of_sig_SNP_miRNA_mus_ALG2_C.txt", header=F)$V1
MusAlg2_C <- as.character(MusAlg2_C)

ALG1_A <- Reduce(union, list(IntAlg1_A, MusAlg1_A, NeuAlg1_A))
ALG1_T <- Reduce(union, list(IntAlg1_T, MusAlg1_T, NeuAlg1_T))
ALG1_G <- Reduce(union, list(IntAlg1_G, MusAlg1_G, NeuAlg1_G))
ALG1_C <- Reduce(union, list(IntAlg1_C, MusAlg1_C, NeuAlg1_C))
ALG2_A <- Reduce(union, list(IntAlg2_A, MusAlg2_A, NeuAlg2_A))
ALG2_T <- Reduce(union, list(IntAlg2_T, MusAlg2_T, NeuAlg2_T))
ALG2_G <- Reduce(union, list(IntAlg2_G, MusAlg2_G, NeuAlg2_G))
ALG2_C <- Reduce(union, list(IntAlg2_C, MusAlg2_C, NeuAlg2_C))

Int_A <- Reduce(union, list(IntAlg1_A, IntAlg2_A))
Int_T <- Reduce(union, list(IntAlg1_T, IntAlg2_T))
Int_G <- Reduce(union, list(IntAlg1_G, IntAlg2_G))
Int_C <- Reduce(union, list(IntAlg1_C, IntAlg2_C))
Neu_A <- Reduce(union, list(NeuAlg1_A, NeuAlg2_A))
Neu_T <- Reduce(union, list(NeuAlg1_T, NeuAlg2_T))
Neu_G <- Reduce(union, list(NeuAlg1_G, NeuAlg2_G))
Neu_C <- Reduce(union, list(NeuAlg1_C, NeuAlg2_C))
Mus_A <- Reduce(union, list(MusAlg1_A, MusAlg2_A))
Mus_T <- Reduce(union, list(MusAlg1_T, MusAlg2_T))
Mus_G <- Reduce(union, list(MusAlg1_G, MusAlg2_G))
Mus_C <- Reduce(union, list(MusAlg1_C, MusAlg2_C))


###########################################################################
#SNPs (A) common to intestine Alg1 (>7nt)

datalist_ALG2_T = list()

for (FileName in IntAlg1_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Intestine_ALG1_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Intestine_ALG1_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to intestine alg1 (>7nt)

datalist_ALG2_T = list()

for (FileName in IntAlg1_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Intestine_ALG1_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Intestine_ALG1_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to intestine Alg1 (>7nt)

datalist_ALG2_G = list()

for (FileName in IntAlg1_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Intestine_ALG1_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Intestine_ALG1_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to intestine alg1 (>7nt)

datalist_ALG2_C = list()

for (FileName in IntAlg1_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,24], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_C[[FileName]] <- sample

}

datalist_ALG2_C <- do.call(rbind.data.frame, datalist_ALG2_C)

SNP_A_ALG2 <- sum(datalist_ALG2_C$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_C$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_C$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_C$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_C$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_C$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_C$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_C$N2_C)

ALG2_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_C) <- c("Total Reads")
ALG2_output_C
write.table(ALG2_output_C, "Intestine_ALG1_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Intestine_ALG1_output_C_7_meta.txt", sep="\t", col.names = NA)

#########################################################################################################################
###########################################################################
#SNPs (A) common to neuron alg1 (>7nt)

datalist_ALG2_A = list()

for (FileName in NeuAlg1_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_A[[FileName]] <- sample

}

datalist_ALG2_A <- do.call(rbind.data.frame, datalist_ALG2_A)

SNP_A_ALG2 <- sum(datalist_ALG2_A$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_A$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_A$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_A$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_A$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_A$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_A$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_A$N2_C)

ALG2_output_A <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_A) <- c("Total Reads")
ALG2_output_A
write.table(ALG2_output_A, "Neuron_ALG1_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_A, "Neuron_ALG1_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to neuron alg1 (>7nt)

datalist_ALG2_T = list()

for (FileName in NeuAlg1_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Neuron_ALG1_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Neuron_ALG1_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to neuron alg1 (>7nt)

datalist_ALG2_G = list()

for (FileName in NeuAlg1_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Neuron_ALG1_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Neuron_ALG1_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to neuron alg1 (>7nt)

datalist_ALG1_C = list()

for (FileName in NeuAlg1_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,26], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG1_C[[FileName]] <- sample

}

datalist_ALG1_C <- do.call(rbind.data.frame, datalist_ALG1_C)

SNP_A_ALG1 <- sum(datalist_ALG1_C$SNP_A)
N2_A_ALG1 <- sum(datalist_ALG1_C$N2_A)
SNP_T_ALG1 <- sum(datalist_ALG1_C$SNP_T)
N2_T_ALG1 <- sum(datalist_ALG1_C$N2_T)
SNP_G_ALG1 <- sum(datalist_ALG1_C$SNP_G)
N2_G_ALG1 <- sum(datalist_ALG1_C$N2_G)
SNP_C_ALG1 <- sum(datalist_ALG1_C$SNP_C)
N2_C_ALG1 <- sum(datalist_ALG1_C$N2_C)

ALG1_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG1_output_C) <- c("Total Reads")
ALG1_output_C
write.table(ALG1_output_C, "Neuron_ALG1_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Neuron_ALG1_output_C_7_meta.txt", sep="\t", col.names = NA)

#########################################################################################################################
###########################################################################
#SNPs (A) common to muscle alg1 (>7nt)

datalist_ALG2_A = list()

for (FileName in MusAlg1_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_A[[FileName]] <- sample

}

datalist_ALG2_A <- do.call(rbind.data.frame, datalist_ALG2_A)

SNP_A_ALG2 <- sum(datalist_ALG2_A$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_A$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_A$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_A$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_A$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_A$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_A$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_A$N2_C)

ALG2_output_A <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_A) <- c("Total Reads")
ALG2_output_A
write.table(ALG2_output_A, "Muscle_ALG1_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_A, "Muscle_ALG1_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to muscle alg1 (>7nt)

datalist_ALG2_T = list()

for (FileName in MusAlg1_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Muscle_ALG1_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Muscle_ALG1_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to muscle alg1 (>7nt)

datalist_ALG2_G = list()

for (FileName in MusAlg1_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Muscle_ALG1_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Muscle_ALG1_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to muscle alg1 (>7nt)

datalist_ALG2_C = list()

for (FileName in MusAlg1_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,28], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_C[[FileName]] <- sample

}

datalist_ALG2_C <- do.call(rbind.data.frame, datalist_ALG2_C)

SNP_A_ALG2 <- sum(datalist_ALG2_C$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_C$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_C$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_C$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_C$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_C$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_C$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_C$N2_C)

ALG2_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_C) <- c("Total Reads")
ALG2_output_C
write.table(ALG2_output_C, "Muscle_ALG1_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Muscle_ALG1_output_C_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (A) common to intestine alg2 (>7nt)

datalist_ALG2_T = list()

for (FileName in IntAlg2_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Intestine_ALG2_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Intestine_ALG2_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to intestine alg2 (>7nt)

datalist_ALG2_T = list()

for (FileName in IntAlg2_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Intestine_ALG2_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Intestine_ALG2_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to intestine alg2 (>7nt)

datalist_ALG2_G = list()

for (FileName in IntAlg2_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Intestine_ALG2_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Intestine_ALG2_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to intestine alg2 (>7nt)

datalist_ALG2_C = list()

for (FileName in IntAlg2_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,25], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_C[[FileName]] <- sample

}

datalist_ALG2_C <- do.call(rbind.data.frame, datalist_ALG2_C)

SNP_A_ALG2 <- sum(datalist_ALG2_C$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_C$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_C$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_C$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_C$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_C$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_C$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_C$N2_C)

ALG2_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_C) <- c("Total Reads")
ALG2_output_C
write.table(ALG2_output_C, "Intestine_ALG2_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Intestine_ALG2_output_C_7_meta.txt", sep="\t", col.names = NA)

#########################################################################################################################
###########################################################################
#SNPs (A) common to neuron alg2 (>7nt)

datalist_ALG2_A = list()

for (FileName in NeuAlg2_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_A[[FileName]] <- sample

}

datalist_ALG2_A <- do.call(rbind.data.frame, datalist_ALG2_A)

SNP_A_ALG2 <- sum(datalist_ALG2_A$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_A$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_A$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_A$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_A$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_A$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_A$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_A$N2_C)

ALG2_output_A <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_A) <- c("Total Reads")
ALG2_output_A
write.table(ALG2_output_A, "Neuron_ALG2_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_A, "Neuron_ALG2_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to neuron alg2 (>7nt)

datalist_ALG2_T = list()

for (FileName in NeuAlg2_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Neuron_ALG2_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Neuron_ALG2_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to neuron alg2 (>7nt)

datalist_ALG2_G = list()

for (FileName in NeuAlg2_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Neuron_ALG2_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Neuron_ALG2_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to neuron alg2 (>7nt)

datalist_ALG2_C = list()

for (FileName in NeuAlg2_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,27], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_C[[FileName]] <- sample

}

datalist_ALG2_C <- do.call(rbind.data.frame, datalist_ALG2_C)

SNP_A_ALG2 <- sum(datalist_ALG2_C$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_C$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_C$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_C$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_C$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_C$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_C$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_C$N2_C)

ALG2_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_C) <- c("Total Reads")
ALG2_output_C
write.table(ALG2_output_C, "Neuron_ALG2_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Neuron_ALG2_output_C_7_meta.txt", sep="\t", col.names = NA)

#########################################################################################################################
###########################################################################
#SNPs (A) common to muscle alg2 (>7nt)

datalist_ALG2_A = list()

for (FileName in MusAlg2_A) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_A[[FileName]] <- sample

}

datalist_ALG2_A <- do.call(rbind.data.frame, datalist_ALG2_A)

SNP_A_ALG2 <- sum(datalist_ALG2_A$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_A$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_A$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_A$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_A$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_A$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_A$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_A$N2_C)

ALG2_output_A <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_A) <- c("Total Reads")
ALG2_output_A
write.table(ALG2_output_A, "Muscle_ALG2_output_A_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_A, "Muscle_ALG2_output_A_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (T) common to muscle alg2 (>7nt)

datalist_ALG2_T = list()

for (FileName in MusAlg2_T) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_T[[FileName]] <- sample

}

datalist_ALG2_T <- do.call(rbind.data.frame, datalist_ALG2_T)

SNP_A_ALG2 <- sum(datalist_ALG2_T$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_T$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_T$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_T$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_T$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_T$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_T$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_T$N2_C)

ALG2_output_T <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_T) <- c("Total Reads")
ALG2_output_T
write.table(ALG2_output_T, "Muscle_ALG2_output_T_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_T, "Muscle_ALG2_output_T_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (G) common to muscle alg2 (>7nt)

datalist_ALG2_G = list()

for (FileName in MusAlg2_G) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_G[[FileName]] <- sample

}

datalist_ALG2_G <- do.call(rbind.data.frame, datalist_ALG2_G)

SNP_A_ALG2 <- sum(datalist_ALG2_G$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_G$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_G$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_G$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_G$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_G$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_G$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_G$N2_C)

ALG2_output_G <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_G) <- c("Total Reads")
ALG2_output_G
write.table(ALG2_output_G, "Muscle_ALG2_output_G_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_G, "Muscle_ALG2_output_G_7_meta.txt", sep="\t", col.names = NA)

###########################################################################
#SNPs (C) common to muscle alg2 (>7nt)

datalist_ALG2_C = list()

for (FileName in MusAlg2_C) {
  sample <- read.delim(FileName)
  sample <- as.data.frame(sample)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,29], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & as.numeric(sample[,23]) > 7, sample[,33], 0)

  sample$FileName <- FileName
  datalist_ALG2_C[[FileName]] <- sample

}

datalist_ALG2_C <- do.call(rbind.data.frame, datalist_ALG2_C)

SNP_A_ALG2 <- sum(datalist_ALG2_C$SNP_A)
N2_A_ALG2 <- sum(datalist_ALG2_C$N2_A)
SNP_T_ALG2 <- sum(datalist_ALG2_C$SNP_T)
N2_T_ALG2 <- sum(datalist_ALG2_C$N2_T)
SNP_G_ALG2 <- sum(datalist_ALG2_C$SNP_G)
N2_G_ALG2 <- sum(datalist_ALG2_C$N2_G)
SNP_C_ALG2 <- sum(datalist_ALG2_C$SNP_C)
N2_C_ALG2 <- sum(datalist_ALG2_C$N2_C)

ALG2_output_C <- rbind(SNP_A_ALG2, SNP_T_ALG2, SNP_G_ALG2, SNP_C_ALG2,
                       N2_A_ALG2, N2_T_ALG2, N2_G_ALG2, N2_C_ALG2)
colnames(ALG2_output_C) <- c("Total Reads")
ALG2_output_C
write.table(ALG2_output_C, "Muscle_ALG2_output_C_7.txt", sep="\t", col.names = NA)
write.table(datalist_ALG2_C, "Muscle_ALG2_output_C_7_meta.txt", sep="\t", col.names = NA)
