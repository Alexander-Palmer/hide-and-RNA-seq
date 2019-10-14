##############################
##############################

##A code to find if I've been wasting my time for the past 7 months

##############################
##############################

sig.genes <- read.csv("Combination vs N2 - Pos.csv", header = FALSE)$V2
sig.genes <- paste0("^", sig.genes, "$")

fileNames <- Sys.glob("*.txt")

fileNames <- as.character(c("Ref_A3p82.txt", "Ref_U3p82.txt", "Ref_G3p82.txt", "Ref_C3p82.txt",
                            "Ref_A5p82.txt", "Ref_U5p82.txt", "Ref_G5p82.txt", "Ref_C5p82.txt"))

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

  sample[paste('comb_av_total')] <- (sample[["av_284"]] + sample[["av_285"]] + sample[["av_286"]]
                                     + sample[["av_287"]] + sample[["av_290"]] + sample[["av_291"]])
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

##############################
##############################

##############################################################################################

##############################
##############################

#SNPs x > A (meta)

datalist_x_A_meta = list()

for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,39], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,39], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,39], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,39], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)

  char <- nchar(as.character(sample[,1]))

  sample$polyA <- ifelse((grepl("A", sample[,22]) & sample[,23] == 1) | (grepl("A", sample[,22]) & sample[,23] == char) , sample[,1], 0)
  sample$polyA_r <- ifelse((grepl("A", sample[,22]) & sample[,23] == 1) | (grepl("A", sample[,22]) & sample[,23] == char) , sample[,39], 0)
  sample$N2_polyA_r <- ifelse((grepl("A", sample[,22]) & sample[,23] == 1) | (grepl("A", sample[,22]) & sample[,23] == char) , sample[,36], 0)
  sample$polyU <- ifelse((grepl("T", sample[,22]) & sample[,23] == 1) | (grepl("T", sample[,22]) & sample[,23] == char) , sample[,1], 0)
  sample$polyU_r <- ifelse((grepl("T", sample[,22]) & sample[,23] == 1) | (grepl("T", sample[,22]) & sample[,23] == char) , sample[,39], 0)
  sample$N2_polyU_r <- ifelse((grepl("T", sample[,22]) & sample[,23] == 1) | (grepl("T", sample[,22]) & sample[,23] == char) , sample[,36], 0)


  sample$FileName <- FileName
  datalist_x_A_meta[[FileName]] <- sample

}

datalist_x_A_meta <- do.call(rbind.data.frame, datalist_x_A_meta)
datalist_x_A_meta <- datalist_x_A_meta[grep(paste(sig.genes, collapse = "|"), datalist_x_A_meta[,1]), ]

SNP_A_A <- sum(datalist_x_A_meta$SNP_A)
N2_A_A <- sum(datalist_x_A_meta$N2_A)
SNP_T_A <- sum(datalist_x_A_meta$SNP_T)
N2_T_A <- sum(datalist_x_A_meta$N2_T)
SNP_G_A <- sum(datalist_x_A_meta$SNP_G)
N2_G_A <- sum(datalist_x_A_meta$N2_G)
SNP_C_A <- sum(datalist_x_A_meta$SNP_C)
N2_C_A <- sum(datalist_x_A_meta$N2_C)
polyA <- sum(datalist_x_A_meta$polyA_r)
N2_polyA <- sum(datalist_x_A_meta$N2_polyA_r)
polyU <- sum(datalist_x_A_meta$polyU_r)
N2_polyU <- sum(datalist_x_A_meta$N2_polyU_r)


##############################################################################################

#SNPs x > T (meta)

datalist_x_T_meta = list()

for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,39], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,39], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,39], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,39], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)

  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample

}

datalist_x_T_meta <- do.call(rbind.data.frame, datalist_x_T_meta)
datalist_x_T_meta <- datalist_x_T_meta[grep(paste(sig.genes, collapse = "|"), datalist_x_T_meta[,1]), ]

SNP_A_T <- sum(datalist_x_T_meta$SNP_A)
N2_A_T <- sum(datalist_x_T_meta$N2_A)
SNP_T_T <- sum(datalist_x_T_meta$SNP_T)
N2_T_T <- sum(datalist_x_T_meta$N2_T)
SNP_G_T <- sum(datalist_x_T_meta$SNP_G)
N2_G_T <- sum(datalist_x_T_meta$N2_G)
SNP_C_T <- sum(datalist_x_T_meta$SNP_C)
N2_C_T <- sum(datalist_x_T_meta$N2_C)

##############################################################################################

#SNPs x > G (meta)

datalist_x_G_meta = list()

for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,39], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,39], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,39], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,39], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)

  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample

}

datalist_x_G_meta <- do.call(rbind.data.frame, datalist_x_G_meta)
datalist_x_G_meta <- datalist_x_G_meta[grep(paste(sig.genes, collapse = "|"), datalist_x_G_meta[,1]), ]

SNP_A_G <- sum(datalist_x_G_meta$SNP_A)
N2_A_G <- sum(datalist_x_G_meta$N2_A)
SNP_T_G <- sum(datalist_x_G_meta$SNP_T)
N2_T_G <- sum(datalist_x_G_meta$N2_T)
SNP_G_G <- sum(datalist_x_G_meta$SNP_G)
N2_G_G <- sum(datalist_x_G_meta$N2_G)
SNP_C_G <- sum(datalist_x_G_meta$SNP_C)
N2_C_G <- sum(datalist_x_G_meta$N2_C)

##############################################################################################

#SNPs x > C (meta)

datalist_x_C_meta = list()

for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,39], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,39], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,39], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,39], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)

  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample

}

datalist_x_C_meta <- do.call(rbind.data.frame, datalist_x_C_meta)
datalist_x_C_meta <- datalist_x_C_meta[grep(paste(sig.genes, collapse = "|"), datalist_x_C_meta[,1]), ]

SNP_A_C <- sum(datalist_x_C_meta$SNP_A)
N2_A_C <- sum(datalist_x_C_meta$N2_A)
SNP_T_C <- sum(datalist_x_C_meta$SNP_T)
N2_T_C <- sum(datalist_x_C_meta$N2_T)
SNP_G_C <- sum(datalist_x_C_meta$SNP_G)
N2_G_C <- sum(datalist_x_C_meta$N2_G)
SNP_C_C <- sum(datalist_x_C_meta$SNP_C)
N2_C_C <- sum(datalist_x_C_meta$N2_C)

##############################
##############################

total_SNP <- sum(SNP_A_A, SNP_T_A, SNP_G_A, SNP_C_A, SNP_A_T, SNP_T_T, SNP_G_T, SNP_C_T,
                 SNP_A_G, SNP_T_G, SNP_G_G, SNP_C_G, SNP_A_C, SNP_T_C, SNP_G_C, SNP_C_C,
                 polyA, polyU)
total_N2 <- sum(N2_A_A, N2_T_A, N2_G_A, N2_C_A, N2_A_T, N2_T_T, N2_G_T, N2_C_T,
                N2_A_G, N2_T_G, N2_G_G, N2_C_G, N2_A_C, N2_T_C, N2_G_C, N2_C_C,
                N2_polyA, N2_polyU)

meta_output_1 <- prettyNum(round(rbind(SNP_A_A, SNP_T_A, SNP_G_A, SNP_C_A, SNP_A_T, SNP_T_T, SNP_G_T, SNP_C_T,
                       SNP_A_G, SNP_T_G, SNP_G_G, SNP_C_G, SNP_A_C, SNP_T_C, SNP_G_C, SNP_C_C,
                       polyA, polyU, total_SNP)),  big.mark = ",")

meta_output_2 <- prettyNum(round(rbind(N2_A_A, N2_T_A, N2_G_A, N2_C_A, N2_A_T, N2_T_T, N2_G_T, N2_C_T,
                       N2_A_G, N2_T_G, N2_G_G, N2_C_G, N2_A_C, N2_T_C, N2_G_C, N2_C_C,
                       N2_polyA, N2_polyU, total_N2)), big.mark = ",")

meta_output_3 <- cbind(meta_output_1, meta_output_2)
colnames(meta_output_3) <- c("SNP reads", "N2 reads")
rownames(meta_output_3) <- c("A>A", "U>A", "G>A", "C>A", "A>U", "U>U", "G>U", "C>U",
                             "A>G", "U>G", "G>C", "C>G", "A>C", "U>C", "G>C", "C>C",
                             "Poly-A", "Poly-U", "Total")
print(meta_output_3, quote=FALSE)

r_total_polyAU <- sum(polyA, polyU)
r_total_SNP2 <- total_SNP - r_total_polyAU

g_SNP_miRNAs <- datalist_x_A_meta$genes
g_SNP_miRNAs <- paste0("^",g_SNP_miRNAs,"$")

#write.table(ALG2_output_T, "Intestine_ALG1_output_A_7.txt", sep="\t", col.names = NA)
#write.table(datalist_ALG2_T, "Intestine_ALG1_output_A_7_meta.txt", sep="\t", col.names = NA)

##############################
##############################

##############################################################################################

##############################
##############################

#Extract read count of canonical miRNA

g_canon_miRNAs <- read.delim("Canonical miRNAs.txt", header = FALSE)$V1
g_canon_miRNAs <- c("TGGAATGTAAAGAAGTATGTA", "CATACTTCCTTACATGCCCATA")
g_canon_miRNAs <- paste0("^",g_canon_miRNAs,"$")

d_all_isomiRs <- read.delim("Reference_no_SNP_meta_2.txt", header = TRUE)
all_isomiRs_1 <- d_all_isomiRs[grep("^3p48$", d_all_isomiRs[,23]), ]
all_isomiRs_2 <- d_all_isomiRs[grep("^5p48$", d_all_isomiRs[,23]), ]
d_all_isomiRs <- rbind(all_isomiRs_1, all_isomiRs_2)


{   d_all_isomiRs[paste('av_284')] <- (d_all_isomiRs[["A284_1"]] + d_all_isomiRs[["A284_2"]])/2
    d_all_isomiRs[paste('av_285')] <- (d_all_isomiRs[["A285_1"]] + d_all_isomiRs[["A285_2"]])/2
    d_all_isomiRs[paste('av_286')] <- (d_all_isomiRs[["A286_1"]] + d_all_isomiRs[["A286_2"]])/2
    d_all_isomiRs[paste('av_287')] <- (d_all_isomiRs[["A287_1"]] + d_all_isomiRs[["A287_2"]])/2
    d_all_isomiRs[paste('av_290')] <- (d_all_isomiRs[["A290_1"]] + d_all_isomiRs[["A290_2"]])/2
    d_all_isomiRs[paste('av_291')] <- (d_all_isomiRs[["A291_1"]] + d_all_isomiRs[["A291_2"]])/2
    d_all_isomiRs[paste('Intcomb')] <- (d_all_isomiRs[["av_284"]] + d_all_isomiRs[["av_285"]])/2
    d_all_isomiRs[paste('Neucomb')] <- (d_all_isomiRs[["av_286"]] + d_all_isomiRs[["av_287"]])/2
    d_all_isomiRs[paste('Muscomb')] <- (d_all_isomiRs[["av_290"]] + d_all_isomiRs[["av_291"]])/2
    d_all_isomiRs[paste('av_N2')] <- (d_all_isomiRs[["N2_1"]] + d_all_isomiRs[["N2_2"]])/2
    d_all_isomiRs[paste('av_ALG1')] <- (d_all_isomiRs[["av_284"]] + d_all_isomiRs[["av_286"]] + d_all_isomiRs[["av_290"]])/3
    d_all_isomiRs[paste('av_ALG2')] <- (d_all_isomiRs[["av_285"]] + d_all_isomiRs[["av_287"]] + d_all_isomiRs[["av_291"]])/3

    d_all_isomiRs[paste('comb_av_total')] <- (d_all_isomiRs[["av_284"]] + d_all_isomiRs[["av_285"]] + d_all_isomiRs[["av_286"]]
                                         + d_all_isomiRs[["av_287"]] + d_all_isomiRs[["av_290"]] + d_all_isomiRs[["av_291"]])
}

d_canon <- d_all_isomiRs[grep(paste(g_canon_miRNAs, collapse = "|"), d_all_isomiRs[,1]), ]
r_canon_miRNAs <- sum(d_canon$comb_av_total)

##############################
##############################

##############################################################################################

##############################
##############################

#Extract read count of frameshift isomiRs (everything else)

g_no_frameshift_miRNAs <- append(g_SNP_miRNAs, g_canon_miRNAs)

outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))
}

g_all_isomiRs_genes <- d_all_isomiRs$genes
g_all_isomiRs_genes <- paste0("^",g_all_isomiRs_genes,"$")
g_frameshift_miRNAs <- outersect(g_no_frameshift_miRNAs, g_all_isomiRs_genes)
d_frameshift <- d_all_isomiRs[grep(paste(g_frameshift_miRNAs, collapse = "|"), d_all_isomiRs[,1]), ]

r_frameshift_miRNAs <- sum(d_frameshift$comb_av_total)

rm(all_isomiRs_1, all_isomiRs_2, sample, char, fileName, FileName, fileNames, meta_output_1, meta_output_2)

##############################
##############################

##############################################################################################

##############################
##############################

#Pie charts

slices1 <- c(SNP_T_A, SNP_G_A, SNP_C_A, SNP_A_T, SNP_G_T, SNP_C_T, SNP_A_G, SNP_T_G,
             SNP_C_G, SNP_A_C, SNP_T_C, SNP_G_C)
slices1 <- round(slices1)

lbls1 <- c("T > A", "G > A", "C > A", "A > T", "G > T", "C > T", "A > G", "T > G",
           "C > G", "A > C", "T > C", "G > C")

pie(slices1, labels = prettyNum(slices1, big.mark = ","), main = "Read distribution of different SNPs", 
    col = rainbow(length(slices1)), cex = 0.7)
legend("bottomright", legend=lbls1, cex=0.7, fill = rainbow(length(slices1)))

##############################################################################################

slices1 <- c(r_frameshift_miRNAs, r_canon_miRNAs, r_total_SNP2)
slices1 <- round(slices1)

lbls1 <- c("Frameshift", "Canonical", "Single nucleotide edits")

pie(slices1, labels = prettyNum(slices1, big.mark = ","), main = "Read distribution of miRNAs and their isomiRs", col = rainbow(length(slices1)))
legend("bottomright", legend=lbls1, cex=0.7, fill = rainbow(length(slices1)))


##############################
##############################

##############################################################################################

##############################
##############################

#Output

write.table(d_canon, file="Canonical miRNAs (miR-48).csv", sep = ",", col.names = NA, qmethod = "double")
write.table(datalist_x_A_meta, file="SNP miRNAs (miR-48).csv", sep = ",", col.names = NA, qmethod = "double")
write.table(d_frameshift, file="Frameshift miRNAs (miR-48).csv", sep = ",", col.names = NA, qmethod = "double")
