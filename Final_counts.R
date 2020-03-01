#Final miRNA and isomiR counts

############
#SNV counts#
############

#Data collection & processing
{
  g_threshold_isomiRs <- paste0("^", read.csv("isomiRs above 50cpm in at least 2 groups.csv", header = TRUE)$genes, "$")
  g_threshold_isomiRs_2 <- setdiff(g_threshold_isomiRs, g_canon_miRNAs)
  g_all_sig_isomiRs <- read.csv("Combination vs N2 - Pos.csv", header = FALSE)$V2
  g_all_sig_isomiRs <- paste0("^", g_all_sig_isomiRs, "$")
  g_all_sig_isomiRs_2 <- setdiff(g_all_sig_isomiRs, g_canon_miRNAs)
  g_sig_isomiRs_Int1 <- read.csv("IntALG1 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Int1 <- paste0("^", g_sig_isomiRs_Int1, "$")
  g_sig_isomiRs_Int1_2 <- setdiff(g_sig_isomiRs_Int1, g_canon_miRNAs)
  g_sig_isomiRs_Int2 <- read.csv("IntALG2 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Int2 <- paste0("^", g_sig_isomiRs_Int2, "$")
  g_sig_isomiRs_Int2_2 <- setdiff(g_sig_isomiRs_Int2, g_canon_miRNAs)
  g_sig_isomiRs_Neu1 <- read.csv("NeuALG1 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Neu1 <- paste0("^", g_sig_isomiRs_Neu1, "$")
  g_sig_isomiRs_Neu1_2 <- setdiff(g_sig_isomiRs_Neu1, g_canon_miRNAs)
  g_sig_isomiRs_Neu2 <- read.csv("NeuALG2 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Neu2 <- paste0("^", g_sig_isomiRs_Neu2, "$")
  g_sig_isomiRs_Neu2_2 <- setdiff(g_sig_isomiRs_Neu2, g_canon_miRNAs)
  g_sig_isomiRs_Mus1 <- read.csv("MusALG1 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Mus1 <- paste0("^", g_sig_isomiRs_Mus1, "$")
  g_sig_isomiRs_Mus1_2 <- setdiff(g_sig_isomiRs_Mus1, g_canon_miRNAs)
  g_sig_isomiRs_Mus2 <- read.csv("MusALG2 vs N2 Top - Pos.csv", header = FALSE)$V2
  g_sig_isomiRs_Mus2 <- paste0("^", g_sig_isomiRs_Mus2, "$")
  g_sig_isomiRs_Mus2_2 <- setdiff(g_sig_isomiRs_Mus2, g_canon_miRNAs)
  g_sig_isomiRs_Int <- union(g_sig_isomiRs_Int1, g_sig_isomiRs_Int2)
  g_sig_isomiRs_Neu <- union(g_sig_isomiRs_Neu1, g_sig_isomiRs_Neu2)
  g_sig_isomiRs_Mus <- union(g_sig_isomiRs_Mus1, g_sig_isomiRs_Mus2)
  g_sig_isomiRs_ALG1 <- union(union(g_sig_isomiRs_Int1, g_sig_isomiRs_Neu1), g_sig_isomiRs_Mus1) 
  g_sig_isomiRs_ALG2 <- union(union(g_sig_isomiRs_Int2, g_sig_isomiRs_Neu2), g_sig_isomiRs_Mus2)
}

fileNames <- Sys.glob("*.txt")

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

datalist_x_A_meta = list()
datalist_x_T_meta = list()
datalist_x_G_meta = list()
datalist_x_C_meta = list()

#SNVs N > A/T/G/C
{
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

  #Int ALG1###
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Int ALG2###
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Neu ALG1###
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Neu ALG2###
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Mus ALG1###
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Mus ALG2#
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), 1, 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), 1, 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Int#
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,33], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,33], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,33], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,33], 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,33], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,33], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,33], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,33], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,33], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,33], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,33], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,33], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,33], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,33], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,33], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,33], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Neu#
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,34], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,34], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,34], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,34], 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,34], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,34], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,34], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,34], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,34], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,34], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,34], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,34], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,34], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,34], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,34], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,34], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #Mus#
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,35], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,35], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,35], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,35], 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,35], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,35], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,35], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,35], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,35], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,35], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,35], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,35], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,35], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,35], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,35], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,35], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #ALG1#
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,37], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,37], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,37], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,37], 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,37], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,37], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,37], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,37], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,37], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,37], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,37], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,37], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,37], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,37], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,37], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,37], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}

  #ALG2
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,38], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,38], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,38], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("A", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("A", sample[,22]), sample[,38], 0)
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
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,38], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,38], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,38], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,38], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("T", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_T_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,38], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,38], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,38], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,38], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("G", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_G_meta[[FileName]] <- sample
  
}
  for (FileName in fileNames) {
  sample <- read.delim(FileName, stringsAsFactors = FALSE)
  sample <- as.data.frame(sample)
  #sample$Pos <- ifelse(as.numeric(sample[,24]) > 10, sample[,23], 0)
  sample$SNP_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,38], 0)
  sample$N2_A <- ifelse(grepl("A", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,38], 0)
  sample$N2_T <- ifelse(grepl("T", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,38], 0)
  sample$N2_G <- ifelse(grepl("G", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  sample$SNP_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,38], 0)
  sample$N2_C <- ifelse(grepl("C", sample[,21]) & grepl("C", sample[,22]), sample[,36], 0)
  
  sample$FileName <- FileName
  datalist_x_C_meta[[FileName]] <- sample
  
}
}

datalist_x_A_meta <- do.call(rbind.data.frame, datalist_x_A_meta)
datalist_x_T_meta <- do.call(rbind.data.frame, datalist_x_T_meta)
datalist_x_G_meta <- do.call(rbind.data.frame, datalist_x_G_meta)
datalist_x_C_meta <- do.call(rbind.data.frame, datalist_x_C_meta)

#Shorten list by significant isomiRs/Change g_sig_isomiRs by what is required
datalist_x_A_meta <- datalist_x_A_meta[grep(paste(g_sig_isomiRs_Mus2, collapse = "|"), datalist_x_A_meta[,1]), ]
datalist_x_T_meta <- datalist_x_T_meta[grep(paste(g_sig_isomiRs_Mus2, collapse = "|"), datalist_x_T_meta[,1]), ]
datalist_x_G_meta <- datalist_x_G_meta[grep(paste(g_sig_isomiRs_Mus2, collapse = "|"), datalist_x_G_meta[,1]), ]
datalist_x_C_meta <- datalist_x_C_meta[grep(paste(g_sig_isomiRs_Mus2, collapse = "|"), datalist_x_C_meta[,1]), ]


#SNV counts
{ SNP_A_A <- sum(datalist_x_A_meta$SNP_A)
  N2_A_A <- sum(datalist_x_A_meta$N2_A)
  SNP_T_A <- sum(datalist_x_A_meta$SNP_T)
  N2_T_A <- sum(datalist_x_A_meta$N2_T)
  SNP_G_A <- sum(datalist_x_A_meta$SNP_G)
  N2_G_A <- sum(datalist_x_A_meta$N2_G)
  SNP_C_A <- sum(datalist_x_A_meta$SNP_C)
  N2_C_A <- sum(datalist_x_A_meta$N2_C)
  SNP_A_T <- sum(datalist_x_T_meta$SNP_A)
  N2_A_T <- sum(datalist_x_T_meta$N2_A)
  SNP_T_T <- sum(datalist_x_T_meta$SNP_T)
  N2_T_T <- sum(datalist_x_T_meta$N2_T)
  SNP_G_T <- sum(datalist_x_T_meta$SNP_G)
  N2_G_T <- sum(datalist_x_T_meta$N2_G)
  SNP_C_T <- sum(datalist_x_T_meta$SNP_C)
  N2_C_T <- sum(datalist_x_T_meta$N2_C)
  SNP_A_G <- sum(datalist_x_G_meta$SNP_A)
  N2_A_G <- sum(datalist_x_G_meta$N2_A)
  SNP_T_G <- sum(datalist_x_G_meta$SNP_T)
  N2_T_G <- sum(datalist_x_G_meta$N2_T)
  SNP_G_G <- sum(datalist_x_G_meta$SNP_G)
  N2_G_G <- sum(datalist_x_G_meta$N2_G)
  SNP_C_G <- sum(datalist_x_G_meta$SNP_C)
  N2_C_G <- sum(datalist_x_G_meta$N2_C)
  SNP_A_C <- sum(datalist_x_C_meta$SNP_A)
  N2_A_C <- sum(datalist_x_C_meta$N2_A)
  SNP_T_C <- sum(datalist_x_C_meta$SNP_T)
  N2_T_C <- sum(datalist_x_C_meta$N2_T)
  SNP_G_C <- sum(datalist_x_C_meta$SNP_G)
  N2_G_C <- sum(datalist_x_C_meta$N2_G)
  SNP_C_C <- sum(datalist_x_C_meta$SNP_C)
  N2_C_C <- sum(datalist_x_C_meta$N2_C)
}

#Data presentation
{ 
  total_SNP <- sum(SNP_A_A, SNP_T_A, SNP_G_A, SNP_C_A, SNP_A_T, SNP_T_T, SNP_G_T, SNP_C_T,
                   SNP_A_G, SNP_T_G, SNP_G_G, SNP_C_G, SNP_A_C, SNP_T_C, SNP_G_C, SNP_C_C)
  total_N2 <- sum(N2_A_A, N2_T_A, N2_G_A, N2_C_A, N2_A_T, N2_T_T, N2_G_T, N2_C_T,
                  N2_A_G, N2_T_G, N2_G_G, N2_C_G, N2_A_C, N2_T_C, N2_G_C, N2_C_C)

  meta_output_1 <- prettyNum(round(rbind(SNP_A_A, SNP_T_A, SNP_G_A, SNP_C_A, SNP_A_T, SNP_T_T, SNP_G_T, SNP_C_T,
                                       SNP_A_G, SNP_T_G, SNP_G_G, SNP_C_G, SNP_A_C, SNP_T_C, SNP_G_C, SNP_C_C,
                                       total_SNP)),  big.mark = ",")

  meta_output_2 <- prettyNum(round(rbind(N2_A_A, N2_T_A, N2_G_A, N2_C_A, N2_A_T, N2_T_T, N2_G_T, N2_C_T,
                                       N2_A_G, N2_T_G, N2_G_G, N2_C_G, N2_A_C, N2_T_C, N2_G_C, N2_C_C,
                                       total_N2)), big.mark = ",")
  
  meta_output_3 <- cbind(meta_output_1, meta_output_2)
  colnames(meta_output_3) <- c("SNP reads", "N2 reads")
  rownames(meta_output_3) <- c("A>A", "U>A", "G>A", "C>A", "A>U", "U>U", "G>U", "C>U",
                             "A>G", "U>G", "G>C", "C>G", "A>C", "U>C", "G>C", "C>C",
                             "Poly-A", "Poly-U", "Total")
  print(meta_output_3, quote=FALSE)
}

########################
#Canonical miRNA counts#
########################

#Data collection and processing
{
  g_canon_miRNAs <- paste0("^", read.csv("Canonical miRNAs.csv", header = TRUE)$genes, "$")
  g_canon_IntALG1 <- paste0("^", read.csv("Intestine ALG1.csv", header = TRUE)$seq, "$")
  g_canon_IntALG2 <- paste0("^", read.csv("Intestine ALG2.csv", header = TRUE)$seq, "$")
  g_canon_NeuALG1 <- paste0("^", read.csv("Neuron ALG1.csv", header = TRUE)$seq, "$")
  g_canon_NeuALG2 <- paste0("^", read.csv("Neuron ALG2.csv", header = TRUE)$seq, "$")
  g_canon_MusALG1 <- paste0("^", read.csv("Muscle ALG1.csv", header = TRUE)$seq, "$")
  g_canon_MusALG2 <- paste0("^", read.csv("Muscle ALG2.csv", header = TRUE)$seq, "$")
  g_canon_Int <- paste0("^", read.csv("Intestine.csv", header = TRUE)$seq, "$")
  g_canon_Neu <- paste0("^", read.csv("Neuron.csv", header = TRUE)$seq, "$")
  g_canon_Mus <- paste0("^", read.csv("Muscle.csv", header = TRUE)$seq, "$")
  g_canon_ALG1 <- paste0("^", read.csv("ALG1.csv", header = TRUE)$seq, "$")
  g_canon_ALG2 <- paste0("^", read.csv("ALG2.csv", header = TRUE)$seq, "$")
  g_canon_significant <- paste0("^", read.csv("All miRNAs.csv", header = TRUE)$seq, "$")
}

d_all_isomiRs <- read.delim("Reference_no_SNP_meta_2.txt", header = TRUE)
{ 
  d_all_isomiRs[paste('av_284')] <- (d_all_isomiRs[["A284_1"]] + d_all_isomiRs[["A284_2"]])/2
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

#Shorten list by significant miRNAs/Change g_canon_miRNAs by what is required
d_canon <- d_all_isomiRs[grep(paste(g_canon_NeuALG2, collapse = "|"), d_all_isomiRs$genes), ]

#miRNA counts
{
  r_canon_miRNAs <- sum(d_canon$comb_av_total)
  r_canon_IntALG1 <- sum(d_canon$av_284)
  r_canon_IntALG2 <- sum(d_canon$av_285)
  r_canon_NeuALG1 <- sum(d_canon$av_286)
  r_canon_NeuALG2 <- sum(d_canon$av_287)
  r_canon_MusALG1 <- sum(d_canon$av_290)
  r_canon_MusALG2 <- sum(d_canon$av_291)
  r_canon_Int <- sum(d_canon$Intcomb)
  r_canon_Neu <- sum(d_canon$Neucomb)
  r_canon_Mus <- sum(d_canon$Muscomb)
  r_canon_ALG1 <- sum(d_canon$av_ALG1)
  r_canon_ALG2 <- sum(d_canon$av_ALG2)
}

##########################
#Frameshift isomiR counts#
##########################

#Data collection and processing
{
  meta.table <- read.delim("Reference_no_SNP_meta_2.txt", header = TRUE)
  canon.list <- read.delim("Canonical miRNAs.txt", header = FALSE)$V1
  canon.list <- paste0("^",canon.list,"$")
  meta.table[paste("comb_reads")] <- ((meta.table[["A284_1"]] + meta.table[["A284_2"]]) / 2) + 
    ((meta.table[["A285_1"]] + meta.table[["A285_2"]]) / 2) +
    ((meta.table[["A286_1"]] + meta.table[["A286_2"]]) / 2) + 
    ((meta.table[["A287_1"]] + meta.table[["A287_2"]]) / 2) +
    ((meta.table[["A290_1"]] + meta.table[["A290_2"]]) / 2) + 
    ((meta.table[["A291_1"]] + meta.table[["A291_2"]]) / 2)
  meta.table[paste("av284")] <- ((meta.table[["A284_1"]] + meta.table[["A284_2"]]) / 2)
  meta.table[paste("av285")] <- ((meta.table[["A285_1"]] + meta.table[["A285_2"]]) / 2)
  meta.table[paste("av286")] <- ((meta.table[["A286_1"]] + meta.table[["A286_2"]]) / 2)
  meta.table[paste("av287")] <- ((meta.table[["A287_1"]] + meta.table[["A287_2"]]) / 2)
  meta.table[paste("av290")] <- ((meta.table[["A290_1"]] + meta.table[["A290_2"]]) / 2)
  meta.table[paste("av291")] <- ((meta.table[["A291_1"]] + meta.table[["A291_2"]]) / 2)
  
  canon_miRNAs <- meta.table[grep(paste(canon.list, collapse = "|"), meta.table[,1]), ]  
  canon_miRNA_miRNA <- canon_miRNAs[,23]
  canon_miRNA_miRNA <- paste0("^",canon_miRNA_miRNA,"$")
}

#Functions
h_count_p <- function(x){
  x_2 <- gsub("\\*", "Z", x)  #x = canon_hairpin
  x_3 <- gsub("Z", "", x_2)
  a <- nchar(x_3)
  b <- strsplit(x_2, "Z*|Z")[[1]][a+1]
  c <- paste0(b, "Z", "*")
  d <- gsub(c, "", x_2)
  str_count(d, "Z")
}
h_count_d <- function(x){
  x_2 <- gsub("\\*", "Z", x)  #x = canon_hairpin
  e <- strsplit(x_2, "Z*|Z")[[1]][2]
  f <- paste0(".", "*", "Z", e, sep = "")
  g <- gsub(f, "", x_2)
  str_count(g, "Z")
}

Ext.Red = list()
for (designation in canon_miRNA_miRNA) {
  canon <- canon_miRNAs[grep(designation, canon_miRNAs[,23]), ]
  canon_hairpin <- canon$h_seq
  canon_proximal <- h_count_p(canon_hairpin)
  canon_distal <- h_count_d(canon_hairpin)
  
  sample <- meta.table[grep(designation, meta.table[,23]), ]
  sample[paste('h_proximal')] <- 0
  sample[paste('h_distal')] <- 0
  sample[paste('5_ext')] <- 0
  sample[paste('3_ext')] <- 0
  sample[paste('5_red')] <- 0
  sample[paste('3_red')] <- 0
  sample[paste('Mixed')] <- 0
  sample[paste('SNP')] <- 0
  sample[paste('Substitution')] <- 0
  sample[paste('Canonical')] <- 0
  
  sample$h_proximal <- sapply(sample$h_seq, h_count_p)
  sample$h_distal <- sapply(sample$h_seq, h_count_d)
  
  sample[,33] <- ifelse(sample[,31] < canon_proximal & sample[,32] == canon_distal, 1, 0)
  sample[,34] <- ifelse(sample[,32] < canon_distal & sample[,31] == canon_proximal, 1, 0)
  sample[,35] <- ifelse(sample[,31] > canon_proximal & sample[,32] == canon_distal, 1, 0)
  sample[,36] <- ifelse(sample[,32] > canon_distal & sample[,31] == canon_proximal, 1, 0)
  sample[,37] <- ifelse(sample[,31] > canon_proximal & sample[,32] > canon_distal |
                          sample[,31] > canon_proximal & sample[,32] < canon_distal |
                          sample[,31] < canon_proximal & sample[,32] > canon_distal |
                          sample[,31] < canon_proximal & sample[,32] < canon_distal, 1, 0)
  sample[,38] <- ifelse(grepl("A", sample[,21]) | grepl("U", sample[,21]) |
                          grepl("G", sample[,21]) | grepl("C", sample[,21]), 1, 0) #do not change from boolean
  sample[,39] <- ifelse(sample[,31] == canon_proximal & sample[,32] == canon_distal & sample[38] == 1, 1, 0)
  sample[,40] <- ifelse(sample[,31] == canon_proximal & sample[,32] == canon_distal & sample[38] == 0, 1, 0) #sample$comb_reads
  
  Ext.Red[[designation]] <- sample
}
Ext.Red <- do.call(rbind.data.frame, Ext.Red)

#Shorten list by significant miRNAs/Change g_threshold_isomiRs by what is required
{
  a <- Ext.Red[grep(paste(g_all_sig_isomiRs, collapse = "|"), Ext.Red$genes), ]
  a2 <- Ext.Red[grep(paste(g_all_sig_isomiRs_2, collapse = "|"), Ext.Red$genes), ]
  b <- Ext.Red[grep(paste(g_threshold_isomiRs, collapse = "|"), Ext.Red$genes), ]
  b2 <- Ext.Red[grep(paste(g_threshold_isomiRs_2, collapse = "|"), Ext.Red$genes), ]
  c <- Ext.Red[grep(paste(g_sig_isomiRs_Int1, collapse = "|"), Ext.Red$genes), ]
  c2 <- Ext.Red[grep(paste(g_sig_isomiRs_Int1_2, collapse = "|"), Ext.Red$genes), ]
  d <- Ext.Red[grep(paste(g_sig_isomiRs_Int2, collapse = "|"), Ext.Red$genes), ]
  d2 <- Ext.Red[grep(paste(g_sig_isomiRs_Int2_2, collapse = "|"), Ext.Red$genes), ]
  e <- Ext.Red[grep(paste(g_sig_isomiRs_Neu1, collapse = "|"), Ext.Red$genes), ]
  e2 <- Ext.Red[grep(paste(g_sig_isomiRs_Neu1_2, collapse = "|"), Ext.Red$genes), ]
  f <- Ext.Red[grep(paste(g_sig_isomiRs_Neu2, collapse = "|"), Ext.Red$genes), ]
  f2 <- Ext.Red[grep(paste(g_sig_isomiRs_Neu2_2, collapse = "|"), Ext.Red$genes), ]
  g <- Ext.Red[grep(paste(g_sig_isomiRs_Mus1, collapse = "|"), Ext.Red$genes), ]
  g2 <- Ext.Red[grep(paste(g_sig_isomiRs_Mus1_2, collapse = "|"), Ext.Red$genes), ]
  h <- Ext.Red[grep(paste(g_sig_isomiRs_Mus2, collapse = "|"), Ext.Red$genes), ]
  h2 <- Ext.Red[grep(paste(g_sig_isomiRs_Mus2_2, collapse = "|"), Ext.Red$genes), ]
}

#isomiR counts
{
  r_ext_5 <- sum(a2$`5_ext`)
  r_ext_3 <- sum(a2$`3_ext`)
  r_red_5 <- sum(a2$`5_red`)
  r_red_3 <- sum(a2$`3_red`)
  r_mix <- sum(a2$Mixed)
  r_sub <- sum(g2$Substitution)
  r_canon <- sum(a2$Canonical)
  r_isomiRs <- sum(a2$av291)
}



###############
#isomiR counts#
###############

g_sig_isomiRs_Int1_2 <- setdiff(g_sig_isomiRs_Int1, g_canon_miRNAs)
g_sig_isomiRs_Int2_2 <- setdiff(g_sig_isomiRs_Int2, g_canon_miRNAs)
g_sig_isomiRs_Neu1_2 <- setdiff(g_sig_isomiRs_Neu1, g_canon_miRNAs)
g_sig_isomiRs_Neu2_2 <- setdiff(g_sig_isomiRs_Neu2, g_canon_miRNAs)
g_sig_isomiRs_Mus1_2 <- setdiff(g_sig_isomiRs_Mus1, g_canon_miRNAs)
g_sig_isomiRs_Mus2_2 <- setdiff(g_sig_isomiRs_Mus2, g_canon_miRNAs)

d_threshold_isomiRs <- d_all_isomiRs[grep(paste(g_threshold_isomiRs_2, collapse = "|"), d_all_isomiRs$genes), ]
d_sig_isomiRs <- d_all_isomiRs[grep(paste(g_all_sig_isomiRs_2, collapse = "|"), d_all_isomiRs$genes), ]
d_IntALG1 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Int1_2, collapse = "|"), d_all_isomiRs$genes), ]
d_IntALG2 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Int2_2, collapse = "|"), d_all_isomiRs$genes), ]
d_NeuALG1 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Neu1_2, collapse = "|"), d_all_isomiRs$genes), ]
d_NeuALG2 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Neu2_2, collapse = "|"), d_all_isomiRs$genes), ]
d_MusALG1 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Mus1_2, collapse = "|"), d_all_isomiRs$genes), ]
d_MusALG2 <- d_all_isomiRs[grep(paste(g_sig_isomiRs_Mus2_2, collapse = "|"), d_all_isomiRs$genes), ]



d_IntALG2 <- d_IntALG2[ ! d_IntALG2$genes %in% g_canon_miRNAs, ]


r_sig_isomiRs <- sum(d_sig_isomiRs$comb_av_total)
r_isomiRs <- sum(d_threshold_isomiRs$comb_av_total)
r_IntALG1 <- sum(d_IntALG1$av_284)
r_IntALG2 <- sum(d_IntALG2$av_285)
r_NeuALG1 <- sum(d_NeuALG1$av_286)
r_NeuALG2 <- sum(d_NeuALG2$av_286)
r_MusALG1 <- sum(d_MusALG1$av_290)
r_MusALG2 <- sum(d_MusALG2$av_291)

