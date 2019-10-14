######################
#Determine frameshift#
#type of isomiRs######
######################

############
#Data input#
############

meta.table <- read.delim("Reference_no_SNP_meta_2.txt", header = TRUE)
canon.list <- read.delim("Canonical miRNAs.txt", header = FALSE)$V1
canon.list <- paste0("^",canon.list,"$")

sig.genes <- read.csv("Combination vs N2 - Pos.csv", header = FALSE)$V2
miR.genes <- meta.table[grep("^5p48$", meta.table[,23]), ]$genes
sig.genes <- intersect(sig.genes, miR.genes)
sig.genes <- paste0("^",sig.genes,"$")

meta.table[paste("comb_reads")] <- ((meta.table[["A284_1"]] + meta.table[["A284_2"]]) / 2) + 
                                   ((meta.table[["A285_1"]] + meta.table[["A285_2"]]) / 2) +
                                   ((meta.table[["A286_1"]] + meta.table[["A286_2"]]) / 2) + 
                                   ((meta.table[["A287_1"]] + meta.table[["A287_2"]]) / 2) +
                                   ((meta.table[["A290_1"]] + meta.table[["A290_2"]]) / 2) + 
                                   ((meta.table[["A291_1"]] + meta.table[["A291_2"]]) / 2)



###########
#Functions#
###########

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

#################
#Data processing#
#################

canon_miRNAs <- meta.table[grep(paste(canon.list, collapse = "|"), meta.table[,1]), ]  
canon_miRNA_miRNA <- canon_miRNAs[,23]
canon_miRNA_miRNA <- paste0("^",canon_miRNA_miRNA,"$")

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
  
  sample[,27] <- ifelse(sample[,25] < canon_proximal & sample[,26] == canon_distal, sample$comb_reads, 0)
  sample[,28] <- ifelse(sample[,26] < canon_distal & sample[,25] == canon_proximal, sample$comb_reads, 0)
  sample[,29] <- ifelse(sample[,25] > canon_proximal & sample[,26] == canon_distal, sample$comb_reads, 0)
  sample[,30] <- ifelse(sample[,26] > canon_distal & sample[,25] == canon_proximal, sample$comb_reads, 0)
  sample[,31] <- ifelse(sample[,25] > canon_proximal & sample[,26] > canon_distal |
                        sample[,25] > canon_proximal & sample[,26] < canon_distal |
                        sample[,25] < canon_proximal & sample[,26] > canon_distal |
                        sample[,25] < canon_proximal & sample[,26] < canon_distal, sample$comb_reads, 0)
  sample[,32] <- ifelse(grepl("A", sample[,21]) | grepl("U", sample[,21]) |
                        grepl("G", sample[,21]) | grepl("C", sample[,21]), 1, 0) #do not change from boolean
  sample[,33] <- ifelse(sample[,25] == canon_proximal & sample[,26] == canon_distal & sample[32] == 1, sample$comb_reads, 0)
  sample[,34] <- ifelse(sample[,25] == canon_proximal & sample[,26] == canon_distal & sample[32] == 0, sample$comb_reads, 0) #sample$comb_reads

  Ext.Red[[designation]] <- sample
}

Ext.Red <- do.call(rbind.data.frame, Ext.Red)
Ext.Red <- Ext.Red[grep(paste(sig.genes, collapse = "|"), Ext.Red[,1]), ] #For previously defined list
#Ext.Red <- Ext.Red[grep("^5p71$", Ext.Red[,23]), ] #For specific miRNAs

s_ext_5 <- sum(Ext.Red$`5_ext`)
s_ext_3 <- sum(Ext.Red$`3_ext`)
s_red_5 <- sum(Ext.Red$`5_red`)
s_red_3 <- sum(Ext.Red$`3_red`)
s_mix <- sum(Ext.Red$Mixed)
s_sub <- sum(Ext.Red$Substitution)
s_canon <- sum(Ext.Red$Canonical)

###################
#Data processing 2#
###################

s_ext_5_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,27] > 0, Ext.Red$comb, 0)) #Ext.Red$comb
s_ext_5_SNP_p <- paste(round((s_ext_5_SNP/s_ext_5) * 100, digits = 2), "%", sep = "")

s_ext_3_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,28] > 0, Ext.Red$comb, 0))
s_ext_3_SNP_p <- paste(round((s_ext_3_SNP/s_ext_3) * 100, digits = 2), "%", sep = "")

s_red_5_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,29] > 0, Ext.Red$comb, 0))
s_red_5_SNP_p <- paste(round((s_red_5_SNP/s_red_5) * 100, digits = 2), "%", sep = "")

s_red_3_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,30] > 0, Ext.Red$comb, 0))
s_red_3_SNP_p <- paste(round((s_red_3_SNP/s_red_3) * 100, digits = 2), "%", sep = "")

s_mix_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,31] > 0, Ext.Red$comb, 0))
s_mix_SNP_p <- paste(round((s_mix_SNP/s_mix) * 100, digits = 2), "%", sep = "")

s_sub_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,33] > 0, Ext.Red$comb, 0))
s_sub_SNP_p <- paste(round((s_sub_SNP/s_sub) * 100, digits = 2), "%", sep = "")

s_canon_SNP <- sum(ifelse(Ext.Red[,32] == 1 & Ext.Red[,34] > 0, Ext.Red$comb, 0))
s_canon_SNP_p <- paste(round((s_canon_SNP/s_canon) * 100, digits = 2), "%", sep = "")

output1 <- rbind(s_ext_5, s_ext_3, s_red_5, s_red_3, s_mix, s_sub, s_canon)
output2 <- rbind(s_ext_5_SNP, s_ext_3_SNP, s_red_5_SNP, s_red_3_SNP, s_mix_SNP, s_sub_SNP, s_canon_SNP)
output3 <- rbind(s_ext_5_SNP_p, s_ext_3_SNP_p, s_red_5_SNP_p, s_red_3_SNP_p, s_mix_SNP_p, s_sub_SNP_p, s_canon_SNP_p)

meta_output <- cbind(output1, output2, output3)
colnames(meta_output) <- c("Total reads", "Reads with SNPs", "% Reads with SNPs")
rownames(meta_output) <- c("5' extension", "3' extention", "5' reduction", "3' reduction", "Mixture", "Substitution", "Canonical")
print(meta_output, quote = FALSE)
write.table(meta_output, "a_miR_48_output_reads_w_cutoff.txt", sep="\t", col.names = NA, quote = FALSE)

###################
#Data processing 3#
###################



####################
#Data visualisation#
####################

slices1 <- c(s_canon, s_ext_5, s_ext_3, s_red_5, s_red_3, s_mix, s_sub)
slices1 <- round(slices1)

lbls1 <- c("Canonical", "5' extension", "3' extension", "5' reduction", "3' reduction", "Mixture", "Substitution")

pie(slices1, labels = prettyNum(slices1, big.mark = ","), cex = 0.8, 
    main = "Read distribution of miRNAs and their isomiRs", col = rainbow(length(slices1)))

legend("bottomright", legend=lbls1, cex=0.7, fill = rainbow(length(slices1)))
