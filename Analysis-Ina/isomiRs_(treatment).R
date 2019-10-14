source("http://bioconductor.org/biocLite.R")
biocLite("matrixTests")
library("matrixTests")

options(max.print=1000000)

fileNames <- Sys.glob("*.txt")

temp = list.files(pattern="*.txt")
myfiles = lapply(temp, read.delim)
for (i in 1:length(temp)) assign(temp[i], read.delim(temp[i]))
drop_col <- c("Av.1", "Av.2", "Rfam", "GtRNAdb", "chr", "start_pos", "end_pos")
#for (i in 1:length(temp)) assign(temp[i], temp[])

data.normalisation <- function(Z) {
  near_output <- (data.matrix(Z))
  output <- prop.table(near_output, 2)
  return(output)
}  

for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample[paste('av_AZ')] <- (sample[["AZ_rep1_S19"]] + sample[["AZ_rep2_S20"]] + sample[["AZ_rep3_S21"]])/3
  sample[paste('av_EB')] <- (sample[["EB_rep1_S22"]] + sample[["EB_rep2_S23"]] + sample[["EB_rep3_S24"]])/3
  sample[paste('av_isp.1')] <- (sample[["isp.1_mutant_rep1_S13"]] + sample[["isp.1_mutant_rep2_S14"]] + sample[["isp.1_mutant_rep3_S15"]])/3
  sample[paste('av_Int')] <- (sample[["SJZ_115.rep1_S10"]] + sample[["SJZ_115.rep2_S11"]] + sample[["SJZ_115.rep3_S12"]])/3
  sample[paste('av_Neu')] <- (sample[["SJZ_114_rep1_S7"]] + sample[["SJZ_114_rep2_S8"]] + sample[["SJZ_114_rep3_S9"]])/3
  sample[paste('av_N2')] <- (sample[["N2_rep1_S1"]] + sample[["N2_rep2_S2"]] + sample[["N2_rep3_S3"]])/3
  sample[paste('av_Mus.inac')] <- (sample[["SJZ_65_rep1_S16"]] + sample[["SJZ_65_rep2_S17"]] + sample[["SJZ_65_rep3_S18"]])/3
  sample[paste('av_Mus')] <- (sample[["SJZ_8_rep1_S4"]] + sample[["SJZ_8_rep2_S5"]] + sample[["SJZ_8_rep3_S6"]])/3
  
  write.table(sample, 
              fileName,
              append = FALSE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)  
}

for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.normalisation(sample)
    write.table(sample, 
              fileName,
              append = FALSE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)  
}

for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- sample[ , !(names(sample) %in% drop_col)]
  write.table(sample, 
              fileName,
              append = FALSE,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)  
}

#########################################

#Mus vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.1 <- row_t_welch(sample[,25:26:27], sample[,22:23:24])
  print(temp.stat.placeholder.1) 
  
  capture.output(fileName, append = TRUE, file = "Mus_Musinac.txt")
  capture.output(temp.stat.placeholder.1, append = TRUE, file = "Mus_Musinac.txt")
  
}

#Neu vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.2 <- row_t_welch(sample[,15:16:17], sample[,22:23:24])
  print(temp.stat.placeholder.2) 
  
  capture.output(fileName, append = TRUE, file = "Neu_Musinac.txt")
  capture.output(temp.stat.placeholder.2, append = TRUE, file = "Neu_Musinac.txt")
  
}

#Int vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.3 <- row_t_welch(sample[,12:13:14], sample[,22:23:24])
  print(temp.stat.placeholder.3) 
  
  capture.output(fileName, append = TRUE, file = "Int_Musinac.txt")
  capture.output(temp.stat.placeholder.3, append = TRUE, file = "Int_Musinac.txt")
  
}

#Isp.1 vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.4 <- row_t_welch(sample[,9:10:11], sample[,22:23:24])
  print(temp.stat.placeholder.4) 
  
  capture.output(fileName, append = TRUE, file = "Isp1_Musinac.txt")
  capture.output(temp.stat.placeholder.4, append = TRUE, file = "Isp1_Musinac.txt")
  
}

#EB vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.5 <- row_t_welch(sample[,6:7:8], sample[,22:23:24])
  print(temp.stat.placeholder.5) 
  
  capture.output(fileName, append = TRUE, file = "EB_Musinac.txt")
  capture.output(temp.stat.placeholder.5, append = TRUE, file = "EB_Musinac.txt")
  
}

#AZ vs Mus.inac
for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.6 <- row_t_welch(sample[,3:4:5], sample[,22:23:24])
  print(temp.stat.placeholder.6) 
  
  capture.output(fileName, append = TRUE, file = "AZ_Musinac.txt")
  capture.output(temp.stat.placeholder.6, append = TRUE, file = "AZ_Musinac.txt")
  
}
#########################################

pre.slice <- readline(prompt = "Enter name of miRNA file: ");
pre.slice <- out2_5p60.txt


pre.cross1 <- readline(prompt = "Set 1_1 (ie. A284_1): ")
pre.cross1 <- 'A284_1'
pre.cross2 <- readline(prompt = "Set 2_1 (ie. A285_1): ")
pre.cross2 <- 'A284_2'
pre.cross3 <- readline(prompt = "Set 1_2 (ie. A284_2): ")
pre.cross3 <- 'A290_1'
pre.cross4 <- readline(prompt = "Set 2_2 (ie. A285_2): ")
pre.cross4 <- 'A290_2'


slices1 <- pre.slice$A284_1
lbls1 <- pre.slice$seq_id
pie(slices1, labels = lbls1, main = "miRNA Variants Set: 284_1 ")

slices2 <- pre.slice$A286_1
lbls2 <- pre.slice$seq_id
pie(slices2, labels = lbls2, main = "miRNA Variants Set: 285_1 ")

slices3 <- pre.slice$A284_2
lbls3 <- pre.slice$seq_id
pie(slices3, labels = lbls3, main = "miRNA Variants Set: 284_2 ")

slices4 <- pre.slice$A286_2
lbls4 <- pre.slice$seq_id
pie(slices4, labels = lbls4, main = "miRNA Variants Set: 285_2 ")
