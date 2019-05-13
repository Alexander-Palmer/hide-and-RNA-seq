source("http://bioconductor.org/biocLite.R")
biocLite("matrixTests")
install.packages("viridis")
library("viridis")
library("matrixTests")

options(max.print=1000000)

fileNames <- Sys.glob("*.txt")

temp = list.files(pattern="*.txt")
myfiles = lapply(temp, read.delim)
for (i in 1:length(temp)) assign(temp[i], read.delim(temp[i]))

data.normalisation <- function(Z) {
  near_output <- (data.matrix(Z))
  output <- prop.table(near_output, 2)
  return(output)
}  

for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample[paste('av_284')] <- (sample[["A284_1"]] + sample[["A284_2"]])/2
  sample[paste('av_285')] <- (sample[["A285_1"]] + sample[["A285_2"]])/2
  sample[paste('av_286')] <- (sample[["A286_1"]] + sample[["A286_2"]])/2
  sample[paste('av_287')] <- (sample[["A287_1"]] + sample[["A287_2"]])/2
  sample[paste('av_290')] <- (sample[["A290_1"]] + sample[["A290_2"]])/2
  sample[paste('av_291')] <- (sample[["A291_1"]] + sample[["A291_2"]])/2
  sample[paste('Intcomb')] <- (sample[["av_284"]] + sample[["av_285"]])
  sample[paste('Neucomb')] <- (sample[["av_286"]] + sample[["av_287"]])
  sample[paste('Muscomb')] <- (sample[["av_290"]] + sample[["av_291"]])
  sample[paste('av_N2')] <- (sample[["N2_1"]] + sample[["N2_2"]])/2
  
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

#########################################

Int1excl <- read.delim("XXList sig miRNA int ALG1.txt", header=F)$V1
Int1excl <- as.character(Int1excl)
Int2excl <- read.delim("XXList sig miRNA int ALG2.txt", header=F)$V1
Int2excl <- as.character(Int2excl)

Mus1excl <- read.delim("XXList sig miRNA mus ALG1.txt", header=F)$V1
Mus1excl <- as.character(Mus1excl)
Mus2excl <- read.delim("XXList sig miRNA mus ALG2.txt", header=F)$V1
Mus2excl <- as.character(Mus2excl)

Neu1excl <- read.delim("XXList sig miRNA neu ALG1.txt", header=F)$V1
Neu1excl <- as.character(Neu1excl)
Neu2excl <- read.delim("XXList sig miRNA neu ALG2.txt", header=F)$V1
Neu2excl <- as.character(Neu2excl)


Int12excl <- union(Int1excl, Int2excl)
Mus12excl <- union(Mus1excl, Mus2excl)
Neu12excl <- union(Neu1excl, Neu2excl)
IntMus <- intersect(Int12excl, Mus12excl)
IntNeu <- intersect(Int12excl, Neu12excl)
MusNeu <- intersect(Mus12excl, Neu12excl)

Int1Mus1excl <- union(Int1excl, Mus1excl)
Int1Neu1excl <- union(Int1excl, Neu1excl)
Int2Mus2excl <- union(Int2excl, Mus2excl)
Int2Neu2excl <- union(Int2excl, Neu2excl)
Mus1Neu1excl <- union(Mus1excl, Neu1excl)
Mus2Neu2excl <- union(Mus2excl, Neu2excl)


#########################################

#284x285 (Int1xInt2)
fileNames1.1 <- intersect(fileNames, Int12excl)
for (fileName in fileNames1.1) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  sample[is.na(sample)] <- 0
  print(fileName)
  
  temp.stat.placeholder.1 <- row_t_welch(sample[,3:4], sample[,5:6])
  print(temp.stat.placeholder.1) 
  
  capture.output(fileName, append = TRUE, file = "out2_284x285.txt")
  capture.output(temp.stat.placeholder.1, append = TRUE, file = "out2_284x285.txt")
  
 }

#284x290 (Int1xMus1)
fileNames1.2 <- intersect(fileNames, Int1Mus1excl)
for (fileName in fileNames1.2) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.2 <- row_t_welch((sample[,3:4], sample[,15:16])
  print(temp.stat.placeholder.2)
  
  capture.output(fileName, append = TRUE, file = "out2_284x290.txt")
  capture.output(temp.stat.placeholder.2, append = TRUE, file = "out2_284x290.txt")

}

#284x291
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.3 <- row_t_welch(sample[,c("A284_1", "A284_2")], sample[,c("A291_1", "A291_2")])
  print(temp.stat.placeholder.3) 
  
  capture.output(fileName, append = TRUE, file = "out2_284x291.txt")
  capture.output(temp.stat.placeholder.3, append = TRUE, file = "out2_284x291.txt")
}

#284x286 (Int1xNeu1)
fileNames1.4 <- intersect(fileNames, Int1Neu1excl)
for (fileName in fileNames1.4) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.4 <- row_t_welch(sample[,3:4], sample[,7:8])
  print(temp.stat.placeholder.4) 
  
  capture.output(fileName, append = TRUE, file = "out2_284x286.txt")
  capture.output(temp.stat.placeholder.4, append = TRUE, file = "out2_284x286.txt")
}

#284x287
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.5 <- row_t_welch(sample[,c("A284_1", "A284_2")], sample[,c("A287_1", "A287_2")])
  print(temp.stat.placeholder.5) 
  
  capture.output(fileName, append = TRUE, file = "out2_284x287.txt")
  capture.output(temp.stat.placeholder.5, append = TRUE, file = "out2_284x287.txt")
}

#285x290
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.6 <- row_t_welch(sample[,c("A285_1", "A285_2")], sample[,c("A290_1", "A290_2")])
  print(temp.stat.placeholder.6) 
  
  capture.output(fileName, append = TRUE, file = "out2_285x290.txt")
  capture.output(temp.stat.placeholder.6, append = TRUE, file = "out2_285x290.txt")
}

#285x291 (Int2xMus2)
fileNames1.7 <- intersect(fileNames, Int2Mus2excl)
for (fileName in fileNames1.7) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.7 <- row_t_welch(sample[,5:6], sample[,17:18])
  print(temp.stat.placeholder.7) 
  
  capture.output(fileName, append = TRUE, file = "out2_285x291.txt")
  capture.output(temp.stat.placeholder.7, append = TRUE, file = "out2_285x291.txt")
}

#285x286
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.8 <- row_t_welch(sample[,c("A285_1", "A285_2")], sample[,c("A286_1", "A286_2")])
  print(temp.stat.placeholder.8) 
  
  capture.output(fileName, append = TRUE, file = "out2_285x286.txt")
  capture.output(temp.stat.placeholder.8, append = TRUE, file = "out2_285x286.txt")
}

#285x287 (Int2xNeu2)
fileNames1.9 <- intersect(fileNames, Int2Neu2excl)
for (fileName in fileNames1.9) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.9 <- row_t_welch(sample[,5:6], sample[,9:10])
  print(temp.stat.placeholder.9) 
  
  capture.output(fileName, append = TRUE, file = "out2_285x287.txt")
  capture.output(temp.stat.placeholder.9, append = TRUE, file = "out2_285x287.txt")
}

#290x291 (Mus1xMus2)
fileNames1.10 <- intersect(fileNames, Mus12excl)
for (fileName in fileNames1.10) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.10 <- row_t_welch(sample[,15:16], sample[,17:18])
  print(temp.stat.placeholder.10)
  
  capture.output(fileName, append = TRUE, file = "out2_290x291.txt")
  capture.output(temp.stat.placeholder.10, append = TRUE, file = "out2_290x291.txt")
}

#290x286 (Mus1xNeu1)
fileNames1.11 <- intersect(fileNames, Mus1Neu1excl)
for (fileName in fileNames1.11) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.11 <- row_t_welch(sample[,15:16], sample[,7:8])
  print(temp.stat.placeholder.11)
  
  capture.output(fileName, append = TRUE, file = "out2_290x286.txt")
  capture.output(temp.stat.placeholder.11, append = TRUE, file = "out2_290x286.txt")
}

#290x287
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.12 <- row_t_welch(sample[,c("A290_1", "A290_2")], sample[,c("A287_1", "A287_2")])
  print(temp.stat.placeholder.12)
  
  capture.output(fileName, append = TRUE, file = "out2_290x287.txt")
  capture.output(temp.stat.placeholder.12, append = TRUE, file = "out2_290x287.txt")
}

#291x286 
#for (fileName in fileNames) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.13 <- row_t_welch(sample[,c("A291_1", "A291_2")], sample[,c("A286_1", "A286_2")])
  print(temp.stat.placeholder.13)
  
  capture.output(fileName, append = TRUE, file = "out2_291x286.txt")
  capture.output(temp.stat.placeholder.13, append = TRUE, file = "out2_291x286.txt")
}

#291x287 (Mus2xNeu2)
fileNames1.14 <- intersect(fileNames, Mus2Neu2excl)
for (fileName in fileNames1.14) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.14 <- row_t_welch(sample[,17:18], sample[,9:10])
  print(temp.stat.placeholder.14)
  
  capture.output(fileName, append = TRUE, file = "out2_291x287.txt")
  capture.output(temp.stat.placeholder.14, append = TRUE, file = "out2_291x287.txt")
}

#286x287 (Neu1xNeu2)
fileNames1.15 <- intersect(fileNames, Neu12excl)
for (fileName in fileNames1.15) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.15 <- row_t_welch(sample[,7:8], sample[,9:10])
  print(temp.stat.placeholder.15)
  
  capture.output(fileName, append = TRUE, file = "out2_286x287.txt")
  capture.output(temp.stat.placeholder.15, append = TRUE, file = "out2_286x287.txt")
}

#IntxMus
fileNames2.1 <- intersect(fileNames, IntMus)
for (fileName in fileNames2.1) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.2.1 <- row_t_welch(sample[,21:22], sample[,25:26])
  print(temp.stat.placeholder.2.1)
  
  capture.output(fileName, append = TRUE, file = "IntxMus.txt")
  capture.output(temp.stat.placeholder.2.1, append = TRUE, file = "IntxMus.txt")
}

#IntxNeu
fileNames2.2 <- intersect(fileNames, IntMus)
for (fileName in fileNames2.2) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.2.2 <- row_t_welch(sample[,21:22], sample[,23:24])
  print(temp.stat.placeholder.2.2)
  
  capture.output(fileName, append = TRUE, file = "IntxNeu.txt")
  capture.output(temp.stat.placeholder.2.2, append = TRUE, file = "IntxNeu.txt")
}

#MusxNeu
fileNames2.3 <- intersect(fileNames, IntMus)
for (fileName in fileNames2.3) {
  sample <- read.delim(fileName)
  sample <- data.matrix(sample)
  print(fileName)
  
  temp.stat.placeholder.2.3 <- row_t_welch(sample[,25:26], sample[,23:24])
  print(temp.stat.placeholder.2.3)
  
  capture.output(fileName, append = TRUE, file = "MusxNeu.txt")
  capture.output(temp.stat.placeholder.2.3, append = TRUE, file = "MusxNeu.txt")
}

#########################################

pre.slice <- readline(prompt = "Enter name of miRNA file: ");
pre.slice <- out2_3p1830.txt

pre.cross1 <- readline(prompt = "Set 1_1 (ie. A284_1): ")
pre.cross1 <- 'A284_1'
pre.cross2 <- readline(prompt = "Set 2_1 (ie. A285_1): ")
pre.cross2 <- 'A284_2'
pre.cross3 <- readline(prompt = "Set 1_2 (ie. A284_2): ")
pre.cross3 <- 'A290_1'
pre.cross4 <- readline(prompt = "Set 2_2 (ie. A285_2): ")
pre.cross4 <- 'A290_2'



slices1 <- pre.slice$A285_1
lbls1 <- pre.slice$A285_1
pie(slices1, labels = lbls1, main = "miRNA Variants Set: 286_1 ", col = rainbow(length(slices1)))
chartlabels1 <- pre.slice$seq_id
legend("bottomright", legend=chartlabels1, cex=0.7, fill = rainbow(length(slices1)))

slices2 <- pre.slice$A287_1
lbls2 <- pre.slice$A287_1
pie(slices2, labels = lbls2, main = "miRNA Variants Set: 287_1 ", col = rainbow(length(slices2)))
chartlabels2 <- pre.slice$seq_id
legend("bottomright", legend=chartlabels2, cex=0.7, fill = rainbow(length(slices1)))

slices3 <- pre.slice$A285_2
lbls3 <- pre.slice$A285_2
pie(slices3, labels = lbls3, main = "miRNA Variants Set: 286_2 ", col = rainbow(length(slices3)))
chartlabels3 <- pre.slice$seq_id
legend("bottomright", legend=chartlabels3, cex=0.7, fill = rainbow(length(slices1)))

slices4 <- pre.slice$A287_2
lbls4 <- pre.slice$A287_2
pie(slices4, labels = lbls4, main = "miRNA Variants Set: 287_2 ", col = rainbow(length(slices4)))
chartlabels4 <- pre.slice$seq_id
legend("bottomright", legend=chartlabels4, cex=0.7, fill = rainbow(length(slices1)))
