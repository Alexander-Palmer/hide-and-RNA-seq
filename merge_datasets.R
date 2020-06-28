#Append isomiR type to isomiR datasets

isomiR_type <- read.delim("META_Ext_Red_type.txt")


isomiR_type$Type <- ifelse(isomiR_type[,8] > 0, "Ext 5'", NA)
isomiR_type$Type <- ifelse(isomiR_type[,9] > 0, "Ext 3'", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,10] > 0, "Red 5'", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,11] > 0, "Red 3'", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,12] > 0, "Mixed Ext/Red", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,13] > 0, "Mixed Ext/Red", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,14] > 0, "Mixed Ext/Red", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,15] > 0, "Mixed Ext/Red", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,16] > 0, "SNP", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,17] > 0, "Substitution", isomiR_type$Type)
isomiR_type$Type <- ifelse(isomiR_type[,18] > 0, "Canonical", isomiR_type$Type)

Int1_Mus1_71 <- read.csv("IntALG1 vs MusALG1 - mir71.csv")
Int1_Mus1_71 <- merge(Int1_Mus1_71, isomiR_type, by="genes", all = FALSE)

Int2_Mus2_71 <- read.csv("IntALG2 vs MusALG2 - mir71.csv")
Int2_Mus2_71 <- merge(Int2_Mus2_71, isomiR_type, by="genes", all = FALSE)

Int1_Neu1_71 <- read.csv("IntALG1 vs NeuALG1 - mir71.csv")
Int1_Neu1_71 <- merge(Int1_Neu1_71, isomiR_type, by="genes", all = FALSE)

Int2_Neu2_71 <- read.csv("IntALG2 vs NeuALG2 - mir71.csv")
Int2_Neu2_71 <- merge(Int2_Neu2_71, isomiR_type, by="genes", all = FALSE)

Neu1_Mus1_71 <- read.csv("NeuALG1 vs MusALG1 - mir71.csv")
Neu1_Mus1_71 <- merge(Neu1_Mus1_71, isomiR_type, by="genes", all = FALSE)

Neu2_Mus2_71 <- read.csv("NeuALG2 vs MusALG2 - mir71.csv")
Neu2_Mus2_71 <- merge(Neu2_Mus2_71, isomiR_type, by="genes", all = FALSE)

#####

Int1_Mus1_90 <- read.csv("IntALG1 vs MusALG1 - mir90.csv")
Int1_Mus1_90 <- merge(Int1_Mus1_90, isomiR_type, by="genes", all = FALSE)

Int2_Mus2_90 <- read.csv("IntALG2 vs MusALG2 - mir90.csv")
Int2_Mus2_90 <- merge(Int2_Mus2_90, isomiR_type, by="genes", all = FALSE)

Int1_Neu1_90 <- read.csv("IntALG1 vs NeuALG1 - mir90.csv")
Int1_Neu1_90 <- merge(Int1_Neu1_90, isomiR_type, by="genes", all = FALSE)

Int2_Neu2_90 <- read.csv("IntALG2 vs NeuALG2 - mir90.csv")
Int2_Neu2_90 <- merge(Int2_Neu2_90, isomiR_type, by="genes", all = FALSE)

Neu1_Mus1_90 <- read.csv("NeuALG1 vs MusALG1 - mir90.csv")
Neu1_Mus1_90 <- merge(Neu1_Mus1_90, isomiR_type, by="genes", all = FALSE)

Neu2_Mus2_90 <- read.csv("NeuALG2 vs MusALG2 - mir90.csv")
Neu2_Mus2_90 <- merge(Neu2_Mus2_90, isomiR_type, by="genes", all = FALSE)
