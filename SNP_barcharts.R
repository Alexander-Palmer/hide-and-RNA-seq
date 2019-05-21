 datalist_Count = list()

  for (Filename in Int_A) {
    aPlcehldr <- a[grep(paste(Filename, collapse = "|"), a[,45]), ]
    aPlcehldr$Count_1 <- NA
    aPlcehldr$Count_1 <- ifelse(as.numeric(aPlcehldr[,24]) == 1, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_2 <- NA
    aPlcehldr$Count_2 <- ifelse(as.numeric(aPlcehldr[,24]) == 2, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_3 <- NA
    aPlcehldr$Count_3 <- ifelse(as.numeric(aPlcehldr[,24]) == 3, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_4 <- NA
    aPlcehldr$Count_4 <- ifelse(as.numeric(aPlcehldr[,24]) == 4, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_5 <- NA
    aPlcehldr$Count_5 <- ifelse(as.numeric(aPlcehldr[,24]) == 5, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_6 <- NA
    aPlcehldr$Count_6 <- ifelse(as.numeric(aPlcehldr[,24]) == 6, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_7 <- NA
    aPlcehldr$Count_7 <- ifelse(as.numeric(aPlcehldr[,24]) == 7, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_8 <- NA
    aPlcehldr$Count_8 <- ifelse(as.numeric(aPlcehldr[,24]) == 8, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_9 <- NA
    aPlcehldr$Count_9 <- ifelse(as.numeric(aPlcehldr[,24]) == 9, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_10 <- NA
    aPlcehldr$Count_10 <- ifelse(as.numeric(aPlcehldr[,24]) == 10, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_11 <- NA
    aPlcehldr$Count_11 <- ifelse(as.numeric(aPlcehldr[,24]) == 11, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_12 <- NA
    aPlcehldr$Count_12 <- ifelse(as.numeric(aPlcehldr[,24]) == 12, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_13 <- NA
    aPlcehldr$Count_13 <- ifelse(as.numeric(aPlcehldr[,24]) == 13, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_14 <- NA
    aPlcehldr$Count_14 <- ifelse(as.numeric(aPlcehldr[,24]) == 14, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_15 <- NA
    aPlcehldr$Count_15 <- ifelse(as.numeric(aPlcehldr[,24]) == 15, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_16 <- NA
    aPlcehldr$Count_16 <- ifelse(as.numeric(aPlcehldr[,24]) == 16, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_17 <- NA
    aPlcehldr$Count_17 <- ifelse(as.numeric(aPlcehldr[,24]) == 17, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_18 <- NA
    aPlcehldr$Count_18 <- ifelse(as.numeric(aPlcehldr[,24]) == 18, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_19 <- NA
    aPlcehldr$Count_19 <- ifelse(as.numeric(aPlcehldr[,24]) == 19, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_20 <- NA
    aPlcehldr$Count_20 <- ifelse(as.numeric(aPlcehldr[,24]) == 20, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_21 <- NA
    aPlcehldr$Count_21 <- ifelse(as.numeric(aPlcehldr[,24]) == 21, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_22 <- NA
    aPlcehldr$Count_22 <- ifelse(as.numeric(aPlcehldr[,24]) == 22, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_23 <- NA
    aPlcehldr$Count_23 <- ifelse(as.numeric(aPlcehldr[,24]) == 23, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_24 <- NA
    aPlcehldr$Count_24 <- ifelse(as.numeric(aPlcehldr[,24]) == 24, 1/NROW(aPlcehldr), 0)
    aPlcehldr$Count_25 <- NA
    aPlcehldr$Count_25 <- ifelse(as.numeric(aPlcehldr[,24]) == 25, 1/NROW(aPlcehldr), 0)

    datalist_Count[[Filename]] <- aPlcehldr
  }

  nt_1 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_1')), na.rm = TRUE)
  nt_2 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_2')), na.rm = TRUE)
  nt_3 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_3')), na.rm = TRUE)
  nt_4 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_4')), na.rm = TRUE)
  nt_5 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_5')), na.rm = TRUE)
  nt_6 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_6')), na.rm = TRUE)
  nt_7 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_7')), na.rm = TRUE)
  nt_8 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_8')), na.rm = TRUE)
  nt_9 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_9')), na.rm = TRUE)
  nt_10 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_10')), na.rm = TRUE)
  nt_11 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_11')), na.rm = TRUE)
  nt_12 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_12')), na.rm = TRUE)
  nt_13 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_13')), na.rm = TRUE)
  nt_14 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_14')), na.rm = TRUE)
  nt_15 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_15')), na.rm = TRUE)
  nt_16 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_16')), na.rm = TRUE)
  nt_17 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_17')), na.rm = TRUE)
  nt_18 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_18')), na.rm = TRUE)
  nt_19 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_19')), na.rm = TRUE)
  nt_20 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_20')), na.rm = TRUE)
  nt_21 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_21')), na.rm = TRUE)
  nt_22 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_22')), na.rm = TRUE)
  nt_23 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_23')), na.rm = TRUE)
  nt_24 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_24')), na.rm = TRUE)
  nt_25 <- as.numeric(do.call(sum, lapply(datalist_Count, '[', 'Count_25')), na.rm = TRUE)

  total_SNP <- c(nt_1, nt_2, nt_3, nt_4, nt_5, nt_6, nt_7, nt_8, nt_9, nt_10, nt_11, nt_12,
                 nt_13, nt_14, nt_15, nt_16, nt_17, nt_18, nt_19, nt_20, nt_21, nt_22, nt_23,
                 nt_24, nt_25)
  total_pos <- c("pos_1", "pos_2", "pos_3", "pos_4", "pos_5", "pos_6", "pos_7", "pos_8", "pos_9",
                 "pos_10", "pos_11", "pos_12", "pos_13", "pos_14", "pos_15", "pos_16", "pos_17",
                 "pos_18", "pos_19", "pos_20", "pos_21", "pos_22", "pos_23", "pos_24", "pos_25")

  xlab <- list(title = "Nucleotide position",
               categoryorder = "array",
               categoryarray = c("pos_1", "pos_2", "pos_3", "pos_4", "pos_5", "pos_6", "pos_7", "pos_8", "pos_9",
                                 "pos_10", "pos_11", "pos_12", "pos_13", "pos_14", "pos_15", "pos_16", "pos_17",
                                 "pos_18", "pos_19", "pos_20", "pos_21", "pos_22", "pos_23", "pos_24", "pos_25"))
  ylab <- list(title = "# reads/isomiRs")

  p <- plot_ly(x = total_pos, y = total_SNP, type = "bar") %>%
    layout(xaxis=xlab, yaxis=ylab)
  p
}
