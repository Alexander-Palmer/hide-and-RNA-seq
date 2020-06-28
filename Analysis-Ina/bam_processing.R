###Ina's sequencing data analysis
###
###Annotate and quantify the presence of transcripts associated with dilated cardiomyopathy in mice treated with novel
###cancer drug doxorubicin. Packages featureCounts and edgeR will be used.

##Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")
BiocManager::install("edgeR")
BiocManager::install("biomaRt")
library(Rsubread)
library(edgeR)
library(biomaRt)

##File import
proc.files <- c("black_6_control_1_left_ventricle_s_dm.bam", "black_6_control_1_right_ventricle_s_dm.bam",
                "black_6_control_2_left_ventricle_s_dm.bam", "black_6_control_2_right_ventricle_s_dm_new7.bam",
                "black_6_control_3_left_ventricle_s_dm.bam", "black_6_control_3_right_ventricle_s_dm_new3.bam",
                "black_6_control_4_left_ventricle_s_dm_new.bam", "black_6_control_4_right_ventricle_s_dm.bam",
                "doxorubicin_treated_non-responder_7_left_ventricle_s_dm.bam", 
                "doxorubicin_treated_non-responder_7_right_ventricle_s_dm.bam",
                "doxorubicin_treated_non-responder_8_left_ventricle_s_dm.bam",
                "doxorubicin_treated_non-responder_8_right_ventricle_s_dm.bam", #This one
                "doxorubicin_treated_non-responder_9_left_ventricle_s_dm.bam",
                "doxorubicin_treated_non-responder_9_right_ventricle_s_dm.bam",
                "doxorubicin_treated_non-responder_10_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_1_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_1_right_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_2_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_2_right_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_3_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_3_right_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_4_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_4_right_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_5_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_6_left_ventricle_s_dm.bam",
                "doxorubicin_treated_responder_6_right_ventricle_s_dm.bam")

proc.files <- c("imb_reid_2014_02_01_FAIRE_control_1_sorted.bam", "imb_reid_2014_02_02_FAIRE_control_2_sorted.bam",
                "imb_reid_2014_02_03_FAIRE_control_3_sorted.bam", "imb_reid_2014_02_04_FAIRE_control_4_sorted.bam",
                "imb_reid_2014_02_05_FAIRE_control_5_sorted.bam", "imb_reid_2014_02_06_FAIRE_control_6_sorted.bam",
                "imb_reid_2014_02_07_FAIRE_control_7_sorted.bam", "imb_reid_2014_02_08_FAIRE_control_8_sorted.bam",
                "imb_reid_2014_02_09_FAIRE_Dox_1_sorted.bam", "imb_reid_2014_02_10_FAIRE_Dox_2_sorted.bam",
                "imb_reid_2014_02_11_FAIRE_Dox_3_sorted.bam", "imb_reid_2014_02_12_FAIRE_Dox_4_sorted.bam",
                "imb_reid_2014_02_13_FAIRE_Dox_5_sorted.bam", "imb_reid_2014_02_14_FAIRE_Dox_6_sorted.bam",
                "imb_reid_2014_02_15_FAIRE_Dox_7_sorted.bam", "imb_reid_2014_02_16_FAIRE_Dox_8_sorted.bam",
                "imb_reid_2014_02_17_FAIRE_Dox_E_1_sorted.bam", "imb_reid_2014_02_18_FAIRE_Dox_E_2_sorted.bam",
                "imb_reid_2014_02_19_FAIRE_Dox_E_3_sorted.bam", "imb_reid_2014_02_20_FAIRE_Dox_E_4_sorted.bam",
                "imb_reid_2014_02_21_FAIRE_Dox_E_5_sorted.bam", "imb_reid_2014_02_22_FAIRE_Dox_E_6_sorted.bam",
                "imb_reid_2014_02_23_FAIRE_Dox_E_7_sorted.bam", "imb_reid_2014_02_24_FAIRE_Dox_E_8_sorted.bam")

##featureCounts analysis
#Annotate reads
meta_results_mm9_3 <-featureCounts(files = proc.files,
                                 annot.ext = "gencode.vM1.annotation.gtf",
                                 isGTFAnnotationFile = TRUE, GTF.featureType = "gene", GTF.attrType = "gene_id",
                                 isPairedEnd = T, nthreads = 18, useMetaFeatures = T,
                                 genome = "Mus_musculus.NCBIM37.67.dna.toplevel.fa",
                                 juncCounts = T, allowMultiOverlap = T)
#Create counts table
write.table(x = data.frame(meta_results_mm9$annotation[,c("GeneID","Length")],
                           meta_results_mm9$counts, stringsAsFactors = FALSE),
                           file = "entpackt_seq_counts_final.txt", quote = FALSE, sep = "\t")

##edgeR analysis
#File import
faire_seq_clean <- read.csv('FAIRE_seq_clean.csv', header = T)
faire_seq_identity <- read.csv('FAIRE_seq_identity2.csv', header = T)

#Processing
faire_seq_clean <- as.data.frame(faire_seq_clean)
faire_seq_clean[is.na(faire_seq_clean)] <- 0
sequences <- faire_seq_identity$GeneID
levels <- letters[1:24]
levels[c(1:8)] <- "Control"
levels[c(9:16)] <- "Dox_non_resp"
levels[c(17:24)] <- "Dox_resp"

#Processing 2
faire_seq_clean <- DGEList(counts=faire_seq_clean, group=levels, genes=sequences)
keep <- rowSums(cpm(faire_seq_clean)>1) >= 1
faire_seq_clean <- faire_seq_clean[keep, , keep.lib.sizes=FALSE]
faire_seq_clean <- calcNormFactors(faire_seq_clean)
faire_seq_clean <- estimateDisp(faire_seq_clean)
design <- model.matrix(~0+levels, data=faire_seq_clean$samples)
colnames(design) <- levels(faire_seq_clean$samples$levels)
fit <- glmQLFit(faire_seq_clean, design)

#Multiple comparisons
#Control_vs_Dox_non_resp
Plcehldr <- glmQLFTest(fit, contrast=c(1,-1,0))
control_dox_nr_a <- topTags(Plcehldr, n=Inf)
control_dox_nr_a <- as.data.frame(control_dox_nr_a)
control_dox_nr_b <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
control_dox_nr_b <- as.data.frame(control_dox_nr_b)

#Control_vs_Dox_resp
Plcehldr <- glmQLFTest(fit, contrast=c(1,0,-1))
control_dox_r_a <- topTags(Plcehldr, n=Inf)
control_dox_r_a <- as.data.frame(control_dox_r_a)
control_dox_r_b <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
control_dox_r_b <- as.data.frame(control_dox_r_b)

#Dox_non_resp_vs_Dox_resp
Plcehldr <- glmQLFTest(fit, contrast=c(0,1,-1))
dox_nr_dox_r_a <- topTags(Plcehldr, n=Inf)
dox_nr_dox_r_a <- as.data.frame(dox_nr_dox_r_a)
dox_nr_dox_r_b <- topTags(Plcehldr, sort.by = "logFC", n=Inf, p.value=0.05)
dox_nr_dox_r_b <- as.data.frame(dox_nr_dox_r_b)

##Annotating results
#Datasets
mRNA_results_mm9 <- read.delim("FAIRE_seq_counts_final.txt", header = TRUE)
mart = useMart("ensembl",host="http://may2012.archive.ensembl.org", dataset = "mmusculus_gene_ensembl")

#All mRNAs
sig_mRNAs <- mRNA_results_mm9$GeneID
genes_placeholder_a <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "gene_biotype",
                                                                         "external_gene_id", "description"),
                                                                        #"goslim_goa_description"),
                                                                         values = sig_mRNAs, mart = mart)
#Control_dox_non_resp
sig_mRNAs <- control_dox_nr_a$genes
genes_placeholder_b <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "gene_biotype",
                                                                         "external_gene_id", "description"),
                                                                        #"goslim_goa_description"),
                                                                         values = sig_mRNAs, mart = mart)
#Control_dox_resp
sig_mRNAs <- control_dox_r_a$genes
genes_placeholder_c <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "gene_biotype",
                                                                         "external_gene_id", "description"),
                                                                        #"goslim_goa_description"),
                                                                         values = sig_mRNAs, mart = mart)
#Dox_non_resp_dox_resp
sig_mRNAs <- dox_nr_dox_r_a$genes
genes_placeholder_d <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "gene_biotype",
                                                                         "external_gene_id", "description"),
                                                                        #"goslim_goa_description"),
                                                                         values = sig_mRNAs, mart = mart)

#Merge datasets and annotations
colnames(genes_placeholder_a) <- c("GeneID", "gene_biotype", "external gene name", "description")
genes_meta <- merge(mRNA_results_mm9, genes_placeholder_a, by="GeneID", all = TRUE)
colnames(genes_placeholder_b) <- c("genes", "gene_biotype", "external gene name", "description")
genes_control_dox_non_resp <- merge(control_dox_nr_a, genes_placeholder_b, by="genes", all = TRUE)
colnames(genes_placeholder_c) <- c("genes", "gene_biotype", "external gene name", "description")
genes_control_dox_resp <- merge(control_dox_r_a, genes_placeholder_c, by="genes", all = TRUE)
colnames(genes_placeholder_d) <- c("genes", "gene_biotype", "external gene name", "description")
genes_dox_non_resp_dox_resp <- merge(dox_nr_dox_r_a, genes_placeholder_d, by="genes", all = TRUE)


##Data output
write.table(genes_meta, "FAIRE-seq all genes.csv", sep=",", col.names = NA, quote = FALSE)
write.table(genes_control_dox_non_resp, "FAIRE-seq control vs dox.csv", sep=",", col.names = NA, quote = FALSE)
write.table(genes_control_dox_resp, "FAIRE-seq control vs dox E.csv", sep=",", col.names = NA, quote = FALSE)
write.table(genes_dox_non_resp_dox_resp, "FAIRE-seq dox vs dox E.csv", sep=",",col.names = NA,quote = FALSE)
