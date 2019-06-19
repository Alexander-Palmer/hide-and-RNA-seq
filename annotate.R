library(BiocManager)
BiocManager::install("biomaRt)
library(biomaRt)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)

alg1 <- read.csv("GSE112753_alg-1_gk214_vs_WT_AD5.csv")
alg2 <- read.csv("GSE112753_alg-2_ok304_vs_WT_AD5.csv")

file <- read.delim("Copy of mir 71 targets.txt", header = FALSE)

mart <- useDataset("celegans_gene_ensembl", useMart("ensembl"))

Eid_file <- file$V4
Eid_alg1 <- alg1$ensembl.ID
Eid_alg2 <- alg2$ensembl.ID
alg1$gene <- NA
alg2$gene <- NA
file$gene <- NA

alg1_genes <- getBM(filters= "chembl", attributes= c("chembl", "gene_biotype", "external_gene_name",
                    "description"), values = Eid_alg1, mart= mart)
alg2_genes <- getBM(filters= "external_transcript_name", attributes= c("ensembl_transcript_id", "gene_biotype", "external_gene_name",
                    "description"), values = Eid_alg2, mart= mart)

file_genes <- getBM(filters= "ensembl_gene_id", attributes= c("wikigene_name", "goslim_goa_description", "family_description"), values = Eid_file, mart= mart)

write.table(file_genes, "gene-ontology_of_select_transcripts.txt", sep="\t", col.names = NA)


algm1 <- merge(x=alg1, y=alg1_genes, by.x="ensembl.ID",by.y="ensembl_transcript_id")
algm2 <- merge(x=alg2, y=alg2_genes, by.x="ensembl.ID",by.y="ensembl_transcript_id")

searchFilters(mart = mart, pattern = "base")
listAttributes()

alg1

alg1 <- alg1[,-6]
