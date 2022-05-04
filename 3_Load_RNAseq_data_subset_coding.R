library(xlsx)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(DESeq2)

#There are 10 result tables from quantification. 

WT_12DPA_R1 <- read.table("Input/AT0375_5019AF_fc_lane_idx.genes.results", sep = "\t", header = T)
WT_12DPA_R2 <- read.table("Input/AT0376_5020AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_12DPA_R1 <- read.table("Input/AT0377_5021AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_12DPA_R2 <- read.table("Input/AT0378_5022AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_12DPA_R3 <- read.table("Input/AT0379_5023AF_fc_lane_idx.genes.results",sep = "\t", header=T)
WT_14DPA_R2 <- read.table("Input/AT0380_5024AF_fc_lane_idx.genes.results",sep = "\t", header=T)
WT_14DPA_R1 <- read.table("Input/AT0381_5025AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_14DPA_R1 <- read.table("Input/AT0382_5026AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_14DPA_R2 <- read.table("Input/AT0383_5027AF_fc_lane_idx.genes.results",sep = "\t", header=T)
MUT_14DPA_R3 <- read.table("Input/AT0384_5028AF_fc_lane_idx.genes.results",sep = "\t", header=T)

identical(WT_12DPA_R1$gene_id, WT_14DPA_R1$gene_id)

identical(WT_12DPA_R1$gene_id, WT_14DPA_R2$gene_id)

counts <- data.frame("gene_id" = MUT_12DPA_R1$gene_id, 
                     "transcript_id" = MUT_12DPA_R1$transcript_id.s., 
                     "MUT_12DPA_R1" = MUT_12DPA_R1$expected_count, 
                     "MUT_12DPA_R2" = MUT_12DPA_R2$expected_count, 
                     "MUT_12DPA_R3" = MUT_12DPA_R3$expected_count,
                     "WT_12DPA_R1" = WT_12DPA_R1$expected_count,
                     "WT_12DPA_R2" = WT_12DPA_R2$expected_count,
                     "MUT_14DPA_R1" = MUT_14DPA_R1$expected_count,
                     "MUT_14DPA_R2" = MUT_14DPA_R2$expected_count,
                     "MUT_14DPA_R3" = MUT_14DPA_R3$expected_count,
                     "WT_14DPA_R1" = WT_14DPA_R1$expected_count,
                     "WT_14DPA_R2" = WT_14DPA_R2$expected_count)

FPKM <- data.frame("gene_id" = MUT_12DPA_R1$gene_id, 
                     "transcript_id" = MUT_12DPA_R1$transcript_id.s., 
                     "MUT_12DPA_R1"= MUT_12DPA_R1$FPKM, 
                     "MUT_12DPA_R2"= MUT_12DPA_R2$FPKM, 
                     "MUT_12DPA_R3"= MUT_12DPA_R3$FPKM,
                     "WT_12DPA_R1" =  WT_12DPA_R1$FPKM,
                     "WT_12DPA_R2" =  WT_12DPA_R2$FPKM,
                     "MUT_14DPA_R1"= MUT_14DPA_R1$FPKM,
                     "MUT_14DPA_R2"= MUT_14DPA_R2$FPKM,
                     "MUT_14DPA_R3"= MUT_14DPA_R3$FPKM,
                     "WT_14DPA_R1" =  WT_14DPA_R1$FPKM,
                     "WT_14DPA_R2" =  WT_14DPA_R2$FPKM)

TPM <- data.frame("gene_id" = MUT_12DPA_R1$gene_id, 
                     "transcript_id" = MUT_12DPA_R1$transcript_id.s., 
                     "MUT_12DPA_R1"= MUT_12DPA_R1$TPM, 
                     "MUT_12DPA_R2"= MUT_12DPA_R2$TPM, 
                     "MUT_12DPA_R3"= MUT_12DPA_R3$TPM,
                     "WT_12DPA_R1" =  WT_12DPA_R1$TPM,
                     "WT_12DPA_R2" =  WT_12DPA_R2$TPM,
                     "MUT_14DPA_R1"= MUT_14DPA_R1$TPM,
                     "MUT_14DPA_R2"= MUT_14DPA_R2$TPM,
                     "MUT_14DPA_R3"= MUT_14DPA_R3$TPM,
                     "WT_14DPA_R1" =  WT_14DPA_R1$TPM,
                     "WT_14DPA_R2" =  WT_14DPA_R2$TPM)

#put names of the genes and restrict to pc

geneIDs2 <- ensembldb::select(EnsDb.Mmusculus.v79, 
                              keys=gsub(pattern = "\\.\\d+$",
                                        x =  counts$gene_id, replacement = ""),
                              keytype = "GENEID", columns = c("GENENAME", "GENEID", "TXBIOTYPE"))

geneIDsprotcod <- unique(geneIDs2[geneIDs2$TXBIOTYPE == "protein_coding",])

geneIDsprotcod <- geneIDsprotcod[,c(1:2)]

write.table(x = geneIDsprotcod, "output/geneIDsprotcod",col.names = T, quote = F)

counts$gene_id <- gsub(pattern = "\\.\\d+$",x =  counts$gene_id, replacement = "")
FPKM$gene_id <- gsub(pattern = "\\.\\d+$",x =  FPKM$gene_id, replacement = "")
TPM$gene_id <- gsub(pattern = "\\.\\d+$",x =  TPM$gene_id, replacement = "")

counts_pc <- counts[counts$gene_id %in% geneIDsprotcod$GENEID,]

counts_pc$gene_name <- ensembldb::select(EnsDb.Mmusculus.v79, keys=counts_pc$gene_id, keytype = "GENEID", columns = c("GENENAME"))

counts_pc <- do.call(data.frame, counts_pc)
counts_pc <- counts_pc[,c(13,14,3,4,5,6,7,8,9,10,11,12)]

write.table(counts_pc, "output/blastema_counts_proteincoding.txt", row.names = F, quote = F, sep = "," )
write.xlsx(counts_pc, "output/blastema_counts_proteincoding.xlsx",row.names = F)

FPKM_pc <- FPKM[FPKM$gene_id %in% geneIDsprotcod$GENEID,]

FPKM_pc$gene_name <- ensembldb::select(EnsDb.Mmusculus.v79, keys=FPKM_pc$gene_id, keytype = "GENEID", columns = c("GENENAME"))

FPKM_pc <- do.call(data.frame, FPKM_pc)
FPKM_pc <- FPKM_pc[,c(13,14,3,4,5,6,7,8,9,10,11,12)]

write.table(FPKM_pc, "output/blastema_fpkm_proteincoding.txt", row.names = F, quote = F, sep = "," )
write.xlsx(FPKM_pc, "output/blastema_fpkm_proteincoding.xlsx",row.names = F)

TPM_pc <- TPM[TPM$gene_id %in% geneIDsprotcod$GENEID,]
TPM_pc$gene_name <- ensembldb::select(EnsDb.Mmusculus.v79, keys=TPM_pc$gene_id, keytype = "GENEID", columns = c("GENENAME"))

TPM_pc <- do.call(data.frame, TPM_pc)
TPM_pc <- TPM_pc[,c(13,14,3,4,5,6,7,8,9,10,11,12)]

write.table(TPM_pc, "output/blastema_tpm_proteincoding.txt", row.names = F, quote = F, sep = "," )
write.xlsx(TPM_pc, "output/blastema_tpm_proteincoding.xlsx",row.names = F)

