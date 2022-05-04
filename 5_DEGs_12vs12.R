#################### DEGs  12 v 12 ################
rm(list = ls())

## load
library(xlsx)
library(pheatmap)
library(rtracklayer)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(ggplot2)
set.seed(123)

counts_pc <- read.table("output/blastema_counts_proteincoding.txt",
                        header = T,sep = ",")

FPKM_pc <- read.table("output/blastema_fpkm_proteincoding.txt",
                      header = T,sep = ",")

TPM_pc <- read.table("output/blastema_tpm_proteincoding.txt",
                     header = T,sep = ",")

geneIDsprotcod <- read.table("output/geneIDsprotcod", header = T)

samplesplan <- read.table("Input/samplesplan.txt", header = T, sep = ",")

rownames(samplesplan) <- samplesplan$sample
rownames(counts_pc)<-counts_pc$gene_name.GENEID

#12 dpa SAMPLES
dds <- DESeqDataSetFromMatrix(countData = round(counts_pc[c(3,4,5,6,7)]),
                              colData = samplesplan[samplesplan$stage=="12DPA",], design = ~genotype)
dds <- DESeq(dds)
res <- results(dds)

saveRDS(dds, "output/dds_12vs12.RDS")
saveRDS(res, "output/res_12vs12.RDS")

dds <- readRDS("output/dds_12vs12.RDS")
res <- readRDS("output/res_12vs12.RDS")

# DEG

res05 <- results(dds, alpha=0.05)
res05ordered <- res05[order(res05$padj),]
summary(res05)

res05ordered_df <- add_rownames(as.data.frame(res05ordered, var = "GENEID"))

res05ordered_df$genename <- geneIDsprotcod$GENENAME[match(res05ordered_df$rowname, geneIDsprotcod$GENEID)]

res05ordered_df <- res05ordered_df[, c(1,8,2,3,4,5,6,7)]

outputtable_005 <- merge(res05ordered_df, TPM_pc, by.x ="rowname", by.y="gene_name.GENEID")

outputtable_005 <- subset(outputtable_005,padj<0.05 & abs(log2FoldChange) > 1.5)

write.xlsx(as.data.frame(outputtable_005), "output/TableS1_DESeq2_12vs12.xlsx", row.names = F)

write.table(as.data.frame(outputtable_005), "output/DESeq2_12vs12.txt", row.names = F,
            quote = F, col.names = T,sep = "\t")


df.sub <- log2(1 + TPM_pc[TPM_pc$gene_name.GENENAME %in% outputtable_005$genename,c(3:7)])
df.sub <- df.sub[rowSums(df.sub[, -(3:7)]) > 0, ]


pdf("Plots/Heatmap_12vs12_significant_scaled.pdf",
    title="Plots/Heatmap_12vs12_significant_scaled.pdf",
    onefile = TRUE, width= 7, height= 15)
breaksListAbs<-c(min(apply(df.sub, 1, scale), -1.5) - 0.0000001,
                 seq(-1.5, 1.5, 0.03),
                 max(apply(df.sub, 1, scale), 1.5) + 0.0000001)
pheatmap(df.sub,
         #labels_row=TPM_pc$gene_name.GENENAME %in% outputtable_005$genename,
         cluster_rows= TRUE,
         cluster_cols= FALSE ,
         breaks = breaksListAbs,
         main="log2(1+TPM) 12dpa WTvsMUT significant genes",
         #clustering_method="correlation",
         clustering_method="ward.D2",
         scale="row",
         color = colorRampPalette(c("blue","white","red"))(103))
dev.off()


wald.test12 <- read.table("output/DESeq2_12vs12.txt", header = T)

wald.test.up12 <- wald.test12[wald.test12$log2FoldChange>0,]
wald.test.down12 <- wald.test12[wald.test12$log2FoldChange<0,]

