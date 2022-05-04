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
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(ggrepel)

counts_pc <- read.table("output/blastema_counts_proteincoding.txt",header = T,sep = ",")

FPKM_pc <- read.table("output/blastema_fpkm_proteincoding.txt",header = T,sep = ",")

TPM_pc <- read.table("output/blastema_tpm_proteincoding.txt",header = T,sep = ",")

geneIDsprotcod <- read.table("output/geneIDsprotcod", header = T) # manually removed row 20823, 21514 they dont have assigned names

#creating the deseq2 object and computing results

samplesplan <- read.table("Input/samplesplan.txt", header = T, sep = ",")
rownames(samplesplan) <- samplesplan$sample
rownames(counts_pc)<-counts_pc$gene_name.GENEID
dds <- DESeqDataSetFromMatrix(countData = round(counts_pc[3:12]),colData = samplesplan, design = ~genotype)
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds <- DESeq(dds)
res <- results(dds)

saveRDS(dds, "output/dds.RDS")
dds <- readRDS("output/dds.RDS")
saveRDS(res, "output/res.RDS")
res <- readRDS("output/res.RDS")
res <- data.frame(res)
resnames <- add_rownames(as.data.frame(res, var = "rowname"))

resnames <- merge(resnames, TPM_pc, by.x ="rowname", by.y="gene_name.GENEID")

writexl::write_xlsx(resnames, "output/WT_VS_MUT_allpcgenes.xlsx")
#vst transformation

vst_counts <- vst(dds)

vst_counts_matrix <- assay(vst_counts)
vst_counts_matrix <- as.data.frame(vst_counts_matrix)

#pca
#pc1-pc2
rv <- rowVars(assay(vst_counts))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]
sample.pca <- prcomp(t(assay(vst_counts)[select,]),
                     center = TRUE,
                     scale. = FALSE)
# Retrieve loadings
new.df <- data.frame(sample.pca$x, "genotype"=samplesplan$genotype, 
                     "stage"=samplesplan$stage, "replicate"=samplesplan$replicate,"sample"=samplesplan$sample, "group"=paste(samplesplan$genotype, ":", samplesplan$stage))
# Evaluate variance explained
var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)
# Plot first 2 PC:
pdf("Plots/FIG_6B_PCAplot_top500_PC1PC2.pdf")
ggplot(new.df, aes(PC1,PC2)) +
    #geom_point( aes(shape=genotype, color=stage), stroke=2, size=1)+
    #geom_text(aes(label=replicate), nudge_x= -2, nudge_y = -2)+
    geom_point( aes(color=group), stroke=2, size=3)+
    ylim(-20,20)+xlim(-50,50)+
    theme_bw() + 
    theme(aspect.ratio = 1)+
    xlab(paste0("PC1: ", var[1],"% variance"))+
    ylab(paste0("PC2: ", var[2],"% variance"))
dev.off()
    
# DEGs are calculated between wt and mut samples to obtain global l2fc for overview

res05 <- results(dds, alpha=0.05)
res05ordered <- res05[order(res05$padj),]
summary(res05)

res05ordered_df <- add_rownames(as.data.frame(res05ordered, var = "GENEID"))

res05ordered_df$genename <- geneIDsprotcod$GENEID[match(res05ordered_df$rowname, geneIDsprotcod$GENENAME)]

res05ordered_df <- res05ordered_df[, c(1,8,2,3,4,5,6,7)]

outputtable_005 <- merge(res05ordered_df, TPM_pc, by.x ="rowname", by.y="gene_name.GENEID")

outputtable_005 <- subset(outputtable_005,padj<0.05 & abs(log2FoldChange)>1.5)

write.table(as.data.frame(outputtable_005), "output/DESeq2_WT_VS_MUT.txt", row.names = F,quote = F, col.names = T,sep = "\t")

openxlsx::write.xlsx(outputtable_005, "output/DESeq2_WT_VS_MUT.xlsx",rowNames=T, overwrite = T)

