###### heatmap inflammation
rm(list = ls())
library(biomaRt)
library(pheatmap)
library(readxl)
library(dplyr)


allpcgenes <- readxl::read_excel("output/WT_VS_MUT_allpcgenes.xlsx")

allpcgenes <- allpcgenes[order(-allpcgenes$log2FoldChange),]

ad_hoc_heatmap_ordered <- function(gene_set, title){
  df.sub <- allpcgenes[allpcgenes$gene_name.GENENAME %in% gene_set,c(8:18)]
  df.sub[,c(2:11)] <- log2(1+df.sub[,c(2:11)])
  df.sub<- df.sub[apply(df.sub[,-1], 1, function(x) !all(x==0)),]
  pdf(paste0("Plots/Heatmap_mine_ordered",title,".pdf"), onefile =T, 
      width= 7.4, height= 5+ 0.22 * length(gene_set))
  breaksListAbs<-c(min(apply(df.sub[,c(2:11)], 1, scale), -1.5) - 0.0001,
                   seq(-1.5, 1.5, 0.03),
                   max(apply(df.sub[,c(2:11)], 1, scale), 1.5) + 0.0001)
  pheatmap(df.sub[,c(2:11)],
           labels_row=pull(df.sub[,1]),
           #labels_row=rep("", nrow(df.sub.il6)),
           cluster_rows= F,
           cluster_cols= FALSE ,
           #kmeans_k = 4,
           breaks = breaksListAbs,
           main=paste0("log2(1+TPM) ", title),
           #clustering_method="complete",
           clustering_method="ward.D2",
           scale="row",
           color = colorRampPalette(c("blue","white","red"))(103),
           fontsize_row = 9,
           cellheight = 16)
  dev.off()
}

genes_inflammation <- c("Cxcl3", "Cxcl2", "S100a9", 
                  "Ccl4", "Lcn2", "Il23a", "Il17f", "Il22ra2", 
                  "Csf2","Il1b","S100a8","Nlrp3","Ccl3","Il17c",
                  "Ccl22","Wdfdc17","Il1f10","Defb3","Tnf", "Zc3h12a","Lgals3")

ad_hoc_heatmap_ordered(genes_inflammation, "Inflammation-related genes")

genes_blastema <- c( "Msx1","Msx2","Mest", "Bmp2",
                                "Bmp4","Arsi","Ltbp2", "Bmp5", "Cxcr4", "Dlx5",
                                "Lhx9", "Serpinf1")

genes_signaling <- c("Rspo4", "Sfrp2", "Mme", "Msx1","Msx2","Crabp1"    )

ad_hoc_heatmap_ordered(genes_blastema, "Blastema genes")

ad_hoc_heatmap_ordered(genes_signaling, "Nail dermis genes")

genes_angiogenesis <- c("Cdh5", "Pecam1","Vcam1","Acta2")

ad_hoc_heatmap_ordered(genes_angiogenesis, "vascular genes")

genes_M2phenotype <- c("Cx3cr1","Mrc1")

genes_antiinflamm <- c("Lgals1", "Siglec1", "Mdk")

ad_hoc_heatmap_ordered(genes_antiinflamm, "Anti-inflammatory")
ad_hoc_heatmap_ordered(genes_M2phenotype, "M2-polarization genes")

