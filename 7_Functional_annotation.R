library(WebGestaltR)

######### 12v14 intersection ######

results12vs12 <- read.table("output/DESeq2_12vs12.txt", header = T)

results12vs12up <- results12vs12[results12vs12$log2FoldChange>0,]
results12vs12dn <- results12vs12[results12vs12$log2FoldChange<0,]

results14vs14 <- read.table("output/DESeq2_14vs14.txt", header = T)

results14vs14up <- results14vs14[results14vs14$log2FoldChange>0,]
results14vs14dn <- results14vs14[results14vs14$log2FoldChange<0,]

resultswtvsmut <- read.table("output/DESeq2_WT_VS_MUT.txt",header = T)

#check no genes are up/down in one and the other way around in the other 
intersect(results12vs12dn$genename, results14vs14up$genename) ##character (0)
intersect(results12vs12up$genename, results14vs14dn$genename) ##character (0)

pdf("Plots/Venn_WT.pdf")
pwt <- plot(eulerr::euler(list("12DPA-DEGs WT"=results12vs12up$genename,
                               "14DPA-DEGs WT"=results14vs14up$genename), 
                          shape="ellipse", legend = T), #labels are manually added with Ps
            fills = list(fill = c( "#00BFC4","#C77CFF")), labels=F)

print(pwt)
dev.off()

pdf("Plots/Venn_MUT.pdf")
pmut <- plot(eulerr::euler(list("12DPA-DEGs MUT"=results12vs12dn$genename,
                                "14DPA-DEGs MUT"=results14vs14dn$genename)),  #labels are manually added with Ps
             fills = list(fill = c( "#F8766D","#7CAE00")),
             legend = F, labels=F)
print(pmut)
dev.off()


#log2fc will be the from the global comparison
res <- readRDS("output/res.RDS")
res <- as.data.frame(res)
res <- tibble::rownames_to_column(res, "geneID")
res <- merge(res, geneIDsprotcod, by.x="geneID", by.y="GENEID")

common12v14 <- res[res$GENENAME %in% intersect(results12vs12$genename,
                                               results14vs14$genename),]

write.xlsx(common12v14, "output/TableS3_common_DEGs_L2fc_From_global.xlsx")

common12v14up <-  common12v14[common12v14$log2FoldChange>0,]
common12v14dn <-  common12v14[common12v14$log2FoldChange<0,]

nrow(common12v14up)+nrow(common12v14dn)

####### webgestalt wrapper #####

ORA_byAleCasIbe <- function(genelist, projectname, genotype){
  set.seed(123)
  genesets <- WebGestaltR::listGeneSet()
  list_a <- genelist$geneID
  dftemp <- WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase = as.vector(genesets[c(2,4,6,7,8,9,10),]$name),
    enrichDatabaseFile = NULL,
    enrichDatabaseType = "ensembl_gene_id",
    enrichDatabaseDescriptionFile = NULL,
    interestGeneFile = NULL,
    interestGene = list_a,
    interestGeneType = "ensembl_gene_id",
    collapseMethod = "mean",
    referenceGeneFile = NULL,
    referenceGene = NULL,
    referenceGeneType = "ensembl_gene_id",
    referenceSet = "genome_protein-coding",
    minNum = 2,maxNum = 500,
    sigMethod = "fdr",
    fdrMethod = "BH",
    fdrThr = 0.05,
    topThr = 10,
    reportNum = 20,
    perNum = 1000,
    gseaP = 1,
    isOutput = TRUE,
    outputDirectory ="output/WebgestaltR",
    projectName = paste0(projectname),
    dagColor = "continuous",
    saveRawGseaResult = FALSE,
    gseaPlotFormat = c("png", "svg"),
    setCoverNum = 10,
    networkConstructionMethod = NULL,
    neighborNum = 10,
    highlightType = "Seeds",
    highlightSeedNum = 10,
    nThreads = 1,
    cache = "output/WebgestaltR/cache1",
    hostName = "http://www.webgestalt.org/")
  for (i in 1:nrow(dftemp)){
    codes <- strsplit(dftemp[i,12], ";")[[1]]
    codes <- ensembldb::select(EnsDb.Mmusculus.v79, keys=codes, keytype = "GENEID", columns = c("GENENAME"))
    dftemp$genesymbols[i] <- toString(codes[,2],)
  }
  writexl::write_xlsx(dftemp, paste0("output/",projectname,".xlsx"))
  if (genotype == "MUT") {
    colors_ <- c(high="orange4", low="orange")
  } else {
    colors_ <- c(high="royalblue4",low="royalblue")
  }
  
  dftemp$log10FDR <- -log10(dftemp$FDR)
  
  pdf(paste0("Plots/Webgestalt_", projectname, "_allterms.pdf"))
  print(ggplot(data = head(dftemp, 20), aes(y=reorder(description, -log10FDR), 
                                            x=log10FDR, fill=log10FDR))+
          geom_bar(stat = "identity")+
          scale_y_discrete(labels = wrap_format(40))+
          theme(axis.text.x=element_text(angle=90,  vjust = 0.3,hjust=1), 
                axis.text.y = element_text(angle=90,  vjust = 0.5,hjust=0.5),
                axis.title.x =   element_blank(), panel.background = element_blank(),
                legend.title = element_blank(),text = element_text(size=10),
                axis.ticks.x = element_blank(), )+
          scale_fill_gradient(high=colors_[[1]], low = colors_[[2]], oob=squish_infinite)+
          xlab("-log10(q-value)")+
          geom_vline(xintercept=1.30103,  col = "red",lty=2)+coord_flip())
  dev.off()
  
}

ORA_byAleCasIbe(common12v14up, "TableS5_12v14commongenes_MUT","MUT")
ORA_byAleCasIbe(common12v14dn, "TableS4_12v14commongenes_WT","WT")

#figure selection 

commonWT <- xlsx::read.xlsx("output/webgestaltR12v14commongenes_WT.xlsx", sheetIndex = 1)

# need 

commonWT[c(5, 8, 18, 19, 27, 30, 33, 35,64, 77), 2]
commonWT$log10FDR <- -log10(commonWT$FDR)

## dotplot

pdf(paste0("Plots/Webgestalt_", "selectionWT_mainfigure", "dots_allterms.pdf"), height = 8, width = 9)
print(ggplot(data = commonWT[c(5, 8, 18, 19, 27, 30, 33, 35,64, 77), ], aes(y=reorder(description, log10FDR), 
                                                                            x=overlap/size))+
        geom_point(stat = "identity", color="black", pch= 21, aes(size=overlap, fill=log10FDR))+
        scale_y_discrete(labels = wrap_format(40))+
        scale_x_continuous(position="top")+
        theme_bw()+
        theme(axis.title.y  = element_blank(),
              axis.text.y = element_text(size = 23))+
        scale_fill_gradient(low="lightblue",high =  "blue")+
        labs(x="Gene Ratio")+
        scale_size_area())

dev.off()

commonMUT <- xlsx::read.xlsx("output/webgestaltR12v14commongenes_MUT.xlsx", sheetIndex = 1)

#need 

commonMUT[c(2,4,8,11,15,25,30,36,46,51), 2]
commonMUT$log10FDR <- -log10(commonMUT$FDR)

## dotplot

pdf(paste0("Plots/Webgestalt_", "selectionMUT_mainfigure", "dots_allterms.pdf"),  height = 8, width = 9)
print(ggplot(data = commonMUT[c(2,4,8,11,15,25,30,36,46,51), ], aes(y=reorder(description, log10FDR), 
                                                                    x=overlap/size))+
        geom_point(stat = "identity", color="black", pch= 21, aes(size=overlap, fill=log10FDR))+
        scale_y_discrete(labels = wrap_format(40))+
        scale_x_continuous(position="top")+
        theme_bw()+
        theme(axis.title.y  = element_blank(),
              axis.text.y = element_text(size = 23))+
        labs(x="Gene Ratio")+
        scale_fill_gradient(low="yellow",high =  "orange")+
        scale_size_area())
dev.off()