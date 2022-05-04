library(ggplot2)
library(xlsx)
library(ggpubr)
library(ggsignif)

av2 <- read.xlsx("quantHCR.xlsx", sheetIndex = 1, header = T)
av2$condition_genotype <- paste(av2$Condition, av2$Genotype, sep="_")

#Fig 2C hcr uninjured

level_order_2 <- c("No_probe_NA","Dorsal_dermis_MUT","Ventral_dermis_WT","Dorsal_dermis_WT")
ggplot(av2[av2$Condition!="12dpa_blastema",],
       aes(y=Value, x= factor(condition_genotype, level=rev(level_order_2)),
           fill = condition_genotype)) + #now its not necessary to separate by genotype
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.1)) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
  scale_x_discrete(labels=c("Dorsal_dermis_WT" = "WT Dorsal dermis",
                              "Ventral_dermis_WT" = "WT Ventral dermis",
                              "Dorsal_dermis_MUT" = "\u0394LARM1/2 Dorsal dermis",
                              "No_probe_NA"="No probe")) +
  # stat_compare_means(method = "anova", angle= 270, label.x = 4, label.y = 120)+
  geom_signif(comparisons = list(c("Dorsal_dermis_WT","No_probe_NA")),
              map_signif_level = T, textsize = 8, test  ="t.test", y_position = c(110,110))+ #2.6e-05
  # geom_signif(comparisons = list(c("Dorsal_dermis_WT","Ventral_dermis_WT")),
  #             map_signif_level = T, textsize = 8, test  ="t.test", y_position = c(90,90))+#0.00043
  geom_signif(comparisons = list(c("Dorsal_dermis_WT","Dorsal_dermis_MUT")),
              map_signif_level = T, textsize = 8, test  ="t.test", y_position = c(100,100))+ #2e-05
  geom_signif(comparisons = list(c("Dorsal_dermis_MUT","No_probe_NA")),
              map_signif_level = T, textsize = 8, test  ="t.test", y_position = c(60,60))+ #0.1
  # geom_signif(comparisons = list(c("Dorsal_dermis_MUT","Ventral_dermis_WT")),
  #             map_signif_level = T, textsize = 8, test  ="t.test", y_position = c(60,60),)+ #0.00031
  ylab("Number of puncta") +
  coord_cartesian(ylim = c(0, 120))

ggsave("Plots/FIG_3D_HCRquant_uninjured_boxplot.jpg", height = 8, width=5)

#### FIG Supp 3 blastema WT by condition_space 

level_order_blast_SPATIAL <- c("Blastema_up","Blastema_bot","Blastema_mid","Blastema_center")
my_yblast_title <- expression(paste(italic("L1/2-/-")))
Okabe_Ito <- c("#D55E00",  "#56B4E9", "#009E73", "#CC79A7")
ggplot(av2[c(18:65),], 
       aes(y=Value, x= factor(Condition_space, level = level_order_blast_SPATIAL)))+  #now its not necessary to separate by genotype
  
  geom_boxplot(outlier.shape = NA, fill=c(Okabe_Ito)) + 
  geom_point(position = position_jitter(0.1)) +
  theme_classic(base_size = 30) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x=element_text(angle=45, vjust=1, hjust = 1))  +
  ylab("Number of puncta") +
  stat_compare_means(method="anova", label.y=125, size = 6)+
  scale_x_discrete(labels=c("Blastema_bot" = "Ventral (v)", 
                            "Blastema_mid" = "Central (c)",
                            "Blastema_center" = "Distal (di)",
                            "Blastema_up"="Dorsal (do)"))+
  coord_cartesian(ylim=(c(0,130)),)
ggsave("Plots/FIG_S3B_HCRquant_blastema_WTlocalization_boxplot.jpg", height = 7, width = 7)

#### blastema

level_order_blast <- c("12dpa_blastema_WT","12dpa_blastema_MUT")

ggplot(av2[av2$Condition=="12dpa_blastema",], 
       aes(y=Value, x= factor(condition_genotype, level=level_order_blast),
           fill = condition_genotype)) + #now its not necessary to separate by genotype
  geom_boxplot(outlier.shape = NA) +  
  geom_point(position = position_jitter(0.1)) +
  theme_classic(base_size = 30) +
  geom_signif(
    comparisons = list(c("12dpa_blastema_WT","12dpa_blastema_MUT")),
    map_signif_level = TRUE, textsize = 10, test  ="t.test") +
  ylab("Number of puncta") +
  scale_y_continuous(breaks=seq(0,160,20),limits = c(0,130))+
  scale_x_discrete(labels=c("12dpa_blastema_WT" = "WT", 
                            "12dpa_blastema_MUT" = "\u0394LARM1/2"))+
  ggtitle("12DPA blastema ")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        # axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(size=22))
  

ggsave("Plots/FIG_S3G_HCRquant_blastema_boxplot.jpg", height = 7, width =5)


                  