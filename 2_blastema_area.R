###### mac2 in blastema
library(xlsx)
library(ggplot2)
library(ggsignif)
library(dplyr)

blastema_area <- read.xlsx("blastema_area.xlsx", header = T, sheetIndex = 1)
level_order <- c("WT","MUT")
exponent <- expression("Blastema area in"~mm^2 )
mean(blastema_area[blastema_area$genotype=="WT",]$Areaµm2)/mean(blastema_area[blastema_area$genotype=="MUT",8],na.rm = T)
#0.5354588
ggplot(blastema_area, aes(y=Areaµm2/1000, x=factor(genotype, level=level_order))) +
  geom_boxplot(aes(fill=genotype)) + 
  geom_point(aes(shape=genotype), position=position_jitter(0.1)) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank())+
  coord_cartesian(ylim=c(0,400))+
  ylab(exponent)+
  geom_signif(comparisons=list(c("WT","MUT")),test = "t.test", map_signif_level = F,textsize = 6,
              y_position = 350)+
  scale_x_discrete(labels=c("WT" = "WT", 
                            "MUT" = "\u0394LARM1/2"))
ggsave("Plots/FIG_5B_Blastema_area_quant.jpg", width = 5, height = 5)