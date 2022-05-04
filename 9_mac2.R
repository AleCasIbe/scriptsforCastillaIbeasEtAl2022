###### mac2 in blastema
library(xlsx)
library(ggplot2)
library(ggsignif)
library(dplyr)

mac2measures <- read.xlsx("mesaurecementswithmore.xlsx", header = T, sheetIndex = 1)

level_order <- c("WT","MUT")
exponent2 <- expression("Positive cells per"~mm^2 )
  

ggplot(mac2measures, aes(y=Num.Positive.per.mm.2, x=factor(genotype, level=level_order))) +
  geom_boxplot(aes(fill=genotype)) + 
  geom_point( position=position_jitter(0.1)) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank())+coord_cartesian(ylim=c(0,900))+
  ylab(exponent2)+
  geom_signif(comparisons=list(c("WT","MUT")),test = "t.test", map_signif_level = F,textsize = 6,
              y_position = 800)+
  scale_x_discrete(labels=c("WT" = "WT", 
                            "MUT" = "\u0394LARM1/2"))
ggsave("Plots/Mac2quant.jpg", width = 6, height = 5)
mean(mac2measures[mac2measures$genotype=="MUT",]$Num.Positive.per.mm.2,na.rm = T)/mean(mac2measures[mac2measures$genotype=="WT",]$Num.Positive.per.mm.2,na.rm = T)
#1.867557
