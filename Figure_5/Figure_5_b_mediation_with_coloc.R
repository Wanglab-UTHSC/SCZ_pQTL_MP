rm(list=ls())
#
library(readxl)
#
library(tidyr)
library(dplyr)
library(ggplot2)
library(corrplot) 
library(ggord)
library(FactoMineR)
library(ggrepel)
library(readr)
##
#
setwd("D:/data/03222021_human_psy_qtl/Figure_code")
#
library(readxl)
#
mediation <- read_excel("rawdata/new_example_mediation.xlsx", skip = 1)
#
names(mediation)[5] <- "p_p" 

names(mediation)[8] <- "p_g" 
#
names(mediation)[11] <- "zscore" 
#
names(mediation)[2] <- "gn" 

library(FactoMineR)
library(ggrepel)
library(qvalue)
library(readr)

merged_lable <- subset(mediation, (zscore < -20 & p_p < 1E-30 ) | gn == "GLRX5"  | gn == "VWA8"  | gn == "PNPO") 
#
library(ggplot2)

ggplot(data = mediation, aes ( x = -log10(p_p), y = -log10(p_g), col = coloc)) + geom_point()+
  theme_classic() +
  xlim(0, 65)+ylim(0,65)+
  geom_text_repel(data = merged_lable,aes ( x = -log10(p_p), y = -log10(p_g), label = gn),
                  size = 3, box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", max.overlaps = 50,
                  show.legend = T) +
  theme_classic()

