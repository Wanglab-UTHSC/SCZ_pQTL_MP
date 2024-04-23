#install.packages('ggtern')
#install.packages("readr")

rm(list=ls())
#
library(ggplot2)
library(tidyverse)
library(ggtern)

setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
###############################################################################################
# coloc plot for pQTL and eQTL
library(readxl)

pqtl_eqtl <- read_excel("rawdata/coloc_results_01312022.xlsx", sheet = "coloc_pQTL")
#
pqtl_eqtl$H0_H1_H2 <- (pqtl_eqtl$PP.H0.abf+pqtl_eqtl$PP.H1.abf+pqtl_eqtl$PP.H2.abf)
#
head(pqtl_eqtl)
#
data1 <- pqtl_eqtl[,c(8:11)]
colnames(data1) <- c("H3","H4","max","H0_H1_H2")
#

#pdf("Ternary_plot_pqtl_eqtl_database_v1.pdf")
#
ggtern(data = data1, 
       aes(x = H0_H1_H2, y = H4, z = H3)) +  
  geom_point(aes(color = 1-H4), size = 2) +
  scale_color_gradient2(low = "#BA0C0C", high = "lightgrey", midpoint = 0.5) + 
  theme_classic() +
  theme_showarrows() + 
  labs(x = "H1+H2", xarrow = "H1+H2",
       y = "H4",     yarrow = "H4",
       z = "H3",       zarrow = "H3") + 
  theme_latex()


