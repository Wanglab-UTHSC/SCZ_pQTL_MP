rm(list=ls())
memory.limit(1000000)

library(stringr)
library(dplyr)
library(tidyr)
library(qvalue)
# READ IN FILES
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
library(eoffice)
library(qvalue)
library(ggplot2)
library(readxl)
library(readr)
#

###############################################################################################################
data4 <- read.delim("rawdata/step_3_forest_rawdata_filter.txt")
#
data4$p_SMR <- ifelse(data4$p_SMR < 1E-10, 1E-10, data4$p_SMR)  #NAֵ??0?滻
#
data4$idnew <- factor(data4$idnew,levels=c(1:23))#

ggplot(data4, aes(b_SMR, idnew, color = group)) +
  
  geom_point(size=3.6) +
  
  geom_errorbarh(aes(xmax = up_smr, xmin = low_smr), height = 0) +
  xlim(-2,2) + 
  
  theme_classic()+   geom_text(aes(label = rownames))



data4$p_GWAS <- ifelse(data4$p_GWAS < 1E-20, 1E-20, data4$p_GWAS)  #NAֵ??0?滻
#
#

ggplot(data4, aes(b_GWAS, idnew,  col = group)) +
  
  geom_point(size=3.6) +
  xlim(-0.3,0.3)+
  
  geom_errorbarh(aes(xmax = up_gwas, xmin = low_gwas), height = 0) +
  
  theme_classic()+   geom_text(aes(label = rownames)) 
  




ggplot(data4, aes(b_eQTL, idnew,  col = group)) +
  
  geom_point(size= 2) + xlim(-2,2)+
  
  geom_errorbarh(aes(xmax = up_qtl, xmin = low_qtl), height = 0.1) +
  
  theme_classic()+   geom_text(aes(label = rownames)) 
