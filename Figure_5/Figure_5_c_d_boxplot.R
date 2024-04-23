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
library(qvalue)
library(readr)
##
#
setwd("D:/data/03222021_human_psy_qtl/Figure_code")
#################################################################################################
mediation <- read_excel("rawdata/new_example_mediation.xlsx", skip = 1)
#
names(mediation)[6] <- "effect_p" 

names(mediation)[9] <- "effect_g" 
#
names(mediation)[13] <- "group" 
#
names(mediation)[3] <- "pid" 


names(mediation)[4] <- "gid" 

mediation$id <- paste(mediation$pid, mediation$gid, sep = "_")
##############################################################################################
# Figure 5D
############################################################################################
gene <- mediation[,c(9,13)]

names(gene) <- c("effect", "group" )


prot <- mediation[,c(6,13)]

names(prot) <- c("effect", "group" )


data3 <- rbind(prot, gene)

ggplot(data= data3, aes( x = as.factor(group), y = abs(effect) ))  +
  geom_boxplot( ) +
  geom_jitter() + ylim(0,0.2)+
  theme_classic()




###############################################################################
#Figure 5C

###################################################################################
library(readr)
#
id <- mediation[,c(3,4,17,13)]


prot <- read_delim("rawdata/QTLtools_PEER_PhenoIn_protein.txt", 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

prot_exp <- prot[,c(4,7:274)]

gene <- read_delim("rawdata/RNA_expression_PEER_PhenoIn.txt", 
                    delim = "\t", escape_double = FALSE,  trim_ws = TRUE)


gene_exp <- gene[,c(4,7:422) ]

names(gene_exp)[1] <- "gid"

prot_sub <- merge(prot_exp, id, by.x = "pid")
#
prot_sub2 <- unique(prot_sub)

prot_sub3 <- prot_sub2[,c(2:269)]
#
rownames(prot_sub3) <- prot_sub2$id
#

gene_sub <- merge(gene_exp, id, by = "gid")
#
gene_sub2 <- unique(gene_sub)
#
gene_sub3 <- gene_sub2[,c(2:417)]
#
rownames(gene_sub3) <- gene_sub2$id
#
sample1 <- names(prot_sub3)
#
sample2 <- names(gene_sub3)
#
sample3 <- intersect(sample1, sample2)

id_list <- rownames(prot_sub3)

gene_sub4 <- gene_sub3[id_list, sample3]

prot_sub4 <- prot_sub3[id_list, sample3]


corr <- as.data.frame(gene_sub3[1,1])
#
for(i in 1:386) { 
  # Head of for-loop
  new <- cor(t(prot_sub4[i,]),t(gene_sub4[i,]))
  corr <- cbind(corr, new)
                                                 
}



#############################################################
corr_sub <- as.data.frame(t(corr[,2:387]))
#
names(corr_sub) <- "corr"
#
corr_sub$id <- rownames(corr_sub)
#
names(corr_sub)
#
merged_cor <- merge(id, corr_sub, by = "id")
#
library(MASS)
library(ggpubr)
#
ggplot(data = merged_cor, aes(y =abs(corr), x = as.factor(group))) + 
  geom_boxplot(width=0.2) +
  ylim(0,1) +
  theme_classic() + stat_compare_means() +
  geom_jitter(height = 0, width = 0.2)
#
table(merged_cor$group)




