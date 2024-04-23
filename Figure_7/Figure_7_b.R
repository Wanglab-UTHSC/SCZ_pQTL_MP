rm(list=ls())
#
library(readxl)
#
library(biomaRt)
library(dplyr)

library(ggplot2)
library(corrplot) 

library(FactoMineR)
library(ggrepel)
library(qvalue)
library(readr)
##
#
setwd("D:/data/03222021_human_psy_qtl/Figure_code")
#
rank <- read.delim("rawdata/step_3_protein_QTL_with_scz_GWAS_rank.txt")
#
rank_sub <- subset(rank, final_rank_SCZ < 61)
#
names(rank_sub)
#
sub1 <- rank_sub[,c(1:3,24)]
#
names(sub1) <- c("gn","p","rank","index")
#
sub1$group <- "pqtl"
#
sub2 <- rank_sub[,c(1,4,5,24)]
#
names(sub2) <- c("gn","p","rank","index")
#
sub2$group <- "eqtl"
#
sub3 <- rank_sub[,c(1,6,7,24)]
#
names(sub3) <- c("gn","p","rank","index")
#
sub3$group <- "coloc"
#
sub4 <- rank_sub[,c(1,10,11,24)]
#
names(sub4) <- c("gn","p","rank","index")
#
sub4$group <- "score"
#
sub5 <- rank_sub[,c(1,17,18,24)]
#
names(sub5) <- c("gn","p","rank","index")
#
sub5$group <- "SCZ"

rownames(sub1) <- sub1$gn
#
sub1$rank <- ifelse(sub1$rank > 500 , 500, sub1$rank) 
sub1$index <- factor(sub1$index, levels = 1:60)
#



##############################################################################
library(pheatmap)
#
mat <- as.data.frame(rank_sub[,c(1,3,5,7,11,18,24)])
#
rownames(mat) <- rank_sub$gn
#
prot <- read_excel("rawdata/step_2_protein_QTL_with_scz_GWAS_v2_624.xlsx",  sheet = "unq")

mat2 <- merge(mat, prot, by.x = "gn", by.y = "gid", all.x = T)

mat3_0 <- as.data.frame(mat2[order(mat2[,7]),])
#
mat3_1 <- unique(mat3_0[,1:7])
#
mat3 <- unique(mat3_1[,2:6])
#
rownames(mat3) <- mat3_1$gn

mat3$rank_list_pQTL <- ifelse(mat3$rank_list_pQTL > 500 , 500, mat3$rank_list_pQTL) 
#
mat3$rank_eQTL <- ifelse(mat3$rank_eQTL > 500 , 500, mat3$rank_eQTL) 
#
mat3$rank_list_coloc <- ifelse(mat3$rank_list_coloc > 500 , 500, mat3$rank_list_coloc) 
#
mat3$rank_list_score <- ifelse(mat3$rank_list_score > 500 , 500, mat3$rank_list_score) 
#
mat3$rank_list_net_SCZ <- ifelse(mat3$rank_list_net_SCZ > 500 , 500, mat3$rank_list_net_SCZ) 

colors <- colorRampPalette(c("#0000FF", "lightgrey"))(100)
#
annotation <- data.frame(mat3_1$final_rank_SCZ)
rownames(annotation) <- rownames(mat3) # check out the row names of annotation

pheatmap(mat3, cluster_rows = F, scale = "none", cluster_cols = F, annotation_row = annotation, 
         color = colors, na_col = "white",cellwidth = 10, cellheight = 8, fontsize_col = 10)
#






