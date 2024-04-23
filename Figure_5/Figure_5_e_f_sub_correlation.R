rm(list=ls())

library(dplyr)
library(tidyr)

# READ IN FILES
setwd("D:/data/03222021_human_psy_qtl/Figure_code")

library(ggplot2)
library(readxl)
library(readr)
###################################################################################
prot <- read_delim("rawdata/QTLtools_PEER_PhenoIn_protein.txt", 
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE)

prot_exp <- prot[,c(4,7:274)]

gene <- read_delim("rawdata/RNA_expression_PEER_PhenoIn.txt", 
                   delim = "\t", escape_double = FALSE,  trim_ws = TRUE)


gene_exp <- gene[,c(4,7:422) ]

names(gene_exp)[1] <- "gid"
##################################################################################################
TRPV2_prot <- subset(prot_exp, pid == "Q9Y5S1")
#
###################
TRPV2_rna <- subset(gene_exp, gid == "ENSG00000187688")
######################################################################
TRPV2_prot2 <- as.data.frame(t(TRPV2_prot[,2:269]))
#
TRPV2_prot2$sample  <- rownames(TRPV2_prot2)
#
TRPV2_rna2 <- as.data.frame(t(TRPV2_rna[,2:417]))
#
TRPV2_rna2$sample  <- rownames(TRPV2_rna2)
#
TRPV2_exp <- merge(TRPV2_prot2, TRPV2_rna2, by = "sample")
#
names(TRPV2_exp) <- c("sample","prot","rna" )
#

ggplot(data = TRPV2_exp, aes(x = prot, y = rna)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black", fill = "lightgray")+
  ylim(-1,1)+ xlim(-0.6,0.6) +
  theme_classic() 
#
#dev.off()


res <- cor.test(TRPV2_exp$prot, TRPV2_exp$rna, method = "pearson")
#
res
#####################################################################
##################################################################################################
GLRX5_prot <- subset(prot_exp, pid == "Q86SX6")
#
###################
GLRX5_rna <- subset(gene_exp, gid == "ENSG00000182512")
######################################################################
GLRX5_prot2 <- as.data.frame(t(GLRX5_prot[,2:269]))
#
GLRX5_prot2$sample  <- rownames(GLRX5_prot2)
#
GLRX5_rna2 <- as.data.frame(t(GLRX5_rna[,2:417]))
#
GLRX5_rna2$sample  <- rownames(GLRX5_rna2)
#
GLRX5_exp <- merge(GLRX5_prot2, GLRX5_rna2, by = "sample")
#
names(GLRX5_exp) <- c("sample","prot","rna" )

#

ggplot(data = GLRX5_exp, aes(x = prot, y = rna)) + 
  geom_point() + 
  xlim(-0.5,0.5) + ylim(-0.3,0.3) +
  geom_smooth(method = "lm", color = "black", fill = "lightgray")+
  theme_classic() 
#
#dev.off()


res <- cor.test(GLRX5_exp$prot, GLRX5_exp$rna, method = "pearson")
#
res
#####################################################################
#############################################