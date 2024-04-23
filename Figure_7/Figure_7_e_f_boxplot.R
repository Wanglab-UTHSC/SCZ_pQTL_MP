rm(list=ls())


###################################################################################
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
PPP2R4_prot <- subset(prot_exp, pid == "Q15257")
#
###################
PPP2R4_rna <- subset(gene_exp, gid == "ENSG00000119383")

######################################################################
PPP2R4_prot2 <- as.data.frame(t(PPP2R4_prot[,2:269]))
#
PPP2R4_prot2$sample  <- rownames(PPP2R4_prot2)
#
PPP2R4_rna2 <- as.data.frame(t(PPP2R4_rna[,2:417]))
#
PPP2R4_rna2$sample  <- rownames(PPP2R4_rna2)
#
PPP2R4_exp <- merge(PPP2R4_prot2, PPP2R4_rna2, by = "sample")
#
names(PPP2R4_exp) <- c("sample","prot","rna" )
#


res <- cor.test(PPP2R4_exp$prot, PPP2R4_exp$rna, method = "pearson")
#
res
#####################################################################
geno_prot <- geno_prot <- read.delim("rawdata/PPP2R4_prot_geno.txt", header=FALSE)
#
head_prot <- read_table2("rawdata/geno_head.txt", skip = 41)
#
names(geno_prot) <- names(head_prot)
#
merge_geno_long <- as.data.frame(t(geno_prot[,10:277]))
#
colnames(merge_geno_long) <- c( "geno")
#
merge_geno_long$sample <- rownames(merge_geno_long)
#
merged_all <- merge(merge_geno_long, PPP2R4_prot2, by = "sample")

merged_all1 <- subset(merged_all, geno != "./.")
#
table(merged_all1$geno)
#
ggplot(merged_all1, aes(x= geno , y= as.numeric(V1))) +
  #geom_violin(trim=FALSE,aes(fill=factor(geno))) +  
  geom_boxplot()+ 
  theme_bw()+ 
  geom_jitter() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
#
table(merged_all1$geno)
#
geno1 <- subset(merged_all1, geno == "0/0")

mean1 <- mean(geno1$V1)
#
geno3 <- subset(merged_all1, geno == "1/1")
#
mean3 <- mean(geno3$V1)
#
lgfc <- mean1-mean3
                   
#################################################################################################
####################################################
#
head_gene <- read_table2("rawdata/head_gene_geno.txt", skip = 41)
#
gene_geno <- read.delim("rawdata/PPP2R4_gene_geno.txt", header=FALSE)
#
names(gene_geno) <- names(head_gene)
#
geno_long <- as.data.frame(t(gene_geno[,10:424]))
names(geno_long) <- "geno"
#
geno_long$sample <- rownames(geno_long)
#
merged_gene <- merge(PPP2R4_rna2, geno_long, by = "sample")


ggplot(merged_gene, aes(x= geno , y= as.numeric(V1))) +
  #geom_violin(trim=FALSE,aes(fill=factor(geno))) +  
  geom_boxplot()+ 
  theme_bw()+ 
  geom_jitter() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))

######
table(merged_gene$geno)
#
geno1 <- subset(merged_gene, geno == "0/0")

mean1 <- mean(geno1$V1)
#
geno3 <- subset(merged_gene, geno == "1/1")
#
mean3 <- mean(geno3$V1)
#
lgfc <- mean1-mean3


