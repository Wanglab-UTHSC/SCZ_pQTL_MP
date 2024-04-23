rm(list=ls())

# READ IN FILES
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")

library(ggplot2)
library(readxl)
library(readr)
###################################################################################
gene_raw <- read.delim("rawdata/RNA_expression_PEER_PhenoIn.txt")
#
############################################################################################
prot_raw <- read.delim("rawdata/QTLtools_PEER_PhenoIn_protein.txt")
#
##################################################################################################
SRR_prot <- subset(prot_raw, pid == "Q9GZT4")
#
###################
SRR_rna <- subset(gene_raw, gene == "ENSG00000167720")
######################################################################
SRR_prot2 <- as.data.frame(t(SRR_prot[,7:274]))
#
names(SRR_prot2) <- "prot_exp"
#
SRR_prot2$sample  <- rownames(SRR_prot2)
#
SRR_rna2 <- as.data.frame(t(SRR_rna[,7:422]))
#
names(SRR_rna2) <- "gene_exp"
#
SRR_rna2$sample  <- rownames(SRR_rna2)
#
SRR_exp <- merge(SRR_prot2, SRR_rna2, by = "sample")
#
names(SRR_exp) <- c("sample","prot","rna" )
#


res <- cor.test(SRR_exp$prot, SRR_exp$rna, method = "pearson")
#
res

SRR_exp$sample <- str_replace(SRR_exp$sample, "X", "")
SRR_exp$sample <- str_replace(SRR_exp$sample, "\\.", "-")


#####################################################################
geno_prot <- geno_prot <- read.delim("rawdata/srr_prot_geno.txt", header=FALSE)
#
head_prot <- read_table2("rawdata/geno_prot_head.txt", skip = 41)
#
names(geno_prot) <- names(head_prot)
#
merge_geno_long <- as.data.frame(t(geno_prot[,10:277]))
#
colnames(merge_geno_long) <- c( "geno")
#
merge_geno_long$sample <- rownames(merge_geno_long)

merged_all <- merge(merge_geno_long, SRR_exp, by = "sample")

merged_all1 <- subset(merged_all, geno != "./.")
#
table(merged_all1$geno)

ggplot(merged_all1, aes(x= geno , y= as.numeric(prot ))) +
  #geom_violin(trim=FALSE,aes(fill=factor(geno))) +  
  geom_boxplot()+ 
  theme_bw()+ 
  geom_jitter() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))


ggplot(merged_all1, aes(x= geno , y= as.numeric(rna ))) +
  #geom_violin(trim=FALSE,aes(fill=factor(geno))) +  
  geom_boxplot()+ 
  theme_bw()+ 
  geom_jitter() +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))

