rm(list=ls())
#
require(dplyr, quietly = T, warn.conflicts = F)
require(ggplot2, quietly = T, warn.conflicts = F)
require(readxl, quietly = T, warn.conflicts = F)
require(tidyr, quietly = T, warn.conflicts = F)
require(RColorBrewer, quietly = T, warn.conflicts = F)
require(tidyverse, quietly = T, warn.conflicts = F)
#
require(limma, quietly = T, warn.conflicts = F)
require(org.hm.eg.db, quietly = T, warn.conflicts = F)
#
require(biomaRt, quietly = T, warn.conflicts = F)
#
library(substr)

setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
## The jump output with all proteins across all batches
library(readr)
###########################################
#input the table from quan analysis(rawdata from JUMP q)
# rawdata 
prot_rawdata  <- read_delim("rawdata/combined_norm_uni_prot.txt", delim = "\t", escape_double = FALSE,   trim_ws = TRUE, skip = 61)
#
colnames(prot_rawdata)[1:4] <- c("Protein_Group","Accession","Description","GN")
################################################
## Uniprot ids mapped to ensembl ids
## RNA expression all samples
rna_all <- read_delim("rawdata/RNA_expression_PEER_PhenoIn.txt", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)[c(4,7:422)]
#
names(rna_all)[1] <- "ensembl"
#



prot_nona2 <- as.data.frame(prot_rawdata[, c(which(grepl(pattern = "sig", x = colnames(prot_rawdata))))])
#
rownames(prot_nona2) <- prot_rawdata$Accession
#
prot_nona2$pid <- rownames(prot_nona2)
#
prot_nona3 <- na.omit(prot_nona2)
#
prot_nona3 <- log2(prot_nona3[,1:319])
#
prot_nona3$pid <- rownames(prot_nona3)
#
prot_all <- as.data.frame(prot_nona2)
#
names(prot_all)[1] <- "pid"
#
prot_all$prot_all <- "prot_all"
#
prot_nona<- as.data.frame(prot_nona3$pid)
#
names(prot_nona) <- "pid"
#
prot_nona$prot_nona <- "prot_nona"
####################################################################################
# A total of 260 samples appeared in all three datasets
### Get ensembl ids for protein datasets

library(readr)

id <- read_delim("rawdata/pid_to_ensembl_unique.txt", 
                 delim = "\t", escape_double = FALSE,  trim_ws = TRUE)


prot_all2 <- prot_all[,c(320:321)]
#
prot_all2 <- prot_all2 %>%
  separate(pid, into = c("id1","uniprot","id2", sep = "\\|"))
#
prot_all3 <- prot_all2[,c(2,5)]

prot_nona <- prot_nona %>%
  separate(pid, into = c("id1","uniprot","id2", sep = "\\|"))
#
prot_nona2 <- prot_nona[,c(2,5)]

prot_all4 <- merge(prot_all3, id,  by.x = "uniprot", by.y = "From")
#
prot_all5 <- unique(prot_all4 )
#
prot_nona3 <- merge( prot_nona2, id,  by.x = "uniprot", by.y = "From" )
#
prot_nona4 <- unique(prot_nona3[,c(2,3)])
#
### Get mean tpm, mark mRNA with protein coverage and plot

g <- rna_all

g$mean <- apply(g[,3:ncol(g)], 1, mean)
#
g2 <- g[,c(1,419)]
#
names(g2)[1] <- "ensembl"
#
g3 <- merge(g2, prot_all5, by.x = "ensembl", by.y = "To" , all.x = T)
#
g4 <- merge(g3, prot_nona4, by.x = "ensembl", by.y = "To" , all.x = T)
#
##############################################################################
names(g4)

ggplot()+
  geom_histogram(data = g4[which(g4$mean > 0),], aes(x = mean, fill = "Transcript",color = "grey" )) +
  geom_histogram(data = g4[which(g4$mean > 0 & g4$prot_all == "prot_all"),], 
                 aes(x = mean, fill = "Proteins"), color = "lightblue", alpha = 0.8) +
  geom_histogram(data = g4[which(g4$mean > 0 & g4$prot_nona == "prot_nona"),], 
                 aes(x = mean, fill = "Proteins (all samples)"), color = "blue")+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  labs(x = "RNA expression Log2(tpm)", y = "Frequency")



