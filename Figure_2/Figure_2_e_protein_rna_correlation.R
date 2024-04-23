#install.packages("ggpointdensity")


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
#
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
## The jump output with all proteins across all batches
library(readr)
library("ggpointdensity")
library(viridis)

###########################################
#input the table from quan analysis(rawdata from JUMP q)
# rawdata 
prot_rawdata  <- read_delim("rawdata/combined_norm_uni_prot.txt", delim = "\t", escape_double = FALSE,   trim_ws = TRUE, skip = 61)
#
colnames(prot_rawdata)[1:4] <- c("Protein_Group","Accession","Description","GN")
################################################
## Uniprot ids mapped to ensembl ids
map <- read.delim("rawdata/pid_to_ensembl_unique.txt")
#
colnames(map) <- c("uniprot","ensembl")
## RNA expression all samples
rna_all <- read_excel("rawdata/Supplementary_Table_v2.xlsx",sheet = "Table S3", skip = 2)
#
names(rna_all)[1] <- "ensembl"
#######################################################################################
prot_nona2 <- as.data.frame(prot_rawdata[, c(which(grepl(pattern = "sig", x = colnames(prot_rawdata))))])
#
rownames(prot_nona2) <- prot_rawdata$Accession

prot_nona3 <- na.omit(prot_nona2)
#
prot_nona3 <- log2(prot_nona3)
#

colors <- RColorBrewer::brewer.pal(n = 7, name = "Blues")

####################################################################################
#PROTEIN COVERAGE OF GENES
## Get mean log2(tpm) for all genes
ra <- rna_all
ra$mean <- apply(ra[, 7:ncol(ra)], 1, mean)
#
ra$rna <- 1
#
## Get ensembl ids mapped to whole protein dataset and reduced protein dataset
## Whole protein dataset
p_all <- prot_rawdata[1:4] %>% separate("Accession", into = c("db", "uniprot","id2"), sep = "\\|") %>% 
  left_join(map, by = "uniprot") 
#
p_all1 <- p_all[c(3,6,7)]
#


#######################p_all1$prot_all <- 1
# Reduced protein dataset
#
prot_nona3$id <- rownames(prot_nona3)
#
p <- prot_nona3[320]  %>% separate("id", into = c("db", "uniprot","id2"), sep = "\\|") %>% 
  left_join(map, by = "uniprot") 
#
p1 <- p[,c(2,4)]
#
p1$prot_nona <- 1
###################################################

#################################################
#PROTEIN VS MRNA EXPRESSION CORRELATION PLOT
## Use reduced protein dataset and all samples RNAseq

prot <-  prot_nona3 %>% separate("id", into = c("db1", "uniprot","id2"), sep = "\\|") %>% 
  left_join(map, by = "uniprot")
## Warning: Expected 2 pieces. Additional pieces discarded in 8231 rows [1, 2, 3,
## 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
rna <- rna_all

## Get mean tpm value and mean protein expression value
rna$rna_mean <- apply(rna[, 3:418], 1, mean)
prot$prot_mean <- apply(prot[, 1:319], 1, mean)

## Join RNA and protein
pr <- na.omit(left_join(prot, rna, by = "ensembl"))

r <- round(cor(x = pr$prot_mean, y = pr$rna_mean), 3)

cor.test(x = pr$prot_mean,y = pr$rna_mean)
# plot

ggplot(data = pr, mapping = aes(x = prot_mean, y = rna_mean)) + 
  geom_pointdensity() + 
  scale_color_viridis() + 
  geom_smooth(method = lm) +  ##ʡ??????????
  scale_x_continuous(breaks = seq(12,28,4)) +
  scale_y_continuous(breaks = seq(0,12,4)) +
  labs(x = "Protein expression(log2(intensity))",  y = "RNA expression(log2(tpm))", title="r = 0.43,  p < 2.2e-16") +
  theme_classic()
 

