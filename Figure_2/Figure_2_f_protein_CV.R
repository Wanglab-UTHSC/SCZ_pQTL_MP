rm(list=ls())
#
require(dplyr, quietly = T, warn.conflicts = F)
require(ggplot2, quietly = T, warn.conflicts = F)
require(readxl, quietly = T, warn.conflicts = F)
require(tidyr, quietly = T, warn.conflicts = F)
require(RColorBrewer, quietly = T, warn.conflicts = F)
require(tidyverse, quietly = T, warn.conflicts = F)
#
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
map <- read.delim("rawdata/pid_to_ensembl_unique.txt")
#
colnames(map) <- c("uniprot","ensembl")
#
prot_nona2 <- as.data.frame(prot_rawdata[, c(which(grepl(pattern = "sig", x = colnames(prot_rawdata))))])
#
rownames(prot_nona2) <- prot_rawdata$Accession

prot_nona3 <- na.omit(prot_nona2)
#
prot_nona3 <- log2(prot_nona3)
#
prot_nona3$id <- rownames(prot_nona3)
#
colors <- RColorBrewer::brewer.pal(n = 7, name = "Blues")

####################################################################################
#PROTEIN COVERAGE OF GENES
## Whole protein dataset
p_all <- prot_nona3 %>% separate("id", into = c("db", "uniprot","id2"), sep = "\\|") %>% 
  left_join(map, by = "uniprot") 
#
##########################################################################
#########################################
#CV PLOT
## Use reduced protein list
protcv <- p_all[, 1:319]
#
rownames(protcv) <- p_all$uniprot
# Get cv for each protein
protcv$mean <- apply(protcv[, which(grepl("sig", colnames(protcv)))], 1, mean)
protcv$cv <- apply(protcv[, which(grepl("sig", colnames(protcv)))], 1, function(x) {
  round(sd(x)/mean(x), 4)
})

pv <- protcv %>% group_by(cv) %>% summarise(n = n())
#
pv <- pv %>% mutate(group = case_when( pv$cv < 0.03 ~ "1", pv$cv >= 0.03 ~ "2"))


## Scatterplot of number of proteins vs cv (similar to model plot)
ggplot(pv, aes( x = as.numeric(cv), y = n, col = group)) + 
  geom_point(size = 0.8) + 
  labs(x = "Number of Proteins", y = "Coefficient of Variance") + 
  #geom_vline(xintercept = (mean(protcv$cv)), lty = 2, color = "red") + 
  #geom_vline(xintercept = 0.005, lty = 2, color = "red") + 
  #geom_vline(xintercept = 0.03, lty = 2, color = "red") + 
  theme_classic()

###########################################
protcv$id <- rownames(protcv)
#
cv_low <- subset(protcv, cv < 0.005)
#
cv_high <- subset(protcv, cv > 0.03)
###########################################################
