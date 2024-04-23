rm(list=ls())
#
library(readxl)
#
library(tidyr)
library(dplyr)
library(eoffice)
library(ggplot2)
library(corrplot) 
library(ggord)
library(FactoMineR)
library(ggrepel)
library(qvalue)

#
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
# CIRCOS QTL PLOT
pqtl <- read.table("rawdata/pQTL_pqtl.nominal.chr.cis_perm_log2_PEER_merge2.WGS", quote="\"", comment.char="")
#
qtl_sig <- read.delim("rawdata/pQTL.cis_perm_log2_PEER_merge2.WGS_qvalue_sig.txt")
#
names(pqtl) <- colnames(qtl_sig[1:19])
#
###############################################################################################
#
pvalue_p <- pqtl$beta_distribution_pvalue
#
qobj <- qvalue(p = pvalue_p)
#
pqtl$qvalue <-qobj$qvalues
#################################################
#
id <- read_delim("rawdata/combined_norm_uni_prot.txt", delim = "\t", escape_double = FALSE,  trim_ws = TRUE, skip = 61)[c(2,4)]
#
names(id) <- c("id","gn")
#
#
id <-  id %>% separate("id", into = c("id1", "uniprot","id2"), sep = "\\|") 
#
id2 <-id[c(2,4)]
#
merged_qtl <- merge(pqtl, id2 , by.x = "Phenotype_ID", by.y = "uniprot")
##########################################################################
#
library(readxl)
#
szgene <- read.delim("rawdata/SCZ_BIP_database_merged_all_v2_11242021.txt")[1]
#
szgene$szgene <- 1
#szgene <- read_excel("E:/data/03222021_human_psy_qtl/GeNets/szgene_group.xlsx")
#
merged_qtl_gz <- merge(merged_qtl, szgene, by.x = "gn", by.y = "Gene", all.x = TRUE)
###################################################################################
#write.csv(merged_qtl_gz, file = "pQTL_results_with_qvalue_gene_name_szgene_02032022.csv", quote = F, row.names = F)
merged_qtl_gz_2 <- subset(merged_qtl_gz, qvalue > 0)


############################################################################################
merged_qtl_gz_2 <- merged_qtl_gz_2  %>% mutate( effect_group  = case_when(Corresponding_regression_slope > 0 ~ "up", 
                                                                Corresponding_regression_slope < 0 ~ "down" )
                            ) 
#
table(merged_qtl_gz_2$effect_group)
#
######################################################################################################
maf <- read.delim("rawdata/SNP_maf.txt")
#
names(merged_qtl_gz_2)
#
merged_maf <- merge(merged_qtl_gz_2, maf , by.x = "ID_of_the_top_variant", by.y  = "ID")
#
merged_lable <- subset(merged_maf, abs(Corresponding_regression_slope) > 0.6 & maf > 0.05)


ggplot(data = merged_maf, aes ( x = maf, y = abs(Corresponding_regression_slope))) + geom_point(color = "grey")+
  theme_classic()+
  #geom_smooth(method = 'gam', span = 0.5, formula = y~s(x, k = 10), se = TRUE, level = 0.95)+
  geom_smooth(method = 'loess', span = 0.5, se = F) +
  xlim(0, 0.5)+ylim(0,2.5)+
  geom_text_repel(data = merged_lable,aes ( x = maf, y = abs(Corresponding_regression_slope), label = gn),
                  size = 3, box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", max.overlaps = 50,
                  show.legend = T) +
  ylim(0,2) + 
  labs(x = "Protein Coverage (%)",  y = "Frequency") + 
  theme_classic()


