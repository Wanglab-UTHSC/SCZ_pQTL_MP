rm(list=ls())
#
library(tidyr)
library(dplyr)
library(ggplot2)
require(biomaRt, quietly = T, warn.conflicts = F)
require(stringr, quietly = T, warn.conflicts = F)
require(tidyr, quietly = T, warn.conflicts = F)
require(qvalue, quietly = T, warn.conflicts = F)
require(tidyverse, quietly = T, warn.conflicts = F)
require(circlize, quietly = T, warn.conflicts = F)

# CIRCOS QTL PLOT
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
eqtl_cis <- read.delim("rawdata/step_4_eQTL_trans_prepare_for_ANNOVA.txt.variant_function", header=FALSE)
#
freq_1 <- as.data.frame(table(eqtl_cis$V1))
#
freq_1$group <- "eqtl_cis"

eqtl_trans <- read.delim("rawdata/step_2_eQTL_trans_prepare_for_ANNOVA.txt.variant_function", header=FALSE)
#
freq_2 <- as.data.frame(table(eqtl_trans$V1))
#
freq_2$group <- "eqtl_trans"
#
pqtl_cis <- read.delim("rawdata/step_2_pQTL_cis_prepare_for_ANNOVA.txt.variant_function", header=FALSE)
#
freq_3 <- as.data.frame(table(pqtl_cis$V1))
#
freq_3$group <- "pqtl_cis"
#
pqtl_trans <- read.delim("rawdata/step_2_pQTL_trans_prepare_for_ANNOVA.txt.variant_function", header=FALSE)
#
freq_4 <- as.data.frame(table(pqtl_trans$V1))
#
freq_4$group <- "pqtl_trans"
#
merged_all <- rbind(freq_1, freq_2, freq_3, freq_4)
#
#write.csv(merged_all, file = "step_3_ANNOVAR_rawdata_04222022.csv", quote =F, row.names = F)
#####################################################################
#
library(plyr)
library('RColorBrewer')
#
all <- read.csv("rawdata/step_3_ANNOVAR_rawdata_merged_with_percetage.csv")
#
all$group <- factor(all$group,levels = c("pqtl_cis","eqtl_cis", "pqtl_trans","eqtl_trans"))
#

color <- brewer.pal(5, "Blues")

#pdf('step_3_all_ANNOVAR_class_11062021.pdf',height = 5, width = 8)
#
ggplot(all) +
  geom_bar(aes(x = group, y = per, fill = Var1),
           stat = "identity") +
  labs(x = "group", y = NULL, fill = "")+
  scale_fill_manual(values = color)+
  theme_classic() +
  coord_flip()
#

#dev.off()

########################################################################
pqtl <- read.delim("rawdata/step_2_pQTL_cis_prepare_for_ANNOVA.txt.exonic_variant_function", header=FALSE)
#
table(pqtl$V2)
#
frq_p <- as.data.frame(table(pqtl$V2))
#
eqtl <- read.delim("rawdata/step_4_eQTL_trans_prepare_for_ANNOVA.txt.exonic_variant_function", header=FALSE)
#
table(eqtl$V2)
#
frq_e <- as.data.frame(table(eqtl$V2))
#
frq_e$group <- "eqtl"
#
frq_p$group <- "pqtl"
#
merged_all <- rbind(frq_e, frq_p)
#
per <- read.csv("rawdata/step_3_ANNOVAR_exonic_rawdata_with_percetage.csv")
#
per$group <- factor(per$group,levels = c("pqtl","vcf", "eqtl"))

#pdf('step_3_all_ANNOVAR_class_exonic_11062021.pdf',height = 5, width = 8)
#
ggplot(per) +
  geom_bar(aes(x = group, y = per, fill = Var1),
           stat = "identity") +
  labs(x = "group", y = NULL, fill = "")+
  scale_fill_manual(values = color)+
  theme_classic() +
  coord_flip()
#


