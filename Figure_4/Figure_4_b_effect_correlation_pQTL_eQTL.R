rm(list=ls())
#
library(biomaRt)
library(dplyr)
library(tidyr)
library(qvalue)
library(ggplot2)
library(readxl)

# READ IN FILES
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
library(readxl)

coloc_all <- read_excel("rawdata/coloc_results_01312022.xlsx", sheet = "coloc_pQTL")
#
data_sub<- subset(coloc_all, PP.H4.abf > 0.8 & PP.H4_max == "Yes")
########################################################
pqtl <- read_excel("rawdata/Supplementary_Table_v2.xlsx", 
                   sheet = "Table S4A", skip = 2)[c(2,3,4,7,10,11)]
names(pqtl) <- c("gn","pid","p_chr","p_pos","q_pqtl","slope_pqtl")
#################################
library(readr)
eqtl <- read_excel("rawdata/Supplementary_Table_v2.xlsx", 
                            sheet = "Table S5A", skip = 2)[c(3,4,7,10,11)]
#
names(eqtl) <- c("gid","g_chr","g_pos","q_eqtl","slope_eqtl")
#
#
merged <- merge(data_sub, pqtl, by.x = "ProteinID", by.y = "pid")
#
merged2 <- merge(merged,eqtl, by.x = "GeneID", by.y = "gid")
#

cor(merged2$slope_pqtl, merged2$slope_eqtl)

res <- cor.test(merged2$slope_pqtl, merged2$slope_eqtl, method = "pearson")
#
res
#
#
count <- merged2  %>% mutate( pqtl = case_when(merged2$slope_pqtl > 0 ~ "1", 
                                             merged2$slope_pqtl < 0 ~ "-1" ),
                            eqtl = case_when(merged2$slope_eqtl > 0 ~ "1", 
                                             merged2$slope_eqtl < 0 ~ "-1"))


count <- count  %>% mutate( group = case_when(count$pqtl ==  -1 & count$eqtl == -1 ~ "neg_neg", 
                                              count$pqtl ==  1 & count$eqtl == -1 ~ "pos_neg", 
                                              count$pqtl ==  -1 & count$eqtl == 1 ~ "neg_pos", 
                                              count$pqtl ==  1 & count$eqtl == 1 ~ "pos_pos", ))

######################################################################
# count the same trend
#
#
count$p_snp<- paste(count$p_chr, count$p_pos, sep = ":")
#
count$e_snp<- paste(count$g_chr, count$g_pos, sep = ":")

#
count <- count  %>% mutate( pos_same = case_when(count$p_snp == count$e_snp  ~ "same", 
                                                 count$p_snp != count$e_snp  ~ "diff"))
table(count$pos_same, count$group)
#
table(count$pos_same)
#
pos_same <- subset(count, pos_same == "same" )

count_sub <- subset(count, pos_same == "same" & (group == "pos_neg" | group == "neg_pos"))
#
p3 <- 75/386
#pdf("dotplot_effect_size_pQTL_eQTL_v2_01122022.pdf")
ggplot(data = pos_same , aes(y = slope_pqtl, x = slope_eqtl, color = group )) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black", fill = "lightgray")+
  theme_classic() +
  geom_hline(aes(yintercept = 0), color ="lightgrey", linetype="dashed")+
  geom_vline(aes(xintercept = 0), color ="lightgrey", linetype="dashed")+
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1))
#
#dev.off()

###############################################################################
gene <- as.data.frame(pos_same[,19] )

names(gene) <- "effect_size"

gene$group <- "gene"


prot <- as.data.frame(pos_same[,19] )

names(prot) <- "effect_size"

prot$group <- "prot"


merged <- rbind(gene, prot)


mean(abs(pos_same$slope_eqtl))

mean(abs(pos_same$slope_pqtl))



median(abs(pos_same$slope_eqtl))

median(abs(pos_same$slope_pqtl))


ggplot(data = merged , aes(y = abs(effect_size) , x = group, color = group )) + 
  geom_boxplot() + 
  theme_classic() +
  scale_y_continuous(limits = c(0,0.4)) +
  geom_jitter()
#
#dev.off()