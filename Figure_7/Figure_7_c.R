rm(list=ls())

require(dplyr, quietly = T, warn.conflicts = F)
require(stringr, quietly = T, warn.conflicts = F)
require(tidyr, quietly = T, warn.conflicts = F)
require(qvalue, quietly = T, warn.conflicts = F)
require(tidyverse, quietly = T, warn.conflicts = F)

library(readxl)
# CIRCOS QTL PLOT
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
data <- read_excel("rawdata/Nature/Supplementary Table 1.xls", sheet = "ST1 624 LD clumps_disc")
####################################################################################
data_sub <- data[,c(2:4,1)]
# 
names(data_sub) <- c("chr","pos","p","snp")
###############################################################################
#prepare the data
sig_data <- data_sub %>% 
  subset(data_sub$p < 0.00000005)
notsig_data <-data_sub %>% 
  subset(data_sub$p >= 0.00000005)  %>%
  group_by(chr) 

gwas_data <- bind_rows(sig_data, notsig_data)

##############################################################################################

chr <- read_excel("rawdata/chr_length_hg19.xlsx")[1:2]
#
names(chr)
#prepare the data
data_cum <- chr %>% 
  mutate(bp_add = lag(cumsum(length*0.0000001), default = 0)) %>% 
  select(chr, bp_add)
#
names(data_cum)[1] <- "chr"
#
gwas_data <- data_sub  %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = pos*0.0000001 + bp_add)

axis_set  <- gwas_data %>% 
  group_by(chr) %>%                         # Specify group indicator
  summarise(center = mean(bp_cum))  

sig <- 5E-8
###############################################################################################
rank <- read.delim("rawdata/step_3_protein_QTL_with_scz_GWAS_rank.txt")
#
gwas_data$p <- ifelse(gwas_data$p < 1E-20, 1E-20, gwas_data$p )
#
gene <- read_excel("rawdata/step_2_protein_QTL_with_scz_GWAS_v2_624.xlsx",sheet = "selected")
#
#
merged <- merge(gwas_data, rank, by = "snp")
#
merged_sig1 <- merge(gwas_data,gene, by.x = "snp", by.y = "snp")
#
#
library(ggrepel)
#

ggplot(gwas_data, aes(x = bp_cum, y = -log10(p), 
                                  color = as.factor(chr))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(3,20)) +
  scale_color_manual(values = rep(c("black", "grey"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y= element_text(angle = 90, size = 8, vjust = 0.5)
  ) +
  geom_label_repel(data= merged_sig1,aes(x = bp_cum, 
      y = -log10(p),label = gid, max.overlaps = 50,show.legend = T, col = "red")) 


#