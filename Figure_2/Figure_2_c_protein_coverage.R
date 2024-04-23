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
require(org.Rn.eg.db, quietly = T, warn.conflicts = F)
#
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
## The jump output with all proteins across all batches
library(readr)
#################################################################################################################
# rawdata 
#input the table from quan analysis(rawdata from JUMP q)
prot_rawdata  <- read_delim("rawdata/combined_norm_uni_prot.txt", delim = "\t", escape_double = FALSE,   trim_ws = TRUE, skip = 61)
#
colnames(prot_rawdata)[1:4] <- c("Protein_Group","Accession","Description","GN")
#
###############################################################################################################
#code for figure 2C
# Protein COVERAGE HISTOGRAM
# Pull out all coverage columns and find mean coverage across all batches for
# each protein
cov <- prot_rawdata[, c(2, which(grepl(pattern = "Coverage", x = colnames(prot_rawdata))))]
cov$mean <- apply(cov[, 2:ncol(cov)], 1, mean)

# Count the number of missing values
cov$zeros <- apply(cov[, 2:30], 1, function(x) {
  sum(x == 0)
})
#
cov <- cov %>% mutate(group = case_when(cov$zeros == 0 ~ "29", 
                                        cov$zeros == 1 ~ "26_28",
                                        cov$zeros == 2 ~ "26_28",
                                        cov$zeros == 3 ~ "26_28",
                                        cov$zeros == 4 ~ "21_25",
                                        cov$zeros == 5 ~ "21_25",
                                        cov$zeros == 6 ~ "21_25",
                                        cov$zeros == 7 ~ "21_25",
                                        cov$zeros == 8 ~ "21_25",
                                        cov$zeros == 9 ~ "16_20",
                                        cov$zeros == 10 ~ "16_20",
                                        cov$zeros == 11 ~ "16_20",
                                        cov$zeros == 12 ~ "16_20",
                                        cov$zeros == 13 ~ "16_20",
                                        cov$zeros == 14 ~ "11_15",
                                        cov$zeros == 15 ~ "11_15",
                                        cov$zeros == 16 ~ "11_15",
                                        cov$zeros == 17 ~ "11_15",
                                        cov$zeros == 18 ~ "11_15",
                                        cov$zeros == 19 ~ "6_10",
                                        cov$zeros == 20 ~ "6_10",
                                        cov$zeros == 21 ~ "6_10",
                                        cov$zeros == 22 ~ "6_10",
                                        cov$zeros == 23 ~ "6_10",
                                        cov$zeros == 24 ~ "1_5",
                                        cov$zeros == 25 ~ "1_5",
                                        cov$zeros == 26 ~ "1_5",
                                        cov$zeros == 27 ~ "1_5",
                                        cov$zeros == 28 ~ "1_5"
                                        ))
############
colors <- RColorBrewer::brewer.pal(n = 7, name = "Blues")
#
ggplot() + geom_histogram(data = cov[which(cov$group == "29"), ], aes(x = mean, fill = "29"), color = "black") + 
  geom_histogram(data = cov[which(cov$group == "26_28"), ], aes(x = mean, fill = "26_28"), color = "black") + 
  geom_histogram(data = cov[which(cov$group == "21_25"), ], aes(x = mean, fill = "21_25"), color = "black") + 
  geom_histogram(data = cov[which(cov$group == "16_20"), ], aes(x = mean, fill = "16_20"), color = "black") + 
  geom_histogram(data = cov[which(cov$group == "11_15"), ], aes(x = mean, fill = "11_15"), color = "black") +
  geom_histogram(data = cov[which(cov$group == "6_10"), ], aes(x = mean, fill = "6_10"), color = "black") + 
  geom_histogram(data = cov[which(cov$group == "1_5"), ], aes(x = mean, fill = "1_5"), color = "black") +
  scale_fill_manual(name = "Batches", breaks = c("29", "26_28", "21_25","16_20", "11_15", "6_10","1_5"), values = colors) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  labs(x = "Protein Coverage (%)",  y = "Frequency") + 
  theme_classic() + 
  theme(legend.position = c(0.9, 0.7))
#
# dev.off()
