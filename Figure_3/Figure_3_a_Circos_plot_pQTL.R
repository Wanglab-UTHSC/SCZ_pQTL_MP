
rm(list=ls())
#

require(dplyr, quietly = T, warn.conflicts = F)
require(stringr, quietly = T, warn.conflicts = F)
require(tidyr, quietly = T, warn.conflicts = F)
require(qvalue, quietly = T, warn.conflicts = F)
require(tidyverse, quietly = T, warn.conflicts = F)
require(circlize, quietly = T, warn.conflicts = F)


# CIRCOS QTL PLOT
setwd("D:/data/03222021_human_psy_qtl/Figure_code/")
#
pqtl_all <- read.table("rawdata/pQTL_pqtl.nominal.chr.cis_perm_log2_PEER_merge2.WGS", quote="\"", comment.char="")
#
pqtl_sig <- read.delim("rawdata/pQTL.cis_perm_log2_PEER_merge2.WGS_qvalue_sig.txt")
#
colnames(pqtl_all) <- colnames(pqtl_sig)[1:19]
#
names(pqtl_all)
####################################################################################
#
## Calculate q-values and get gene names for cis in protein
### Protein qtl
#input the table from quan analysis(rawdata from JUMP q)
# rawdata 
prot_gn  <- read_delim("rawdata/combined_norm_uni_prot.txt", delim = "\t", escape_double = FALSE,   trim_ws = TRUE, skip = 61)[1:4]
#
colnames(prot_gn)[1:4] <- c("id","Accession","Description","GN")
#
prot_gn <- prot_gn %>% separate("Accession", into = c("db", "uniprot","id2"), sep = "\\|") 
#
prot_gn2 <- prot_gn[,c(3,6)]
#
prot_merged <- merge(pqtl_all, prot_gn2, by.x = "Phenotype_ID", by.y = "uniprot")
###################################################################
# count qvalue
pvalue_p <- prot_merged$beta_distribution_pvalue
#
qobj <- qvalue(p = pvalue_p)
#
prot_merged$qvalue <-qobj$qvalues
#
#write.table(prot_merged, file = "pQTL_all_with_qvalue_gene_name.txt", quote = F, row.names = F, sep = "\t")
## Generate bed-ish file for plot
prot_merged$qvalue2 <- prot_merged$nominal_pvalue_of_association_between_Phenotype_and_top_variant_in_cis
#
prot_merged$qvalue2  <- ifelse(prot_merged$qvalue2 < 1e-40, 1e-40, prot_merged$qvalue2 )
#
pbed <- prot_merged[which(prot_merged$Total_number_of_Variants_in_cis != 0 & prot_merged$nominal_pvalue_of_association_between_Phenotype_and_top_variant_in_cis > 1e-60),c(9:11,22,20)]

pbed$ChrID_of_top_variant <- as.character(pbed$ChrID_of_top_variant)
pbed$ChrID_of_top_variant <- str_replace(pbed$ChrID_of_top_variant, "23", "X")
pbed$ChrID_of_top_variant <- str_replace(pbed$ChrID_of_top_variant, "24", "Y")
pbed$ChrID_of_top_variant <- str_replace(pbed$ChrID_of_top_variant, "^", "chr")
#
pbed$qvalue2 <- -log10(pbed$qvalue2)
#
pbed_sig <- pbed[which(pbed$qvalue2 > 7.3),]
########################################################
pbed1 <- subset(pbed, qvalue2 > 10)
#
pbed1_sig <- pbed1[which(pbed1$qvalue2 > 15),]
#
beds_high <- list(pbed1, pbed1_sig)
#####################################
pbed2 <- subset(pbed, qvalue2 <=  10 & qvalue2 > 0)
#
pbed2_sig <- pbed2[which(pbed2$qvalue2 > 2),]
#
beds_low <- list(pbed2, pbed2_sig)
#########################################

pbed_sig <- pbed[which(pbed$qvalue2 > 7.3),]
#

beds <- list(pbed, pbed_sig)
####################################################################################
#
## Get trans files ready
### pQTL trans
# added the sig protein cis name
tp_name <- subset(prot_merged, nominal_pvalue_of_association_between_Phenotype_and_top_variant_in_cis < 1E-40 & GN != "NA")

tp_name$GN <- str_replace(tp_name$GN,"C","c")
tp_name$ChrID_of_Phenotype <- str_replace(tp_name$ChrID_of_Phenotype, "23", "X")
tp_name$ChrID_of_Phenotype <- str_replace(tp_name$ChrID_of_Phenotype, "24", "Y")
tp_name$ChrID_of_Phenotype <- str_replace(tp_name$ChrID_of_Phenotype, "^", "chr")


labels_p <- unique(na.omit(data.frame("chrom" = tp_name$ChrID_of_Phenotype, "start" = as.numeric(tp_name$Start_position_of_Phenotype), "end" = as.numeric(tp_name$End_position_of_Phenotype), "name" = tp_name$GN)))


# Circos plot with links
circos.clear()

#names(CELL_META)

#pdf("circos_plot_with_sig_pQTL_sig_eQTL_lable_top_trans_protein_top_protein_name_version3.pdf",width=8,height=8)

circos.par("gap.degree" = c(rep(c(4, 4), 11),4,15), "start.degree" = 0)
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicLabels(labels_p, labels.column = 4, side = "outside",col = "black",facing = "clockwise")
#circos.initializeWithIdeogram(plotType = c("ideogram"), chromosome.index = chroms)
#
circos.genomicTrackPlotRegion(beds, panel.fun = function(region, value, ...){
     i = getI(...)
  circos.genomicPoints(region, value, cex = 0.2, pch = 19, col = i, ...)
}, bg.border = NA)
#

#circos.track(ylim = c(0,1), track.height = 0.05)

#
circos.genomicIdeogram(track.height = 0.02)
circos.track( ylim = c(0, 0.1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  #xlim = CELL_META$xlim
  #ylim = CELL_META$ylim
  #circos.rect(xlim[1], 0, xlim[2], 1, col = "white")
  circos.text(0, 0.01, chr, cex = 0.6, col = "black",
              facing = "downward", niceFacing = T)
}, track.height = 0.05, bg.border = NA)


