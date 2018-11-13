# Nature communication resubmission
# 
# Function of the script: plotting individual genes across limb
# (1) create a 1:1:1:1 orthology list for the five species considered: human, anolis, chicken, alligator, mouse
# (2) calculate TPMs for each species for this gene list.
# (3) plot these values across the limb (d1-d5).
# 
# June 12 2018
# TAS

# initialize packages
library("reshape")

# make orthology list for all species studied ####
# mouse, anolis, human -- v85
mm_ac_hs_id <- read.delim("raw_data/ensemble_v85_orthology_symbol.txt", na.strings = "")
mm_ac_hs_id <- mm_ac_hs_id[c("gene_symbol", "mm_gene_id", "ac_gene_id", "hs_gene_id")]
mm_ac_hs_id <- na.omit(mm_ac_hs_id)

# alligator
am_gene_id <- read.delim("raw_data/am_gene_annotations.txt", na.strings = "")
am_gene_id <- am_gene_id[c("am_gene_id", "gene_symbol")]

# chicken
gg_gene_id <- read.delim("raw_data/GalGal5.0_gene_length_symbol.txt", na.strings="", col.names=c("gg_gene_id","length_gg", "gene_symbol"))
gg_gene_id <- gg_gene_id[c("gg_gene_id", "gene_symbol")]
gg_gene_id <- na.omit(gg_gene_id)
gg_gene_id <- unique(gg_gene_id)


# combining alligator and chicken to other list
hs_ac_am_mm <- merge(am_gene_id, mm_ac_hs_id, "gene_symbol")
hs_ac_am_mm_gg <- merge(hs_ac_am_mm, gg_gene_id, "gene_symbol")
hs_ac_am_mm_gg_one2one <- hs_ac_am_mm_gg[!(duplicated(hs_ac_am_mm_gg$gene_symbol) | duplicated(hs_ac_am_mm_gg$gene_symbol, fromLast = TRUE)), ] # 1049 now.


# reduce count data to genes of interest ####
# human ####
HTSeq_hs <- read.delim("raw_data/hs_HTSeq.txt")
gene_length_hs <- read.delim("raw_data/hs_GRCh37.82_gene_length_median.txt")

merged_hs <-merge(gene_length_hs, HTSeq_hs, by="hs_gene_id") 
merged_hs <-merge(merged_hs, hs_ac_am_mm_gg_one2one, by="hs_gene_id")
merged_hs <- merged_hs[order(merged_hs$gene_symbol),] 


hs_ac_am_mm_gg_one2one <- merged_hs[c("gene_symbol", "hs_gene_id", "am_gene_id", "mm_gene_id", "ac_gene_id", "gg_gene_id")]
write.table(hs_ac_am_mm_gg_one2one, file = "output_files/acamhsmmgg_one2one_8045.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gene_list <- hs_ac_am_mm_gg_one2one


# chicken ####
HTSeq_gg <- read.delim ("raw_data/gg_HTseq.txt")
gg_gene_id <- read.delim ("raw_data/gal_names_length_median.txt", na.strings = "", col.names = c("gg_gene_id", "length_gg", "gene_symbol"))

merged_gg <-merge(gg_gene_id, HTSeq_gg, by="gg_gene_id")
merged_gg <-merge(merged_gg, gene_list, by="gene_symbol")
merged_gg <- merged_gg[order(merged_gg$gene_symbol),] 


# calculating TPMs ####
fn.TPM_gg <-function(sample) {
  (sample*34*10^6)/(merged_gg$length_gg*sum((sample*34)/merged_gg$length_gg))
}
TPM_gg <- as.data.frame(sapply(merged_gg[4:(length(merged_gg)-5)], fn.TPM_gg))


# build data table for plotting ####




attach(TPM_gg)
gene_symbol <- gene_list$gene_symbol

# chicken forelimb - early ####
gg_FL_early <- data.frame(gene_symbol, digit_1=c(rep(NA, length(gene_symbol))), gg_FL_D1_E, gg_FL_D2_E, gg_FL_D3_E)
gg_FL_early <- melt(gg_FL_early, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken forelimb - early", nrow(gg_FL_early)))
digit <- c(rep("D1", 1*length(gg_FL_D1_E)), 
           rep("D2", 1*length(gg_FL_D1_E)), 
           rep("D3", 1*length(gg_FL_D1_E)),
           rep("D4", 1*length(gg_FL_D1_E)))
gg_FL_early <- data.frame (limb, digit, gg_FL_early)

# chicken forelimb - late ####
gg_FL_late <- data.frame(gene_symbol, digit_1=c(rep(NA, length(gene_symbol))), gg_FL_D1_L, gg_FL_D2_L, gg_FL_D3_L)
gg_FL_late <- melt(gg_FL_late, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken forelimb - late", nrow(gg_FL_late)))
digit <- c(rep("D1", 1*length(gg_FL_D1_E)), 
           rep("D2", 1*length(gg_FL_D1_E)), 
           rep("D3", 1*length(gg_FL_D1_E)),
           rep("D4", 1*length(gg_FL_D1_E)))
gg_FL_late <- data.frame (limb, digit, gg_FL_late)


# chicken hindlimb - early ####
gg_HL_early <- data.frame(gene_symbol, gg_HL_D1_E, gg_HL_D2_E, gg_HL_D3_E, gg_HL_D4_E)
gg_HL_early <- melt(gg_HL_early, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken hindlimb - early", nrow(gg_HL_early)))
digit <- c(rep("D1", 1*length(gg_FL_D1_E)), 
           rep("D2", 1*length(gg_FL_D1_E)), 
           rep("D3", 1*length(gg_FL_D1_E)),
           rep("D4", 1*length(gg_FL_D1_E)))
gg_HL_early <- data.frame (limb, digit, gg_HL_early)

# chicken hindlimb - late ####
gg_HL_late <- data.frame(gene_symbol, gg_HL_D1_L, gg_HL_D2_L, gg_HL_D3_L, gg_HL_D4_L)
gg_HL_late <- melt(gg_HL_late, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken hindlimb - late", nrow(gg_HL_late)))
digit <- c(rep("D1", 1*length(gg_HL_D1_L)), 
           rep("D2", 1*length(gg_HL_D1_L)), 
           rep("D3", 1*length(gg_HL_D1_L)),
           rep("D4", 1*length(gg_HL_D1_L)))
gg_HL_late <- data.frame (limb, digit, gg_HL_late)

detach(TPM_gg)

# Aggregate tables
data<- rbind (gg_FL_early, gg_FL_late,gg_HL_early, gg_HL_late)
names(data)[names(data)=="value"] <- "TPM"

remove(digit)
remove(limb)
remove(gene_symbol)

attach(data)
data2 <- data.frame(gene_symbol, TPM, digit, limb)
detach(data)

data=data2

TPMfigure <- function(gene){
  data_specific <- data[data$gene_symbol==gene,]
  q <- qplot(digit, sqrt(TPM), data=data_specific, group=limb,  ylim=c(0,NA)) +
    geom_line(aes(colour = limb)) +
    theme_bw() + ggtitle(gene) + 
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
    theme(legend.position="none") # If you want grid lines, turn off the previous line
  pdf(paste("output_figures/plotting_genes_sqrt_tpm_gg","_", gene, ".pdf", sep = ''), width=3, height=3,onefile=T)
  print(q)
  dev.off()
}

TPMfigure("SALL1")
x <- c("DLX6", "HOXD11", "HOXD12", "LHX9", "TFAP2B", "TBX2", "TBX3", "SALL1")
lapply(x,TPMfigure)
