# Nature Resubmission
# To make Figure 4 b -- heatmap of relative TPM changes between D1 and D2
# June 6 2018
# TAS


# initialize packages ####
library(reshape2); library(ggplot2)


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

# mouse ####
HTSeq_mm <- read.table("raw_data/mm_HTSeq.txt", sep="\t", header=T)
gene_length_mm <- read.delim ("raw_data/mm_gene_length_median.txt")

merged_mm <- merge(gene_length_mm, HTSeq_mm, by="mm_gene_id")
merged_mm <- merge(merged_mm, gene_list, by="mm_gene_id")
merged_mm <- merged_mm[order(merged_mm$gene_symbol),] 

# alligator ####
HTSeq_am <- read.delim("raw_data/am_HTSeq.txt")
gene_length_am <- read.delim("raw_data/am_gene_length_median.txt")
am_gene_id <- read.delim("raw_data/am_gene_annotations.txt", na.strings = "")
am_gene_id <- am_gene_id[c("am_gene_id", "gene_symbol")]

merged_am <-merge(am_gene_id, HTSeq_am, by="am_gene_id")
merged_am <-merge(gene_length_am, merged_am, by="am_gene_id")
merged_am <-merge(gene_list, merged_am, by="gene_symbol") 
merged_am <- merged_am[order(merged_am$gene_symbol),] 


# anolis ####
HTSeq_ac <- read.delim("raw_data/ac_HTSeq.txt")
gene_length_ac <- read.delim("raw_data/ac_gene_length_median.txt")

merged_ac <- merge(gene_length_ac, HTSeq_ac, by="ac_gene_id")
merged_ac <- merge(merged_ac, gene_list, by="ac_gene_id")
merged_ac <- merged_ac[order(merged_ac$gene_symbol),]


# chicken ####
HTSeq_gg <- read.delim ("raw_data/gg_HTseq.txt")
gg_gene_id <- read.delim ("raw_data/gal_names_length_median.txt", na.strings = "", col.names = c("gg_gene_id", "length_gg", "gene_symbol"))

merged_gg <-merge(gg_gene_id, HTSeq_gg, by="gg_gene_id")
merged_gg <-merge(merged_gg, gene_list, by="gene_symbol")
merged_gg <- merged_gg[order(merged_gg$gene_symbol),] 


# calculating TPMs ####
# human ####
fn.TPM_hs <-function(sample) {
  (sample*75*10^6)/(merged_hs$length_hs*sum((sample*75)/merged_hs$length_hs))
}
TPM_hs <- as.data.frame(sapply(merged_hs[3:(length(merged_hs)-5)], fn.TPM_hs))

# mouse ####
fn.TPM_mm <-function(sample) {
  (sample*75*10^6)/(merged_mm$length_mm*sum((sample*75)/merged_mm$length_mm))
}
TPM_mm <- as.data.frame(sapply(merged_mm[3:(length(merged_mm)-5)], fn.TPM_mm))

# alligator ####
fn.TPM_am <-function(sample) {
  (sample*75*10^6)/(merged_am$length_am*sum((sample*75)/merged_am$length_am))
}
TPM_am <- as.data.frame(sapply(merged_am[9:length(merged_am)], fn.TPM_am))

# anolis ####
fn.TPM_ac <-function(sample) {
  (sample*75*10^6)/(merged_ac$length_ac*sum((sample*75)/merged_ac$length_ac))
}
TPM_ac <- as.data.frame(sapply(merged_ac[3:(length(merged_ac)-5)], fn.TPM_ac))

# chicken ####
fn.TPM_gg <-function(sample) {
  (sample*34*10^6)/(merged_gg$length_gg*sum((sample*34)/merged_gg$length_gg))
}
TPM_gg <- as.data.frame(sapply(merged_gg[4:(length(merged_gg)-5)], fn.TPM_gg))


# calculate difference between D1 and D2 ####
# chicken #### 
attach(TPM_gg)

# chicken hindlimb
tpm <- (gg_HL_D2_E-gg_HL_D1_E) / ((gg_HL_D2_E+gg_HL_D1_E)/2)
limb <- c(rep("chicken h early", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
gg_HL_early <- data.frame (limb, comparison, gene_symbol, tpm)

tpm <- (gg_HL_D2_L-gg_HL_D1_L) / ((gg_HL_D2_L+gg_HL_D1_L) /2)
limb <- c(rep("chicken h late", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
gg_HL_late <- data.frame (limb, comparison, gene_symbol, tpm)

detach(TPM_gg)


# alligator ####
attach(TPM_am)

# alligator forelimb - early
tpm <- (rowMeans(cbind(FL_E_D2_1, FL_E_D2_2, FL_E_D2_3))-rowMeans(cbind(FL_E_D1_1, FL_E_D1_2, FL_E_D1_3))) /rowMeans(cbind(FL_E_D2_1, FL_E_D2_2, FL_E_D2_3, FL_E_D1_1, FL_E_D1_2, FL_E_D1_3)) 
limb <- c(rep("alligator forelimb early", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
am_FL_early <- data.frame (limb, comparison, gene_symbol, tpm)

# alligator forelimb - late 
tpm <- (rowMeans(cbind(FL_L_D2_1, FL_L_D2_2, FL_L_D2_3))-rowMeans(cbind(FL_L_D1_1, FL_L_D1_2, FL_L_D1_3))) / rowMeans(cbind(FL_L_D2_1, FL_L_D2_2, FL_L_D2_3, FL_L_D1_1, FL_L_D1_2, FL_L_D1_3))
limb <- c(rep("alligator forelimb late", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
am_FL_late <- data.frame (limb, comparison, gene_symbol, tpm)

detach(TPM_am)


# human ####
attach(TPM_hs)

# Human forelimb 88
tpm <- (Hs_FL_D25_11088-Hs_FL_D1_11088) / ((Hs_FL_D25_11088+Hs_FL_D1_11088)/2)
limb <- c(rep("human forelimb 88", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
hs_FL_88 <- data.frame (limb, comparison, gene_symbol, tpm)

# Human forelimb 17
tpm <- (Hs_FL_D25_11117-Hs_FL_D1_11117) / ((Hs_FL_D25_11117+Hs_FL_D1_11117)/2)
limb <- c(rep("human forelimb 17", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
hs_FL_17 <- data.frame (limb, comparison, gene_symbol, tpm)

# Human hindlimb 88
tpm <- (Hs_HL_D25_11088-Hs_HL_D1_11088) / ((Hs_HL_D25_11088+Hs_HL_D1_11088)/2)
limb <- c(rep("human hindlimb 88", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
hs_HL_88 <- data.frame (limb, comparison, gene_symbol, tpm)

# Human hindlimb 17
tpm <- (Hs_HL_D25_11117-Hs_HL_D1_11117) / ((Hs_HL_D25_11117+Hs_HL_D1_11117)/2)
limb <- c(rep("human hindlimb 17", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
hs_HL_17 <- data.frame (limb, comparison, gene_symbol, tpm)

detach(TPM_hs)


# anolis ####
attach(TPM_ac)

tpm <- (rowMeans(cbind( ac_FL_D2_1, ac_FL_D2_2, ac_FL_D2_3))-rowMeans(cbind(ac_FL_D1_1, ac_FL_D1_2, ac_FL_D1_3))) / rowMeans(cbind( ac_FL_D2_1, ac_FL_D2_2, ac_FL_D2_3, ac_FL_D1_1, ac_FL_D1_2, ac_FL_D1_3))
limb <- c(rep("anolis forelimb", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
ac_FL <- data.frame (limb, comparison, gene_symbol, tpm)

detach(TPM_ac)


# mouse ####
attach(TPM_mm)

tpm <- (rowMeans(cbind( mm_FL_D2_1, mm_FL_D2_2, mm_FL_D2_3))-rowMeans(cbind(mm_FL_D1_1, mm_FL_D1_2, mm_FL_D1_3))) / rowMeans(cbind( mm_FL_D2_1, mm_FL_D2_2, mm_FL_D2_3, mm_FL_D1_1, mm_FL_D1_2, mm_FL_D1_3))
limb <- c(rep("mouse forelimb", length(tpm)))
comparison <- c(rep("1v2", length(tpm)))
mm_FL <- data.frame (limb, comparison, gene_symbol, tpm)

detach(TPM_mm)

# plot heatmap for D1 D2 comparison ####
data <- as.data.frame(rbind(gg_HL_early, gg_HL_late,
                            am_FL_early, am_FL_late,
                            ac_FL, 
                            mm_FL, 
                            hs_FL_88, hs_FL_17,
                            hs_HL_88, hs_HL_17))

genelist <- c("ALX1", "DLX6", "DMRT2", "HAND2", "HOXD11", "HOXD12", "PAX9", "PITX1", "SALL1", "TFAP2B")

test <- data[data$gene_symbol %in% genelist,]
test <- test[order(test$gene_symbol),]

q <- qplot (x=limb, y=gene_symbol, data= test, fill=tpm, geom="tile")
q <- q + scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-2,2))
q <- q + theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_4_b_D1D2_heat.pdf", width = 8, height = 8, units=c("in"))
