# Nature communication resubmission
# 
# Function of the script: plotting individual genes across limb
# (1) create a 1:1:1:1 orthology list for the five species considered: human, anolis, chicken, alligator, mouse
# (2) calculate TPMs for each species for this gene list.
# (3) plot these values across the limb (d1-d5).
# 
# June 12 2018
# TAS


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


# plotting multiple species on one plot ####
TPM_all <- data.frame(gene_list$gene_symbol, 
                      TPM_gg, TPM_hs, 
                      TPM_am, TPM_ac, TPM_mm)
names(TPM_all)[names(TPM_all)=="gene_list.gene_symbol"] <- "gene_symbol"


# average TPMs across replcates ####
attach(TPM_all)

TPM_mm_1 <- rowSums(TPM_mm[1:3])
TPM_mm_2 <- rowSums(TPM_mm[4:6])
TPM_mm_3 <- rowSums(TPM_mm[7:9])
TPM_mm_4 <- rowSums(TPM_mm[10:12])
TPM_mm_5 <- rowSums(TPM_mm[13:15])
TPM_mm_avg <- data.frame(TPM_mm_1,TPM_mm_2,TPM_mm_3,TPM_mm_4,TPM_mm_5)

TPM_ac_1 <- rowSums(TPM_ac[1:3])
TPM_ac_2 <- rowSums(TPM_ac[4:6])
TPM_ac_3 <- rowSums(TPM_ac[7:9])
TPM_ac_4 <- rowSums(TPM_ac[10:12])
TPM_ac_5 <- rowSums(TPM_ac[13:15])
TPM_ac_avg <- data.frame(TPM_ac_1,TPM_ac_2,TPM_ac_3,TPM_ac_4,TPM_ac_5)

TPM_am_1_e <- rowSums(TPM_am[1:3])
TPM_am_2_e <- rowSums(TPM_am[4:6])
TPM_am_3_e <- rowSums(TPM_am[7:8])
TPM_am_4_e <- rowSums(TPM_am[9:11])
TPM_am_5_e <- rowSums(TPM_am[12:14])
TPM_am_e_avg <- data.frame(TPM_am_1_e, TPM_am_2_e, TPM_am_3_e, TPM_am_4_e, TPM_am_5_e)

TPM_am_1_l <- rowSums(TPM_am[15:17])
TPM_am_2_l <- rowSums(TPM_am[18:20])
TPM_am_3_l <- data.frame(FL_L_D3_1, FL_L_D3_2, FL_L_D3_3)
TPM_am_3_l <- rowSums(TPM_am_3_l)
TPM_am_4_l <- data.frame(FL_L_D4_1, FL_L_D4_2, FL_L_D4_3)
TPM_am_4_l <- rowSums(TPM_am_4_l)
TPM_am_5_l <- data.frame(FL_L_D5_1, FL_L_D5_2, FL_L_D5_3)
TPM_am_5_l <- rowSums(TPM_am_5_l)
TPM_am_l_avg <- data.frame(TPM_am_1_l, TPM_am_2_l, TPM_am_3_l, TPM_am_4_l, TPM_am_5_l)

# Figure 4 C ####
# build data table for plotting ####
# alligator forelimb - early
am_FL_early <- data.frame(gene_symbol,TPM_am_e_avg)
am_FL_early <- melt(am_FL_early)

limb <- c(rep("alligator forelimb - early", nrow(am_FL_early)))
digit <- c(rep("I", 1*length(FL_E_D1_1)), 
           rep("II", 1*length(FL_E_D1_1)), 
           rep("III", 1*length(FL_E_D1_1)), 
           rep("IV", 1*length(FL_E_D1_1)), 
           rep("V", 1*length(FL_E_D1_1)) )
am_FL_early <- data.frame (limb, digit, am_FL_early)


# alligator forelimb - late
am_FL_late <- data.frame(gene_symbol, TPM_am_l_avg)
am_FL_late <- melt(am_FL_late)

limb <- c(rep("alligator forelimb - late", nrow(am_FL_late)))
digit <- c(rep("I", 1*length(FL_E_D1_1)), 
           rep("II", 1*length(FL_E_D1_1)), 
           rep("III", 1*length(FL_E_D1_1)), 
           rep("IV", 1*length(FL_E_D1_1)), 
           rep("V", 1*length(FL_E_D1_1)) )
am_FL_late <- data.frame (limb, digit, am_FL_late)

# anolis forelimb
ac_FL <- data.frame(gene_symbol, TPM_ac_avg)
ac_FL <- melt(ac_FL)

limb <- c(rep("anolis forelimb", nrow(ac_FL)))
digit <- c(rep("I", 1*length(ac_FL_D1_1)), 
           rep("II", 1*length(ac_FL_D1_1)), 
           rep("III", 1*length(ac_FL_D1_1)), 
           rep("IV", 1*length(ac_FL_D1_1)), 
           rep("V", 1*length(ac_FL_D1_1)) )
ac_FL <- data.frame (limb, digit, ac_FL)


# mouse
mm_FL <- data.frame(gene_symbol,TPM_mm_avg)
mm_FL <- melt(mm_FL)

limb <- c(rep("mouse forelimb", nrow(mm_FL)))
digit <- c(rep("I", 1*length(mm_FL_D1_1)), 
           rep("II", 1*length(mm_FL_D1_1)), 
           rep("III", 1*length(mm_FL_D1_1)), 
           rep("IV", 1*length(mm_FL_D1_1)), 
           rep("V", 1*length(mm_FL_D1_1)) )
mm_FL <- data.frame (limb, digit, mm_FL)


# chicken hindlimb - early ####
gg_HL_early <- data.frame(gene_symbol, gg_HL_D1_E, gg_HL_D2_E, gg_HL_D3_E, gg_HL_D4_E, digit_V=c(rep(NA, length(gene_symbol)))) 
gg_HL_early <- melt(gg_HL_early, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken hindlimb - early", nrow(gg_HL_early)))
digit <- c(rep("I", 1*length(gg_FL_D1_E)), 
           rep("II", 1*length(gg_FL_D1_E)), 
           rep("III", 1*length(gg_FL_D1_E)),
           rep("IV", 1*length(gg_FL_D1_E)),
           rep("V", 1*length(gg_FL_D1_E)) )
gg_HL_early <- data.frame (limb, digit, gg_HL_early)

# chicken hindlimb - late ####
gg_HL_late <- data.frame(gene_symbol, gg_HL_D1_L, gg_HL_D2_L, gg_HL_D3_L, gg_HL_D4_L, 
                         digit_V=c(rep(NA, length(gene_symbol)))) 
gg_HL_late <- melt(gg_HL_late, na.rm = FALSE, id.vars = "gene_symbol")

limb <- c(rep("chicken hindlimb - late", nrow(gg_HL_late)))
digit <- c(rep("I", 1*length(gg_HL_D1_L)), 
           rep("II", 1*length(gg_HL_D1_L)), 
           rep("III", 1*length(gg_HL_D1_L)),
           rep("IV", 1*length(gg_HL_D1_L)),
           rep("V", 1*length(gg_HL_D1_L)))
gg_HL_late <- data.frame (limb, digit, gg_HL_late)

detach(TPM_all)


# Aggregate tables
data<- rbind (am_FL_early, am_FL_late, 
              gg_HL_early, gg_HL_late,
              ac_FL, mm_FL)
names(data)[names(data)=="value"] <- "TPM"

remove(digit)
remove(limb)

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
  pdf(paste("output_figures/plotting_genes_sqrt_tpm_multispecies","_", gene, ".pdf", sep = ''), width=3, height=3,onefile=T)
  print(q)
  dev.off()
}

TPMfigure("SALL1")

x <- c("TBX2", "TBX3", "TBX15", "SALL1")
lapply(x,TPMfigure)
