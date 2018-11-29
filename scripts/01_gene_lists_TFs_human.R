# Script made for Nat.Commun. resubmission, specifically as part of the species-specific list of transcription factors.
# last modified: May 27, 2018
# TAS

# load data ####
TFs_human <- read.delim("raw_data/orthology_only_entrez_id.txt", na.strings = "")
TFs_human <- TFs_human[c("hs_gene_id", "gene_symbol")]

# format for the combined table ####
colnames(TFs_human) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("transcription_factors", nrow(TFs_human)))
species <- c(rep("human", nrow(TFs_human)))
gene_list_TFs_human <- cbind (gene_list, species,TFs_human)
gene_list_TFs_human <- unique(gene_list_TFs_human)

write.table(gene_list_TFs_human, file = "output_files/gene_list_TFs_human.txt", sep = "\t", row.names = FALSE, col.names = TRUE)