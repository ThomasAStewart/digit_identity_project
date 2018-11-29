# This is used to generate a list of 1:1:1 orthologous genes by gene symbol for the Nature Communications submission.
# May 27 2018
# TAS

# load data ####
all_tfs <- read.delim("raw_data/orthology_only_entrez_id.txt", na.strings="")

# mouse
mm_tfs <- all_tfs[c("mm_gene_id", "gene_symbol")]
mm_tfs <- na.omit(mm_tfs)
mm_tfs <- unique(mm_tfs)

# alligator
am_gene_id <- read.delim("raw_data/am_gene_annotations.txt", na.strings = "")
am_gene_id <- am_gene_id[c("am_gene_id", "gene_symbol")]

# chicken
gg_gene_id <- read.delim("raw_data/GalGal5.0_gene_length_symbol.txt", na.strings="", col.names=c("gg_gene_id","length_gg", "gene_symbol"))
gg_gene_id <- gg_gene_id[c("gg_gene_id", "gene_symbol")]
gg_gene_id <- na.omit(gg_gene_id)
gg_gene_id <- unique(gg_gene_id)

# anolis
ac_tfs <- all_tfs[c("ac_gene_id", "gene_symbol")]
ac_tfs <- na.omit(ac_tfs)
ac_tfs <- unique(ac_tfs)

# alligator, chicken, mouse -- combine across species and reduce duplicates
tfs_ammm <- merge(am_gene_id, mm_tfs, "gene_symbol") # 1459 # mouse and alligator.
tfs_ammmgg <- merge(tfs_ammm, gg_gene_id, "gene_symbol") #1135 #  mouse+alligator and chicken
tfs_ammmgg_one2one <- tfs_ammmgg[!(duplicated(tfs_ammmgg$gene_symbol) | duplicated(tfs_ammmgg$gene_symbol, fromLast = TRUE)), ] # 1049 now.

write.table(tfs_ammmgg_one2one, file = "output_files/amggmm_TFs_1049.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# pentadactyl species --  combine across species and reduce duplicates
test <- merge(ac_tfs, mm_tfs, "gene_symbol")
test <- merge(am_gene_id, test, "gene_symbol") # 1459 # mouse and alligator.
test <- test[!(duplicated(test$gene_symbol) | duplicated(test$gene_symbol, fromLast = TRUE)), ] # 1346.
test <- test[!(duplicated(test$ac_gene_id) | duplicated(test$ac_gene_id, fromLast = TRUE)), ] # 1346.
write.table(test, file = "output_files/acammm_TFs_1121.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
