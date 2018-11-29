# This script was written for the resubmission to Nature Communications.
#
# This shows BH method for FDR correction.
# script functions to:
# (1) create all alligator, Anolis, chicken venn diagrams for pairwise digit and also ANOVA comparisons
# (2) create all pentadactyl venn diagrams for pairwise digit comparisons
#
# May 27th, 2018
# TAS

library(qvalue)

# load data ####
pval_mouse_12 <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_lrt.1vs2_all_ONgene.txt")
pval_mouse_23 <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_lrt.2vs3_all_ONgene.txt")
pval_mouse_34 <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_lrt.3vs4_all_ONgene.txt")
pval_mouse_45 <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_lrt.4vs5_all_ONgene.txt")

pval_alligator_12 <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_lrt.1vs2_all_ONgene.txt")
pval_alligator_23 <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_lrt.2vs3_all_ONgene.txt")
pval_alligator_34 <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_lrt.3vs4_all_ONgene.txt")
pval_alligator_45 <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_lrt.4vs5_all_ONgene.txt")

pval_anolis_12 <- read.delim("raw_data/edgeR_outputs/anolis_pairwise_PC1_all_ONgene_pvals/anolis_edgeR_lrt.pc.1vs2_all_ONgene.txt")
pval_anolis_23 <- read.delim("raw_data/edgeR_outputs/anolis_pairwise_PC1_all_ONgene_pvals/anolis_edgeR_lrt.pc.2vs3_all_ONgene.txt")
pval_anolis_34 <- read.delim("raw_data/edgeR_outputs/anolis_pairwise_PC1_all_ONgene_pvals/anolis_edgeR_lrt.pc.3vs4_all_ONgene.txt")
pval_anolis_45 <- read.delim("raw_data/edgeR_outputs/anolis_pairwise_PC1_all_ONgene_pvals/anolis_edgeR_lrt.pc.4vs5_all_ONgene.txt")

pval_chicken_12 <- read.delim("raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.1vs2_all_ONgene.txt")
pval_chicken_23 <- read.delim("raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.2vs3_all_ONgene.txt")
pval_chicken_34 <- read.delim("raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.3vs4_all_ONgene.txt")

pval_anova_alligator  <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_anova_NULL_all_ONgene.txt")
pval_anova_mouse <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_anova_NULL_all_ONgene.txt")
pval_anova_anolis <- read.delim("raw_data/edgeR_outputs/anolis_pairwise_PC1_all_ONgene_pvals/anolis_edgeR_anova_PC1_all_ONgene.txt")
pval_anova_chicken <- read.delim("raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_anova_NULL_all_ONgene.txt")


# calculate q values with BH ####
# These q values are printed in figures 2, 5 a, and Extended Data fig 4 a
# note: the original output FDR value has that method; I also did it explicitly for the ANOVA data sets.
qval_mouse_12 <- pval_mouse_12[(pval_mouse_12$FDR < 0.05),]
qval_mouse_23 <- pval_mouse_23[(pval_mouse_23$FDR < 0.05),]
qval_mouse_34 <- pval_mouse_34[(pval_mouse_34$FDR < 0.05),]
qval_mouse_45 <- pval_mouse_45[(pval_mouse_45$FDR < 0.05),]

qval_alligator_12 <- pval_alligator_12[(pval_alligator_12$FDR < 0.05),]
qval_alligator_23 <- pval_alligator_23[(pval_alligator_23$FDR < 0.05),]
qval_alligator_34 <- pval_alligator_34[(pval_alligator_34$FDR < 0.05),]
qval_alligator_45 <- pval_alligator_45[(pval_alligator_45$FDR < 0.05),]

qval_anolis_12 <- pval_anolis_12[(pval_anolis_12$FDR < 0.05),]
qval_anolis_23 <- pval_anolis_23[(pval_anolis_23$FDR < 0.05),]
qval_anolis_34 <- pval_anolis_34[(pval_anolis_34$FDR < 0.05),]
qval_anolis_45 <- pval_anolis_45[(pval_anolis_45$FDR < 0.05),]

qval_chicken_12 <- pval_chicken_12[p.adjust(pval_chicken_12$PValue, method = "BH") < 0.05,] # still needs to be calculatd for these files
qval_chicken_23 <- pval_chicken_23[p.adjust(pval_chicken_23$PValue, method = "BH") < 0.05,]
qval_chicken_34 <- pval_chicken_34[p.adjust(pval_chicken_34$PValue, method = "BH") < 0.05,]

qval_mouse_anova <- pval_anova_mouse[p.adjust(pval_anova_mouse$PValue, method = "BH") < 0.05,]
qval_alligator_anova <- pval_anova_alligator[p.adjust(pval_anova_alligator$PValue, method = "BH") < 0.05,]
qval_anolis_anova <- pval_anova_anolis[p.adjust(pval_anova_anolis$PValue, method = "BH") < 0.05,]
qval_chicken_anova <- pval_anova_chicken[p.adjust(pval_anova_chicken$PValue, method = "BH") < 0.05,]


# Count transcription factors ####
# These TF countsare printed in figures 2, 5 a, and Extended Data fig 4 a
TFs_mouse <- read.delim("output_files/gene_list_TFs_mouse.txt")
TFs_mouse <- as.vector(unique(TFs_mouse$gene_id))
TFs_mouse_12 <- qval_mouse_12[rownames(qval_mouse_12) %in% TFs_mouse,]
TFs_mouse_23 <- qval_mouse_23[rownames(qval_mouse_23) %in% TFs_mouse,]
TFs_mouse_34 <- qval_mouse_34[rownames(qval_mouse_34) %in% TFs_mouse,]
TFs_mouse_45 <- qval_mouse_45[rownames(qval_mouse_45) %in% TFs_mouse,]

TFs_alligator <- read.delim("output_files/gene_list_TFs_alligator.txt")
TFs_alligator <- as.vector(unique(TFs_alligator$gene_id))
TFs_alligator_12 <- qval_alligator_12[rownames(qval_alligator_12) %in% TFs_alligator,]
TFs_alligator_23 <- qval_alligator_23[rownames(qval_alligator_23) %in% TFs_alligator,]
TFs_alligator_34 <- qval_alligator_34[rownames(qval_alligator_34) %in% TFs_alligator,]
TFs_alligator_45 <- qval_alligator_45[rownames(qval_alligator_45) %in% TFs_alligator,]

TFs_anolis <- read.delim("output_files/gene_list_TFs_anolis.txt")
TFs_anolis <- as.vector(unique(TFs_anolis$gene_id))
TFs_anolis_12 <- qval_anolis_12[rownames(qval_anolis_12) %in% TFs_anolis,]
TFs_anolis_23 <- qval_anolis_23[rownames(qval_anolis_23) %in% TFs_anolis,]
TFs_anolis_34 <- qval_anolis_34[rownames(qval_anolis_34) %in% TFs_anolis,]
TFs_anolis_45 <- qval_anolis_45[rownames(qval_anolis_45) %in% TFs_anolis,]

TFs_chicken <- read.delim("output_files/gene_list_TFs_chicken.txt")
TFs_chicken <- as.vector(unique(TFs_chicken$gene_id))
TFs_chicken_12 <- qval_chicken_12[rownames(qval_chicken_12) %in% TFs_chicken,]
TFs_chicken_23 <- qval_chicken_23[rownames(qval_chicken_23) %in% TFs_chicken,]
TFs_chicken_34 <- qval_chicken_34[rownames(qval_chicken_34) %in% TFs_chicken,]

TFs_mouse_anova <- qval_mouse_anova[rownames(qval_mouse_anova) %in% TFs_mouse,]
TFs_alligator_anova <- qval_alligator_anova[rownames(qval_alligator_anova) %in% TFs_alligator,]
TFs_chicken_anova <- qval_chicken_anova[rownames(qval_chicken_anova) %in% TFs_chicken,]
TFs_anolis_anova <- qval_anolis_anova[rownames(qval_anolis_anova) %in% TFs_anolis,]


# chicken, alligator, mouse -- intersection with 1:1 list ####
one2one_list <- read.delim("output_files/amggmm_TFs_1049.txt") # this is generated by the script "mklist_one_to_one_TFS_ammmgg.R"

# Figure 4 and extended data figure 4
# first reduce to that list.
one_to_one_gg_12 <- one2one_list[one2one_list$gg_gene_id %in% rownames(qval_chicken_12),1]
one_to_one_gg_23 <- one2one_list[one2one_list$gg_gene_id %in% rownames(qval_chicken_23),1]
one_to_one_gg_34 <- one2one_list[one2one_list$gg_gene_id %in% rownames(qval_chicken_34),1]

one_to_one_mm_12 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_12),1]
one_to_one_mm_23 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_23),1]
one_to_one_mm_34 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_34),1]
one_to_one_mm_45 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_45),1]

one_to_one_am_12 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_12),1]
one_to_one_am_23 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_23),1]
one_to_one_am_34 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_34),1]
one_to_one_am_45 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_45),1]

# then count DE genes per species for the ven diagrams
c(length(one_to_one_gg_12), length(one_to_one_gg_23), length(one_to_one_gg_34))
c(length(one_to_one_mm_12), length(one_to_one_mm_23), length(one_to_one_mm_34), length(one_to_one_mm_45))
c(length(one_to_one_am_12), length(one_to_one_am_23), length(one_to_one_am_34), length(one_to_one_am_45))

# venn diagram -- center of circles
intersect(intersect(one_to_one_gg_12, one_to_one_mm_12), one_to_one_am_12) # 10 same ones as before.
intersect(intersect(one_to_one_gg_23, one_to_one_mm_23), one_to_one_am_23) # 0
intersect(intersect(one_to_one_gg_34, one_to_one_mm_34), one_to_one_am_34) # 1
intersect(one_to_one_am_45, one_to_one_mm_45)

# venn diagram -- chicken and mouse
intersect(one_to_one_gg_12, one_to_one_mm_12) # 35 - 10
intersect(one_to_one_gg_23, one_to_one_mm_23) # 0
intersect(one_to_one_gg_34, one_to_one_mm_34) # 3 - 1

# venn diagram -- chicken and alligator
intersect(one_to_one_gg_12, one_to_one_am_12) 
intersect(one_to_one_gg_23, one_to_one_am_23) 
intersect(one_to_one_gg_34, one_to_one_am_34)

# venn diagram -- mouse and alligator
intersect(one_to_one_mm_12, one_to_one_am_12) 
intersect(one_to_one_mm_23, one_to_one_am_23) 
intersect(one_to_one_mm_34, one_to_one_am_34)
intersect(one_to_one_mm_45, one_to_one_am_45)


# ANOVA differential expression ####
# intersect with 1:1:1 list
one_to_one_mm_anova <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_anova),1]
one_to_one_am_anova <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_anova),1]
one_to_one_gg_anova <- one2one_list[one2one_list$gg_gene_id %in% rownames(qval_chicken_anova),1]
# center of venn diagram
intersect(intersect(one_to_one_gg_anova, one_to_one_mm_anova), one_to_one_am_anova) # n = 49
length(intersect(one_to_one_gg_anova, one_to_one_mm_anova)) - 49
length(intersect(one_to_one_am_anova, one_to_one_mm_anova)) - 49
length(intersect(one_to_one_gg_anova, one_to_one_am_anova)) - 49




# pentadactyl species -- intersection with 1:1 list ####
# should run lines 1-61 first.
one2one_list <- read.delim("output_files/acammm_TFs_1133.txt") # this is generated by the script "mklist_one_to_one_TFS_ammmgg.R"

# Figure 4 and extended data figure 4
# first reduce to that list.
one_to_one_ac_12 <- one2one_list[one2one_list$ac_gene_id %in% rownames(qval_anolis_12),1]
one_to_one_ac_23 <- one2one_list[one2one_list$ac_gene_id %in% rownames(qval_anolis_23),1]
one_to_one_ac_34 <- one2one_list[one2one_list$ac_gene_id %in% rownames(qval_anolis_34),1]
one_to_one_ac_45 <- one2one_list[one2one_list$ac_gene_id %in% rownames(qval_anolis_45),1]

one_to_one_mm_12 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_12),1]
one_to_one_mm_23 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_23),1]
one_to_one_mm_34 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_34),1]
one_to_one_mm_45 <- one2one_list[one2one_list$mm_gene_id %in% rownames(qval_mouse_45),1]

one_to_one_am_12 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_12),1]
one_to_one_am_23 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_23),1]
one_to_one_am_34 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_34),1]
one_to_one_am_45 <- one2one_list[one2one_list$am_gene_id %in% rownames(qval_alligator_45),1]

# then count DE genes per species for the ven diagrams
c(length(one_to_one_ac_12), length(one_to_one_ac_23), length(one_to_one_ac_34), length(one_to_one_ac_45))
c(length(one_to_one_mm_12), length(one_to_one_mm_23), length(one_to_one_mm_34), length(one_to_one_mm_45))
c(length(one_to_one_am_12), length(one_to_one_am_23), length(one_to_one_am_34), length(one_to_one_am_45))

# venn diagram -- center of circles
intersect(intersect(one_to_one_ac_12, one_to_one_mm_12), one_to_one_am_12) # 3
intersect(intersect(one_to_one_ac_23, one_to_one_mm_23), one_to_one_am_23) # 0
intersect(intersect(one_to_one_ac_34, one_to_one_mm_34), one_to_one_am_34) # 0
intersect(intersect(one_to_one_ac_45, one_to_one_am_45), one_to_one_mm_45) # 1

# venn diagram -- anolis and mouse
intersect(one_to_one_mm_12, one_to_one_ac_12) 
intersect(one_to_one_mm_23, one_to_one_ac_23) 
intersect(one_to_one_mm_34, one_to_one_ac_34)
intersect(one_to_one_mm_45, one_to_one_ac_45)

# venn diagram -- anolis and alligator
intersect(one_to_one_am_12, one_to_one_ac_12) 
intersect(one_to_one_am_23, one_to_one_ac_23) 
intersect(one_to_one_am_34, one_to_one_ac_34)
intersect(one_to_one_am_45, one_to_one_ac_45)

# venn diagram -- mouse and alligator
intersect(one_to_one_mm_12, one_to_one_am_12) 
intersect(one_to_one_mm_23, one_to_one_am_23) 
intersect(one_to_one_mm_34, one_to_one_am_34)
intersect(one_to_one_mm_45, one_to_one_am_45)