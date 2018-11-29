# For Nature Communication Resubmission
#
# script functions to
# (1) compile all the PCA and bootsrap values as a supplementary table
# (2) compile all orthology and gene lists to a supplementary table.
#
# May 27 2018.
# TAS

# combining all PCA bootstrap values ####
fig1b <- read.delim("output_files/fig_1_b_mouse_limb.txt", head=TRUE)
fig1c <- read.delim("output_files/fig_1_c_alligator_limb.txt", head=TRUE)
fig1d <- read.delim("output_files/fig_1_d_anolis_limb.txt", head=TRUE)

fig5c_mouse <- read.delim("output_files/fig_5_c_mouse.txt", head=TRUE)
fig5c_alligator  <- read.delim("output_files/fig_5_c_alligator.txt", head=TRUE)
fig5c_Anolis  <- read.delim("output_files/fig_5_c_anolis.txt", head=TRUE)
fig5c_chicken  <- read.delim("output_files/fig_5_c_chicken.txt", head=TRUE)

fig6a  <- read.delim("output_files/fig_6_a_chicken.txt", head=TRUE)

extend_fig1a <- read.delim("output_files/extended_data_1_a_mouse.txt", head=TRUE)
extend_fig1b <- read.delim("output_files/extended_data_1_b_alligator.txt", head=TRUE)
extend_fig1c <- read.delim("output_files/extended_data_1_c_anolis.txt", head=TRUE)

extend_fig2a <- read.delim("output_files/extended_data_2_a_mouse.txt", head=TRUE)
extend_fig2b <- read.delim("output_files/extended_data_2_b_alligator.txt", head=TRUE)
extend_fig2c <- read.delim("output_files/extended_data_2_c_anolis.txt", head=TRUE)

extend_fig3a <- read.delim("output_files/extended_data_3_a_chicken_all.txt", head=TRUE)
extend_fig3b <- read.delim("output_files/extended_data_3_b_chicken_TFs.txt", head=TRUE)
extend_fig3c <- read.delim("output_files/extended_data_3_c_chicken_limb.txt", head=TRUE)

extend_fig5a <- read.delim("output_files/extended_data_5_a_chicken_all.txt", head=TRUE)
extend_fig5b <- read.delim("output_files/extended_data_5_b_chicken_TFs.txt", head=TRUE)
extend_fig5c <- read.delim("output_files/extended_data_5_c_chicken_limb.txt", head=TRUE)

all_bootstrap <- rbind(fig1b, fig1c, fig1d, 
                       fig5c_mouse, fig5c_alligator, fig5c_Anolis, 
                       fig5c_chicken, 
                       fig6a, 
                       extend_fig1a, extend_fig1b, extend_fig1c, 
                       extend_fig2a, extend_fig2b, extend_fig2c, 
                       extend_fig3a, extend_fig3b, extend_fig3c, 
                       extend_fig5a, extend_fig5b, extend_fig5c)
write.table(all_bootstrap, file = "output_files/PCA_bootstrap_all_plots.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# combining all gene lists ####
TFs_am <- read.delim("output_files/gene_list_TFs_alligator.txt")
TFs_ac <- read.delim("output_files/gene_list_TFs_anolis.txt")
TFs_gg <- read.delim("output_files/gene_list_TFs_chicken.txt")
TFs_hs <- read.delim("output_files/gene_list_TFs_human.txt")
TFs_mm <- read.delim("output_files/gene_list_TFs_mouse.txt")

limb_am <- read.delim("output_files/gene_list_limb_alligator.txt")
limb_ac <- read.delim("output_files/gene_list_limb_anolis.txt")
limb_gg <- read.delim("output_files/gene_list_limb_chicken.txt")
limb_mm <- read.delim("output_files/gene_list_limb_mouse.txt")

CDEG_am <- read.delim("output_files/gene_list_CDEG_alligator.txt")
CDEG_ac <- read.delim("output_files/gene_list_CDEG_anolis.txt")
CDEG_gg <- read.delim("output_files/gene_list_CDEG_chicken.txt")
CDEG_mm <- read.delim("output_files/gene_list_CDEG_mouse.txt")


one2one_TFs <- read.delim("output_files/amggmm_TFs_1049.txt")
one2one_TFs_1049_am <- data.frame(gene_list = c(rep("one2one_TFs_amggmm1049", 1049)),
                                  species = c(rep("alligator", 1049)),
                                  gene_id = one2one_TFs$am_gene_id,
                                  gene_symbol = one2one_TFs$gene_symbol)
one2one_TFs_1049_mm <- data.frame(gene_list = c(rep("one2one_TFs_amggmm1049", 1049)),
                                  species = c(rep("mouse", 1049)),
                                  gene_id = one2one_TFs$mm_gene_id,
                                  gene_symbol = one2one_TFs$gene_symbol)
one2one_TFs_1049_gg <- data.frame(gene_list = c(rep("one2one_TFs_amggmm1049", 1049)),
                                  species = c(rep("chicken", 1049)),
                                  gene_id = one2one_TFs$gg_gene_id,
                                  gene_symbol = one2one_TFs$gene_symbol)

# and again for the anolis one to one?


all_genelists <- rbind(TFs_am, TFs_ac, TFs_gg, TFs_mm, TFs_hs,
                       limb_am, limb_ac, limb_gg, limb_mm,
                       CDEG_am, CDEG_ac, CDEG_gg, CDEG_mm,
                       one2one_TFs_1049_am, one2one_TFs_1049_mm, one2one_TFs_1049_gg)
  
write.table(all_bootstrap, file = "output_files/PCA_bootstrap_all_plots.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
