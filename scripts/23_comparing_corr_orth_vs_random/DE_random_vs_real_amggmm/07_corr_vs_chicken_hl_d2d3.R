# This script is to be run after 01_chicken_mouse_alligator_vs_random.R
# It will recover both the correlation between orthologs among the 3 species for genes differentially expressed at:
# D2/D3 of chicken 
# and also between the differentially expressed genes and random genes of similar expression level.
# This is repeated for both early and late alligator samples, as well as chicken early and late samples.
# Correlation values are saved to a table in the 'output_files' folder.
# T.A.S.

limb <- c("mm",  "gg_e", "gg_l",  "am_e", "am_l")
simulated_gene_matrix <- vector("list",length(limb))

for(k in 1:length(limb)){
  taxonid <- c("mm_gene_id", "gg_gene_id","gg_gene_id", "am_gene_id", "am_gene_id")
  TPM_lists <- data.frame(TPM_mm_2,
                          "TPM_gg_HL_2_e"= TPM_gg$gg_HL_D2_E,
                          "TPM_gg_HL_2_l"= TPM_gg$gg_HL_D2_L,
                          TPM_am_2_e, TPM_am_2_l)
  id_lists <- list(merged_mm, merged_gg, merged_gg, merged_am, merged_am)
  foldchange_lists <- list(data.frame(tpm_fold_change_mm[,1:2], tpm_fold_change_mm["change2v3"]),
                           data.frame(tpm_fold_change_gg[,1:2], tpm_fold_change_gg["change2v3_e"]),
                           data.frame(tpm_fold_change_gg[,1:2], tpm_fold_change_gg["change2v3_l"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change2v3_e"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change2v3_l"]))
  
  # generate the list of TPMs from chicken HL E D1 to generate randoms
  gg_gene_id <- merged_gg$gg_gene_id
  TPM_gg_withID <- cbind(gg_gene_id, TPM_gg)
  target_list <- merge(TPM_gg_withID, one_to_one_gg_23, by="gg_gene_id")
  target_list <- target_list[order(target_list$gene_symbol),]
  
  target_list_TPM <- target_list$gg_HL_D2_E
  target_list_fold <- foldchange_lists[[2]]$gg_gene_id %in% one_to_one_gg_23$gg_gene_id
  target_list_fold <- foldchange_lists[[2]][target_list_fold,]

  target_list_fold <- target_list_fold[order(target_list_fold$gene_symbol),3]
  
  
  limb <- limb[k]
  taxonid <- taxonid[k]
  id_lists <- id_lists[[k]]
  
  # real
  real_fold_list <- merge(one_to_one_gg_23, foldchange_lists[[k]], by=taxonid, sort=F)
  
  # random
  find_nearest_TPM <- function(i){
    which.min( as.matrix( abs( TPM_lists[k] - target_list_TPM[i])))
  }
  
  simulated_gene_list <- as.vector(unlist(lapply(c(1:nrow(target_list)), find_nearest_TPM)))
  simulated_gene_list <- as.data.frame(id_lists[simulated_gene_list,taxonid,1])
  names(simulated_gene_list) = c(taxonid)
  simulated_gene_list <- merge(simulated_gene_list, foldchange_lists[[k]], by=taxonid, sort=F)
  
  cor(simulated_gene_list[3], target_list_fold)
  
  simulated_gene_matrix[[k]] <- c(cor(simulated_gene_list[3], target_list_fold), cor(real_fold_list[6], target_list_fold))
}
corr_table <- matrix(unlist(simulated_gene_matrix), nrow=2, ncol=5)
colnames(corr_table) <- c("vs_mm", "vs_gg_e", "vs_gg_l", "vs_am_e", "vs_am_l")
corr_table <- cbind (position="d2_d3", DE_source=rep("gg_e", 2), comparison=c("random", "orthologs"), corr_table)
write.table(corr_table,file="output_files/corr_table_orth_vs_random_amggmm.txt", append=TRUE, row.names=FALSE, col.names=FALSE)

# repeat for late state chicken ####
limb <- c("mm",  "gg_e", "gg_l",  "am_e", "am_l")
simulated_gene_matrix <- vector("list",length(limb))

for(k in 1:length(limb)){
  taxonid <- c("mm_gene_id", "gg_gene_id","gg_gene_id", "am_gene_id", "am_gene_id")
  TPM_lists <- data.frame(TPM_mm_2,
                          "TPM_gg_HL_2_e"= TPM_gg$gg_HL_D2_E,
                          "TPM_gg_HL_2_l"= TPM_gg$gg_HL_D2_L,
                          TPM_am_2_e, TPM_am_2_l)
  id_lists <- list(merged_mm, merged_gg, merged_gg, merged_am, merged_am)
  foldchange_lists <- list(data.frame(tpm_fold_change_mm[,1:2], tpm_fold_change_mm["change2v3"]),
                           data.frame(tpm_fold_change_gg[,1:2], tpm_fold_change_gg["change2v3_e"]),
                           data.frame(tpm_fold_change_gg[,1:2], tpm_fold_change_gg["change2v3_l"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change2v3_e"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change2v3_l"]))
  
  # generate the list of TPMs from chicken HL E D1 to generate randoms
  gg_gene_id <- merged_gg$gg_gene_id
  TPM_gg_withID <- cbind(gg_gene_id, TPM_gg)
  target_list <- merge(TPM_gg_withID, one_to_one_gg_23, by="gg_gene_id")
  target_list <- target_list[order(target_list$gene_symbol),]
  
  target_list_TPM <- target_list$gg_HL_D2_L
  target_list_fold <- foldchange_lists[[3]]$gg_gene_id %in% one_to_one_gg_23$gg_gene_id
  target_list_fold <- foldchange_lists[[3]][target_list_fold,]
  
  target_list_fold <- target_list_fold[order(target_list_fold$gene_symbol),3]
  
  
  limb <- limb[k]
  taxonid <- taxonid[k]
  id_lists <- id_lists[[k]]
  
  # real
  real_fold_list <- merge(one_to_one_gg_23, foldchange_lists[[k]], by=taxonid, sort=F)
  
  # random
  find_nearest_TPM <- function(i){
    which.min( as.matrix( abs( TPM_lists[k] - target_list_TPM[i])))
  }
  
  simulated_gene_list <- as.vector(unlist(lapply(c(1:nrow(target_list)), find_nearest_TPM)))
  simulated_gene_list <- as.data.frame(id_lists[simulated_gene_list,taxonid,1])
  names(simulated_gene_list) = c(taxonid)
  simulated_gene_list <- merge(simulated_gene_list, foldchange_lists[[k]], by=taxonid, sort=F)
  
  cor(simulated_gene_list[3], target_list_fold)
  
  simulated_gene_matrix[[k]] <- c(cor(simulated_gene_list[3], target_list_fold), cor(real_fold_list[6], target_list_fold))
}
corr_table <- matrix(unlist(simulated_gene_matrix), nrow=2, ncol=5)
colnames(corr_table) <- c("vs_mm", "vs_gg_e", "vs_gg_l", "vs_am_e", "vs_am_l")
corr_table <- cbind (position="d2_d3", DE_source=rep("gg_l", 2), comparison=c("random", "orthologs"), corr_table)
write.table(corr_table,file="output_files/corr_table_orth_vs_random_amggmm.txt", append=TRUE, row.names=FALSE, col.names=FALSE)

