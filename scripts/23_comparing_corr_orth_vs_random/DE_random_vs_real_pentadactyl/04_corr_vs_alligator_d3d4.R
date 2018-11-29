# This script is to be run after 01_anolis_mouse_alligator_vs_random.
# It will recover both the correlation between orthologs among the 3 species for genes differentially expressed at D2/D3 of alligator and also between the differentially expressed genes and random genes of similar expression level.
# This is repeated for both early and late alligator samples, and the correlation values are saved to a table in the 'output_files' folder.
# T.A.S.

limb <- c("mm", "ac", "am_e", "am_l")
simulated_gene_matrix <- vector("list",length(limb))

for(k in 1:length(limb)){
  taxonid <- c("mm_gene_id", "ac_gene_id", "am_gene_id", "am_gene_id")
  TPM_lists <- data.frame(TPM_mm_3,
                          TPM_ac_3,
                          TPM_am_3_e, TPM_am_3_l)
  id_lists <- list(merged_mm, merged_ac, merged_am, merged_am)
  foldchange_lists <- list(data.frame(tpm_fold_change_mm[,1:2], tpm_fold_change_mm["change3v4"]),
                           data.frame(tpm_fold_change_ac[,1:2], tpm_fold_change_ac["change3v4"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change3v4_e"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change3v4_l"]))
  
  TPM_withID <- TPM_am_avg
  
  target_list <- merge(TPM_withID, one_to_one_am_34, by="am_gene_id")
  target_list <- target_list[order(target_list$gene_symbol),]
  
  target_list_TPM <- target_list$TPM_am_3_e
  target_list_fold <- foldchange_lists[[3]]$am_gene_id %in% one_to_one_am_34$am_gene_id
  target_list_fold <- foldchange_lists[[3]][target_list_fold,]
  
  target_list_fold <- target_list_fold[order(target_list_fold$gene_symbol),3]
  
  limb <- limb[k]
  taxonid <- taxonid[k]
  id_lists <- id_lists[[k]]
  
  # real
  real_fold_list <- merge(one_to_one_am_34, foldchange_lists[[k]], by=taxonid, sort=F)
  
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

corr_table <- matrix(unlist(simulated_gene_matrix),nrow=2, ncol=4)
colnames(corr_table) <- c("vs_mm", "vs_ac", "vs_am_e", "vs_am_l")
corr_table <- cbind (position="d3_d4", DE_source=rep("am_e", 2), comparison=c("random", "orthologs"), corr_table)

write.table(corr_table,file="output_files/corr_table_orth_vs_random_pentadactyl.txt", append=TRUE, row.names=FALSE, col.names=FALSE)


# repeat for late stage alligator
limb <- c("mm", "ac", "am_e", "am_l")
simulated_gene_matrix <- vector("list",length(limb))

for(k in 1:length(limb)){
  taxonid <- c("mm_gene_id", "ac_gene_id", "am_gene_id", "am_gene_id")
  TPM_lists <- data.frame(TPM_mm_3,
                          TPM_ac_3,
                          TPM_am_3_e, 
                          TPM_am_3_l)
  id_lists <- list(merged_mm, merged_ac, merged_am, merged_am)
  foldchange_lists <- list(data.frame(tpm_fold_change_mm[,1:2], tpm_fold_change_mm["change3v4"]),
                           data.frame(tpm_fold_change_ac[,1:2], tpm_fold_change_ac["change3v4"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change3v4_e"]),
                           data.frame(tpm_fold_change_am[,1:2], tpm_fold_change_am["change3v4_l"]))
  
  TPM_withID <- TPM_am_avg
  
  target_list <- merge(TPM_withID, one_to_one_am_34, by="am_gene_id")
  target_list <- target_list[order(target_list$gene_symbol),]
  
  target_list_TPM <- target_list$TPM_am_3_l
  target_list_fold <- foldchange_lists[[4]]$am_gene_id %in% one_to_one_am_34$am_gene_id
  target_list_fold <- foldchange_lists[[4]][target_list_fold,]
  
  target_list_fold <- target_list_fold[order(target_list_fold$gene_symbol),3]
  
  limb <- limb[k]
  taxonid <- taxonid[k]
  id_lists <- id_lists[[k]]
  
  # real
  real_fold_list <- merge(one_to_one_am_34, foldchange_lists[[k]], by=taxonid, sort=F)
  
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

corr_table <- matrix(unlist(simulated_gene_matrix),nrow=2, ncol=4)
colnames(corr_table) <- c("vs_mm", "vs_ac", "vs_am_e", "vs_am_l")
corr_table <- cbind (position="d3_d4", DE_source=rep("am_l", 2), comparison=c("random", "orthologs"), corr_table)

write.table(corr_table,file="output_files/corr_table_orth_vs_random_pentadactyl.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
