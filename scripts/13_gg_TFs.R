# This script corresponds to the Nature Communication resubmission.
# This code generates Extended data figure 3 b, Extended data figure 5 b.
#
# May 27 2018
# TAS


# initialize packages
library(pvclust); library(ggplot2); library(gplots); library(reshape2); library(bootSVD)


# load data ####
gene_list <- read.delim("raw_data/orthology_only_entrez_id.txt")
gene_list <- gene_list[c("gene_symbol")]

HTSeq_gg <- read.delim ("raw_data/gg_HTseq.txt")

gg_gene_id <- read.delim ("raw_data/gal_names_length_median.txt", na.strings="", col.names=c("gg_gene_id","length_gg", "gene_symbol"))

merge_gg <-merge(gg_gene_id, HTSeq_gg, by="gg_gene_id")
merge_gg <-merge (merge_gg, gene_list, by="gene_symbol")


# Save TF list ####
gene_list_TFs <-merge_gg[c("gg_gene_id", "gene_symbol")]
colnames(gene_list_TFs) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("transcription_factors", nrow(gene_list_TFs)))
species <- c(rep("chicken", nrow(gene_list_TFs)))
gene_list_TFs <- cbind (gene_list, species, gene_list_TFs)
gene_list_TFs <- unique(gene_list_TFs)
write.table(gene_list_TFs, file = "output_files/gene_list_TFs_chicken.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

merge_gg <- unique(merge_gg[-1]) # this reduces it to all unique ensemble genes; the above table will include non- 1:orthologous genes.

# calc TPM ####
fn.TPM_gg <-function(sample) {
  (sample*34*10^6)/(merge_gg$length_gg*sum((sample*34)/merge_gg$length_gg))
}
TPM_gg <- sapply(merge_gg[, 3:length(merge_gg)], fn.TPM_gg)
TPM_gg<- as.data.frame(sqrt(TPM_gg))
TPM_gg <- TPM_gg[rowSums(TPM_gg[-1])>0, ] 

# batch correction
attach(TPM_gg)
gg_FLE <- data.frame(gg_FL_D1_E, gg_FL_D2_E, gg_FL_D3_E) 
gg_FLL <- data.frame(gg_FL_D1_L, gg_FL_D2_L, gg_FL_D3_L) 
gg_HLE <- data.frame(gg_HL_D1_E, gg_HL_D2_E, gg_HL_D3_E, gg_HL_D4_E) 
gg_HLL <- data.frame(gg_HL_D1_L, gg_HL_D2_L, gg_HL_D3_L, gg_HL_D4_L) 
gg_adj.FLE <- (gg_FLE-rowMeans(gg_FLE))
gg_adj.FLL <- (gg_FLL-rowMeans(gg_FLL))
gg_adj.HLE <- (gg_HLE-rowMeans(gg_HLE))
gg_adj.HLL <- (gg_HLL-rowMeans(gg_HLL))
detach(TPM_gg)


# Extended Figure 3 b -- hindlimb analyses only ####
data <- data.frame(gg_adj.HLE, gg_adj.HLL)

# PCA ####
digit.number <- c('D1', 'D2', 'D3', 'D4', 'D1', 'D2', 'D3', 'D4')
stage <- c(rep('early', 4), rep('late', 4))
limb <- c(rep('chicken hindlimb', 8))
pheno.digit <- data.frame(sample=1:length(digit.number), digit.number, stage, limb, row.names=names(data))

pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

g<- qplot(x=PC1, y=PC2, data=scores, 
          colour=factor(pheno.digit$digit.number),
          shape=factor(pheno.digit$limb)) + 
  geom_point(size=2) + 
  theme_bw() + theme(legend.position="none") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio=1:1,
    text = element_text(size=8))
print(g)
ggsave("output_figures/extended_fig_3_b_gg_TFs_PCA.pdf", width = 6, height = 6, units=c("cm"))


# PCA boostrap ####
data2 <- as.matrix(data)
pca<-prcomp(t(data))
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_3_b_chicken", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extend_data_3_b_chicken_TFs.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA ####
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="chicken TFs"),
       file="output_figures/extended_fig_3_b_gg_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="chicken TFs") # for reference
x <- as.vector( noquote(hc$hclust$order))
x <- as.vector(c(1,5, 2,6, 3,7, 8,4)) # ordered for readability


# heatmap ####
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_fig_3_b_gg_heat_ordered.pdf", width = 8, height = 8, units=c("in"))



# Extended Figure 5 b -- forelimb and hindlimb analyses ####
data<- data.frame(gg_adj.FLE, gg_adj.HLE, gg_adj.FLL, gg_adj.HLL)

# PCA ####
digit.number <- c('D2', 'D3', 'D4', 'D1', 'D2', 'D3', 'D4',
                  'D2', 'D3', 'D4', 'D1', 'D2', 'D3', 'D4')
stage <- c(rep('early', 7), rep('late', 7))
limb <- c(rep('chicken forelimb', 3), 
          rep('chicken hindlimb', 4),
          rep('chicken forelimb', 3), 
          rep('chicken hindlimb', 4))
pheno.digit <- data.frame(sample = 1:length(digit.number),
                          digit.number, stage, limb,
                          row.names = names(data))

pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

g<- qplot(x=PC1, y=PC2, data=scores, 
          colour=factor(pheno.digit$digit.number),
          shape=factor(pheno.digit$limb)) + 
  geom_point(size=2) + 
  theme_bw() + theme(legend.position="none") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio=1:1,
    text = element_text(size=8))
print(g)
ggsave("output_figures/extended_fig_5_b_gg_TFs_PCA.pdf", width = 6, height = 6, units=c("cm"))


# PCA bootstraps ####
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_5_b_chicken", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extend_data_5_b_chicken_TFs.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA ####
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="chicken TFs"),
       file="output_figures/extended_fig_5_b_gg_TFs_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="chicken TFs") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # the figure for pub is saved in the original order.


# HEATMAP
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_fig_5_b_gg_TFs_heat.pdf", width = 8, height = 8, units=c("in"))
