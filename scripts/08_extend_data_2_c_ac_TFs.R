# This script corresponds to the Nature Communication resubmission.
# Specifically, it has been modified to extract standard error values for the PCA via bootstrapping. 
#
# May 27 2018
# TAS


# initialize packages ####
library(pvclust); library(ggplot2); library(gplots); library(reshape2); library(pvclust)


# load data ####
gene_list <- read.delim("raw_data/orthology_only_entrez_id.txt")
gene_list <- unique(gene_list[c("gene_symbol", "ac_gene_id")])

HTSeq_ac <- read.delim("raw_data/ac_HTSeq.txt")
gene_length_ac <- read.delim("raw_data/ac_gene_length_median.txt")

merged_ac <- merge(gene_length_ac, HTSeq_ac, by="ac_gene_id")
merged_ac <- merge(merged_ac, gene_list, by="ac_gene_id")
merged_ac <- unique(merged_ac)


# save TF list ####
gene_list_TFs <-merged_ac[c("ac_gene_id", "gene_symbol")]
colnames(gene_list_TFs) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("transcription_factors", nrow(gene_list_TFs)))
species <- c(rep("anolis", nrow(gene_list_TFs)))
gene_list_TFs <- cbind (gene_list, species, gene_list_TFs)
write.table(gene_list_TFs, file = "output_files/gene_list_TFs_anolis.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

merged_ac <- unique(merged_ac[-length(merged_ac)]) # remove gene symbol and just keep the orthology list.


# calc TPM ####
fn.TPM_ac <-function(sample) {
  (sample*75*10^6)/(merged_ac$length_ac*sum((sample*75)/merged_ac$length_ac))
}
TPM_ac <- as.data.frame(sapply(merged_ac[3:length(merged_ac)], fn.TPM_ac))
TPM_ac<- as.data.frame(sqrt(TPM_ac))
TPM_ac <- TPM_ac[rowSums(TPM_ac[-1])>0, ]

data <- data.frame(TPM_ac)


# PCA ####
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('anolis forelimb', length(data)))
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
ggsave("output_figures/extended_fig_2_c_ac_TFs_pca.pdf", width = 6, height = 6, units=c("cm"))

# quick plot with legend
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$digit.number)) + geom_point(size=5)


# PCA bootstraps ####
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_2_c_anolis", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extend_data_2_c_anolis.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA ####
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="anolis TFs"),
       file="output_figures/extended_fig_2_c_ac_TFs_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="alligator TFs") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(1, 
                 2, 5, 7, 3, 6, 14, 15, 9, 4, 11, 10, 13,
                 8, 12))

# heatmap ####
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") + 
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(.9749,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_fig_2_c_ac_TFs_heat_ordered.pdf", width = 8, height = 8, units=c("in"))