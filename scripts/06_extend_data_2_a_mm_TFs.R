# This script corresponds to the Nature Communication resubmission.
# Specifically, it has been modified to extract standard error values for the PCA via bootstrapping. 
#
# May 27 2018
# TAS


# initialize packages ####
library(ggplot2); library(gplots); library(reshape2); library(pvclust); library(bootSVD)


# load data ####
gene_list <- read.delim("raw_data/orthology_only_entrez_id.txt")
gene_list <- gene_list[,(c("mm_gene_id", "gene_symbol"))]

HTSeq_mm <- read.delim("raw_data/mm_HTSeq.txt")
gene_length_mm <- read.delim("raw_data/mm_gene_length_median.txt")

merged_mm <- merge(gene_length_mm, HTSeq_mm, by="mm_gene_id")
merged_mm <- merge(merged_mm, gene_list, by="mm_gene_id") 
merged_mm <- unique(merged_mm)


# Save TF list ####
gene_list_TFs_mouse <-merged_mm[c("mm_gene_id","gene_symbol")]
colnames(gene_list_TFs_mouse) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("transcription_factors", nrow(gene_list_TFs_mouse)))
species <- c(rep("mouse", nrow(gene_list_TFs_mouse)))
gene_list_TFs_mouse <- cbind (gene_list, species, gene_list_TFs_mouse)
gene_list_TFs_mouse <- unique(gene_list_TFs_mouse)
dim(gene_list_TFs_mouse)
write.table(gene_list_TFs_mouse, file = "output_files/gene_list_TFs_mouse.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# calc TPM ####
fn.TPM_mm <-function(sample) {
  (sample*75*10^6)/(merged_mm$length_mm*sum((sample*75)/merged_mm$length_mm))
}
TPM_mm <- as.data.frame(sapply(merged_mm[3:(length(merged_mm)-1)], fn.TPM_mm))
TPM_mm<- as.data.frame(sqrt(TPM_mm))
TPM_mm <- TPM_mm[rowSums(TPM_mm[-1])>0, ]

data <- data.frame(TPM_mm)


# PCA ####
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('mouse forelimb', length(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number), digit.number, stage, limb, row.names = names(data))

pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

g<- qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$digit.number), shape=factor(pheno.digit$limb)) + 
  geom_point(size=2) + 
  theme_bw() + theme(legend.position="none") +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio=1:1,
    text = element_text(size=8))
print(g)
ggsave("output_figures/extended_fig_2_a_mm_TFs_PCA.pdf", width = 6, height = 6, units=c("cm"))


# PCA bootstraps ####
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_2_a_mouse", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extended_data_2_a_mouse.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# clustering ####
pdf("output_figures/extended_fig_2_a_mm_TFs_clust.pdf")
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
plot (hc, print.pv=TRUE, print.num=FALSE, main="mouse - tfs")
dev.off() # for some reason not saving properly?

x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(3, 1, 2,
                 5, 4, 6, # 2s
                 7, 9, 8, 11, #3s,4
                 12, #4
                 15, #5
                 10, 13, 14 # 455
                 ))

# heatmap ####
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(.9749,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_data_2_a_mm_TFs_heat_ordered.pdf", width = 8, height = 8, units=c("in"))