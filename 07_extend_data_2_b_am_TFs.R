# This script corresponds to the Nature Communication resubmission.
# Specifically, it has been modified to extract standard error values for the PCA via bootstrapping. 
#
# May 27 2018
# TAS


# initialize packages
library(ggplot2); library(gplots); library(reshape2); library(pvclust); library(bootSVD)


# load data ####
gene_list <- read.delim("raw_data/orthology_only_entrez_id.txt")
gene_list <- unique(gene_list["gene_symbol"])

am_gene_id <- read.delim("raw_data/am_gene_annotations.txt")
am_gene_id <- am_gene_id[c("am_gene_id", "gene_symbol")]

alligator_tfs <- merge(gene_list, am_gene_id, "gene_symbol")

HTSeq_am <- read.delim("raw_data/am_HTSeq.txt")

gene_length_am <- read.delim("raw_data/am_gene_length_median.txt")

merged <-merge(am_gene_id, HTSeq_am, by="am_gene_id") # Combine counts & ids
merged <-merge(gene_length_am, merged, by="am_gene_id")
merged <-merge(merged, alligator_tfs, by="am_gene_id")


# Save TF list ####
gene_list_TFs <- alligator_tfs[c("am_gene_id","gene_symbol")]
colnames(gene_list_TFs) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("transcription_factors", nrow(gene_list_TFs)))
species <- c(rep("alligator", nrow(gene_list_TFs)))
gene_list_TFs <- cbind (gene_list, species, gene_list_TFs)
gene_list_TFs <- unique(gene_list_TFs)
write.table(gene_list_TFs, file = "output_files/gene_list_TFs_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

merged <- merged[-length(merged)]

# calc TPM ####
fn.TPM_am <-function(sample) {
  (sample*75*10^6)/(merged$length*sum((sample*75)/merged$length))
}
TPM_am <- sapply(merged[, 4:length(merged)], fn.TPM_am)
TPM_am<- as.data.frame(sqrt(TPM_am))
TPM_am <- TPM_am[rowSums(TPM_am[-1])>0, ]

# Batch corrections
attach(TPM_am)
am_FLE <- data.frame(FL_E_D1_1, FL_E_D1_2, FL_E_D1_3,
                     FL_E_D2_1, FL_E_D2_2, FL_E_D2_3,
                     FL_E_D3_1, FL_E_D3_2,
                     FL_E_D4_1, FL_E_D4_2, FL_E_D4_3,
                     FL_E_D5_1, FL_E_D5_2, FL_E_D5_3)
am_FLL <- data.frame(FL_L_D1_1, FL_L_D1_2, FL_L_D1_3,
                     FL_L_D2_1, FL_L_D2_2, FL_L_D2_3,
                     FL_L_D3_1, FL_L_D3_2, FL_L_D3_3,
                     FL_L_D4_1, FL_L_D4_2, FL_L_D4_3,
                     FL_L_D5_1,FL_L_D5_2, FL_L_D5_3)
am_adj.FLE <- (am_FLE-rowMeans(am_FLE))
am_adj.FLL <- (am_FLL-rowMeans(am_FLL))
detach(TPM_am)

data <- data.frame(am_adj.FLE, am_adj.FLL)


# PCA ####
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',2), rep('D4',3), rep('D5',3),
                  rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3), rep('D5',3))
stage <- c(rep('early', 14), rep('late', 15))
limb <- c(rep('alligator forelimb', length(data)))
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
ggsave("output_figures/extended_fig_2_b_am_TFs_PCA.pdf", width = 6, height = 6, units=c("cm"))

# quick plot for legend.
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage), shape=factor(pheno.digit$stage)) + geom_point(size=4)


# PCA bootstraps ####
data2 <- as.matrix(data)
pca<-prcomp(t(data))
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_2_b_alligator", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extended_data_2_b_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# clustering ####
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="alligator TFs"),
       file="output_figures/extended_fig_2_b_am_TFs_clust.pdf")


plot(hc, print.pv=TRUE, print.num=FALSE, main="alligator TFs") # for reference

x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(2,17,15,16,# 1,1,1,1 
                 1,3,19,# 1,1,2,
                 5,4,6,# 2,2,2.
                 18,# 2
                 8,21,23,# all 3s
                 20,22, 7,12,# 2,3,3,5 
                 25,9,10,11,24,26,# fours
                 27,28,29,13,14)) # fives

# heatmap ####
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_fig_2_b_am_TFs_heat_ordered.pdf", width = 8, height = 8, units=c("in"))

