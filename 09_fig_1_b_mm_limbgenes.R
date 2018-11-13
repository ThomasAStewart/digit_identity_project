# Modified for Nat. Commun. Resubmission
# May 27 2018
# TAS

# initialize packages ####
library(pvclust); library(ggplot2); library(gplots); library(reshape2); library(bootSVD)


# load data ####
mm_gene_symbols <- read.delim("raw_data/ensemble_v85_orthology_symbol.txt", na.strings = "")
mm_gene_symbols <- mm_gene_symbols[c("mm_gene_id", "gene_symbol")]
mm_gene_symbols <- unique(mm_gene_symbols)

patterning_genes <- read.table ("raw_data/limb_patterning_genes.txt", sep="\t", header=T)

gene_list <- merge(mm_gene_symbols, patterning_genes, by="gene_symbol")
gene_list <- unique(na.omit(gene_list))

HTSeq_mm <- read.table("raw_data/mm_HTSeq.txt", sep="\t", header=T)

gene_length_mm <- read.delim ("raw_data/mm_gene_length_median.txt")

merged_mm <- merge(gene_length_mm, HTSeq_mm, by="mm_gene_id")
merged_mm <- merge(merged_mm, gene_list, by="mm_gene_id") # 1765


# save limb patterning gene list ####
gene_list_limb <- na.omit(gene_list)
colnames(gene_list_limb) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("limb_genes", nrow(gene_list_limb)))
species <- c(rep("mouse", nrow(gene_list_limb)))
gene_list_limb <- cbind (gene_list, species, gene_list_limb)
gene_list_limb <- unique(gene_list_limb)
write.table(gene_list_limb, file = "output_files/gene_list_limb_mouse.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

merged_mm <- merged_mm[-length(merged_mm)]
  
# calc TPM ####
fn.TPM_mm <-function(sample) {
  (sample*75*10^6)/(merged_mm$length_mm*sum((sample*75)/merged_mm$length_mm))
}
TPM_mm <- as.data.frame(sapply(merged_mm[3:(length(merged_mm))], fn.TPM_mm))
TPM_mm<- as.data.frame(sqrt(TPM_mm))

data <- data.frame(TPM_mm)


# PCA ####
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('mouse forelimb', length(data)))
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
ggsave("output_figures/fig_1_b_mm_limb_pca.pdf", width = 6, height = 6, units=c("cm"))


# mouse bootstraps ####
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_1_b_mouse", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_1_b_mouse_limb.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA ####
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="mouse limb"),
       file="output_figures/fig_1_b_mm_limb_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="mouse limb") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x<- as.vector(c( 3, 1, 2, 
                 5, 4, 6,
                 7, 8, 9,
                 11, 12,
                 10, 
                 13, 14, 15 ))


# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") + 
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(.887,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_1_b_mm_limb_heat_ordered.pdf", width = 8, height = 8, units=c("in"))