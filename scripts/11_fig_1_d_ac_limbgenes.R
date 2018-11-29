# Modified for Nat. Commun. Resubmission, specifically to get bootstrap support values
# May 27 2018
# TAS

# initialize packaages ####
library(pvclust); library(ggplot2); library(gplots); library(reshape2); library(bootSVD)


# load data ####
ac_gene_symbols <- read.delim("raw_data/ensemble_v85_orthology_symbol.txt", na.strings = "")
ac_gene_symbols <- ac_gene_symbols[c("ac_gene_id", "gene_symbol")]
ac_gene_symbols <- unique(na.omit(ac_gene_symbols))

patterning_genes <- read.delim("raw_data/limb_patterning_genes.txt")
gene_list <- merge(patterning_genes, ac_gene_symbols, by="gene_symbol")

HTSeq_ac <- read.delim("raw_data/ac_HTSeq.txt")

gene_length_ac <- read.delim("raw_data/ac_gene_length_median.txt")

merged_ac <- merge(gene_length_ac, HTSeq_ac, by="ac_gene_id")
merged_ac <- merge(merged_ac, gene_list, by="ac_gene_id")
merged_ac <- unique(merged_ac[-length(merged_ac)])


# save limb patterning gene list ####
gene_list_limb <- na.omit(gene_list[c(2,1)])
colnames(gene_list_limb) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("limb_genes", nrow(gene_list_limb)))
species <- c(rep("anolis", nrow(gene_list_limb)))
gene_list_limb <- cbind (gene_list, species, gene_list_limb)
gene_list_limb <- unique(gene_list_limb)
dim(gene_list_limb)
write.table(gene_list_limb, file = "output_files/gene_list_limb_anolis.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# calc TPM ####
fn.TPM_ac <-function(sample) {
  (sample*75*10^6)/(merged_ac$length_ac*sum((sample*75)/merged_ac$length_ac))
}
TPM_ac <- as.data.frame(sapply(merged_ac[3:length(merged_ac)], fn.TPM_ac))
TPM_ac<- as.data.frame(sqrt(TPM_ac))
TPM_ac <- TPM_ac[rowSums(TPM_ac[-1])>0, ]

data <- data.frame(TPM_ac)


# PCA
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('anolis forelimb', length(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number),
                          digit.number, stage, limb,
                          row.names = names(data))
pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

g <- qplot(x=PC1, y=PC2, data=scores, 
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
ggsave("output_figures/fig_1_d_ac_limb_pca_12.pdf", width = 6, height = 6, units=c("cm"))


# PCA bootstrap ####
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_1_d_anolis", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_1_d_anolis_limb.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# clustering
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="mouse limb"),
       file="output_figures/fig_1_d_ac_limb_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="mouse limb") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x<- as.vector(c(1,
                2,3,5, 
                4,6,
                13,14,15,
                7,10,11,9,8,12))


# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") + 
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(.887,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_1_d_ac_limb_heat_ordered.pdf", width = 8, height = 8, units=c("in"))

