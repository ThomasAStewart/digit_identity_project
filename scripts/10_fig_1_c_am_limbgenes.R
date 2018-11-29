# Modified for Nat. Commun. Resubmission, specifically to get bootstrap support values and also to save the alligator orthology list for each data type.
# May 20, 2018
# TAS


# initialize packages
library(pvclust); library(ggplot2); library(gplots); library(reshape2); library(bootSVD)


# load data ####
HTSeq_am <- read.delim("raw_data/am_HTSeq.txt")

gene_length_am <- read.delim("raw_data/am_gene_length_median.txt")

am_gene_id <- read.delim("raw_data/am_gene_annotations.txt", na.strings = "")
am_gene_id <- am_gene_id[c("am_gene_id", "gene_symbol")]

gene_list <- read.table("raw_data/limb_patterning_genes.txt", sep = "\t", header = TRUE)

merged <-merge(am_gene_id, HTSeq_am, by="am_gene_id")
merged <-merge(gene_length_am, merged, by="am_gene_id")
merged_am <-merge(gene_list, merged, by="gene_symbol") 
merged_am<- merged_am[!duplicated(merged_am[,1]),]


# save limb patterning gene list ####
gene_list_limb <- merged_am[c("am_gene_id", "gene_symbol")]
colnames(gene_list_limb) <- c("gene_id", "gene_symbol")
gene_list <- c(rep("limb_genes", nrow(gene_list_limb)))
species <- c(rep("alligator", nrow(gene_list_limb)))
gene_list_limb <- cbind (gene_list, species, gene_list_limb)
gene_list_limb <- unique(gene_list_limb)
write.table(gene_list_limb, file = "output_files/gene_list_limb_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# calc TPM
fn.TPM_am <-function(sample) {
  (sample*75*10^6)/(merged_am$length_am*sum((sample*75)/merged_am$length_am))
}
TPM_am <- sapply(merged_am[4:(length(merged_am))], fn.TPM_am)
TPM_am<- as.data.frame(sqrt(TPM_am)) # Stabilize variance
TPM_am <- TPM_am[rowSums(TPM_am[-1])>0, ] 

# bulk correction for the 2 different stages
TPM_am_early <- TPM_am[1:14]
TPM_am_late <- TPM_am[15:length(TPM_am)]
adj.FLE <- (TPM_am_early-rowMeans(TPM_am_early))
adj.FLL <- (TPM_am_late-rowMeans(TPM_am_late))

data <- as.matrix(cbind(adj.FLE, adj.FLL))

# PCA
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',2), rep('D4',3), rep('D5',3), 
                  rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3), rep('D5',3))
stage <- c(rep('early', 14), rep('late', 15))
limb <- c(rep('alligator forelimb', ncol(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number), digit.number, stage, limb, row.names = names(data))

pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

g <- qplot(x=PC1, y=PC2, data=scores, 
           colour=factor(pheno.digit$digit.number),
           shape=factor(pheno.digit$limb),
           alpha=factor(pheno.digit$stage)) + 
  geom_point(size=2) + 
  theme_bw() + theme(legend.position="none") +
  theme( plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         aspect.ratio=1:1,
         text = element_text(size=8))
print(g)
ggsave("output_figures/fig_1_c_alligator_limb_genes_PCA.pdf", width = 6, height = 6, units=c("cm"))

# quick plots with legends.
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage)) + geom_point(size=5)


# bootstrap PCA #### 
boots <-bootPCA(centerSamples = TRUE, Y=data, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_1_c", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_1_c_alligator_limb.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="alligator limb"),
       file="output_figures/fig_1_c_am_limb_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="alligator limb") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(2, 17, 15, 16,  1,  3, #1s
                 5, 19, 18, 20, #2s
                 4, 6, #2s
                 22, 7, 8, 21, 23, #3s
                 25, 24, 26, 11, 9, 10, #4s 
                 12, 27, 28, 29, 13, 14)) #5s))

# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[,x]), method=c("pearson"))), fill=value, geom="tile") + 
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_1_c_am_limb_heat.pdf", width = 8, height = 8, units=c("in"))

