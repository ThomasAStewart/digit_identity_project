# This script corresponds to the Nature Communication resubmission.
# Specifically, it has been modified to extract standard error values for the PCA via bootstrapping. 
#
# May 27 2018
# TAS


# initialize packages ####
library(pvclust); library(ggplot2); library(reshape2); library(bootSVD)


# load data ####
HTSeq_am <- read.delim("raw_data/am_HTSeq.txt")
gene_length_am <- read.delim("raw_data/am_gene_length_median.txt")

merged_am <-merge(gene_length_am, HTSeq_am, by="am_gene_id")


# calc TPM ####
fn.TPM_am <-function(sample) {
  (sample*75*10^6)/(merged_am$length_am*sum((sample*75)/merged_am$length_am))
}
TPM_am <- sapply(merged_am[3:(length(merged_am))], fn.TPM_am)
TPM_am<- as.data.frame(sqrt(TPM_am)) # Stabilize variance
TPM_am <- TPM_am[rowSums(TPM_am[-1])>0, ] 


# bulk correction for the 2 different stages
TPM_am_early <- TPM_am[1:14]
TPM_am_late <- TPM_am[15:length(TPM_am)]
adj.FLE <- (TPM_am_early-rowMeans(TPM_am_early))
adj.FLL <- (TPM_am_late-rowMeans(TPM_am_late))

data <- as.matrix(cbind(adj.FLE, adj.FLL))


# PCA #### 
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
ggsave("output_figures/extended_fig_1_b_am_PCA.pdf", width = 6, height = 6, units=c("cm"))

# quick plot, to have legend for stages.
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage)) + geom_point(size=5)


# PCA bootstraps ### 
boots <-bootPCA(centerSamples = TRUE, Y=data, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("extend_data_1_b", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/extended_data_1_b_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# HCA ####
pdf("output_figures/extended_fig_1_b_am_allgenes_clust.pdf")
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
plot (hc, print.pv=TRUE, print.num=FALSE, main="alligator - all genes")
dev.off()

x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(1, 3, 17, 15, 16, # 1s.
                  18, 19,  #2s
                  2, 20, 5, 22, # 1222
                  4, 6,  #22
                  21, 8, 23, # 333 
                  7, 12, # 3, 5
                  9, 10, 11, 25, 24, 26, # 4s
                  13, 28, 29, 14, 27 #5s
                  )) # edited for readability in final figure, once HCA rotated around some nodes

# heatmap ####
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[,x]), method=c("pearson"))), fill=value, geom="tile") + 
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/extended_data_1_b_am_heatmap_ordered.pdf", width = 8, height = 8, units=c("in"))
