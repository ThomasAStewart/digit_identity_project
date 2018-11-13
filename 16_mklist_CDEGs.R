# For Nature Communications revision.
# This script will
# (1)  generate the CDEG list - Figure 5 B
# (2) Run all of the associated analyses/plots -- Figure 5 C, 6A, 6D, 
#
# May 27, 2019
# TAS

library(bootSVD); library(pvclust); library(ggplot2)

# Generate CDEG list by comparing ANOVAs for mouse chicken alligator ####
am_anova <- read.delim("raw_data/edgeR_outputs/alligator_all_ONgene_pvals/alligator_edgeR_anova_NULL_all_ONgene.txt")
am_anova<-cbind(row.names(am_anova), am_anova)
row.names(am_anova) <- NULL
names(am_anova)[names(am_anova)=="row.names(am_anova)"] <- "am_gene_id"
am_anova <- am_anova[which(am_anova$FDR<0.05),] # 4677

mm_anova <- read.delim("raw_data/edgeR_outputs/mouse_all_ONgene_pvals/mouse_edgeR_anova_NULL_all_ONgene.txt")
mm_anova<-cbind(row.names(mm_anova), mm_anova)
row.names(mm_anova) <- NULL
names(mm_anova)[names(mm_anova)=="row.names(mm_anova)"] <- "mm_gene_id"
mm_anova <- mm_anova[which(mm_anova$FDR<0.05),] # 2950

gg_anova <- read.delim("raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_anova_NULL_all_ONgene.txt")
gg_anova<-cbind(row.names(gg_anova), gg_anova)
row.names(gg_anova) <- NULL
names(gg_anova)[names(gg_anova)=="row.names(gg_anova)"] <- "gg_gene_id"
gg_anova <- gg_anova[(p.adjust(gg_anova$PValue, method = "BH") < 0.05),] # 1505

one2one_list <- read.delim("output_files/amggmm_TFs_1049.txt")


# these are the ANOVA DE genes for each species that are found in the 1:1:1 list
tf_common_mm <- one2one_list[(one2one_list$mm_gene_id %in% mm_anova$mm_gene_id),1]
tf_common_am <- one2one_list[(one2one_list$am_gene_id %in% am_anova$am_gene_id),1]
tf_common_gg <- one2one_list[(one2one_list$gg_gene_id %in% gg_anova$gg_gene_id),1]

# figure 5 b, venn diagram total values for each spp.
length(tf_common_mm) # 242
length(tf_common_am) # 282
length(tf_common_gg) # 128

# generate the CDEG list by intersecting the species-specific lists for the 1048
tf_common <- one2one_list[(one2one_list$mm_gene_id %in% mm_anova$mm_gene_id),] # the mouse genes in 1:1:1
tf_common <- tf_common[(tf_common$am_gene_id %in% am_anova$am_gene_id),] # mouse and alligator genes in 1:1:1
tf_common <- tf_common[(tf_common$gg_gene_id %in% gg_anova$gg_gene_id),]

write.table(tf_common, file = "output_files/CDEGs_49.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Figure 5 b ####
tf_common_symbol <- sort(tf_common[,1]) # 49 CDEGs in alphabetical order
length(unique(intersect(tf_common_mm, tf_common_gg)))-49  # shared between mouse and chicken, excluding center
length(unique(intersect(tf_common_mm, tf_common_am)))-49  # shared between mouse and alligator, excluding center
length(unique(intersect(tf_common_am, tf_common_gg)))-49  # shared between alligator and chicken, excluding center


# Below are the analyses of CDEGs for all species
# (1) saving the gene lists, (2) PCA, including bootstraps; (3) HCA; (4) heatmaps; (5) projection of chicken samples into the PCA of other species.


# CDEG mouse ####
HTSeq_mm <- read.delim("raw_data/mm_HTSeq.txt")
gene_length_mm <- read.delim("raw_data/mm_gene_length_median.txt")

merged_mm <- merge(gene_length_mm, HTSeq_mm, by="mm_gene_id")
merged_mm <- merge(tf_common, merged_mm, by="mm_gene_id")
merged_mm <- merged_mm[match(tf_common$mm_gene_id, merged_mm$mm_gene_id),] # This sorts the rows alphabetically by row name, so that the TPM data sets can be combined later, if desired.

# calc TPM
fn.TPM_mm <-function(sample) {
  (sample*75*10^6)/(merged_mm$length_mm*sum((sample*75)/merged_mm$length_mm))
}
TPM_mm <- as.data.frame(sapply(merged_mm[6:length(merged_mm)], fn.TPM_mm))
TPM_mm<- as.data.frame(sqrt(TPM_mm)) # normalize var.

data <- data.frame(TPM_mm)

# PCA
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('mouse forelimb', length(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number), digit.number, stage, limb, row.names = names(data))
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
ggsave("output_figures/fig_5_c_mouse_49_CDEGs_PCA.pdf", width = 6, height = 6, units=c("cm"))

# bootstrap PCA 
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_5_c_mouse.txt", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_5_c_mouse.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# HCA 
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="anolis TFs"),
       file="output_figures/fig_5_c_mouse_49_CDEGs_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="chicken all genes") # for reference
x <- as.vector( noquote(hc$hclust$order) ) # so you can see original order for heatmap
x <- as.vector (c(3, 1, 2, 
                  5, 4, 6, 
                  9, 7, 8, 
                  12, 11, 10, 13, 14, 15)) # re-ordered so that it reads 1 > 5 on the final figure.

# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=(c(0.75,1.0))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_5_c_mouse_49_CDEGs_heat.pdf", width = 8, height = 8, units=c("in"))


# CDEG alligator ####
HTSeq_am <- read.delim("raw_data/am_HTSeq.txt")
gene_length_am <- read.delim("raw_data/am_gene_length_median.txt")

merged_am <- merge(gene_length_am, HTSeq_am, by="am_gene_id")
merged_am <- merge(tf_common, merged_am, by="am_gene_id")
merged_am <- merged_am[match(tf_common$am_gene_id, merged_am$am_gene_id),] # This sorts the rows alphabetically by row name, so that the TPM data sets can be combined later, if desired.

# calc TPM
fn.TPM_am <-function(sample) {
  (sample*75*10^6)/(merged_am$length_am*sum((sample*75)/merged_am$length_am))
}
TPM_am <- as.data.frame(sapply(merged_am[6:(length(merged_am))], fn.TPM_am))
TPM_am<- as.data.frame(sqrt(TPM_am)) # normalize variance

# bulk correction
TPM_am_early <- TPM_am[1:14]
TPM_am_late <- TPM_am[15:length(TPM_am)]
adj.FLE <- (TPM_am_early-rowMeans(TPM_am_early))
adj.FLL <- (TPM_am_late-rowMeans(TPM_am_late))

data <- data.frame(adj.FLE, adj.FLL)

# PCA
digit.number <- c(rep('D1',3), rep('D2',3),  rep('D3',2), rep('D4',3), rep('D5',3),
                  rep('D1',3), rep('D2',3),  rep('D3',3), rep('D4',3), rep('D5',3))
stage <- c(rep('early', 14), rep('late', 15))
limb <- c(rep('alligator forelimb', length(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number),
                          digit.number, stage, limb,
                          row.names = names(data))
pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage)) + geom_point(size=4)
g <- qplot(x=PC1, y=PC2, data=scores, 
          colour=factor(pheno.digit$digit.number),
          shape=factor(pheno.digit$stage)) + 
  geom_point(size=2) + 
  theme_bw() + theme(legend.position="none") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio=1:1,
    text = element_text(size=8))
print(g)
ggsave("output_figures/fig_5_c_alligator_49_CDEGs_PCA.pdf", width = 6, height = 6, units=c("cm"))

# PCA bootstraps
data2 <- as.matrix(data)
pca <- prcomp(t(data))
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_5_c_alligator", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_5_c_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# HCA
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="alligator CDEGs"),
       file="output_figures/fig_5_c_alligator_49_CDEGs_clust.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE, main="alligator CDEGs")
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(1, 3, 15, 16, 2, 17,
                 4, 6, 5, 19, 18, 20, 
                 8, 21, 23, 7, 22,
                 11, 9, 10, 26, 24, 25,
                 12, 27, 28, 29, 13, 14)) # reordered for final figure.

# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_5_c_alligator_49_CDEGs_heat.pdf", width = 8, height = 8, units=c("in"))


# CDEG chicken ####
HTSeq_gg <- read.delim ("raw_data/gg_HTseq.txt")
gg_gene_id <- read.delim ("raw_data/gal_names_length_median.txt", col.names=c("gg_gene_id","length_gg", "gene_symbol"))

merged <- merge(gg_gene_id, HTSeq_gg, by="gg_gene_id")
merged_gg <- merge(tf_common, merged, by="gg_gene_id")
merged_gg <- merged_gg[match(tf_common$gg_gene_id, merged_gg$gg_gene_id),] # This sorts the rows alphabetically by row name, so that the TPM data sets can be combined later, if desired.

 
# Calc TPM
fn.TPM_gg <-function(sample) {
  (sample*34*10^6)/(merged_gg$length_gg*sum((sample*34)/merged_gg$length_gg))
}
TPM_gg <- sapply(merged_gg[, 7:length(merged_gg)], fn.TPM_gg) # Apply function
TPM_gg<- as.data.frame(sqrt(TPM_gg)) # Variance stabilization

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

data <- data.frame( gg_adj.FLE, gg_adj.HLE, gg_adj.FLL, gg_adj.HLL)

# PCA
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

qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage)) + geom_point(size=4)

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
ggsave("output_figures/fig_6_a_chicken_49_CDEGs_PCA_all_limbs.pdf", width = 6, height = 6, units=c("cm"))

# PCA bootstraps
data2 <- as.matrix(data)
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_6_a_chicken", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_6_a_chicken.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# HCA
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="chicken CDEGs"),
       file="output_figures/fig_6_a_chicken_49_CDEGs_clust_all_limbs.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE)
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(1, 8, 4, 11, 5, 12, 6, 13, 2, 9, 7, 14,3, 10))

# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_6_a_chicken_49_CDEGs_heat_all_limbs.pdf", width = 8, height = 8, units=c("in"))


# CDEG chicken hindlimb ####
data <- data.frame (gg_adj.HLE, gg_adj.HLL)

# PCA
digit.number <- c('D1', 'D2', 'D3', 'D4',
                  'D1', 'D2', 'D3', 'D4')
stage <- c(rep('early', 4), rep('late', 4))
limb <- c(rep('chicken hindlimb', 8))
pheno.digit <- data.frame(sample = 1:length(digit.number),
                          digit.number, stage, limb,
                          row.names = names(data))
pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)

qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$stage)) + geom_point(size=4)

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
ggsave("output_figures/fig_5_c_chicken_49_CDEGs_PCA_hindlimb.pdf", width = 6, height = 6, units=c("cm"))

# PCA bootstrap
data2 <- as.matrix(data)
pca<-prcomp(t(data))
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_5_c_chicken", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_5_c_chicken.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# HCA
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="chicken CDEGs"),
       file="output_figures/fig_5_c_chicken_49_CDEGs_clust_hindlimb.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE)
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x <- as.vector(c(1, 5, 2, 6, 3, 7, 4, 8))

# HEATMAP
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_5_c_chicken_49_CDEGs_heat_hindlimbs.pdf", width = 8, height = 8, units=c("in"))


# CDEG anolis ####
HTSeq_ac <- read.delim("raw_data/ac_HTSeq.txt")
gene_length_ac <- read.delim("raw_data/ac_gene_length_median.txt")

one2one_list <- read.delim("raw_data/orthology_only_entrez_id.txt", na.strings="")
one2one_list <- one2one_list[c("ac_gene_id", "gene_symbol")]
anolis_CDEGs <- merge(tf_common, one2one_list, by="gene_symbol") 
anolis_CDEGs <- na.omit(anolis_CDEGs[c("gene_symbol", "ac_gene_id")])
anolis_CDEGs <- anolis_CDEGs[!(duplicated(anolis_CDEGs$gene_symbol) | duplicated(anolis_CDEGs$gene_symbol, fromLast = TRUE)), ] # This removes all genes that had more than one gene symbol (ie., not 1:1)

merged_ac <- merge(gene_length_ac, HTSeq_ac, by="ac_gene_id")
merged_ac <- merge(anolis_CDEGs, merged_ac, by="ac_gene_id")
merged_ac <- merged_ac[order(merged_ac$gene_symbol),] # alphebetizes the list by gene symbol.

# calc TPM
fn.TPM_ac <-function(sample) {
  (sample*75*10^6)/(merged_ac$length_ac*sum((sample*75)/merged_ac$length_ac))
}
TPM_ac <- as.data.frame(sapply(merged_ac[4:(length(merged_ac))], fn.TPM_ac))
TPM_ac<- as.data.frame(sqrt(TPM_ac)) # Stabilize variance

data <- data.frame(TPM_ac)

# PCA
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3), rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('anolis forelimb', length(data)))
pheno.digit <- data.frame(sample = 1:length(digit.number),
                          digit.number, stage, limb,
                          row.names = names(data))
pca<-prcomp(t(data))
scores <- data.frame(pheno.digit, pca$x[,1:3])
summary(pca)
qplot(x=PC1, y=PC2, data=scores, colour=factor(pheno.digit$digit.number)) + geom_point(size=4)
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
ggsave("output_figures/fig_5_c_anolis_42_CDEGs_PCA.pdf", width = 6, height = 6, units=c("cm"))

# PCA bootstrap
data2 <- as.matrix(data)
pca<-prcomp(t(data))
boots <-bootPCA(centerSamples = TRUE, Y=data2, K=2, B=1000)
boots <- matrix(unlist(boots$LD_moments$sdPCs), ncol = 2, byrow = TRUE)
supplementary_table <- data.frame(rep("fig_5_c_anolis", length(pca$x[,1])), row.names(pca$x), pca$x[,1], boots[,1], pca$x[,2], boots[,2])
colnames(supplementary_table) <- c("plot", "sample", "PC1_loading", "PC1_SE", "PC2_loading", "PC2_SE")
write.table(supplementary_table, file = "output_files/fig_5_c_anolis.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# HCA
hc <- pvclust(data, method.hclust="average", method.dist="correlation", nboot=1000)
ggsave(plot (hc, print.pv=TRUE, print.num=FALSE, main="anolis CDEGs"),
       file="output_figures/fig_5_c_anolis_42_CDEGs_clust_hindlimb.pdf")

plot(hc, print.pv=TRUE, print.num=FALSE)
x <- as.vector( noquote(hc$hclust$order) ) # Saves order for heatmap
x<- as.vector( c(1,  2,  3, 
                 5,4, 6, 15, 13, 14, # right to left 2+5 with one two moved.
                 9,8,12,11,7,10))
# heatmap
q <- qplot (x=Var1, y=Var2, data= (melt(cor(data[x]), method=c("pearson"))), fill=value, geom="tile") +
  scale_fill_gradientn(colours = c("steelblue4","white", "firebrick3"), limits=(c(0.75,1.0))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(q)
ggsave("output_figures/fig_5_c_anolis_42_CDEGs_heat.pdf", width = 8, height = 8, units=c("in"))


# Figure 6 D -- projecting chicken onto other PCAs. ####
ac_adj <- (TPM_ac-rowMeans(TPM_ac))
mm_adj <- (TPM_mm-rowMeans(TPM_mm))
am_adj <- cbind(adj.FLE, adj.FLL)
gg_adj <- cbind(gg_adj.FLE, gg_adj.FLL, gg_adj.HLE, gg_adj.HLL) 

digit.number <- c('D2', 'D3', 'D4', 
                  'D2', 'D3', 'D4', 
                  'D1', 'D2', 'D3', 'D4',  
                  'D1', 'D2', 'D3', 'D4')
stage <- c(rep('early', 3), rep('late', 3), rep('early', 4), rep('late', 4))
limb <- c(rep('chicken forelimb', 6),  rep('chicken hindlimb', 8))
pheno.digit_gg <- data.frame(digit.number, stage, limb,
                             row.names = names(gg_adj))


# Mouse, chicken on top ####
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('mouse forelimb', length(mm_adj)))
pheno.digit_mm <- data.frame(digit.number, stage, limb, row.names = names(mm_adj))

pca<-prcomp(t(mm_adj))
scores <- data.frame(pheno.digit_mm, pca$x[,1:3])
scale_mm <- scale(t(gg_adj), pca$center, pca$scale) %*% pca$rotation 
scores_ggam <- data.frame(pheno.digit_gg, scale_mm[,1:3])
scores_combined <- rbind(scores,scores_ggam)

g <- qplot(x=PC1, y=PC2, data=scores_combined, shape=factor(scores_combined$limb), colour=factor(scores_combined$digit.number)) + geom_point(size=4)
print(g)
ggsave("output_figures/fig_6_d_chicken_on_mouse_PCA.pdf", width = 6, height = 4, units=c("in"))


# alligator, chicken on top ####
digit.number <- c(rep('D1',3), rep('D2',3),  rep('D3',2), rep('D4',3), rep('D5',3),
                  rep('D1',3), rep('D2',3),  rep('D3',3), rep('D4',3), rep('D5',3))
stage <- c(rep('early', 14), rep('late', 15))
limb <- c(rep('alligator forelimb', length(am_adj)))
pheno.digit_am <- data.frame(digit.number, stage, limb,
                             row.names = names(am_adj))

pca<-prcomp(t(am_adj))
scores <- data.frame(pheno.digit_am, pca$x[,1:3])
scale_am <- scale(t(gg_adj), pca$center, pca$scale) %*% pca$rotation 
scores_ggam <- data.frame(pheno.digit_gg, scale_am[,1:3])
scores_combined <- rbind(scores,scores_ggam)

g <- qplot(x=PC1, y=PC2, data=scores_combined, shape=factor(scores_combined$limb), colour=factor(scores_combined$digit.number)) + geom_point(size=4)
print(g)
ggsave("output_figures/fig_6_d_chicken_on_alligator_PCA.pdf", width = 6, height = 4, units=c("in"))


# Anolis, chicken on top ####
# Calculate the new TPM for chicken with 41 gene list.
CDEG41 <- merged_ac[c("gene_symbol")] 
colnames(CDEG41) <- c("gene_symbol.x")
CDEG41_gg <- merge(merged_gg, CDEG41, by="gene_symbol.x")

CDEG41_gg_fn <-function(sample) {
  (sample*34*10^6)/(CDEG41_gg$length_gg*sum((sample*34)/CDEG41_gg$length_gg))
}
TPM_CDEG4_gg <- sapply(CDEG41_gg[, 7:length(CDEG41_gg)], CDEG41_gg_fn) # Apply function
TPM_CDEG4_gg<- as.data.frame(sqrt(TPM_CDEG4_gg)) # Variance stabilization

# batch correction
attach(TPM_CDEG4_gg)
gg_FLE <- data.frame(gg_FL_D1_E, gg_FL_D2_E, gg_FL_D3_E) 
gg_FLL <- data.frame(gg_FL_D1_L, gg_FL_D2_L, gg_FL_D3_L) 
gg_HLE <- data.frame(gg_HL_D1_E, gg_HL_D2_E, gg_HL_D3_E, gg_HL_D4_E) 
gg_HLL <- data.frame(gg_HL_D1_L, gg_HL_D2_L, gg_HL_D3_L, gg_HL_D4_L) 
gg_adj.FLE <- (gg_FLE-rowMeans(gg_FLE))
gg_adj.FLL <- (gg_FLL-rowMeans(gg_FLL))
gg_adj.HLE <- (gg_HLE-rowMeans(gg_HLE))
gg_adj.HLL <- (gg_HLL-rowMeans(gg_HLL))
detach(TPM_CDEG4_gg)
gg_adj <- cbind(gg_adj.FLE, gg_adj.FLL, gg_adj.HLE, gg_adj.HLL) 

# run PCA on anolis adjusted
digit.number <- c(rep('D1',3), rep('D2',3), rep('D3',3), rep('D4',3),  rep('D5',3))
stage <- c(rep('early', 15))
limb <- c(rep('anolis forelimb', length(ac_adj)))
pheno.digit_ac <- data.frame(digit.number, stage, limb, row.names = names(ac_adj))

pca<-prcomp(t(ac_adj))
scores <- data.frame(pheno.digit_ac, pca$x[,1:3])
scale_ac <- scale(t(gg_adj), pca$center, pca$scale) %*% pca$rotation 
scores_ggac <- data.frame(pheno.digit_gg, scale_ac[,1:3])
scores_combined <- rbind(scores,scores_ggac)

qplot(x=PC1, y=PC2, data=scores_combined, shape=factor(scores_combined$limb), colour=factor(scores_combined$digit.number)) + geom_point(size=4)
ggsave("output_figures/fig_6_d_chicken_on_anolis_PCA.pdf", width = 6, height = 4, units=c("in"))


# Saving CDEG list ####
tf_common_mm <- tf_common[c("mm_gene_id", "gene_symbol")]
tf_common_am <- tf_common[c("am_gene_id", "gene_symbol")]
tf_common_gg <- tf_common[c("gg_gene_id", "gene_symbol")]
tf_common_ac <- merged_ac[c("ac_gene_id", "gene_symbol")]

gene_list_cdeg <- tf_common_mm
names(gene_list_cdeg) <- c("gene_id", "gene_symbol") 
gene_list <- c(rep("CDEG", nrow(gene_list_cdeg)))
species <- c(rep("mouse", nrow(gene_list_cdeg)))
gene_list_cdeg <- cbind (gene_list, species, gene_list_cdeg)
gene_list_cdeg <- unique(gene_list_cdeg)
write.table(gene_list_cdeg, file = "output_files/gene_list_CDEG_mouse.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gene_list_cdeg <- tf_common_am
names(gene_list_cdeg) <- c("gene_id", "gene_symbol") 
gene_list <- c(rep("CDEG", nrow(gene_list_cdeg)))
species <- c(rep("alligator", nrow(gene_list_cdeg)))
gene_list_cdeg <- cbind (gene_list, species, gene_list_cdeg)
gene_list_cdeg <- unique(gene_list_cdeg)
write.table(gene_list_cdeg, file = "output_files/gene_list_CDEG_alligator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gene_list_cdeg <- tf_common_ac
names(gene_list_cdeg) <- c("gene_id", "gene_symbol") 
gene_list <- c(rep("CDEG", nrow(gene_list_cdeg)))
species <- c(rep("anolis", nrow(gene_list_cdeg)))
gene_list_cdeg <- cbind (gene_list, species, gene_list_cdeg)
gene_list_cdeg <- unique(gene_list_cdeg)
write.table(gene_list_cdeg, file = "output_files/gene_list_CDEG_anolis.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

gene_list_cdeg <- tf_common_gg
names(gene_list_cdeg) <- c("gene_id", "gene_symbol") 
gene_list <- c(rep("CDEG", nrow(gene_list_cdeg)))
species <- c(rep("chicken", nrow(gene_list_cdeg)))
gene_list_cdeg <- cbind (gene_list, species, gene_list_cdeg)
gene_list_cdeg <- unique(gene_list_cdeg)
write.table(gene_list_cdeg, file = "output_files/gene_list_CDEG_chicken.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

