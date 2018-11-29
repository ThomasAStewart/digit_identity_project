# For nature communications resubmission -- tidied the chicken DE analyses and saved outputs into a distinct files.
#
# This makes Figure 5 a, Extended data figure 4 a.
#
# May 27, 2018
# TAS


# initialize packages
library(edgeR); library(ggplot2); library(qvalue); library(gridExtra); library(locfit); library(hexbin); library(pvclust); library(reshape2)


# load data ####
counts <- read.table(file="raw_data/gg_HTseq.txt", sep='\t', as.is=T, header=T, row.names=1)
counts <- counts[1:(dim(counts)[1]-5), ] # clean up bottom few lines
counts <- counts[c("gg_HL_D1_E", "gg_HL_D2_E", "gg_HL_D3_E", "gg_HL_D4_E", 
                   "gg_HL_D1_L", "gg_HL_D2_L", "gg_HL_D3_L", "gg_HL_D4_L")] # Partition to just hindlimbs.


group <- c("D1", "D2", "D3", "D4", "D1", "D2", "D3", "D4")
stage <- c(rep("early",4), rep("late",4))


# specify model for edgeR ####
y.stage <- DGEList(counts = counts, group = group)
keep <- rowSums(cpm(y.stage) > 1) >= 2 # filter low count genes
y.stage <- y.stage[keep, , keep.lib.sizes = F]

y.stage <- calcNormFactors(y.stage) # calculate normalization factor -- by median of reads, while the factor of whole library size is taken internally in the package
design.stage <- model.matrix(~group + stage) # build design matrix removing Stage effect.
y.stage <- estimateDisp(y.stage, design.stage) # estimate dispersion

fit.stage <- glmFit(y.stage, design.stage)


# run the pairwise tests ####
lrt.stage.1vs2 <- glmLRT(fit.stage, contrast = c(0,1,0,0,0)) # compare digit 1 and 2
lrt.stage.2vs3 <- glmLRT(fit.stage, contrast = c(0,-1,1,0,0)) # compare digit 2 and 3
lrt.stage.3vs4 <- glmLRT(fit.stage, contrast = c(0,0,-1,1,0)) # compare digit 3 and 4


# pairwise: count genes for extended data figure 4 ####
# total number of DE genes
gg_12 <- row.names(lrt.stage.1vs2$table[p.adjust(lrt.stage.1vs2$table$PValue, method = "BH") < 0.05,])
gg_23 <- row.names(lrt.stage.2vs3$table[p.adjust(lrt.stage.2vs3$table$PValue, method = "BH") < 0.05,])
gg_34 <- row.names(lrt.stage.3vs4$table[p.adjust(lrt.stage.3vs4$table$PValue, method = "BH") < 0.05,])

# number of TFs
gg_tfs <- read.delim("output_files/gene_list_TFs_chicken.txt")
gg_tfs <- as.vector(gg_tfs["gene_id"])

sum(gg_tfs$gene_id %in% gg_12) # 102
sum(gg_tfs$gene_id %in% gg_23) # 17
sum(gg_tfs$gene_id %in% gg_34) # 22


# save table ####
write.table(lrt.stage.1vs2$table, file="raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.1vs2_all_ONgene.txt", sep= "\t", row.names=TRUE, col.names = TRUE)
write.table(lrt.stage.2vs3$table, file="raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.2vs3_all_ONgene.txt", sep= "\t", row.names=TRUE, col.names = TRUE)
write.table(lrt.stage.3vs4$table, file="raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_lrt.stage.3vs4_all_ONgene.txt", sep= "\t", row.names=TRUE, col.names = TRUE)


# Plot pval distribution ####
# function for combining plots on a page.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot12 <- qplot(lrt.stage.1vs2$table$PValue, geom="histogram", binwidth=.01, xlab="p value") + theme_bw() + labs(title = "D1/D2") + ylim(0,1500)
plot23 <- qplot(lrt.stage.2vs3$table$PValue, geom="histogram", binwidth=.01, xlab="p value") + theme_bw() + labs(title = "D2/D3") + ylim(0,1500)
plot34 <- qplot(lrt.stage.3vs4$table$PValue, geom="histogram", binwidth=.01, xlab="p value") + theme_bw() + labs(title = "D3/D4") + ylim(0,1500)

# This needs to be saved by hand. Still not sure how to explort automatically.
extended_fig_4_a <- multiplot(plot12, plot23, plot34, cols=1)


# ANOVA for chicken hindlimb ####
fit.stage <- glmFit(y.stage, design.stage)
lrt.stage.anova <- glmLRT(fit.stage, coef=2:4)

write.table(lrt.stage.anova$table, file="raw_data/edgeR_outputs/chicken_pairwise_all_ONgene_pvals/chicken_edgeR_anova_NULL_all_ONgene.txt", sep= "\t", row.names=TRUE, col.names = TRUE)

# ANOVA : count genes for figure 5 a ####
# total gene number
gg_anova_count <- row.names(lrt.stage.anova$table[p.adjust(lrt.stage.anova$table$PValue, method = "BH") < 0.05,]) # 1505
# number TFs
sum(gg_tfs$gene_id %in% gg_anova_count) # 102


# ANOVA: save p value distribution plot ####
q <- qplot(lrt.stage.anova$table$PValue, geom="histogram", binwidth=.02, xlab="p value") + theme_bw() + labs(title = "chicken hindlimb ANOVA") + ylim(0,4000)
print(q)
ggsave("output_figures/fig_5_a_gg_HL_pval_anova.pdf", width = 3, height = 3, units=c("in"))
