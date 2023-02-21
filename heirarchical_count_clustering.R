library(tidyverse)


gct_deg <- t(norm_counts$norm.dat) %>% as.data.frame(.) %>% select(top_found_genes) %>% t(.)
gene_dist <- dist(gct_deg)

gene_hclust <- hclust(gene_dist, method = "ward.D2")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 1600, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram


#gene_cluster <- cutree(gene_hclust, k = 2)

gene_cluster <- cutree(gene_hclust, k = 2) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

head(gene_cluster)

test <- gene_cluster[which(gene_cluster$cluster==2),"gene"]

for (i in 1:nrow(test)) {
  cat(test$gene[i], "\n")
}
branch_genes <- test[1:nrow(test),"gene"]
write.table(branch_genes, file = "branch_genes_Clust2_12k.txt", row.names = F, col.names = F, quote = F)

########################################
#Using Pheatmap

library(limma)

# Subset SNM count matrix to most significant genes
sig_norm_counts <- norm_counts$norm.dat[row.names(gct) %in% top_found_genes,]
sig_norm_counts <- gct[row.names(gct) %in% top_found_genes,]

bio.var.mod = model.matrix(~  steadyState.vs.voc , data=pheno)

# Fit a linear model to the data
#fit <- lmFit(gct, bio.var.mod)
#fit <- lmFit(norm_counts$norm.dat, bio.var.mod)
fit <- lmFit(sig_norm_counts, bio.var.mod)

# Perform empirical Bayes moderation to shrink the standard errors of the coefficients
ebfit <- eBayes(fit)

# Extract the log fold changes and associated p-values
lfc <- topTable(ebfit, coef = 2, number = nrow(sig_norm_counts), adjust.method = "BH")

# View the log fold changes and associated p-values
head(lfc)

library(pheatmap)

# Z-score normalize the gene count matrix
# Manually scaling instead of using pheatmaps scaling feature results in a better image because
# pheatmap will make the color scale too large to easily discriminate differences
#gct_zscore <- t(scale(t(gct)))
#gct_zscore <- t(scale(t(norm_counts$norm.dat)))
gct_zscore <- t(scale(t(sig_norm_counts)))

# Load the list of differentially expressed genes
differentially_expressed_genes <- rownames(lfc[which(sqrt((lfc$logFC)^2) > 1 & lfc$adj.P.Val < 0.0001), ])
#differentially_expressed_genes <- rownames(lfc[which(lfc$adj.P.Val < 0.0001), ])
#differentially_expressed_genes <- rownames(lfc[which(sqrt((lfc$logFC)^2) > 1), ])

# Subset the Z-score normalized gene count matrix to include only the differentially expressed genes

gct_zscore_filtered <- gct_zscore[rownames(gct_zscore) %in% differentially_expressed_genes, ]

######## Plot the heat map ##########
# Set up column annotations
anno <- data.frame(pheno[,c('steadyState.vs.voc', "Chronic.Pain.")])
rownames(anno) <- row.names(pheno)

pheatmap(gct_zscore_filtered, clustering_method = "ward.D2", fontsize = 7, annotation_col = anno, 
         show_rownames = F, show_colnames = F,
)



# Using SNM
gct_snm <- t(norm_counts$norm.dat) %>% as.data.frame(.) %>% select(top_found_genes) %>% t(.)
gct_zcore_snm <- t(scale(t(gct_snm)))
pheatmap(gct_zcore_snm, clustering_method = "ward.D2", fontsize = 7, annotation_col = anno, 
         show_rownames = F, show_colnames = F,)
