
# Calculate log fc between groups

# Subset to the top significant genes from SNM
sig_norm_counts <- norm_counts$norm.dat[row.names(gct) %in% top_found_genes,]

# Add the absolute value of the minimum entry in each row to all entries in that row to remove negative numbers.
sig_norm_counts <- sig_norm_counts + abs(apply(sig_norm_counts, 1, min, na.rm = TRUE))

# Get the average expression for each gene across samples in each group
inpatient_avg <- apply(sig_norm_counts[, pheno$steadyState.vs.voc == "Inpatient VOC"], 1, mean)
steady_state_avg <- apply(sig_norm_counts[, pheno$steadyState.vs.voc == "Steady State"], 1, mean)

# Create a new data frame with the averages
gene_averages <- data.frame(
  Gene = rownames(sig_norm_counts),
  "Inpatient VOC" = inpatient_avg,
  "Steady State" = steady_state_avg,
  "Log2FC_Inpatient VOC" = log2(inpatient_avg/steady_state_avg)
)

# Z-score normalize
sig_norm_counts_zscore <- t(scale(t(sig_norm_counts)))

# Subset z-score normalized matrix to genes with |log2FC| > 1
sig_norm_counts_zscore <- sig_norm_counts_zscore[row.names(sig_norm_counts_zscore) %in% 
                              gene_averages[which(abs(gene_averages$Log2FC_Inpatient.VOC)>1),"Gene"],]

# Inverse quantile normalization
library(preprocessCore)
sig_norm_counts_iqn <- t(normalize.quantiles(t(sig_norm_counts)))
colnames(sig_norm_counts_iqn) <- colnames(sig_norm_counts)
row.names(sig_norm_counts_iqn) <- row.names(sig_norm_counts)
sig_norm_counts_iqn <- sig_norm_counts_iqn[row.names(sig_norm_counts_iqn) %in% 
                              gene_averages[which(abs(gene_averages$Log2FC_Inpatient.VOC)>1),"Gene"],]

######## Plot the heat map ##########
library(pheatmap)
# Set up column annotations
anno <- data.frame(pheno[,c('steadyState.vs.voc', "Chronic.Pain.")])
rownames(anno) <- row.names(pheno)

pheatmap(sig_norm_counts_iqn, clustering_method = "ward.D2", fontsize = 7, annotation_col = anno, 
         show_rownames = F, show_colnames = F,
)
