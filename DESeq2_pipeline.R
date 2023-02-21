library(DESeq2)

pheno$steadyState.vs.voc <- sub(' ', '_', pheno$steadyState.vs.voc)
########### Prep gene count matrix ##########################
gct <- read.table(gzfile("../../pharmhu/aggregate_files/pharmhu_topmed_to5_rnaseq_gene_reads.gct.gz"), 
                  header = T, check.names = F, skip = 2)
#gct2 <- gct
#gct <- gct2


# Remove rows where gene name is a duplicate (why are they present?) and remove gene name column
gct <- gct[!duplicated(gct$Description),]
rownames(gct) <- gct$Description
gct <- gct[,-c(1,2)]

# Remove ERCC transcripts
gct <- gct[!startsWith(row.names(gct), "ERCC-"), ]

# Remove sample not found in metadata and order the count matrix with respect to the metadata
cols_to_keep <- colnames(gct) %in% rownames(pheno)
gct <- gct[, cols_to_keep]
gct <- gct[, order(match(colnames(gct), rownames(pheno)))]

# Remove genes that are not expressed in at least 10% of samples
threshold <- 0.05
gct <- gct[which(rowSums(gct > 0)/ncol(gct) >= threshold), ]

# Remove genes where counts are zero for all samples
#gct <- gct[which(rowSums(gct) != 0),]

# Remove genes where count average is <= 
gct <- gct[which(rowSums(gct)/ncol(gct) > 10), ]
#gct <- gct[which(rowSums(gct)/ncol(gct) > 1), ]


### Perfrom DEG analysis with DESeq2
dds <- DESeqDataSetFromMatrix(countData = gct,
                              colData = pheno,
                              design= ~ CD71_libprep_batch + steadyState.vs.voc + sv1)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="steadyState.vs.voc_Steady_State_vs_Inpatient_VOC")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="steadyState.vs.voc_Steady.State_vs_Inpatient.VOC", type="apeglm")

test <- data.frame("Gene" = row.names(res),
                   log2FC = res$log2FoldChange,
                   pval_adj = res$padj)

sig <- test[which(sqrt((test$log2FC)^2) > 1 & test$pval_adj < 0.05),]
#sig <- test[which(test$pval_adj < 0.05),]

for(i in 1:nrow(sig)){
  cat(sig$Gene[i], '\n')
}

### Volcano plot
library(ggplot2)

ggplot(test, aes(x = log2FC, y = -log10(pval_adj))) +
  geom_point(size = 1, color = ifelse(test$log2FC > 0 & test$pval_adj <= 0.05, "red", 
                                      ifelse(test$log2FC <= 0 & test$pval_adj <= 0.05, "blue", "black"))) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("VOC vs Steady State") +
  theme(plot.title = element_text(hjust = 0.5))


library(pheatmap)

# Log2(1+X) transform the count matrix
gct <- apply(gct, 2, function(x) log2(x+1))

# Z-score normalize the gene count matrix
# Manually scaling instead of using pheatmaps scaling feature results in a better image because
# pheatmap will make the color scale too large to easily discriminate differences
gct_zscore <- t(scale(t(gct)))

# Load the list of differentially expressed genes
differentially_expressed_genes <- sig$Gene
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
