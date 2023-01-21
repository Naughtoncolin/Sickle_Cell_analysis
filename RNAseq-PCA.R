# Load the necessary libraries
library(ggplot2)
library(prcomp)
library(stats)
#gct2<-gct
#gct <- gct2

gct_pca_matrix <- t(gct)
gct_pca_matrix <- t(norm_counts$norm.dat)

# Perform PCA on the RNA-seq count matrix
pca_results <- prcomp(gct_pca_matrix)

# Extract the first two principal components
pc1 <- pca_results$x[,1]
pc2 <- pca_results$x[,2]

# Create a data frame with the principal component scores and the sample categories
pc_data <- data.frame(pc1, pc2, baseline.vs.voc = pheno$CD71_libprep_batch)
pc_data <- data.frame(pc1, pc2, baseline.vs.voc = pheno$baseline.vs.voc)
pc_data <- data.frame(pc1, pc2, baseline.vs.voc = pheno$Timepoint..at.RNA.collection.)
pc_data <- data.frame(pc1, pc2, baseline.vs.voc = pheno$Sex)
pc_data <- data.frame(pc1, pc2, baseline.vs.voc = pheno$Chronic.Pain.)

# Plot the first two principal components with different colors for each category
ggplot(pc_data, aes(x = pc1, y = pc2, color = baseline.vs.voc)) + 
  geom_point() + 
  xlab("PC1") + 
  ylab("PC2")

