# Load the necessary libraries
library(ggplot2)
library(prcomp)
library(stats)
#gct2<-gct
#gct <- gct2

gct_pca_matrix <- t(gct)
gct_pca_matrix2 <- t(norm_counts$norm.dat)

# Perform PCA on the RNA-seq count matrix
pca_results <- prcomp(gct_pca_matrix)
pca_results2 <- prcomp(gct_pca_matrix2)


# Extract the first two principal components
pc1.1 <- pca_results$x[,1]
pc2.1 <- pca_results$x[,2]
pc1.2 <- pca_results2$x[,1]
pc2.2 <- pca_results2$x[,2]

# Create a data frame with the principal component scores and the sample categories
# Include basic timepoint
pc_data <- data.frame(pc1.1, pc2.1, Batch = pheno$CD71_libprep_batch, baseline.vs.voc = pheno$baseline.vs.voc,
                       Sex = pheno$Sex, Chronic.Pain. = pheno$Chronic.Pain.)
pc_data2 <- data.frame(pc1.2, pc2.2, Batch = pheno$CD71_libprep_batch, baseline.vs.voc = pheno$baseline.vs.voc,
                      Sex = pheno$Sex, Chronic.Pain. = pheno$Chronic.Pain.)

# Plot the first two principal components with different colors for each category
library(gridExtra)
for(i in 3:ncol(pc_data)){
  p1 <- ggplot(pc_data, aes(x= pc1.1, y = pc2.1, color = pc_data[,i])) + 
    geom_point() + 
    xlab("PC1") + 
    ylab("PC2") +
    ggtitle("Before supervised normalization") +
    scale_color_discrete(name=colnames(pc_data)[i])
  
  p2 <- ggplot(pc_data2, aes(x = pc1.2, y = pc2.2, color = pc_data2[,i])) + 
    geom_point() + 
    xlab("PC1") + 
    ylab("PC2") +
    ggtitle("After supervised normalization") +
    scale_color_discrete(name=colnames(pc_data2)[i])
  
  # Plot the first two principal components with different colors for each category
  grid.arrange(p1,p2, ncol=2)
}
p1 <- ggplot(pc_data, aes(x= pc1.1, y = pc2.1, color = baseline.vs.voc)) + 
  geom_point() + 
  xlab("PC1") + 
  ylab("PC2") +
  ggtitle("Before supervised normalization")

p2 <- ggplot(pc_data2, aes(x = pc1.2, y = pc2.2, color = baseline.vs.voc)) + 
  geom_point() + 
  xlab("PC1") + 
  ylab("PC2") +
  ggtitle("After supervised normalization")

# Plot the first two principal components with different colors for each category
grid.arrange(p1,p2, ncol=2)
