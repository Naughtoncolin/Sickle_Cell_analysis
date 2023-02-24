setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")

library(dplyr)
library(sva)
library(snm)
library(pvca)
library(Biobase)
library(ExpressionNormalizationWorkflow)

########### Prep metadata ####################################
pheno <- readxl::read_xlsx("../pain-omics-phenotype/Pain omics Phenotype_230122.xlsx", 
                           sheet = "TORIDs Passing QC with measurem")
pheno <- as.data.frame(pheno)

# Subset down to samples with NWGC IDs associated with CD71 TORIDs
pheno <- pheno[!is.na(pheno$CD71_NWGC_ID),]

#Select columns for downstream analysis
feature_col <- c("CD71_NWGC_ID", "Subject_ID", "CD71_libprep_batch", "Sex", "Timepoint..at.RNA.collection.", "Chronic.Pain.")
pheno <- pheno %>% select(any_of(feature_col))

#Drop rows with NA
# IS THIS NECESSARY?
pheno <- na.omit(pheno)

# Make metadata column grouping points\
base <- c("Baseline", "HU", "HU MTD", "Metformin Baseline", "Metformin", "Metformin MTD")
pheno$steadyState.vs.voc <- ifelse(pheno$Timepoint..at.RNA.collection. %in% base, "Steady State", pheno$Timepoint..at.RNA.collection.)
pheno <- pheno[which(pheno$steadyState.vs.voc!="Inpatient F/U"),] # Get id of VOC F/U


#Subset to unique individuals with baseline & VOC
# There are still duplicates accross the two categories
pheno <- pheno %>% 
  select(CD71_NWGC_ID, Subject_ID, steadyState.vs.voc, Timepoint..at.RNA.collection., CD71_libprep_batch, Sex, Chronic.Pain.) %>%
  group_by(steadyState.vs.voc) %>% 
  distinct(Subject_ID, .keep_all = T) %>%
  ungroup() %>%
  as.data.frame()
row.names(pheno) <- pheno$CD71_NWGC_ID
pheno$CD71_NWGC_ID <- NULL

#Change metadata to numeric for numeric features
#pheno$HGB <- as.numeric(sub(" .*", "", pheno$HGB))

# Select 10 random rows for each "Baseline" and "Inpatient VOC" group
#pheno <- pheno %>%
# select(CD71_NWGC_ID, Subject_ID, baseline.vs.voc, CD71_libprep_batch, Sex) %>%
#   group_by(baseline.vs.voc) %>% 
#   distinct(Subject_ID, .keep_all = T) %>%
#   sample_n(10) %>%
#   ungroup() %>%
#   as.data.frame()

# Change metadata character columns to factors
# This may need to be edited once more metadata features are used
pheno[] <- lapply(pheno, function(x) {
  if (is.character(x)) {
    as.factor(x)
  } else {
    x
  }
})
pheno$Subject_ID <- as.factor(pheno$Subject_ID)


########### Prep gene count matrix ##########################
#gct <- read.table("cd71_tpm_counts.gct", header = T, check.names = F)
gct <- read.table(gzfile("../../pharmhu/aggregate_files/pharmhu_topmed_to5_rnaseq_gene_tpm.gct.gz"), 
                  header = T, check.names = F, skip = 2)
#gct2 <- gct
#gct <- gct2


# Remove rows where gene name is a duplicate (why are they present?) and remove gene name column
gct <- gct[!duplicated(gct$Description),]
rownames(gct) <- gct$Description
gct <- gct[,-c(1,2)]
#gct <- gct[,-1] #Keeping for testing purposes so no need to unzip count matrix


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

# Remove genes where count average is <= 1 TPM
gct <- gct[which(rowSums(gct)/ncol(gct) > .1), ]
#gct <- gct[which(rowSums(gct)/ncol(gct) > 1), ]


# Log2(1+X) transform the count matrix
gct <- apply(gct, 2, function(x) log2(x+1))

##################### PVCA with original covariates #################
covariates <- c("CD71_libprep_batch", "Sex", "steadyState.vs.voc", "Chronic.Pain.")
pct_threshold <- 0.75
pd <- new("AnnotatedDataFrame", data = pheno)
inpData <- ExpressionSet(assayData = gct, phenoData = pd)

# Make PVCA plot
pvcaPlot1 <- pvcAnaly(inpData, pct_threshold, covariates) #Need an "Expression Set" Object as input

##################### SVA Analysis ########################################

### Make formulas for modeling matrix in next step. Useful for when using large numbers of features
###########!!!!!!!!!!!!!!Issue occuring with the "%" sign in column names
#regform <- as.formula(paste("~ ", paste(names(pheno),collapse="+")))
#regform2 <- as.formula(paste("~ ", paste(names(pheno)[-6],collapse="+")))
#Example output: ~CD71 + CD71_Batch + SEX + RACE + ETHNICITY + HGB

### Model matrices
mod = model.matrix(~ CD71_libprep_batch + Sex + Chronic.Pain. + steadyState.vs.voc , data=pheno)
mod0 = model.matrix(~ CD71_libprep_batch + Sex + Chronic.Pain.,data=pheno)

### Determine number of surrogate variables to calculate manually
n.sv <- num.sv(gct,mod,method="leek")
#num.sv(gct,mod,method="be")

### Run SVA; automatically determine number of surrogate variables to calculate
sva.out <- sva(gct,mod,mod0, n.sv = n.sv)

##############  Determine relationship of surrogate variables to metadata ##########
# Add surrogate variables to metadata
pheno$sv1 <- sva.out$sv[,1]
#pheno$sv2 <- sva.out$sv[,2]
#pheno$sv3 <- sva.out$sv[,2]
## Fit a generalized linear model for sv1
glm.sv1 <- glm(sv1 ~ CD71_libprep_batch + Sex + steadyState.vs.voc, data = pheno) 
summary(glm.sv1)
glm.sv2 <- glm(sv2 ~ CD71_libprep_batch + Sex + steadyState.vs.voc, data = pheno) 
summary(glm.sv2)
glm.sv3 <- glm(sv3 ~ CD71_libprep_batch + Sex + steadyState.vs.voc, data = pheno) 
summary(glm.sv3)
coef(summary(glm.sv1))[,4]

######### PVCA including surrogate variables
#pheno2 <- pheno
#pheno <- pheno2

# conTocat requires two things in the list.
# Tried to use it on inpData, but wouldn't allow duplicate columns.
# Instead used on pheno, which was then used to make inpData variable
#var_names <- c("sv1", "sv2", "sv3")
var_names <- c("sv1", "sv1")
pheno<-conTocat(pheno, var_names)

#Only use when 1 SV
pheno <- pheno[,-ncol(pheno)]
pd2 <- new("AnnotatedDataFrame", data = pheno)
inpData <- ExpressionSet(assayData = gct, phenoData = pd2)

# Make PVCA plot including surrogate variables
pvcaPlot2 <- pvcAnaly(inpData, pct_threshold, covariates)

# Make PVCA plot but include interactions between covariates
gibPlot2 <- pvcAnaly(inpData, pct_threshold, covariates) # Residual: 0.752 > 0.739

# Graph variance explained
# bp <- barplot(pvcaObj$dat,
#               ylab = "Variance explained",
#               ylim= c(0,1.1),col = c("blue"), las=2,
#               main="PVCA estimation bar chart")
# axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
# values = pvcaObj$dat
# new_values = round(values , 3)
# text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)

#pheno <- pheno
# pheno$sv1 <- sva.out$sv[,1]
# inpData <- expSetobj(gct, pheno)
# var_names <- c("sv1") 
# pData(inpData)<-conTocat(pData(inpData), var_names) 
# cvrts_eff_var <- c("Subject_ID","CD71_libprep_batch","Sex","baseline.vs.voc", "sv1")
# pvcAnaly(inpData, pct_threshold, cvrts_eff_var)


############## SNM Analysis #############################
# Extract the surrogate variables from the 'sv' object
sv_vars = sva.out$sv

# Use the 'snm()' function to normalize the count matrix using the surrogate variables
#norm_counts = snm(gct, bio.var=pheno ,adj.var=sv_vars)

# Test
# Create model matrices for the biological variable and adjustment variable
#bio.var.mod = model.matrix(~ CD71_libprep_batch , data=pheno)
#adj.var.mod = model.matrix(~sv_vars)

#bio.var.mod = model.matrix(~ Sex + steadyState.vs.voc + Chronic.Pain., data=pheno)
#adj.var.mod = model.matrix(~ CD71_libprep_batch + sv_vars, data = pheno)

bio.var.mod = model.matrix(~  steadyState.vs.voc , data=pheno)
adj.var.mod = model.matrix(~ Sex + CD71_libprep_batch + sv1 + Chronic.Pain., data = pheno)

# Normalize the count matrix using the biological variable and adjustment variable model matrices
# NOTE: This leave negative values in the matrix!
norm_counts = snm(gct, bio.var=bio.var.mod, adj.var=adj.var.mod, rm.adj=T)
head(summary(norm_counts))
#ks.test(norm_counts$pval, "punif")

######### Enrichment Analysis ############################
library(qvalue)


# Calculate q-values and order from smallest to largest
q <- qvalue(norm_counts$pval)
deg <- data.frame("gene" = row.names(gct), "qval" = q$qvalues, "pval" = q$pvalues)
deg <- deg[order(deg$qval),]

number_of_DEG <- 3778

### Print top DEGs
# Below bonferroni threshold
# bonferroni_cutoff <- .05/nrow(gct)
# counter <- 0 
# for (i in 1:nrow(deg)) {
#   if(deg$pval[i]<bonferroni_cutoff){
#     counter <- counter +1
#     cat(deg$gene[i], "\n")
#   }
# }

# Based on q-value
# for (i in 1:number_of_DEG) {
#   cat(deg$gene[i], "\n")
# }
top_found_genes <- deg[1:number_of_DEG,"gene"]
write.table(top_found_genes, file = "found_genes_12k_1e-04.txt", row.names = F, col.names = F, quote = F)
  