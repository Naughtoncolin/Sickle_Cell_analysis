setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")

library(sva)

########### Prep gene count matrix ##########################
gct <- read.table("CD71_tpm_counts_log2_subset.gct", header = T, check.names = F )

# Remove rows where gene name is a duplicate (Is this necessary?)
gct <- gct[which(!duplicated(gct$Description)),]
#non_unique_values <- gct$Description[duplicated(gct$Description)]

# Remove gene column from count matrix
gct <- gct[,-1]

# Remove genes where counts are zero for all samples
gct <- gct[which(rowSums(gct) != 0),]

# SVA doesn't allow data frames, must convert to matrix
gct <- as.matrix(gct)

########### Prep metadata ####################################
pheno <- read.table("CD71_metadata_subset2.txt", header = T, check.names = F, sep = "\t")

# Subset metadata down to desired features
features <- c("CD71", "CD71_Batch", "SEX", "RACE", "ETHNICITY", "HGB")
pheno <- select(pheno, features)

#Change metadata to numeric for numeric features
pheno$HGB <- as.numeric(sub(" .*", "", pheno$HGB))
#pheno$`LYMPH%` <- as.numeric(sub(" .*", "", pheno$`LYMPH%`))

# Change metadata character columns to factors
pheno[] <- lapply(pheno, function(x) {
  if (is.character(x)) {
    as.factor(x)
  } else {
    x
  }
})
pheno$CD71 <- as.factor(pheno$CD71)
pheno$CD71_Batch <- as.factor(pheno$CD71_Batch)

##################### Run SVA ########################################

# Make formulas for modeling matrix in next step. Useful for when using large numbers of features
###########!!!!!!!!!!!!!!Issue occuring with the "%" sign in column names
regform <- as.formula(paste("~ ", paste(names(pheno),collapse="+")))
regform2 <- as.formula(paste("~ ", paste(names(pheno)[-6],collapse="+")))
# ~CD71 + CD71_Batch + SEX + RACE + ETHNICITY + HGB

# Model matrices
# NWGC ID and RACE were causing issues with sva due to error: "Lapack routine dgesv: system is exactly singular: U[197,197] = 0"
mod = model.matrix(~ CD71 +CD71_Batch + SEX + ETHNICITY + HGB, data=pheno)
mod0 = model.matrix(~ CD71 +CD71_Batch + SEX  + ETHNICITY  ,data=pheno)

#mod = model.matrix(~ CD71_Batch +SEX+ ETHNICITY + HGB, data=pheno)
#mod0 = model.matrix(~ CD71_Batch +SEX+ETHNICITY ,data=pheno)

# Determine number of surrogate variables to calculate manually
#num.sv(gct,mod,method="leek")

# Run SVA; automatically determine number of surrogate variables to calculate
sva.out <- sva(gct,mod,mod0)
