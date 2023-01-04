
library(dplyr)

# Subset down to samples with NWGC IDs associated with CD71 TORIDs
test <- results[!is.na(results$`CD71`),]

# Get list of NWGC IDs
sample_ids <- unique(test$CD71)



# Subset data to have only 1 row per NWGC ID
# Since most NWGC have multiple rows (observation dates), take the row with the most measurements.
# Should take the row that has the date that the RNA was collected #############################!!!
singlet__result <- data.frame()
for (id in sample_ids) {
  # Get all rows associated with a given NWGC ID
  rows <- test[test$CD71 == id, ]
  # find the row with the most non-NA values
  best_row <- which.max(rowSums(!is.na(rows)))
  # append the row to the result data frame
  singlet__result <- rbind(singlet__result, rows[best_row,])
}

#Remove unneeded columns e.g. not features, batch effect, sample ID for relevent cell type
nonfeature_col <- c("CD45", "CD45_Batch", "TORID_1_Date_of_collection", "TORID_2_Date_of_collection", 
                    "TORID_3_Date_of_collection", "TOR ID_1", "TOR ID_2", "TOR ID_3", "Sample ID_1", "Sample ID_2", "Sample_3", "RESULT_DATE")
singlet__result <- select(singlet__result, -nonfeature_col)

# Identify columns with all missing values & drop them
na_columns <- which(colSums(is.na(singlet__result)) == nrow(singlet__result))
singlet__result <- select(singlet__result, -na_columns)

# Get non-NA counts from each column and place in dataframe
nonna_counts <- sapply(singlet__result[,1:length(colnames(singlet__result))], function(x) sum(!is.na(x)))
nonna_df <- data.frame(Measurement = names(nonna_counts), "non NA" = nonna_counts)
nonna_df <- nonna_df[rev(order(nonna_df$non.NA)), ]

# Get the features that have the most observations
relevent_columns <- head(row.names(nonna_df), 19)
foo <- singlet__result[,relevent_columns]

#foo <- rev(foo)

# Write subsetted metadata to file
write.table(foo, file="CD71_metadata_subset2.txt", row.names = F, sep = "\t", quote = F)

# Write CD71 NWGC IDs to file
cd71_nwgc_id <- foo$CD71
#write.table(cd71_nwgc_id, file="cd71_nwgc_id.txt", row.names = F, col.names = F)

############ Edit the count matrix
# Create count matrix with only genes and NWGC IDs
# zcat ../../pharmhu/aggregate_files/pharmhu_topmed_to5_rnaseq_gene_tpm.gct | tail -n +3 | cut -f2- > cd71_tpm_counts.gct #Edit this; the output doesn't have only CD71
gct <- read.table("cd71_tpm_counts.gct", header = T, check.names = F)

# Subset count matrix to CD71 NWGC IDs
# The following line must have had "Description" manually added to the beginning when I made it
ids <- readLines("../../pharmhu/aggregate_files/test/cd71_nwgc_id.txt" )
gct <- gct %>% select(ids)

# Replace values in every column except the first with the log2(x+1) transformation of the TPM
gct[,-1] <- apply(gct[,-1], 2, function(x) log2(x+1))
write.table(gct, file="CD71_tpm_counts_log2_subset.gct", row.names = F, sep = "\t", quote = F)


