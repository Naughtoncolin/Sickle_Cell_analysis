setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")

library(dplyr)
library(tidyr)
library(purrr)

# Read in file of lab measurements
metadata <- readxl::read_xlsx("deident_Copy of KM_SC_GENOTYPE_LAB_VALUES - from past to 020618 on 020718.xlsx")

# For testing
#metadata2 <- metadata
#metadata <- metadata2

# Identify columns with all missing values & drop them
na_columns <- which(colSums(is.na(metadata)) == nrow(metadata))
metadata <- select(metadata, -na_columns)

# Replace character "NA" values in the TOR ID_1 column with the NULL NA value in R
metadata$`TOR ID_1` <- ifelse(metadata$`TOR ID_1` == "NA", NA, metadata$`TOR ID_1`)

# Identify rows with an NA value in columns 2, 3, and 4 (TORID Columns) & drop them
na_rows <- which(is.na(metadata[,2]) & is.na(metadata[,3]) & is.na(metadata[,4]))
metadata <- filter(metadata, !(row_number() %in% na_rows))

# Identify rows with an NA value for whatever value was trying to be measured & drop them.
na_values <- which(is.na(metadata$ORD_VALUE))
metadata <- filter(metadata, !(row_number() %in% na_values))

#metadata <- unique(metadata[,c("TOR ID_1", "COMPONENT_NAME", "ORD_VALUE", "ORDERING_DATE")])
#Remove "CULTURE SOURCE" measurement due to repeat values.
metadata <- metadata[which(metadata$COMPONENT_NAME!="CULTURE SOURCE"),]

# Make new column with measured value with its respective units.
metadata$COMPONENT_VALUE <- paste(metadata$ORD_VALUE, metadata$REFERENCE_UNIT, ' ')

# Remove unneccessary columns interfering with spread.
unneeded_col <- c("PAT_ID", "MRN", "PAT_ENC_CSN_ID", "PROC_CODE", "PROC_NAME", "ORDERING_DATE", "ORD_VALUE", "ORD_NUM_VALUE", "REFERENCE_UNIT", "COMPONENT_ID", "COMPONENT_ABBR")
metadata <- select(metadata, -unneeded_col)

# Remove duplicate rows
metadata <- unique(metadata)

# Remove rows that have "CANCELLED" within the measured value
metadata <- filter(metadata, !grepl("CANCELLED", COMPONENT_VALUE))

# Find the rows with duplicated values in the specified columns implicated in spread error, for example if a patient on the same date had 2 different blood type results
metadata <- metadata[order(metadata$`TOR ID_1`, metadata$RESULT_DATE, metadata$COMPONENT_NAME),]


# Trying to use the spread function threw errors for certain rows, for example if a patient on the same date had 2 different blood type results
# Find the rows with duplicated values in the specified columns implicated in spread error, for example if a patient on the same date had 2 different blood type results
# Copied rows from output to txt file and ran the bash command: cat rows3.txt | sed 's/\*//' | tr '\n' ',' | sed 's/, /\n/g;s/,/\n/' > rows3-long.txt
duplicate.rows <- readLines("rows3-long.txt")
duplicate.rows <- as.numeric(duplicate.rows)

# Write excel file with duplicates
#duplicates <- metadata[duplicate.rows,]
#write.xlsx(duplicates, file="duplicates.xlsx")


# Remove duplicate rows
metadata <- metadata[-duplicate.rows,]

# Transform long data to wide data: Make each measurement type a column with measurements as values in the column e.g. CREATININE & 0.59 MG/DL
results <- metadata %>% 
  spread("COMPONENT_NAME", "COMPONENT_VALUE")
#write.csv(results, file="TOR_lab-values.csv", row.names = F )

############## Summarize number of measurements for each feature and place in bins ##################################
# Get non-NA counts from  feature column
nonna_counts <- sapply(results[,12:length(colnames(results))], function(x) sum(!is.na(x)))
nonna_df <- data.frame(Measurement = names(nonna_counts), "non NA" = nonna_counts)
nonna_df <- nonna_df[rev(order(nonna_df$non.NA)), ]
#xlsx::write.xlsx(nonna_df, file="TOR_lab-value-counts.xlsx", row.names = F)

# Define the breaks for the bins
breaks <- c(0, 500, 1000, 5000, 10000)

# Use the cut function to bin the counts
bins <- cut(nonna_counts, breaks = breaks)

# Use the split function to split the counts into a list of vectors
# with each vector containing the counts that fall within a specific bin
binned_counts <- split(nonna_counts, bins)

# Create a data frame for each bin
binned_dfs <- lapply(binned_counts, function(x) data.frame(Measurement = names(x), "non NA" = x))

# Write the data frames to the Excel file using write.xlsx()
#xlsx::write.xlsx(binned_dfs[[1]], file="binned_lab-measurement_counts.xlsx", sheetName = "0-500", row.names = F)
#xlsx::write.xlsx(binned_dfs[[2]], file="binned_lab-measurement_counts.xlsx", sheetName = "501-1000", row.names = F, append = T)
#xlsx::write.xlsx(binned_dfs[[3]], file="binned_lab-measurement_counts.xlsx", sheetName = "1001-5000", row.names = F, append = T)
#xlsx::write.xlsx(binned_dfs[[4]], file="binned_lab-measurement_counts.xlsx", sheetName = "5001-10000", row.names = F, append = T)



################# Associate NWGC IDs with TORIDs in metadata ###########################
torids <- readxl::read_xlsx("lookup_pharmhu_topmed_to5_rnaseq_1_final.xlsx")
#torids <- torids[which(is.na(torids$Notes)),"TORID"] # Anything that had info in "Notes" was a bad sample.
torids <- torids[which(is.na(torids$Notes)),c("NWGC Sample ID" , "TORID", "Tissue Type")] # Anything that had info in "Notes" was a bad sample.

# See if TORIDs from lab measurements are found in TORIDs with RNA-seq data
matching_rows <- which(results$`TOR ID_1` %in% torids$TORID | results$`TOR ID_2` %in% torids$TORID | results$`TOR ID_3` %in% torids$TORID)
results <- results[matching_rows,]

# Determine which TORIDs in the lab measurements metadata aren't found to the pharmhu_topmed lookup sheet
#missing_torids <- results[-matching_rows,] 

# Determine which TORIDs in the pharmhu_topmed lookup table can be found in the lab measurement metadata.
present_TORID_rows <- which(torids$TORID %in% results$`TOR ID_1` | torids$TORID %in% results$`TOR ID_2` | torids$TORID %in% results$`TOR ID_3`)
torids <- torids[present_TORID_rows,]

# For testing code
#results2 <- results
#results <- results2

# Add new columns that will contain NWGC IDs
results$"CD71" <- NA
results$"CD45" <- NA

# Iterate through lab measurement metadata rows and add NWGC IDs
# This is probably where the issue is happening due to multiple CD45 or CD71 NWGS IDs
for (i in 1:nrow(results)) {
  
  #find rows in "torids" dataframe that match any of the values in the "TOR ID_1", "TOR ID_2", or "TOR ID_3" columns of the "results" dataframe
  matching_rows <- torids[torids$TORID %in% results[i, c("TOR ID_1", "TOR ID_2", "TOR ID_3")],]

  #loop through each matching row in "torids" dataframe
  for (j in 1:nrow(matching_rows)) {
    # check if "Tissue Type" column contains "CD71"
    if (grepl("CD71", matching_rows[j, "Tissue Type"])) {
      # fill in value from "NWGC Sample ID" column in "results" dataframe
      results[i, "CD71"] <- matching_rows[j, "NWGC Sample ID"]
    }
    # check if "Tissue Type" column contains "CD45"
    if (grepl("CD45", matching_rows[j, "Tissue Type"])) {
      # fill in value from "NWGC Sample ID" column in "results" dataframe
      results[i, "CD45"] <- matching_rows[j, "NWGC Sample ID"]
    }
  }
}

# Make the CD45+ & CD71+ NWGC IDs appear first.
results <- results %>%
  select("CD45", "CD71", everything())
#write.csv(results, file="TOR_lab-values.csv", row.names = F )

# Count number of unique CD45+ & CD71+ NWGC IDs
length(unique(results$`CD45`)) # 8
length(unique(results$`CD71`)) # 188

