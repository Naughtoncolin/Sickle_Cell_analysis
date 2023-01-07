setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")

master <- readxl::read_xlsx("RNA Masterlist with names.xlsx", sheet = "CD71+")

pain <- readxl::read_xlsx("../pain-omics-phenotype/Pain omics Phenotype_221117.xlsx")

###################### Add blood collection and MRN info to "Pain omics Phenotype" ##########################3
pain$Blood_collection_date <- NA
pain$Blood_processing_date <- NA
pain$RNA_extraction_date <- NA
pain$MRN <- NA

# Look through file with date info
for (excell_sheet in c("CD45+", "CD71+")){
  master <- readxl::read_xlsx("RNA Masterlist with names.xlsx", sheet = excell_sheet)
  
  # Loop through the rows and TORIDs within "Pain omics phenotype"
  for (i in 1:nrow(pain)) {
    for (col in c("TOPMed.RNA.ID..CD45..", "TOPMed.RNA.ID..CD71..")) {
      if (!is.na(pain[i, col])) {
        # Find the row in the results data frame where the TOR.ID_1, TOR.ID_2, or TOR.ID_3 column matches the current entry
        match_row <- match(pain[i, col], master$`TOR ID`)
      
      # Check if a matching row was found
      if (!is.na(match_row)) {
        # Copy the PAT_ID and MRN entries from the matching row in the results data frame to the corresponding columns in the pain data frame
        if(is.na(pain$Blood_collection_date[i])){pain$Blood_collection_date[i] <- as.character(master$`Date of collection`[match_row])}
        if(is.na(pain$Blood_processing_date[i])){pain$Blood_processing_date[i] <- as.character(master$`Date of processing`[match_row])}
        if(is.na(pain$RNA_extraction_date[i])){pain$RNA_extraction_date[i] <- as.character(master$`date of extraction`[match_row])}
        if(is.na(pain$MRN[i])){pain$MRN[i] <- master$MRN[match_row]}
        if(is.na(pain$Timepoint..at.RNA.collection.[i])){pain$MRN[i] <- master$Group[match_row]} #Didn't yield
        if(is.na(pain$Therapy..at.RNA.collection.[i])){pain$MRN[i] <- master$Treatment[match_row]} #Didn't yield
        }
      }
    }
  }
}

# Replace character NA with NULL NA
pain$MRN <- sub("N/A", NA, pain$MRN)

# Reorder the columns
pain <- pain %>%
  select("TOPMed.RNA.ID..CD45..", "TOPMed.RNA.ID..CD71..", "MRN", "Blood_collection_date", 
         "Blood_processing_date","RNA_extraction_date", "Chronic.Pain.", "Timepoint..at.RNA.collection.", 
         "Therapy..at.RNA.collection.", everything())

##################### Add NWGC IDs, batch, and notes from sequencing center####################################################
pain$CD45_NWGC_ID <- NA
pain$CD71_NWGC_ID <- NA
pain$CD45_libprep_batch <- NA
pain$CD71_libprep_batch <- NA
pain$CD45_notes <- NA
pain$CD71_notes <- NA

torids <- readxl::read_xlsx("lookup_pharmhu_topmed_to5_rnaseq_1_final.xlsx")
torids <- torids[,c("Notes", "NWGC Sample ID" , "TORID", "Pick Plate ID(s) - Indicates which samples were prepped together")]

for (i in 1:nrow(torids)){
  matching_row <- which(pain$TOPMed.RNA.ID..CD45.. %in% torids[i,"TORID"])
  if(length(matching_row) == 1){
    pain[matching_row,"CD45_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    pain[matching_row,"CD45_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    pain[matching_row,"CD45_notes"] <- torids[i,"Notes"]
    matching_row <- NA
  }
  matching_row <- which(pain$TOPMed.RNA.ID..CD71.. %in% torids[i,"TORID"])
  if(length(matching_row) == 1){
    pain[matching_row,"CD71_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    pain[matching_row,"CD71_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    pain[matching_row,"CD71_notes"] <- torids[i,"Notes"]
    matching_row <- NA
  }
  matching_row <- NA
}

###### Add Missing MRN & DNA IDs ########################
ids <- readxl::read_xlsx("Matched IDs_deident.xlsx")

#pain2 <- pain
pain <- pain2

for (i in 1:nrow(pain)) {
  # Check if the value in the "Investigator.RNA.Plasma.ID" column in "pain"
  # can be found in any of the "Sample" columns in "ids"
  match_result <- match(pain$`Investigator.RNA.Plasma.ID`[i], unlist(ids[, grep("Sample", colnames(ids))]))
  if (!is.na(match_result)) {
    # Check if MRN or DNA ID exists, if not then fill it in
    if (is.na(pain$MRN[i])) {
      pain$MRN[i] <- ids$MRN[match_result]
    }
    if(is.na(pain$Investigator.DNA.ID[i]))
      pain$Investigator.DNA.ID[i] <- ids$`DNA ID`[match_result]
  }
}


#################### Write xslx file ###################################################

# Reorder the columns
pain <- pain %>%
  select("CD45_notes", "TOPMed.RNA.ID..CD45..","CD45_NWGC_ID", "CD71_notes", "TOPMed.RNA.ID..CD71..", "CD71_NWGC_ID", "CD45_libprep_batch", 
         "CD71_libprep_batch", "MRN", "Blood_collection_date", "Blood_processing_date","RNA_extraction_date", 
         "Chronic.Pain.", "Timepoint..at.RNA.collection.", "Therapy..at.RNA.collection.", everything())


all.TORID <- pain[which(!is.na(pain$TOPMed.RNA.ID..CD45..)|!is.na(pain$TOPMed.RNA.ID..CD71..)),] #Get rows that have TORID
mrn.TORID <- all.TORID[which(!is.na(all.TORID$MRN)),]


library(xlsx)
write.xlsx(pain, file = "../pain-omics-phenotype/Pain omics Phenotype_230106.xlsx", sheetName = "All",
           col.names = TRUE, row.names = FALSE, append = FALSE, showNA=FALSE)
write.xlsx(all.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230106.xlsx", sheetName = "All TORIDs",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
