setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")
library(dplyr)
library(tidyr)

master <- readxl::read_xlsx("Omics_Master list(2).xlsx", sheet = "Omics Master List_wide")

# Remove excess empty rows; 964 not empty upon visual inspection
master <- master[1:964,]

##################### Add NWGC IDs, batch, and notes from sequencing center####################################################
master$CD45_NWGC_ID <- NA
master$CD71_NWGC_ID <- NA
master$CD45_libprep_batch <- NA
master$CD71_libprep_batch <- NA
master$CD45_notes <- NA
master$CD71_notes <- NA
master$Investigator.Sex <- NA

torids <- readxl::read_xlsx("lookup_pharmhu_topmed_to5_rnaseq_1_final.xlsx")
# Some of the sexes were missing, but the sequencing center determined them, so
# create a unified list of sex; if the investigator sex is NA, then replace with
# the NWGC assay-determined sex.
torids$`Investigator Sex` <- sub("Unknown", NA, torids$`Investigator Sex`)
torids$`Assay Sex` <- sub("N/A", NA, torids$`Assay Sex`)
#torids$`Investigator Sex` <- ifelse(is.na(torids$`Investigator Sex`), torids$`Assay Sex`, torids$`Investigator Sex`)

torids[!is.na(torids$`Family ID`),] <- torids[!is.na(torids$`Family ID`),] %>%
  group_by(`Family ID`) %>%
  fill(c(`Investigator Sex`, `Assay Sex`), .direction = "updown") %>%
  ungroup()
torids$`Investigator Sex` <- ifelse(is.na(torids$`Investigator Sex`), torids$`Assay Sex`, torids$`Investigator Sex`)

#Include missing sex mismatch
torids[which(torids$`NWGC Sample ID` == "517120"),"Notes"] <- "Sex mismatch"

# 2294B, 2544B, 2562B, and 2711B all have sex mismatches, but have a complimentary ID which confirms the sex via the assay sex
# Thus I'm filling in all Investigator Sex IDs from the assay sex; the IDs that had a mismatchs sex that couldn't be rectified 
# Still have "Sex mismatch" in the notes. 2294B, 2544B, 2562B all came from the same MRN, which makes it weird that they all had
# mismatched sex.
torids$`Assay Sex` <- ifelse(is.na(torids$`Assay Sex`), torids$`Investigator Sex`, torids$`Assay Sex`)
# Now assay sex has the correct IDs, except for those with "Sex mismatch" in note, so the IDs can be transferred over.
torids$`Investigator Sex` <- torids$`Assay Sex`

# Check if any of the TORIDs can be found in the metadata provided by the NWGC and fill in info.
# TODO: This code code have 1 if statement if you regex out the "CD45" or "CD71" from the columns in "master".
for (i in 1:nrow(torids)){
  matching_row <- which(master$`TOR ID_CD45` %in% torids[i,"TORID"])
  if(length(matching_row)==1){
    master[matching_row,"CD45_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    master[matching_row,"CD45_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    master[matching_row,"CD45_notes"] <- torids[i,"Notes"]
    master[matching_row,"Investigator.Sex"] <- torids[i,"Investigator Sex"]
    matching_row <- NULL
  }
  matching_row <- NULL
  matching_row <- which(master$`TOR ID_CD71` %in% torids[i,"TORID"])
  if(length(matching_row)==1){
    master[matching_row,"CD71_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    master[matching_row,"CD71_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    master[matching_row,"CD71_notes"] <- torids[i,"Notes"]
    master[matching_row,"Investigator.Sex"] <- torids[i,"Investigator Sex"]
    matching_row <- NULL
  }
  matching_row <- NULL
}


###################### Add blood collection dates and MRN info to "Pain omics Phenotype" ##########################3
# Add columns for data that will be added later
master[c("CD45.Blood_collection_date","CD71.Blood_collection_date", "RNA_extraction_date", 
       "Birthday", "Race", "Ethnicity", "Sex")] <- NA

###Try with lab measurement MRN and Birthday
lab.measurements <- read.csv("TOR_lab-values.csv")
correct_mets <- read.csv("Met_corrections.csv", header = T)

# # Replace incorrect MET IDs in Sample.ID_1, Sample.ID_2, and Sample_3 in the lab measurements with correct versions
# lab.measurements$Sample.ID_1 <- ifelse(match(lab.measurements$Sample.ID_1, correct_mets$Wrong) > 0, correct_mets$Correct[match(lab.measurements$Sample.ID_1, correct_mets$Wrong)], lab.measurements$Sample.ID_1)
# lab.measurements$Sample.ID_2 <- ifelse(match(lab.measurements$Sample.ID_2, correct_mets$Wrong) > 0, correct_mets$Correct[match(lab.measurements$Sample.ID_2, correct_mets$Wrong)], lab.measurements$Sample.ID_2)
# lab.measurements$Sample_3 <- ifelse(match(lab.measurements$Sample_3, correct_mets$Wrong) > 0, correct_mets$Correct[match(lab.measurements$Sample_3, correct_mets$Wrong)], lab.measurements$Sample_3)

for(j in 1:nrow(lab.measurements)){
  #print(lab.measurements$Sample.ID_1[j] %in% correct_mets$Wrong)
  # lab.measurements$Sample.ID_1[j] <- ifelse(match(lab.measurements$Sample.ID_1[j], correct_mets$Wrong) > 0, 
  #                                           correct_mets$Correct[match(lab.measurements$Sample.ID_1[j], correct_mets$Wrong)], 
  #                                           lab.measurements$Sample.ID_1[j])
  # lab.measurements$Sample.ID_2[j] <- ifelse(match(lab.measurements$Sample.ID_2[j], correct_mets$Wrong) > 0, 
  #                                           correct_mets$Correct[match(lab.measurements$Sample.ID_2[j], correct_mets$Wrong)], 
  #                                           lab.measurements$Sample.ID_2[j])
  # lab.measurements$Sample_3[j] <- ifelse(match(lab.measurements$Sample_3[j], correct_mets$Wrong) > 0, 
  #                                        correct_mets$Correct[match(lab.measurements$Sample_3[j], correct_mets$Wrong)], 
  #                                        lab.measurements$Sample_3[j])
  if(lab.measurements$Sample.ID_1[j] %in% correct_mets$Wrong){
    lab.measurements$Sample.ID_1[j] <- correct_mets$Correct[match(lab.measurements$Sample.ID_1[j], correct_mets$Wrong)]
  }
  else if(lab.measurements$Sample.ID_1[j] %in% correct_mets$Correct){
    lab.measurements$Sample.ID_1[j] <- NA
  }
  if(lab.measurements$Sample.ID_2[j] %in% correct_mets$Wrong){
    lab.measurements$Sample.ID_2[j] <- correct_mets$Correct[match(lab.measurements$Sample.ID_2[j], correct_mets$Wrong)]
  }
  else if(lab.measurements$Sample.ID_2[j] %in% correct_mets$Correct){
    lab.measurements$Sample.ID_2[j] <- NA
  }
  if(lab.measurements$Sample_3[j] %in% correct_mets$Wrong){
    lab.measurements$Sample_3[j] <- correct_mets$Correct[match(lab.measurements$Sample_3[j], correct_mets$Wrong)]
  }
  else if(lab.measurements$Sample_3[j] %in% correct_mets$Correct){
    lab.measurements$Sample_3[j] <- NA
  }
  
  # Format the dates to "YYYY-MM-DD"
  if(!is.na(lab.measurements$RESULT_DATE[j])){
    date_string <- lab.measurements$RESULT_DATE[j]
    date_object <- as.Date(date_string, format = "%d-%b-%y")
    lab.measurements$RESULT_DATE[j] <- format(date_object, format = "%Y-%m-%d")
  }
  if(!is.na(lab.measurements$DOB[j])){
    date_string <- lab.measurements$DOB[j]
    date_object <- as.Date(date_string, format = "%d-%b-%y")
    lab.measurements$DOB[j] <- format(date_object, format = "%Y-%m-%d")
  }
}


#master2 <- master
#master <- master2
# Fill in missing MRN
set.seed(123)
master$MRN[is.na(master$MRN)] <- paste0(sample(1e6, sum(is.na(master$MRN)), replace=TRUE), "_fake")


for (i in master$MRN){
  if (i %in% lab.measurements$MRN){
    master[master$MRN == i, "Birthday"] <- lab.measurements[lab.measurements$MRN == i, "DOB"][1]
    master[master$MRN == i, "Race"] <- lab.measurements[lab.measurements$MRN == i, "RACE"][1]
    master[master$MRN == i, "Ethnicity"] <- lab.measurements[lab.measurements$MRN == i, "ETHNICITY"][1]
    master[master$MRN == i, "Sex"] <- lab.measurements[lab.measurements$MRN == i, "SEX"][1]
    #master[master$MRN == i, "Birthday"] <- lab.measurements[lab.measurements$MRN == i, "DOB"][1]
  }
}


for(i in 1:nrow(master)){
  if(!is.na(master$`TOR ID_CD45`[i])){
    for(j in 1:nrow(lab.measurements)){
      if(master$`TOR ID_CD45`[i] == lab.measurements$TOR.ID_1[j]){
        master$CD45.Blood_collection_date[i] <- lab.measurements$TORID_1_Date_of_collection[j]
      } else if(!is.na(lab.measurements$TOR.ID_2[j])){
        if(master$`TOR ID_CD45`[i] == lab.measurements$TOR.ID_2[j]){
          master$CD45.Blood_collection_date[i] <- lab.measurements$TORID_2_Date_of_collection[j]
        }
      } else if(!is.na(lab.measurements$TOR.ID_3[j])){
        if(master$`TOR ID_CD45`[i] == lab.measurements$TOR.ID_3[j]){
          master$CD45.Blood_collection_date[i] <- lab.measurements$TORID_3_Date_of_collection[j]
        }
      }
    }
  }
}
for(i in 1:nrow(master)){
  if(!is.na(master$`TOR ID_CD71`[i])){
    for(j in 1:nrow(lab.measurements)){
      if(master$`TOR ID_CD71`[i] == lab.measurements$TOR.ID_1[j]){
        master$CD71.Blood_collection_date[i] <- lab.measurements$TORID_1_Date_of_collection[j]
      } else if(!is.na(lab.measurements$TOR.ID_2[j])){
        if(master$`TOR ID_CD71`[i] == lab.measurements$TOR.ID_2[j]){
          master$CD71.Blood_collection_date[i] <- lab.measurements$TORID_2_Date_of_collection[j]
        }
      } else if(!is.na(lab.measurements$TOR.ID_3[j])){
        if(master$`TOR ID_CD71`[i] == lab.measurements$TOR.ID_3[j]){
          master$CD71.Blood_collection_date[i] <- lab.measurements$TORID_3_Date_of_collection[j]
        }
      }
    }
  }
}
master$Race <- sub("Unable to Obtain", NA, master$Race)

master[!is.na(master$`Unique ID`),] <- master[!is.na(master$`Unique ID`),] %>%
  group_by(`Unique ID`) %>%
  fill(c("Investigator.Sex", 
         "Birthday", "Race", "Ethnicity", "Sex"), .direction = "updown") %>%
  ungroup()

master <- master %>%
  select(`Unique ID`, "Investigator.Sex", "Sex", "Birthday", "DOB", "CD71.Blood_collection_date", 
         "CD45.Blood_collection_date", "Date of collection_CD71", "Date of collection_cd45", everything())
# NWGC ID 517065 doesn't have matching Sex after adding from the pharmhu and lab measurements file, but
# the lab measurements seem more trustworthy because they're more comprehensive. Searching for NWGC ID 517065
# after adding from both indicates the lab measurements file is correct.

#TOR472994 has two different MRN in lab measurements
# MET043 has two different MRN; it's also the ID with the paired sex mismatch
# MRN 3001485403 has different DOB; can have "Sex" transferred over?


master$Investigator.Sex <- sub("Male", "M", master$Investigator.Sex)
master$Investigator.Sex <- sub("Female", "F", master$Investigator.Sex)
master$Sex <- ifelse(is.na(master$Sex), master$Investigator.Sex, master$Sex)

# Find unmatching sex
# test_sex2 <- data.frame()
# for(i in 1:nrow(master)){
#   if(!is.na(master$Sex[i])){
#     if(master$Sex[i]!=master$Investigator.Sex[i]){
#       test_sex2 <- rbind(test_sex2, master[i,])
#     }
#   }
# }
# Unique IDs 140 and 338 had mismatched sexes, but had additional libraries to confirm values in "Sex" was correct
# as opposed to "Investigator.Sex". Thus, all the values in "Sex" shold be trusted and the "Investigator.Sex"
# column can be dropped.
master$Investigator.Sex <- NULL


#master2 <- master
#Find unmatching DOB
# master$DOB<- format(master$DOB, format = "%Y-%m-%d")
# test_dob2 <- data.frame()
# for(i in 1:nrow(master)){
#   if(!is.na(master$Birthday[i])){
#     if(master$DOB[i]!=master$Birthday[i]){
#       test_dob2 <- rbind(test_dob2, master[i,])
#     }
#     
#   }
# }
#All matching; can remove the "Birthday" column since it's redundant.
master$Birthday <- NULL

# test_cd45 <- data.frame()
# master$`Date of collection_cd45`  <- format(master$`Date of collection_cd45`, format = "%Y-%m-%d")
# for(i in 1:nrow(master)){
#   if(!is.na(master$CD45.Blood_collection_date[i])){
#     if(master$CD45.Blood_collection_date[i]!=master$`Date of collection_cd45`[i]){
#       test_cd45 <- rbind(test_cd45, master[i,])
#     }
#   }
# }
#All matching; can remove the "CD45.Blood_collection_date" column since it's redundant.
master$CD45.Blood_collection_date <- NULL

# test_cd71 <- data.frame()
# master$`Date of collection_CD71`  <- format(master$`Date of collection_CD71`, format = "%Y-%m-%d")
# for(i in 1:nrow(master)){
#   if(!is.na(master$CD71.Blood_collection_date[i])){
#     if(master$CD71.Blood_collection_date[i]!=master$`Date of collection_CD71`[i]){
#       test_cd71 <- rbind(test_cd71, master[i,])
#     }
#   }
# }
#All matching; can remove the "CD71.Blood_collection_date" column since it's redundant.
master$CD71.Blood_collection_date <- NULL

# Add timepoint, therapy at timepoint, and chronic pain status from pain omics phenotype
pain <- readxl::read_xlsx("../pain-omics-phenotype/Pain omics Phenotype_221030.xlsx")

master$Timepoint..at.RNA.collection. <- NA
master$Therapy..at.RNA.collection. <- NA
master$Chronic.Pain. <- NA

for (i in master$`TOR ID_CD71`){
  if(!is.na(i)){
    if (i %in% pain$TOPMed.RNA.ID..CD71..){
      master[which(master$`TOR ID_CD71` == i), "Timepoint..at.RNA.collection."] <- pain[which(pain$TOPMed.RNA.ID..CD71.. == i), "Timepoint..at.RNA.collection."][1]
      master[which(master$`TOR ID_CD71` == i), "Therapy..at.RNA.collection."] <- pain[which(pain$TOPMed.RNA.ID..CD71.. == i), "Therapy..at.RNA.collection."][1]
      master[which(master$`TOR ID_CD71` == i), "Chronic.Pain."] <- pain[which(pain$TOPMed.RNA.ID..CD71.. == i), "Chronic.Pain."][1]
    }
  }
}
for (i in master$`TOR ID_CD45`){
  if(!is.na(i)){
    if (i %in% pain$TOPMed.RNA.ID..CD45..){
      master[which(master$`TOR ID_CD45` == i), "Timepoint..at.RNA.collection."] <- pain[which(pain$TOPMed.RNA.ID..CD45.. == i), "Timepoint..at.RNA.collection."][1]
      master[which(master$`TOR ID_CD45` == i), "Therapy..at.RNA.collection."] <- pain[which(pain$TOPMed.RNA.ID..CD45.. == i), "Therapy..at.RNA.collection."][1]
      master[which(master$`TOR ID_CD45` == i), "Chronic.Pain."] <- pain[which(pain$TOPMed.RNA.ID..CD45.. == i), "Chronic.Pain."][1]
    }
  }
}

# Change metadata character columns to factors
master[] <- lapply(master, function(x) {
  if (lubridate::is.POSIXct(x)) {
    format(x, format = "%Y-%m-%d")
  } else {
    x
  }
})

#Add Additional subject IDs and sort by subject IDs
na_rows <- which(is.na(master$`Unique ID`))
master$`Unique ID`[na_rows] <- seq(701, 700 + length(na_rows))
master <- master %>% arrange(`Unique ID`)

#################### Combine IDs that can be found in a VCF into single column ######################
master$Combined_VCF_IDs <- NA
master$Combined_VCF_IDs <- ifelse(is.na(master$`St. Jude IDs`) | is.na(master$NWD_ID), NA, paste0(master$`St. Jude IDs`, "_", pain$NWD_ID))
master$Combined_VCF_IDs <- ifelse(is.na(master$Combined_VCF_IDs), master$`St. Jude IDs`, master$Combined_VCF_IDs)
master$Combined_VCF_IDs <- ifelse(is.na(master$Combined_VCF_IDs), master$NWD_ID, master$Combined_VCF_IDs)

col_order <- c("Unique ID", "DOB", "Sex", "Race", "Ethnicity", "Chronic.Pain.", "Timepoint..at.RNA.collection.", "Therapy..at.RNA.collection.", "TOR ID_CD45", "Sample ID_CD45", "CD45_NWGC_ID", "CD45_libprep_batch", "CD45_notes", "Date of collection_cd45", "date of extraction_CD45",
               "Age(yr) at collection_CD45", "TOR ID_CD71", "Sample ID_CD71", "CD71_NWGC_ID", "CD71_libprep_batch", "CD71_notes", "Date of collection_CD71", "date of extraction_CD71",
               "Age (years) at collection_CD71", "NWD_ID", "St. Jude IDs", "Combined_VCF_IDs")
master <- master %>% select(col_order, everything())

#master2 <- master
master_updated <- merge(master, lab.measurements[,c("MRN", "RESULT_DATE", colnames(lab.measurements)[20:ncol(lab.measurements)])], 
                        by.x =c("MRN", "Date of collection_CD71"), by.y=c("MRN", "RESULT_DATE"), all.x=TRUE)
master_updated <- master_updated %>% arrange(`Unique ID`)

# # Create an empty data frame to store the merged data
# master_updated <- data.frame()
# 
# # Loop through each row in the master data frame
# for (i in 1:nrow(master)) {
#   if(!is.na(master$MRN[i])){
#     # Subset the lab.measurements data frame for the current MRN in the master data frame
#     current_mrn_data <- subset(lab.measurements, lab.measurements$MRN == master[i, "MRN"])
#     #current_mrn_data <- lab.measurements[which(lab.measurements$MRN == master[i, "MRN"])]
#     #if(nrow(current_mrn_data)>1){print(nrow(current_mrn_data)) & print(master[i, "MRN"])}
#     #print(which(lab.measurements$MRN == master[i, "MRN"]))
#     
#     # Check if the RESULT_DATE in lab.measurements matches either the Date of collection_CD71 or Date of collection_cd45 in the master data frame
#     matching_dates <- which(current_mrn_data$RESULT_DATE == master[i, "Date of collection_CD71"] | current_mrn_data$RESULT_DATE == master[i, "Date of collection_cd45"])
#     
#     # If a match is found, copy the columns 20 and on from lab.measurements to the master data frame
#     if (length(matching_dates) > 0) {
#       current_row <- cbind(master[i,], current_mrn_data[matching_dates, 20:ncol(current_mrn_data)])
#       master_updated <- rbind(master_updated, current_row)
#     }
#   }
# }
# Arrange lab measurements by those that are most common
test <- data.frame(column=colnames(master_updated), count = colSums(!is.na(master_updated)))
test <- test %>% arrange(desc(count))
master_updated <- master_updated %>%
  select(test$column) %>%
  select(colnames(master), everything())

# Replace the original master data frame with the merged data frame
#master <- master_updated


############# 338 sex mismatch no longer, change note
# Summarize
master_updated$MRN <- NULL
master_updated$`Name 2` <- NULL
#pain <- as.data.frame(pain)
all.TORID <- master_updated[which(!is.na(master_updated$`TOR ID_CD45`)|!is.na(master_updated$`TOR ID_CD71`)),] #Get rows that have TORID
#all.TORID <- master[which(!is.na(master$`TOR ID_CD45`)|!is.na(master$`TOR ID_CD71`)),] #Get rows that have TORID
#passqc.TORID <- all.TORID[is.na(all.TORID$CD45_notes)&is.na(all.TORID$CD71_notes),]
passqc.TORID <- all.TORID[is.na(all.TORID$CD45_notes)&is.na(all.TORID$CD71_notes)&(!is.na(all.TORID$CD45_NWGC_ID)|!is.na(all.TORID$CD71_NWGC_ID)),]




passqc.TORID$Timepoint..at.RNA.collection. <- factor(passqc.TORID$Timepoint..at.RNA.collection.)
passqc.TORID$Chronic.Pain. <- factor(passqc.TORID$Chronic.Pain.)

summary.table <- passqc.TORID %>%
  group_by(`Timepoint..at.RNA.collection.`, `Chronic.Pain.`, .drop = F) %>%
  summarize(`CD45 TORID #` = sum(!is.na(`TOR ID_CD45`)),
            `CD45 Subject #` = sum(!is.na(`TOR ID_CD45`) & !is.na(`Unique ID`)),
            `CD45 Collection Date #` = sum(!is.na(`TOR ID_CD45`) & !is.na(`Date of collection_cd45`)),
            `CD45 VCF ID #` = sum(!is.na(`TOR ID_CD45`) & !is.na(Combined_VCF_IDs)),
            `CD71 TORID #` = sum(!is.na(`TOR ID_CD71`)),
            `CD71 Subject #` = sum(!is.na(`TOR ID_CD71`) & !is.na(`Unique ID`)),
            `CD71 Collection Date #` = sum(!is.na(`TOR ID_CD71`) & !is.na(`Date of collection_CD71`)),
            `CD71 VCF ID #` = sum(!is.na(`TOR ID_CD71`) & !is.na(Combined_VCF_IDs))
            )

summary.table <- summary.table %>%
  mutate(`sum_TORID` = `CD45 TORID #` + `CD71 TORID #`)

summary.table <- summary.table %>%
  select("Chronic.Pain.", "Timepoint..at.RNA.collection.", everything())

#################### Write xslx file ###################################################
#rm(pain, master, lab.measurements, current_mrn_data, correct_mets, test)
master_updated <- as.data.frame(master_updated)
passqc.TORID <- as.data.frame(passqc.TORID)
summary.table <- as.data.frame(summary.table)
options(java.parameters = "-Xmx8000m")
library(xlsx)
write.xlsx(master_updated, file = "../pain-omics-phenotype/Pain omics Phenotype_230203.xlsx", sheetName = "All",
           col.names = TRUE, row.names = FALSE, append = FALSE, showNA=FALSE)
write.xlsx(all.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230203.xlsx", sheetName = "All TORIDs",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
write.xlsx(passqc.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230203.xlsx", sheetName = "TORIDs Passing QC",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
# write.xlsx(passqc.labmeasurements.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230203.xlsx", sheetName = "TORIDs Passing QC with measurements",
#            col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
write.xlsx(summary.table, file = "../pain-omics-phenotype/Pain omics Phenotype_230203.xlsx", sheetName = "Summary",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
