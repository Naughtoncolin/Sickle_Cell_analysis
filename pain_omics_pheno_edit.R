setwd("C:/Users/Naugh/Dropbox (GaTech)/Gibson/Working/SickleCell/metadata/RNA/")
library(dplyr)
library(tidyr)

pain <- readxl::read_xlsx("../pain-omics-phenotype/Pain omics Phenotype_221117.xlsx")

###################### Add blood collection dates and MRN info to "Pain omics Phenotype" ##########################3
# Add columns for data that will be added later
pain[c("Blood_collection_date", "Blood_processing_date", "RNA_extraction_date", "MRN", 
       "Birthday", "Race", "Ethnicity", "Sex")] <- NA

###Try with lab measurement MRN and Birthday
lab.measurements <- read.csv("TOR_lab-values.csv")

# Format the dates to "YYYY-MM-DD"
for(j in 1:nrow(lab.measurements)){
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
for(i in 1:nrow(pain)){
  if(!is.na(pain$Investigator.RNA.Plasma.ID[i])){
    for(j in 1:nrow(lab.measurements)){
      if(pain$Investigator.RNA.Plasma.ID[i] == lab.measurements$Sample.ID_1[j]){
        pain$MRN[i] <- lab.measurements$MRN[j]
        pain$Birthday[i] <- lab.measurements$DOB[j]
        pain$Race[i] <- lab.measurements$RACE[j]
        pain$Ethnicity[i] <- lab.measurements$ETHNICITY[j]
        pain$Sex[i] <- lab.measurements$SEX[j]
      } else if(!is.na(lab.measurements$Sample.ID_2[j])){
        if(pain$Investigator.RNA.Plasma.ID[i] == lab.measurements$Sample.ID_2[j]){
          pain$MRN[i] <- lab.measurements$MRN[j]
          pain$Birthday[i] <- lab.measurements$DOB[j]
          pain$Race[i] <- lab.measurements$RACE[j]
          pain$Ethnicity[i] <- lab.measurements$ETHNICITY[j]
          pain$Sex[i] <- lab.measurements$SEX[j]
        }
      } else if(!is.na(lab.measurements$Sample_3[j])){
        if(pain$Investigator.RNA.Plasma.ID[i] == lab.measurements$Sample_3[j]){
          pain$MRN[i] <- lab.measurements$MRN[j]
          pain$Birthday[i] <- lab.measurements$DOB[j]
          pain$Race[i] <- lab.measurements$RACE[j]
          pain$Ethnicity[i] <- lab.measurements$ETHNICITY[j]
          pain$Sex[i] <- lab.measurements$SEX[j]
        }
      }        
    }
  }
}

# Look through file with date and MRN info
for (excell_sheet in c("CD45+", "CD71+")){
  master <- readxl::read_xlsx("RNA Masterlist with names.xlsx", sheet = excell_sheet)
  # Edit the file to correct where "sex" and "race" were mixed up
  # Within the CD45+ sheet, "Gender" was manually changed to "Sex"(CD71+ had "Sex" as column name) and the 2nd "Race"
  # column was removed. 
  # This currently doesn't fill in anything on "pain omics phenotype" because the issues with mixed-up#########################
  # race & sex wasn't fixed. EDIT IN THE FUTURE ###############################################################################################################
  for (i in 1:nrow(master)) {
    # Check if the value in Sex is not Male or Female and the value in Race is not NA
    if (!master$Sex[i] %in% c("Male", "Female") & !is.na(master$Race[i])) {
      # Swap the values in the Sex and Race columns
      temp <- master$Sex[i]
      master$Sex[i] <- master$Race[i]
      master$Race[i] <- temp
    }
  }
  # Loop through TORIDs within "Pain omics phenotype"
  for (i in 1:nrow(pain)) {
    for (col in c("TOPMed.RNA.ID..CD45..", "TOPMed.RNA.ID..CD71..")) {
      if (!is.na(pain[i, col])) {
        # Find the row in "RNA Masterlist with names.xlsx" where the TOR.ID_1, TOR.ID_2, or TOR.ID_3 column matches the current TORID
        # and copy over the information from the matching row.
        match_row <- match(pain[i, "Investigator.RNA.Plasma.ID"], master$`Sample ID`)
        if (!is.na(match_row)) {
          if(is.na(pain$Blood_collection_date[i])){pain$Blood_collection_date[i] <- as.character(master$`Date of collection`[match_row])}
          if(is.na(pain$Blood_processing_date[i])){pain$Blood_processing_date[i] <- as.character(master$`Date of processing`[match_row])}
          if(is.na(pain$RNA_extraction_date[i])){pain$RNA_extraction_date[i] <- as.character(master$`date of extraction`[match_row])}
          if(is.na(pain$MRN[i])){pain$MRN[i] <- master$MRN[match_row]}
          if(is.na(pain$Timepoint..at.RNA.collection.[i])){pain$MRN[i] <- master$Group[match_row]} #Didn't yield
          if(is.na(pain$Therapy..at.RNA.collection.[i])){pain$MRN[i] <- master$Treatment[match_row]} #Didn't yield
          if(is.na(pain$Investigator.DNA.ID[i])){pain$Investigator.DNA.ID[i] <- master$`Sample ID`[match_row]}
          
          # Calculations for age at blood collection.
          # Note: The DOB doesn't seem to be correct as there are some individuals with blood collection dates also as DOB
          if(is.na(pain$Birthday[i])){pain$Birthday[i] <- as.character(master$DOB[match_row])} #Changed from DOB to Birthday because of conflicts during merge further below
          #start_date <- as.Date(pain$DOB[i], "%Y-%m-%d")
          #print(start_date)
          #end_date <- as.Date(pain$Blood_collection_date[i], "%Y-%m-%d")
          #print(end_date)
          #time_diff <- difftime(end_date,start_date, units = "days")
          #time_diff <- as.numeric(time_diff/365.25636) # 365.25636 days per year via WolframAlpha
          #pain$Age_at_collection[i] <- time_diff
        }
        match_row <- NA
      }
    }
    #Can this be removed???##################################################################
    if (is.na(pain[i, "MRN"])) {
      match_row <- match(pain[i, "Investigator.RNA.Plasma.ID"], master$`Sample ID`)
      if (!is.na(match_row)) {
        if(is.na(pain$MRN[i])){pain$MRN[i] <- master$MRN[match_row]}
      }
    }
  }
}

# Replace character NA with NULL NA
pain$MRN <- sub("N/A", NA, pain$MRN)

##################### Add NWGC IDs, batch, and notes from sequencing center####################################################
pain$CD45_NWGC_ID <- NA
pain$CD71_NWGC_ID <- NA
pain$CD45_libprep_batch <- NA
pain$CD71_libprep_batch <- NA
pain$CD45_notes <- NA
pain$CD71_notes <- NA
pain$Investigator.Sex <- NA

torids <- readxl::read_xlsx("lookup_pharmhu_topmed_to5_rnaseq_1_final.xlsx")
# Some of the sexes were missing, but the sequencing center determined them, so
# create a unified list of sex; if the investigator sex is NA, then replace with
# the NWGC assay-determined sex.
torids$`Investigator Sex` <- sub("Unknown", NA, torids$`Investigator Sex`)
torids$`Assay Sex` <- sub("N/A", NA, torids$`Assay Sex`)
torids$`Investigator Sex` <- ifelse(is.na(torids$`Investigator Sex`), torids$`Assay Sex`, torids$`Investigator Sex`)

# Check if any of the TORIDs can be found in the metadata provided by the NWGC and fill in info.
# TODO: This code code have 1 if statement if you regex out the "CD45" or "CD71" from the columns in "pain".
for (i in 1:nrow(torids)){
  matching_row <- which(pain$TOPMed.RNA.ID..CD45.. %in% torids[i,"TORID"])
  if(length(matching_row) == 1){
    pain[matching_row,"CD45_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    pain[matching_row,"CD45_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    pain[matching_row,"CD45_notes"] <- torids[i,"Notes"]
    pain[matching_row,"Investigator.Sex"] <- torids[i,"Investigator Sex"]
    matching_row <- NA
  }
  matching_row <- which(pain$TOPMed.RNA.ID..CD71.. %in% torids[i,"TORID"])
  if(length(matching_row) == 1){
    pain[matching_row,"CD71_NWGC_ID"] <- torids[i,"NWGC Sample ID"]
    pain[matching_row,"CD71_libprep_batch"] <- torids[i,"Pick Plate ID(s) - Indicates which samples were prepped together"]
    pain[matching_row,"CD71_notes"] <- torids[i,"Notes"]
    pain[matching_row,"Investigator.Sex"] <- torids[i,"Investigator Sex"]
    matching_row <- NA
  }
  matching_row <- NA
}

##################### Add Missing MRN & DNA IDs ########################
ids <- readxl::read_xlsx("Matched IDs_deident.xlsx")

for (i in 1:nrow(pain)) {
  if(!is.na(pain$`Investigator.RNA.Plasma.ID`[i])){
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
    match_result <- NA
  }
}

########## Fill in missing values with values from other rows where they should be the same###################

# Create a new data frame with the unique values in the "Investigator.DNA.ID" column
# and their corresponding "MRN" values
unique_ids <- unique(pain[which(!is.na(pain$Investigator.DNA.ID)&!is.na(pain$MRN)), c("Investigator.DNA.ID", "MRN", "Investigator.Sex")])
# MRN 3001485403 had duplicate entries with DNA IDs "GM014" & "MET039". Workaround with the following code:
# Need to remember to manually add other entry ##############################!!!!!!!!!!!!!!!!!!!!!!!!
# Fill in missing values within the unique ID manifest
## COuld this be used to get rid of the code after?
unique_ids <- unique_ids %>%
  group_by(MRN) %>%
  fill(-MRN, .direction = "updown") %>%
  slice(1) %>%
  ungroup()

# Fill in MRN, DNA ID, and Sex
# Should this be done after filling in with the lab values sheet since it also contains more comprehensive Race/Ethnicity Info?
for (i in 1:nrow(pain)) {
  # Check if the "MRN" value is NA
  if (is.na(pain$MRN[i])) {
    if (pain$`Investigator.DNA.ID`[i] %in% unique_ids$`Investigator.DNA.ID`) {
      # If the "MRN" value is NA, look up the corresponding "MRN" value
      # in the "unique_ids" data frame and set it in the "pain" data frame
      pain$MRN[i] <- unique_ids$MRN[unique_ids$`Investigator.DNA.ID` == pain$`Investigator.DNA.ID`[i]]
    }
  }
  if (is.na(pain$`Investigator.DNA.ID`[i])) {
    if (pain$MRN[i] %in% unique_ids$MRN) {
      # If the "MRN" value is NA, look up the corresponding "MRN" value
      # in the "unique_ids" data frame and set it in the "pain" data frame
      pain$`Investigator.DNA.ID`[i] <- unique_ids$`Investigator.DNA.ID`[unique_ids$MRN == pain$MRN[i]]
    }
  }
  if (is.na(pain$`Investigator.Sex`[i])) {
    if (pain$MRN[i] %in% unique_ids$MRN) {
      # If the "MRN" value is NA, look up the corresponding "MRN" value
      # in the "unique_ids" data frame and set it in the "pain" data frame
      pain$`Investigator.Sex`[i] <- unique_ids$`Investigator.Sex`[unique_ids$MRN == pain$MRN[i]]
    }
  }
}


# Visual inspection indicated no more missing IDs could be filled in based on rows where they should
# be duplicates e.g. the same MRN
sum(!is.na(pain$MRN)) #585
sum(!is.na(pain$Investigator.DNA.ID)) #811

####### Make subject ID column#################
pain$Subject_ID <- pain$MRN
for(i in 1:nrow(pain)){
  if(is.na(pain$Subject_ID[i])){
    pain$Subject_ID[i] <- pain$Investigator.DNA.ID[i]
  }
}

#################### Combine IDs that can be found in a VCF into single column ######################

pain$Combined_VCF_IDs <- NA
pain$Combined_VCF_IDs <- ifelse(is.na(pain$St.Jude.DNA.ID) | is.na(pain$NWD_ID), NA, paste0(pain$St.Jude.DNA.ID, "_", pain$NWD_ID))
pain$Combined_VCF_IDs <- ifelse(is.na(pain$Combined_VCF_IDs), pain$St.Jude.DNA.ID, pain$Combined_VCF_IDs)
pain$Combined_VCF_IDs <- ifelse(is.na(pain$Combined_VCF_IDs), pain$NWD_ID, pain$Combined_VCF_IDs)

#
sum(!is.na(pain$TOPMed.RNA.ID..CD45..) & !is.na(pain$Combined_VCF_IDs)) #306
sum(!is.na(pain$TOPMed.RNA.ID..CD71..) & !is.na(pain$Combined_VCF_IDs)) #330

pain <- pain %>%
  select("MRN","Investigator.DNA.ID", "St.Jude.DNA.ID", "NWD_ID", "Chronic.Pain.", "TOPMed.RNA.ID..CD45..","CD45_NWGC_ID", "CD71_notes", "TOPMed.RNA.ID..CD71..", "CD71_NWGC_ID", "CD45_libprep_batch", 
         "CD71_libprep_batch", "Blood_collection_date", "Blood_processing_date","RNA_extraction_date", 
         "Chronic.Pain.", "Timepoint..at.RNA.collection.", "Therapy..at.RNA.collection.", everything())

##################### Include info lab measurements file #################################3
# Also has demographic and sex info that could be used to fill in missing info
lab.measurements <- read.csv("TOR_lab-values.csv")

# Format the dates to "YYYY-MM-DD"
for(j in 1:nrow(lab.measurements)){
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

###Testing
# Reorder the columns
pain <- pain %>%
  select("MRN","Investigator.DNA.ID", "Subject_ID", "Chronic.Pain.", 
         "Birthday", "Blood_collection_date", "Blood_processing_date","RNA_extraction_date", 
         "Chronic.Pain.", "Timepoint..at.RNA.collection.", "Therapy..at.RNA.collection.", "Race",
         "Ethnicity", "TOPMed.RNA.ID..CD45..", "TOPMed.RNA.ID..CD71..", everything())
pain <- pain %>%
  select("MRN","Investigator.DNA.ID", "Subject_ID", "Chronic.Pain.", 
         "Birthday", "Blood_collection_date", "Race", "Ethnicity", "Combined_VCF_IDs", 
         "Chronic.Pain.", "Timepoint..at.RNA.collection.", "Therapy..at.RNA.collection.", "Race",
         "Ethnicity", "TOPMed.RNA.ID..CD45..", "TOPMed.RNA.ID..CD71..", everything())
### Add Met IDs 
#met<- readxl::read_xlsx("Edited Metformin IDs.xlsx")
#pain$Met_ID <- NA
#pain$Met_ID <- met$"Investigator ID1 (DNA ID1)"[match(pain$Investigator.RNA.Plasma.ID, met$`Sample ID_CD71 (RNA ID)`)]
#There are inconsistencies between the MET IDs in this excell and what is documented in the metadata file already

# Remove patient info and make anonymized subject ID
mrn <- data.frame("mrn"=unique(pain$Subject_ID), "id"=1:length(unique(pain$Subject_ID)))
pain$Subject_ID <- mrn$id[match(pain$Subject_ID, mrn$mrn)]

#pain$DOB <- NULL # DOB column was interrupting merge
pain <- as.data.frame(pain)
all.TORID <- pain[which(!is.na(pain$TOPMed.RNA.ID..CD45..)|!is.na(pain$TOPMed.RNA.ID..CD71..)),] #Get rows that have TORID
passqc.TORID <- all.TORID[is.na(all.TORID$CD45_notes)&is.na(all.TORID$CD71_notes),]


# Make merged data frame with lab measurements based on collection date, then make another dataframe
# based on whether the plasma IDs match, then integrate that with the main metadata.
# THIS NEEDS TO BE REDONE IN A BETTER WAY
#merged_df <- merge(pain, lab.measurements[,c("Sample.ID_1", "Sample.ID_2", "Sample_3", "RESULT_DATE")], by.x = "Blood_collection_date", by.y = "RESULT_DATE")
merged_df <- merge(passqc.TORID, lab.measurements, by.x = c("MRN","Blood_collection_date"), 
                   by.y = c("MRN","RESULT_DATE"))

filtered_df <- data.frame()
for(i in 1:nrow(merged_df)){
  if(is.na(merged_df$Investigator.RNA.Plasma.ID[i]) || is.na(merged_df$Sample.ID_1[i])){
    next
    # Bind the current row to filtered_df in plasma IDs match
  } else if(merged_df$Investigator.RNA.Plasma.ID[i] == merged_df$Sample.ID_1[i]){ 
    filtered_df <- rbind(filtered_df, merged_df[i,])
  } else if(is.na(merged_df$Investigator.RNA.Plasma.ID[i]) || is.na(merged_df$Sample.ID_2[i])){
    next
  } else if(merged_df$Investigator.RNA.Plasma.ID[i] == merged_df$Sample.ID_2[i]){ 
    filtered_df <- rbind(filtered_df, merged_df[i,])
  } else if(is.na(merged_df$Investigator.RNA.Plasma.ID[i]) || is.na(merged_df$Sample_3[i])){
    next
  } else if(merged_df$Investigator.RNA.Plasma.ID[i] == merged_df$Sample_3[i]){ 
    filtered_df <- rbind(filtered_df, merged_df[i,])
  }
}

# Merge lab measurements with TORIDs that passed QC; this makes the other higher level dataframes
# easier to read.
passqc.labmeasurements.TORID <- merge(passqc.TORID, filtered_df, by = colnames(passqc.TORID), all.x = T)

# Fill in missing sex for 1 individual
passqc.labmeasurements.TORID$Sex <- ifelse(is.na(passqc.labmeasurements.TORID$Sex), passqc.labmeasurements.TORID$SEX, passqc.labmeasurements.TORID$Sex)
passqc.labmeasurements.TORID$Sex <- ifelse(!is.na(passqc.labmeasurements.TORID$Sex), passqc.labmeasurements.TORID$Investigator.Sex, passqc.labmeasurements.TORID$Sex)
passqc.labmeasurements.TORID$Sex <- sub("Male", "M", passqc.labmeasurements.TORID$Sex)
passqc.labmeasurements.TORID$Sex <- sub("Female", "F", passqc.labmeasurements.TORID$Sex)
passqc.labmeasurements.TORID$Investigator.Sex <- NULL

#Remove columns that had been in lab measurements that are no longer useful.
col_to_remove <- c("TORID_1_Date_of_collection", "TORID_2_Date_of_collection", "TORID_2_Date_of_collection", "TOR.ID_1", "TOR.ID_2", 
                   "TOR.ID_3", "Sample.ID_1", "Sample.ID_2", "Sample_3", "SEX", "DOB", "RACE", "ETHNICITY")
passqc.labmeasurements.TORID <- passqc.labmeasurements.TORID %>% select(-one_of(col_to_remove))

#### Make data frame with number of occurrences of different lab measurements
# Create empty data frame with correct number of rows
results_df <- data.frame(matrix(nrow = length(passqc.labmeasurements.TORID)))
# Count number of non-NA values in each column
results_df$non_na_count <- colSums(!is.na(passqc.labmeasurements.TORID))
# Add column names as row names
row.names(results_df) <- colnames(passqc.labmeasurements.TORID)

# Remove any columns that have more than a certain number of NA
threshold <- 100
passqc.labmeasurements.TORID <- passqc.labmeasurements.TORID[, colSums(is.na(passqc.labmeasurements.TORID)) 
                                                             < nrow(passqc.labmeasurements.TORID) - threshold]

#### Fill in missing things again...
#delete <- passqc.labmeasurements.TORID
#passqc.labmeasurements.TORID <- delete
sum(!is.na(passqc.labmeasurements.TORID$Birthday))
sum(!is.na(passqc.labmeasurements.TORID$Race))
passqc.labmeasurements.TORID <- passqc.labmeasurements.TORID %>%
  group_by(MRN) %>%
  fill(c(Birthday, Race, Ethnicity, Sex, Combine_VCF_IDs, St.Jude.DNA.ID, NWD_ID,
         Avg.pain.admissions.per.yr), .direction = "updown") %>%
  ungroup()


############################# Summarize #########################
#pain <- as.data.frame(pain)
#all.TORID <- pain[which(!is.na(pain$TOPMed.RNA.ID..CD45..)|!is.na(pain$TOPMed.RNA.ID..CD71..)),] #Get rows that have TORID
#passqc.TORID <- all.TORID[is.na(all.TORID$CD45_notes)&is.na(all.TORID$CD71_notes),]

#test <- all.TORID[which(!is.na(all.TORID$TOPMed.RNA.ID..CD45..)& all.TORID$Timepoint..at.RNA.collection.=="Baseline"),]

# Keeping rows with all 0 while grouping requires the grouping variables to be factors
passqc.labmeasurements.TORID$Timepoint..at.RNA.collection. <- factor(passqc.TORID$Timepoint..at.RNA.collection.)
passqc.labmeasurements.TORID$Chronic.Pain. <- factor(passqc.TORID$Chronic.Pain.)

summary.table <- passqc.labmeasurements.TORID %>%
  group_by(`Timepoint..at.RNA.collection.`, `Chronic.Pain.`, .drop = F) %>%
  summarize(`CD45 TORID #` = sum(!is.na(`TOPMed.RNA.ID..CD45..`)),
            `CD45 Subject #` = sum(!is.na(`TOPMed.RNA.ID..CD45..`) & !is.na(`Subject_ID`)),
            `CD45 Collection Date #` = sum(!is.na(TOPMed.RNA.ID..CD45..) & !is.na(Blood_collection_date)),
            `CD45 VCF ID #` = sum(!is.na(TOPMed.RNA.ID..CD45..) & !is.na(Combined_VCF_IDs)),
            `CD71 TORID #` = sum(!is.na(`TOPMed.RNA.ID..CD71..`)),
            `CD71 Subject #` = sum(!is.na(`TOPMed.RNA.ID..CD71..`) & !is.na(`Subject_ID`)),
            `CD71 Collection Date #` = sum(!is.na(TOPMed.RNA.ID..CD71..) & !is.na(Blood_collection_date)),
            `CD71 VCF ID #` = sum(!is.na(TOPMed.RNA.ID..CD71..) & !is.na(Combined_VCF_IDs)))

summary.table <- summary.table %>%
  mutate(`sum_TORID` = `CD45 TORID #` + `CD71 TORID #`)

summary.table <- summary.table %>%
  select("Chronic.Pain.", "Timepoint..at.RNA.collection.", everything())

# Sex Count
sum(!is.na(all.TORID$Investigator.Sex) & all.TORID$Investigator.Sex == "Male")

# Ids removed based on NWGS QC
sum(!is.na(all.TORID$TOPMed.RNA.ID..CD45..) & !is.na(all.TORID$CD45_notes)) #29
sum(!is.na(all.TORID$TOPMed.RNA.ID..CD71..) & !is.na(all.TORID$CD71_notes)) #17

# Unique VCF IDs
sum(!is.na(unique(passqc.labmeasurements.TORID$Combined_VCF_IDs)))



#################### Write xslx file ###################################################

pain <- as.data.frame(pain)
summary.table <- as.data.frame(summary.table)

library(xlsx)
write.xlsx(pain, file = "../pain-omics-phenotype/Pain omics Phenotype_230111.xlsx", sheetName = "All",
           col.names = TRUE, row.names = FALSE, append = FALSE, showNA=FALSE)
write.xlsx(all.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230111.xlsx", sheetName = "All TORIDs",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
write.xlsx(passqc.TORID, file = "../pain-omics-phenotype/Pain omics Phenotype_230111.xlsx", sheetName = "TORIDs Passing QC",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
write.xlsx(summary.table, file = "../pain-omics-phenotype/Pain omics Phenotype_230111.xlsx", sheetName = "Summary",
           col.names = TRUE, row.names = FALSE, append = TRUE, showNA=FALSE)
