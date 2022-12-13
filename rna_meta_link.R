library(dplyr)
library(tidyr)
library(xlsx)
library(purrr)

metadata <- readxl::read_xlsx("deident_Copy of KM_SC_GENOTYPE_LAB_VALUES - from past to 020618 on 020718.xlsx")
#metadata2 <- metadata
metadata <- metadata2
#c("TOR ID_1", "COMPONENT_NAME", "ORD_VALUE")

# Identify columns with all missing values
na_columns <- which(colSums(is.na(metadata)) == nrow(metadata))

# Drop columns with all missing values
metadata <- select(metadata, -na_columns)

# Replace character "NA" values in the TOR ID_1 column with the NULL NA value in R
metadata$`TOR ID_1` <- ifelse(metadata$`TOR ID_1` == "NA", NA, metadata$`TOR ID_1`)

# Identify rows with an NA value in columns 2, 3, and 4
na_rows <- which(is.na(metadata[,2]) & is.na(metadata[,3]) & is.na(metadata[,4]))

# Remove rows with an NA value in these columns
metadata <- filter(metadata, !(row_number() %in% na_rows))

metadata<- unique(metadata)
#metadata <- unique(metadata[,c("TOR ID_1", "COMPONENT_NAME", "ORD_VALUE", "ORDERING_DATE")])
metadata <- metadata[which(metadata$COMPONENT_NAME!="CULTURE SOURCE"),]

# Trying to use the spread function threw errors for these rows, for example if a patient on the same date had 2 different blood type results
duplicate.rows <- c(109541, 109564, 128672, 128673, 180823, 180837, 180822, 180834, 180821, 180833, 195919, 195921, 195923, 195925, 195926, 195928, 195932, 195934, 195935, 209319, 209321, 209323, 209325, 209335, 209326, 209328, 209330, 209332, 209333, 109540, 109563, 109549, 109562, 9657, 9662, 19929, 19940, 14776, 14783, 14831, 14835, 55883, 55891, 60153, 60158, 60161, 93419, 93426, 110441, 110447, 125194, 125200, 125254, 125259, 125380, 125385, 176921, 176927, 177094, 177100, 177117, 177119, 177126, 177130, 177382, 177388, 201637, 201645, 201647, 215334, 215338, 215339, 215341, 215342, 215398, 215404, 215407, 215526, 215532, 259617, 259619, 14727, 14732, 14778, 14785, 14829, 14834, 14884, 14886, 45517, 45525, 45573, 45575, 45622, 45625, 55800, 55804, 55885, 55894, 93416, 93425, 93547, 93555, 108568, 108574, 108653, 108656, 110443, 110448, 110457, 110460, 125199, 125202, 125256, 125258, 201422, 201424, 201639, 201644, 215254, 215258, 215400, 215403, 215469, 215473, 259621, 259626, 60157, 60160, 60162, 110452, 110458, 176923, 176929, 177099, 177101, 177114, 177120, 177122, 177129, 177384, 177390, 201633, 201641, 201646, 215330, 215336, 215340, 215394, 215402, 215406, 93409, 93411, 180820, 180832, 19926, 19938, 19896, 19911, 14723, 14730, 55797, 55802, 55881, 55889, 14726, 14731, 55793, 55801, 55886, 55892, 14722, 14729, 14775, 14782, 14832, 14836, 14883, 14885, 45522, 45523, 45526, 45527, 45571, 45576, 45624, 45626, 55798, 55803, 55882, 55890, 93420, 93427, 93553, 93556, 93605, 93610, 108573, 108575, 108650, 108655, 110440, 110446, 110456, 110459, 110467, 110472, 125195, 125201, 125257, 125260, 125325, 125326, 125379, 125384, 201423, 201425, 201636, 201642, 215251, 215257, 215333, 215337, 215397, 215408, 215409, 215466, 215472, 259618, 259625, 55880, 55888)
duplicates <- metadata[duplicate.rows,]
write.xlsx(duplicates, file="duplicates.xlsx")
metadata <- metadata[-duplicate.rows,]

# Identify rows with an NA value for whatever value was trying to be measured.
na_values <- which(is.na(metadata$ORD_VALUE))

# Remove rows with an NA value in these columns
metadata <- filter(metadata, !(row_number() %in% na_values))


test <- metadata
test$COMPONENT_VALUE <- paste(test$ORD_VALUE, test$REFERENCE_UNIT, ' ')
unneeded_col <- c("PAT_ID", "MRN", "PAT_ENC_CSN_ID", "PROC_CODE", "PROC_NAME", "ORDERING_DATE", "ORD_VALUE", "ORD_NUM_VALUE", "REFERENCE_UNIT", "COMPONENT_ID", "COMPONENT_ABBR")
test <- select(test, -unneeded_col)
test <- filter(test, !grepl("CANCELLED", COMPONENT_VALUE))
test <- unique(test)
# Trying to use the spread function threw errors for these rows, for example if a patient on the same date had 2 different blood type results
duplicate.rows <- readLines("rows2-long.txt")
duplicate.rows <- as.numeric(duplicate.rows)
duplicates2 <- test[duplicate.rows,]
duplicates2 <- duplicates2[order(duplicates2$`TOR ID_1`, duplicates2$RESULT_DATE, duplicates2$COMPONENT_NAME),]
write.xlsx(duplicates2, file="duplicates2.xlsx")
test <- test[-duplicate.rows,]
test <- test[order(test$`TOR ID_1`, test$RESULT_DATE, test$COMPONENT_NAME),]


results_test <- test %>% 
  spread("COMPONENT_NAME", "COMPONENT_VALUE")
write.xlsx(results_test, file="TOR_lab-values.xlsx")

results <- metadata %>% 
  spread("COMPONENT_NAME", "ORD_VALUE")

# write.xlsx(results, file="TOR_lab-values.xlsx")

colnames(results_test)
summary(results_test)
