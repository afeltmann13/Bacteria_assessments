### This script is intended to be run after the Markdown file is run. The objects
### created by that file are used to isolate raw data to QC by hand.

# Read in data for the 2020 listings and the 2022 assessment results

bact_2020_listings <- readxl::read_xlsx("E:/2022 Integrated Report/Data Assessments/BACTERIA ASSESSMENTS/data QCs/Bact_Listing_2020.xlsx")
bact_2022_results <- read.csv("E:/2022 Integrated Report/Data Assessments/BACTERIA ASSESSMENTS/output/BACTERIA ASSESSMENTS.csv")

# Create the objects for 2022 listings and the 2020 listings that were also in
# the 2022 data set.

failures_2022 <- bact_2022_results %>%
        filter(ASSESSMENT == "FAIL") %>% 
        mutate(QC_Type = "2022 Failure")

bact_carryover <- bact_2022_results %>%
        filter(AU %in% bact_2020_listings$AU) %>% 
        mutate(QC_Type = "AUs In 2020 List")

# The previous mutates designated whether the object is intended for a failure
# or previous listing check. Those can then be consolidated into one object.

QC_Listing <- bind_rows(failures_2022,
                        bact_carryover)

# Use the QC listings as a filter to isolate the raw data we need to export to
# assess by hand. This is done for individual sample and geometric mean methods.

is_qc_data <- bact_data_is_filter %>%
        filter(AU %in% QC_Listing$AU) %>% 
        arrange(AU, CollectionDate)

gm_qc_data <- lake_ORW_raw %>% 
        filter(AU %in% QC_Listing$AU) %>%
        arrange(CollectionDate)

# Write the data to xlsx files that we can use to run the assessment QCs by hand

writexl::write_xlsx(is_qc_data, "E:/2022 Integrated Report/Data Assessments/BACTERIA ASSESSMENTS/data QCs/QC_Data_Files/IS_QC_Data.xlsx")
writexl::write_xlsx(gm_qc_data, "E:/2022 Integrated Report/Data Assessments/BACTERIA ASSESSMENTS/data QCs/QC_Data_Files/GM_QC_Data.xlsx")
