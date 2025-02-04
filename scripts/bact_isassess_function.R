#start_date <- "2021-01-01"
#end_date <- "2021-12-31"


Individual_Sample_Assessment <- function(assessment_data, start_date, end_date){
        ### NOTE: this script is a source file for the "Bacteria Assessment" markdown.
        ### This means that it cannot be run until a specific portion of the markdown
        ### has been run first.
        
        bact_time <- data.frame(Date = seq(as.Date(start_date),
                                           as.Date(end_date),
                                           "day")) 
        
        # Column explanations:
        # 1. The month number is extracted from each date
        # 2. The extracted month number is used to index the "month.abb" object which
        # contains character values for month abbreviations. That is then pasted with
        # the extracted day number from each date
        # 3. An empty vector for season that will be added to later
        
        bact_time <- bact_time %>%
                mutate(bact_month = month(Date),
                       bact_total_time = paste(month.abb[bact_month],
                                               day(Date),
                                               sep = "-"),
                       season = rep(NA, 365))
        
        #designating sequential primary and secondary season dates throughout the POR
        
        Secondary1 <- seq.Date(as.Date("2018-4-1"), as.Date("2018-4-30"), by = "days")
        Secondary2 <- seq.Date(as.Date("2018-10-1"), as.Date("2019-4-30"), by = "days")
        Secondary3 <- seq.Date(as.Date("2019-10-1"), as.Date("2020-4-30"), by = "days")
        Secondary4 <- seq.Date(as.Date("2020-10-1"), as.Date("2021-4-30"), by = "days")
        Secondary5 <- seq.Date(as.Date("2021-10-1"), as.Date("2022-4-30"), by = "days")
        Secondary6 <- seq.Date(as.Date("2022-10-1"), as.Date("2023-3-31"), by = "days")
        
        Primary1 <- seq.Date(as.Date("2018-5-1"), as.Date("2018-9-30"), by = "days")
        Primary2 <- seq.Date(as.Date("2019-5-1"), as.Date("2019-9-30"), by = "days")
        Primary3 <- seq.Date(as.Date("2020-5-1"), as.Date("2020-9-30"), by = "days")
        Primary4 <- seq.Date(as.Date("2021-5-1"), as.Date("2021-9-30"), by = "days")
        Primary5 <- seq.Date(as.Date("2022-5-1"), as.Date("2022-9-30"), by = "days")
        
        # Because the row numbers will be consistent for every date value across most
        # years (except for leap year), hard-coding the index values for season will
        # be safe to then add in character data for whether that value is during the
        # primary or secondary season.
        
        bact_time$season[121:273] <- "PRIMARY"
        bact_time$season[c(1:120, 274:365)] <- "SECONDARY"
        
        # # We can then use the season designations to index the total date sequence. This
        # # will then allow us to index our raw data date values by these distinct objects.
        # 
        primary_season <- bact_time %>%
                filter(season == "PRIMARY")
        
        secondary_season <- bact_time %>%
                filter(season == "SECONDARY")
        
        # We then filter to our bacteria data and add in season qualifiers. The `mutate()`
        # is fairly similar to the "bact_time" data, but the `case_when()` will assign
        # designations depending on which raw data value is in each respective object's
        # date series. This works because those specific objects only contain dates that
        # occur during a given season. Because there are more nuances involved with these
        # groupings, we can create another pairing code specific to bacteria data. This
        # will be done to ease some processing later.
        
        # bact_data_is <- primary_data %>%
        #         mutate(month_num = month(CollectionDate),
        #                month_day = (paste(month.abb[month_num], day(CollectionDate), sep = "-")),
        #                bact_season = case_when(
        #                        month_day %in% primary_season$bact_total_time ~ "PRIMARY",
        #                        month_day %in% secondary_season$bact_total_time ~ "SECONDARY"),
        #                bact_pairing = paste(AU, bact_season, SiteType, ORW, year, sep = "-"))
        
        bact_data_is <- primary_data %>%
                mutate(bact_spec_season = case_when(
                        CollectionDate %in% Primary1 ~ "Primary1",
                        CollectionDate %in% Primary2 ~ "Primary2",
                        CollectionDate %in% Primary3 ~ "Primary3",
                        CollectionDate %in% Primary4 ~ "Primary4",
                        CollectionDate %in% Primary5 ~ "Primary5",
                        CollectionDate %in% Secondary1 ~ "Secondary1",
                        CollectionDate %in% Secondary2 ~ "Secondary2",
                        CollectionDate %in% Secondary3 ~ "Secondary3",
                        CollectionDate %in% Secondary4 ~ "Secondary4",
                        CollectionDate %in% Secondary5 ~ "Secondary5",
                        CollectionDate %in% Secondary6 ~ "Secondary6"),
                       month_num = month(CollectionDate),
                       month_day = (paste(month.abb[month_num], day(CollectionDate), sep = "-")),
                       bact_season = case_when(
                               month_day %in% primary_season$bact_total_time  ~ "PRIMARY",
                               month_day %in% secondary_season$bact_total_time  ~ "SECONDARY")) %>%
                mutate(bact_season = case_when(month_day %in% primary_season$bact_total_time & WS_SQMI < 10 & SiteType == "River/Stream" ~ "SECONDARY_WS",
                                               TRUE ~ bact_season),
                       bact_pairing = paste(AU, bact_season, SiteType, ORW, sep = "-"))
        
        
        
        # We can then create a summary with the primary filter of less than 8 samples.
        # This will filter out lakes and ORW's, but those will be filtered to next.
        
        bact_allwater <- bact_data_is %>%
                group_by(AU, bact_pairing, bact_spec_season, SiteType, ORW ) %>%
                summarise(sample_per_season = n(),
                          .groups = "drop") %>%
                arrange(AU) %>%
                mutate(bact_seasid = paste(AU, bact_spec_season, sep = "-"))
        
        # Filtering data that does not meet the minimum sample requirement
        
        bact_allwater <- bact_allwater %>%
                filter(sample_per_season >= 8)
        
        # This will be the data that can be assessed under the individual sample methodology
        
        bact_data_is_filter <- bact_data_is %>%
                mutate(bact_seasid = paste(AU, bact_spec_season, sep = "-")) %>%
                filter(bact_seasid %in% bact_allwater$bact_seasid) %>%
                select(-bact_seasid)
        
        #calculating excrescences and number of samples
        
        bact_assessment_is_setup <- bact_data_is_filter %>%
                group_by(bact_pairing) %>%
                summarise(n_sample = n(),
                          n_exceed = case_when(all(bact_season == "PRIMARY" & SiteType =="River/Stream" & ORW == "N") ~ sum(ResultValue > 410),
                                               all(bact_season == "SECONDARY" & SiteType =="River/Stream" & ORW == "N") ~ sum(ResultValue > 2050),
                                               all(bact_season == "PRIMARY" & SiteType =="River/Stream" & ORW == "Y") ~ sum(ResultValue > 298),
                                               all(bact_season == "SECONDARY" & SiteType =="River/Stream" & ORW == "Y") ~ sum(ResultValue > 1490),
                                               all(bact_season == "PRIMARY" & SiteType =="Lake/Reservoir") ~ sum(ResultValue > 298),
                                               all(bact_season == "SECONDARY" & SiteType =="Lake/Reservoir") ~ sum(ResultValue > 1490),
                                               all(bact_season == "SECONDAR_WS") ~ sum(ResultValue > 2050)),
                          .groups = "drop")
        
        # calculating exceedance rates and designating au's where the 25% exceedances rate is violated 
        
        bact_assessment_is_sum <- bact_assessment_is_setup %>%
                mutate(exceed_rate = n_exceed/n_sample) %>%
                mutate(ASSESSMENT = case_when(exceed_rate >= .25 ~ "FAIL",
                                              T ~ "ATTAIN"))
        
        bact_assessment_season_sum <- bact_data_is_filter %>%
                group_by(bact_spec_season, bact_pairing) %>%
                summarise(n_sample = n())
        
        #Setting up the final export for the main code, filtering out some unnecessary parameters
        
        bact_assessment_is <- bact_assessment_is_sum %>%
                select(-n_sample, -n_exceed, -exceed_rate)
        
        bact_assessment_is$method <- "IS"
        
        return(bact_assessment_is)
}
