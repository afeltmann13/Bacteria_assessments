
Geometric_Mean_Assessment<- function(assessment_data, start_date, end_date){
        
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
        
        # Because the row numbers will be consistent for every date value across most
        # years (except for leap year), hard-coding the index values for season will
        # be safe to then add in character data for whether that value is during the
        # primary or secondary season.
        
        bact_time$season[121:273] <- "PRIMARY"
        bact_time$season[c(1:120, 274:365)] <- "SECONDARY"
        
        # We can then use the season designations to index the total date sequence. This
        # will then allow us to index our raw data date values by these distinct objects.
        
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
        
        bact_data_gm <- primary_data %>%
            mutate(month_num = month(CollectionDate),
                   month_day = (paste(month.abb[month_num], day(CollectionDate), sep = "-")),
                   bact_season = case_when(
                       month_day %in% primary_season$bact_total_time ~ "PRIMARY",
                       month_day %in% secondary_season$bact_total_time ~ "SECONDARY"),
                   bact_pairing = paste(AU, bact_season, SiteType, ORW, sep = "-"))
        
        # We can then create a summary with the primary filter of less than 8 samples.
        # This will filter out lakes and ORW's, but those will be filtered to next.
        
        # bact_primary_filter <- bact_data %>%
        #     group_by(AU, bact_season, SiteType, ORW) %>%
        #     summarise(sample_per_season = n(),
        #               .groups = "drop") %>%
        #     arrange(AU) %>%
        #     filter(sample_per_season < 5) %>%
        #     filter(SiteType != "LAKES") %>%
        #     filter(ORW != "Y")
        
        # # Filter down to just lakes and ORWs. We will then need to filter further to
        # data not having 5 samples in a 30-day period. 
        
        bact_lake_ORW <- bact_data_gm %>%
            group_by(AU, bact_season, SiteType, ORW) %>%
            summarise(sample_per_season = n(),
                      .groups = "drop") %>%
            arrange(AU) %>%
            filter(SiteType == "LAKES" | ORW == "Y") %>%
            mutate(bact_pairing = paste(AU, bact_season, SiteType, ORW, sep = "-"))
        
        bact_lake_ORW <- bact_lake_ORW %>%
            filter(sample_per_season >= 5)
        
        # Filter to the raw data values needed for the date difference calculation.
        
        lake_ORW_raw <- bact_data_gm %>%
            filter(bact_pairing %in% bact_lake_ORW$bact_pairing) %>%
            arrange(AU, CollectionDate)
        
        # Calculate the rolling date difference using `diff()` and setting the lag accordingly.
        
        lag_data <- lake_ORW_raw %>%
            group_by(bact_pairing) %>%
            mutate(lag = c(rep(NA, 4), diff(CollectionDate, lag = 4, na.pad = T)))
        
        #convert to days
        
        #lag_data$lag_data <- lag_data$lag_data/86400
        
        #lag_data$lag_data <- gsub('secs', '', lag_data$lag_data)
        
        # Filter to the instances passing the 30-day period threshold. We can then use 
        # this data to filter the "bact_lake_ORW" object to every pairing failing the
        # threshold. 
        
        lag_filter <- lag_data %>%
            filter(lag <= 30) %>%
            distinct(bact_pairing)
        
        lake_ORW_raw <- lake_ORW_raw %>%
            filter(bact_pairing %in% lag_filter$bact_pairing)
        
        # We can then add that lag data back in with the raw data.
        
        lake_ORW_raw <- left_join(lake_ORW_raw,
                                  lag_data)
        
        # Calculate Geometric Means 
        
        # If its stupid and it works, it isn't stupid for each bact pairing we calculate 
        # the rolling geometric mean for 5 sample cycles that are taken within 30 days of
        # each other NAs are filled where a 5 day geomean cannot be calculated
        
        lake_ORW_raw <- lake_ORW_raw %>%
                group_by(bact_pairing) %>%
                mutate(window =interval(CollectionDate -30, CollectionDate),
                       samples = purrr::map_int(window, function(x) sum(.$CollectionDate %within% x)),
                       # roll_geo_og = case_when(lag <= 30 ~rollapply(ResultValue, 5, geometric.mean, fill =  NA, align = "right")), FOR QC ONLY
                       roll_geo = case_when(lag <= 30~ rollapply(ResultValue, samples, geometric.mean, fill =  NA, align = "right")))%>%
                select(-c(window, samples))
        # we're gonna want to retain the NAs as it provides valuable QC insights
        
        # summarizing
        
        # Applying same case_when logic as previous assessments, the code will summarize each sample as failing or attaining
        
        lake_orw_sum <- lake_ORW_raw %>% 
            group_by(bact_pairing) %>%
            summarise(ASSESSMENT = case_when(bact_season == "PRIMARY" & roll_geo >= 126 ~ "FAIL",
                                             bact_season == "PRIMARY" & roll_geo < 126 ~ "ATTAIN",
                                             bact_season == "Secondary" & roll_geo >= 630 ~ "FAIL",
                                             bact_season == "Secondary" & roll_geo < 630 ~ "ATTAIN"))
        
        # Here the Aus that have non attainments are collected for further processing and 
        # the next block will colate all the aus that attain by filtering out the aus that fail
        
        # Failures
        
        lake_orw_sum_fail <- lake_orw_sum %>%
            group_by(bact_pairing) %>%
            filter(ASSESSMENT == "FAIL") %>%
            distinct()
        
        # attainments
        
        lake_orw_sum_attain <- lake_orw_sum %>%
            group_by(bact_pairing) %>%
            filter(ASSESSMENT == "ATTAIN") %>%
            filter(!bact_pairing %in% lake_orw_sum_fail$bact_pairing) %>%
            distinct()
        
        #Combine and setup and identifier and there ya go
        
        bact_assessment_gm <- rbind(lake_orw_sum_fail,
                                    lake_orw_sum_attain)
        
        bact_assessment_gm$method <- "GM"
        
        return(bact_assessment_gm)

}
