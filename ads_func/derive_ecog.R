# df is a dataframe containing the patients in the cohort of interest and index date (plus the index_date_year) 
# index_date is the date which defines the ECOG setting of interest 
# systemic_therapy_df is the name of the dataframe which contains cleaned systemic therapies data for the cohort of interest
## function assumes this dataframe has variable "start_date" and "end_date" for the therapy start/end date
# radiation_df is the name of the dataframe which contains cleaned radiation therapy data for the cohort of interest
## function assumes this dataframe has variable "start_date" "end_date" for the radiation start/end dates
# surgery_df is the name of the dataframe which contains cleaned primary site surgery data for the cohort of interest
## function assumes this dataframe has variable "surgery_date" for the surgery date
#
# This function [will eventually] additionally handles imputed dates to prevent exclusion of dates due to imputation.

derive_ecog <- function(df, index_date, performance_df, systemic_therapy_df, radiation_df, surgery_df){

  index_date <- enquo(index_date)
  
  # Clean ECOG data
  ecog <- performance_df %>%
    # Combine all performance statuses into a single column and all dates into a single column
    rename(performancedate = cr_datetherapyecog,
           performancestatus = cr_therapyecog) %>%
    select(patientid, performancedate, performancestatus) %>%
    bind_rows(performance_df %>%
                rename(performancedate = datetherapyecog,
                       performancestatus = therapyecog) %>%
                select(patientid, performancedate, performancestatus)) %>%
    bind_rows(performance_df %>%
                rename(performancedate = datemetstherapyecog,
                       performancestatus = therapymetsecog) %>%
                select(patientid, performancedate, performancestatus)) %>%
    # we have to have a date to determine ECOG at specific time points
    filter(!is.na(performancedate)) %>%
    map_oc(performancestatus, "pf_status_values") %>%
    # NAs will be replaced with "Unknown" anyways
    filter(performancestatus != "Unknown") %>%
    mutate(performancemethod = str_extract(performancestatus, "ECOG|KPS"),
           performancescore = str_extract(performancestatus, "\\d.*")) %>%
    # convert KPS to ECOG: https://oncologypro.esmo.org/oncology-in-practice/practice-tools/performance-scales
    mutate(ecog = case_when(is.na(performancemethod) ~ performancescore,
                            performancemethod == "ECOG" ~ performancescore,
                            suppressWarnings(as.numeric(performancescore)) >= 90 ~ "0",
                            suppressWarnings(as.numeric(performancescore)) >= 70 ~ "1",
                            suppressWarnings(as.numeric(performancescore)) >= 50 ~ "2",
                            suppressWarnings(as.numeric(performancescore)) >= 30 ~ "3",
                            suppressWarnings(as.numeric(performancescore)) >= 10 ~ "4",
                            TRUE ~ performancescore)) %>%
    rowwise() %>% # necessary for str_extract_all
    # some ECOG values are given in a range (e.g., "0-1") and we extract the highest value in the range
    mutate(ecog_date = ymd(performancedate),
           # is.na case is for SQA and empty df
           ecog = ifelse(is.na(ecog),
                         NA_real_,
                         max(as.numeric(str_extract_all(ecog, "\\d")[[1]])))) %>%
    ungroup() %>%
    mutate(ecog = as.numeric(ecog)) %>% # for SQA empty df case
    select(patientid, ecog_date, ecog)
  
  ### use already defined MEDS, SURGERY, RADIATION DATA to check ECOG date with 30 days of first treatment
  
  # First treatment post-dx for ECOG determination
  first_treatment_post_index <- bind_rows(systemic_therapy_df$after%>% 
                                            select(patientid, contains("start_date"), contains("end_date")) %>% 
                                            rename_with(~ str_remove_all(., "therapy_")),
                                          radiation_df %>% 
                                            select(patientid, contains("start_date"), contains("end_date")) %>% 
                                            rename_with(~ str_remove_all(., "radiation_")),
                                          surgery_df %>% 
                                            rename(start_date = surgery_date,
                                                   start_date_granularity = surgery_date_granularity)%>% 
                                            select(patientid, contains("start_date"))
                                          ) %>%
    mutate(start_date = pmin(start_date, end_date, na.rm = TRUE)) %>%
    inner_join(df %>%
                 select(patientid, !!index_date),
               by = "patientid") %>%
    filter(start_date >= !!index_date) %>%
    group_by(patientid) %>%
    summarize(first_tx_post_index_date = min(start_date, na.rm = TRUE)) %>% 
    ungroup()
  
  ## EXPANSION
  # check against imputed dates to 
  # IF ecog_date <= first therapy date AND ecog_date >= first therapy date - days(30), THEN INCLUDE
  # IF first therapy date granularity == "MONTH" AND month(ecog_date) == month(first therapy date) AND 
  # year(ecog_date) == year(first therapy date), THEN INCLUDE
  # ELSE EXCLUDE
  
  # ECOG after index and in 30 days prior to first treatment
  ecog %>%
    inner_join(df %>%
                 select(patientid, !!index_date),
               by = "patientid") %>%
    inner_join(first_treatment_post_index,
               by = "patientid") %>%
    # after diagnosis
    filter(ecog_date >= !!index_date) %>%
    # in the 30 day window prior to first treatment
    filter(ecog_date <= first_tx_post_index_date &
             ecog_date >= first_tx_post_index_date - days(30)) %>%
    mutate(distance = day_difference(ecog_date, first_tx_post_index_date)) %>%
    arrange(patientid, ecog_date) %>%
    group_by(patientid) %>%
    # scores closest to the first_tx_post_index_date
    filter(distance == min(distance)) %>%
    # if ties remain, take the highest ECOG
    slice_max(ecog) %>%
    ungroup() %>%
    select(patientid, ecog, ecog_date) %>%
    # this shouldn't be possible (ECOG = 5 = Dead)
    mutate(ecog = ifelse(ecog == 5, "Unknown", as.character(ecog)),
           ecog = as.character(ecog) # for SQA emty df case
    ) %>% 
    distinct()
  
}
