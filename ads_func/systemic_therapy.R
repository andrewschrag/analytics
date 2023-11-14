#
# NOTE(ian) 04/24/2023:
# This function is used to retreive systemic therapies
# for patients, which is used to created a nested data frame
# within the ADS.
#
# This function also logs a few data quality checks, for example
# counting null records, etc.
#
get_systemic_therapies <- function(cohort_data, cohort_name) {
    message(glue::glue("* * * Systemic Therapy * * *"))

    # Pull in OpenClinica systemic therapies
    sys_tx_query = glue::glue("SELECT
                            patientid,
                            formdata::json -> 'systemic_therapy' -> 'therapy' AS sys_tx
                          FROM openclinica.formdata
                          WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cohort_name}'
                          AND formdata::json -> 'systemic_therapy' -> 'therapy' IS NOT NULL")
    meds_oc <- run_query(sys_tx_query,
                         connection = spmd) %>%
      mutate(sys_tx = purrr::map(sys_tx, fromJSON_na)) %>%
      unnest(sys_tx, keep_empty = TRUE) %>%
      rename_oc() %>%
      filter(therapyrow != "") %>% # removes empty rows of data
      # only doing this for enriched cohort and some registry patients have abstracted systemic therapies
      semi_join(cohort_data$ads %>%
                  filter(enriched_cohort_flag),
                by = "patientid")
    # Check for missing therapy names needing resolution by CTRs
    meds_oc_na_check <- meds_oc %>%
      filter(is.na(therapyname))
    if (dim(meds_oc_na_check)[1] > 0) {
      message(glue::glue("{dim(meds_oc_na_check)[1]} systemic therapy records missing therapy names - these will be set to unknown and processed as such until resolved!"))
      meds_oc_na_check %>%
        left_join(tbl(spmd, in_schema("ca", "map_spmd_participantid")) %>%
                    filter(site == cohort_name) %>%
                    collect(),
                  by = "patientid") %>%
        select(participant_id, site_name, abstract_final_values, therapyrow, therapyname, therapynameother) %>%
        print(n = 50)
    }
    # Wrangle OpenClinica data
    meds_oc <- meds_oc %>%
      map_oc(routeofadministration, "routeOfAdministration (r)") %>%
      map_oc(startreason, "startReason (r)") %>%
      map_oc(stopreason, "stopReason (r)") %>%
      map_oc(therapyname, "therapy_values") %>%
      mutate(therapyname = tolower(therapyname)) %>%
      replace_na(list(therapyname = "unknown")) %>% # these cases need to be resolved by a CTR
      convert_oc_yes_no(clinicaltrial) %>%
      # standardizes trial drug name and fixes disagreement between clinicaltrial and therapyname = trial drug
      mutate(therapyname = ifelse(grepl("trial", therapyname), "trial drug", therapyname),
             clinicaltrial = ifelse(therapyname == "trial drug", "Yes", clinicaltrial)) %>%
      convert_oc_yes_no(therapyongoing) %>%
      convert_empty_string(therapystartdate) %>%
      convert_empty_string(therapystartdateym) %>%
      convert_empty_string(therapystartdatey) %>%
      convert_empty_string(therapyenddate) %>%
      convert_empty_string(therapyenddateym) %>%
      convert_empty_string(therapyenddatey) %>%
      impute_dates_oc(keep_partial_date_field = TRUE, impute_year = FALSE) %>%
      # Recode NAs
      replace_na(list(routeofadministration = "Unknown",
                      startreason = "Unknown",
                      therapyongoing = "Unknown")) %>%
      # Stop reason must be handled differently because then therapy is ongoing, there won't be a stop reason
      # If stopreason is NA and therapy is ongoing, leave as is
      # otherwise recode NAs to Unknown (ongoing therapies will have NAs)
      mutate(stopreason = if_else(therapyongoing == "Yes",
                                  stopreason,
                                  coalesce(stopreason, "Unknown"),
                                  coalesce(stopreason, "Unknown"))) %>%
      # When therapy is ongoing we fill in with date of last abstraction since abstracted date of last contact
      # cannot be relied upon to reflect the date through which we know all data has been abstracted.
      left_join(followup_ma, by = "patientid") %>%
      mutate(therapyenddate = if_else(therapyongoing == "Yes",
                                      last_abstraction_date,
                                      therapyenddate,
                                      therapyenddate),
             therapyenddate_granularity = if_else(therapyongoing == "Yes" &
                                                    (!is.na(last_abstraction_date)),
                                                  "DAY",
                                                  therapyenddate_granularity,
                                                  therapyenddate_granularity)) %>%
      extract_year_oc(therapystartdate) %>%
      extract_year_oc(therapyenddate)
    # Check for illogical date sequences not due to imputation (start after end date)
    meds_oc_start_after_end <- meds_oc %>%
      semi_join(cohort_data$ads %>% filter(sourceschema == "ma"), by = "patientid") %>%
      filter(therapystartdate > therapyenddate) %>%
      # if they are in different years, different months, or the granularity is day for both, then imputation is
      # not the culprit
      filter(therapystartdate_year != therapyenddate_year |
               month(therapystartdate) != month(therapyenddate) |
               (therapystartdate_granularity == "DAY" & therapyenddate_granularity == "DAY")) %>%
      filter(therapyongoing != "Yes") %>%
      mutate(therapystartdate = if_else(therapystartdate_granularity == "MONTH", NA_Date_, therapystartdate),
             therapyenddate = if_else(therapyenddate_granularity == "MONTH", NA_Date_, therapyenddate))
    if (dim(meds_oc_start_after_end)[1] > 0) {
      message(glue::glue("{dim(meds_oc_start_after_end)[1]} systemic therapy records with start date after end date - these will exclude patients from LoTs until resolved!"))
      meds_oc_start_after_end %>%
        left_join(tbl(spmd, in_schema("ca", "map_spmd_participantid")) %>%
                    filter(site == cohort_name) %>%
                    collect(),
                  by = "patientid") %>%
        select(site_name, participant_id, patientid, therapyname, therapystartdate, therapystartdateym, therapyenddate, therapyenddateym) %>%
        arrange(site_name, participant_id, therapyname, therapystartdate, therapyenddate) %>%
        print(n = 50)
    }
    # Check for illogical dates of last abstraction/contact when filling therapy end date - therapies that start
    # after date of last abstraction will be removed, so these won't all cause patient exclusions.
    meds_oc_ongoing_before_start <- meds_oc %>%
      semi_join(cohort_data$ads %>% filter(sourceschema == "ma"), by = "patientid") %>%
      filter(therapystartdate > therapyenddate) %>%
      # if they are in different years, different months, or the granularity is day for both, then imputation is
      # not the culprit
      filter(therapystartdate_year != therapyenddate_year |
               month(therapystartdate) != month(therapyenddate) |
               (therapystartdate_granularity == "DAY" & therapyenddate_granularity == "DAY")) %>%
      filter(therapyongoing == "Yes")
    if (dim(meds_oc_ongoing_before_start)[1] > 0) {
      message(glue::glue("{dim(meds_oc_ongoing_before_start)[1]} systemic therapy records with therapy ongoing and start date after date of last abstraction - these records will not be included!"))
      meds_oc_ongoing_before_start %>%
        left_join(tbl(spmd, in_schema("ca", "map_spmd_participantid")) %>%
                    filter(site == cohort_name) %>%
                    select(-last_abstraction_date) %>%
                    collect(),
                  by = "patientid") %>%
        select(site_name, participant_id, patientid, therapyname,
               therapystartdate, therapystartdateym, therapyenddate, therapyenddateym,
               therapyongoing, last_abstraction_date) %>%
        arrange(site_name, participant_id, therapyname, therapystartdate, therapyenddate) %>%
        print(n = 50)
    }
    # Correct date imputations when imputed date is in the same month as un-imputed date and causes start to be
    # before the end date. Subset to columns for LoTs.
    meds_oc <- meds_oc %>%
      # we impute to the first of the month when imputed start date comes after the end date in the same month
      mutate(therapystartdate = case_when(therapystartdate_granularity != "MONTH" ~ therapystartdate,
                                          therapystartdate > therapyenddate &
                                            therapystartdate_year == therapyenddate_year &
                                            month(therapystartdate) == month(therapyenddate) ~ ymd(therapystartdateym), # by definition first of the month
                                          TRUE ~ therapystartdate),
             # we impute to the end of the month when imputed end date comes before the start date in the same month
             therapyenddate = case_when(therapyenddate_granularity != "MONTH" ~ therapyenddate,
                                        therapystartdate > therapyenddate &
                                          therapystartdate_year == therapyenddate_year &
                                          month(therapystartdate) == month(therapyenddate) ~ impute_date_last_of_month(therapyenddateym),
                                        TRUE ~ therapyenddate))

    # Variables added to the ADS
    # Systemic therapies data frame in data frame
    systemic_therapy_oc <- meds_oc %>%
      rename_with(~ str_replace(., "therapystartdate", "start_date")) %>%
      rename_with(~ str_replace(., "therapyenddate", "end_date")) %>%
      rename(start_reason = startreason,
             stop_reason = stopreason,
             therapy_ongoing_flag = therapyongoing,
             route_of_administration = routeofadministration) %>%
      select(patientid,
             therapyname,
             start_date,
             start_date_granularity,
             start_date_year,
             start_reason,
             end_date,
             end_date_granularity,
             end_date_year,
             stop_reason,
             therapy_ongoing_flag,
             clinicaltrial) %>%
      mutate(clinicaltrial = ifelse(is.na(clinicaltrial), "Unknown", clinicaltrial)) %>%
      distinct()  %>%
      semi_join(cohort_data$patients_enriched, by = "patientid")

    systemic_therapy_oc_ads <- systemic_therapy_oc %>%
      select(patientid,
             therapyname,
             start_date,
             start_date_granularity,
             start_date_year,
             start_reason,
             end_date,
             end_date_granularity,
             end_date_year,
             stop_reason,
             therapy_ongoing_flag) %>%
      nest(systemic_therapy = c(!patientid))

      #
      # NOTE(ian) 04/24/2023:
      # I'm returning both of these objects in a list, b/c
      # they are used separately in building the ADS.
      #
      # As of now, systemic_therapy_oc is used in deriving
      # Neoadjuvant chemo flags, while systemic_therapy_oc_ads
      # is used in building the final ADS. We can revisit this
      # and possibly pair it down a bit more later on.
      #
      return(list(systemic_therapy_oc = systemic_therapy_oc,
                  systemic_therapy_oc_ads = systemic_therapy_oc_ads))
}
