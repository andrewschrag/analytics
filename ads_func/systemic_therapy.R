# get_systemic_therapy ----
get_systemic_therapy_oc <- function(patients, cohort_name, spmd = spmd_con()){
  run_query(glue::glue("SELECT
                        patientid,
                        participantid,
                        formdata::json -> 'systemic_therapy' -> 'therapy' AS sys_tx
                      FROM openclinica.formdata
                      WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cohort_name}'
                      AND formdata::json -> 'systemic_therapy' -> 'therapy' IS NOT NULL"),
            connection = spmd) %>%
    mutate(sys_tx = purrr::map(sys_tx, fromJSON_na)) %>%
    unnest(sys_tx, keep_empty = TRUE) %>%
    rename_oc() %>%
    # only doing this for enriched cohort and some registry patients have abstracted systemic therapies
    semi_join(patients,
              by = c("patientid", "participantid"))
}

# clean_systemic_therapy ----
clean_systemic_therapy_oc <- function(.data, cohort_name, spmd = spmd_con()){

  meds_oc <- .data %>%
    filter(!therapyname %in% c("NULL", "")) %>%
    {if (cohort_name %in% c("lung", "breast")) mutate(., therapyname = coalesce(therapyname, therapynameother)) else .}

  # Check for missing therapy names needing resolution by CTRs
  meds_oc_na_check <- meds_oc %>%
    filter(is.na(therapyname),
           (startdate != "" | enddate != "" | therapyongoing!= "" | clinicaltrial != "" | therapyrow > "1"))
  if (dim(meds_oc_na_check)[1] > 0) {
    message(glue::glue("{dim(meds_oc_na_check)[1]} systemic therapy records missing therapy names - these will be set to unknown and processed as such until resolved!"))
    meds_oc_na_check %>%
      left_join(tbl(spmd, in_schema("ca", "map_spmd_participantid")) %>%
                  filter(site == cohort_name) %>%
                  collect(),
                by = "patientid") %>%
      select(participant_id, site_name, abstract_final_values, therapyrow, therapyname, therapynameother, startdate, enddate,
             therapyongoing, clinicaltrial) %>%
      print(n = 50)
  }

  meds_oc %>%
    filter(!is.na(therapyname)) %>%
    map_oc(routeofadministration, "routeOfAdministration (r)") %>%
    map_oc(startreason, "startReason (r)") %>%
    map_oc(stopreason, "stopReason (r)") %>%
    map_oc(startsource, "startSource (r)") %>%
    map_oc(stopsource, "stopSource (r)") %>%
    map_oc(therapyname, "therapy_values") %>%
    # mutate(therapyname = if_else(therapyname=='x',
    #                              tolower(str_remove_all(therapynameother, pattern = '\\"')),
    #                              tolower(get_oc_values(therapyname, "therapy_values")))) %>%
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
    # Stop reason must be handled differently because when therapy is ongoing, there won't be a stop reason
    # If therapy is ongoing, set to NA (almost all are NA, but there are a couple edge cases)
    # otherwise recode NAs to Unknown
    mutate(stopreason = if_else(therapyongoing == "Yes",
                                NA_character_,
                                coalesce(stopreason, "Unknown"),
                                coalesce(stopreason, "Unknown"))) %>%
    # Stop source is also handled differently because we have Medrio migrated data that differs from the OC-native
    # data. In Medrio, it was possible to have a stopreason = Death/Unknown and still have a stopsource. This
    # option is now hidden for stopreason = Death/Unknown patients in OC so these will be NA in OC. So, it is more
    # appropriate to leave stopsource = NA if it is missing and stopreason = Death/Unknown. All other cases
    # should be back-filled with Unknown.
    mutate(stopsource = na_if(stopsource, "NULL")) %>%
    mutate(stopsource = if_else(stopreason %in% c("Death", "Unknown"),
                                stopsource, # whether NA or not
                                coalesce(stopsource, "Unknown"),
                                coalesce(stopsource, "Unknown"))) %>%
    # start source is also handled differently because startsource is only collected if we have startreason =
    # Yes. Thus, for startreason = No/Unknown leave as NA and for startreason = Yes back-fill with Unknown.
    mutate(startsource = if_else(startreason == "Yes",
                                 coalesce(startsource, "Unknown"),
                                 NA_character_,
                                 NA_character_))

}

# wrangle_systemic_therapy
wrangle_systemic_therapy <- function(.data, cohort_name, followup, morphology, spmd = spmd_con()){

  meds_oc <- .data %>%
    # When therapy is ongoing we fill in with date of last abstraction since abstracted date of last contact
    # cannot be relied upon to reflect the date through which we know all data has been abstracted.
    left_join(followup, by = c("patientid", "participantid")) %>%
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

  # Correct date imputations when imputed date is in the same month as un-imputed diagnosis date and causes start to be
  # before the diagnosis date.
  meds_oc <- meds_oc %>%
    left_join(morphology %>% select(patientid, diagnosis_date, diagnosis_date_year)) %>%
    # we impute to the diagnosis date when imputed start date comes after diagnosis date in the same month
    mutate(therapystartdate = case_when(therapystartdate_granularity != "MONTH" ~ therapystartdate,
                                        therapyenddate_granularity == "DAY" &
                                          therapyenddate <= diagnosis_date ~ therapystartdate,
                                        therapystartdate_year == diagnosis_date_year &
                                          month(therapystartdate) == month(diagnosis_date) ~ diagnosis_date,
                                        TRUE ~ therapystartdate),
           # we impute to the end of the month when imputed end date comes before diagnosis date in the same month
           therapyenddate = case_when(therapyenddate_granularity != "MONTH" ~ therapyenddate,
                                      diagnosis_date > therapyenddate &
                                        therapyenddate_year == diagnosis_date_year &
                                        month(therapyenddate) == month(diagnosis_date) ~ impute_date_last_of_month(therapyenddateym),
                                      TRUE ~ therapyenddate))

  # Correct date imputations when imputed date is in the same month as abstraction date and causes start to be
  # after the abstraction date.
  meds_oc <- meds_oc %>%
    # we impute to the last abstraction date when imputed start date comes after the last abstraction date
    mutate(therapystartdate = case_when(therapystartdate_granularity != "MONTH" ~ therapystartdate,
                                        therapystartdate > last_abstraction_date &
                                          therapystartdate_year == year(last_abstraction_date) &
                                          month(therapystartdate) == month(diagnosis_date) ~ last_abstraction_date,
                                        TRUE ~ therapystartdate),
           # we impute to the last abstraction date when imputed end date comes after the last abstracted date
           therapyenddate = case_when(therapyenddate_granularity != "MONTH" ~ therapyenddate,
                                      last_abstraction_date < therapyenddate &
                                        therapyenddate_year == year(last_abstraction_date) &
                                        month(therapyenddate) == month(last_abstraction_date) ~ last_abstraction_date,
                                      TRUE ~ therapyenddate))

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

  # Filter out therapies that start after the date of last abstraction - we are not confident that dates
  # beyond date of last abstraction indicate the patient was fully refreshed and hence will ignore these
  # until the patient is refreshed and the date of last abstraction is updated. Do this after date imputations
  # are corrected in case that resolves any of these problems.
  meds_after_abstraction <- meds_oc %>%
    filter(therapystartdate > last_abstraction_date)
  meds_oc <- meds_oc %>%
    anti_join(meds_after_abstraction, colnames(meds_oc))

  # Variables added to the ADS
  #### Systemic therapies data frame in data frame ----
  systemic_therapy_oc <- meds_oc %>%
    rename_with(~ str_replace(., "therapystartdate", "start_date")) %>%
    rename_with(~ str_replace(., "therapyenddate", "end_date")) %>%
    select(patientid,
           therapyname,
           start_date,
           start_date_granularity,
           start_date_year,
           startreason,
           startsource,
           end_date,
           end_date_granularity,
           end_date_year,
           stopreason,
           stopsource,
           therapyongoing,
           clinicaltrial,
           route_of_administration = routeofadministration # 100% NA as of 5/24
    ) %>%
    mutate(clinicaltrial = ifelse(is.na(clinicaltrial), "Unknown", clinicaltrial),
           ## fix some obvious spelling mistakes/ issues
           therapyname = case_when(therapyname %in% c("sacituzimab","sacitizumab/trodelvy",
                                                      "sacituzumab govitecan-hziy",
                                                      "sacituzumab") ~ "sacituzumab govitecan",
                                   therapyname %in% c("avelumab ")~ "avelumab",
                                   therapyname %in% c("utomilumab ")~ "utomilumab",
                                   therapyname %in% c("velaparib") ~ "veliparib",
                                   therapyname %in% c("trastuzumab-pkrb", "trastuzumab-anns") ~ "trastuzumab",
                                   therapyname %in% c("keytruda/pembrolizumab")~ "pembrolizumab",
                                   therapyname %in% c("nab-paclitaxel") ~ "paclitaxel",
                                   therapyname %in% c("leuprolide acetate") ~ "leuprolide",
                                   therapyname %in% c("doxorubuicin liposome injection", "doxorubicin liposome injection") ~ "doxorubicin liposome",
                                   therapyname %in% c("trastuzumab and hyaluronidase-oysk") ~ "trastuzumab/hyaluronidase-oysk",
                                   therapyname == "11202" ~ "vincristine",
                                   TRUE ~ as.character(therapyname))
    ) %>%
    distinct()

  cohort_name_sym <- sym(cohort_name)
  class_sym <- sym(glue::glue("{cohort_name}_class"))
  subclass_sym <- sym(glue::glue("{cohort_name}_subclass"))

  therapy_class_subclass_map <- tbl(spmd, in_schema("ca", "map_therapy_class_subclass")) %>% collect()

  systemic_therapy_oc <- systemic_therapy_oc %>%
    left_join(therapy_class_subclass_map %>%
                filter({{cohort_name_sym}} == "X") %>%
                select(therapyname, class = {{class_sym}}, subclass = {{subclass_sym}}) %>%
                mutate(therapyname = tolower(therapyname)),
              by = "therapyname")

  class_not_mapped <- systemic_therapy_oc %>%
    filter(is.na(class) | is.na(subclass)) %>%
    distinct(therapyname)

  if(class_not_mapped %>% nrow() > 0){
    rlang::warn(glue::glue('{class_not_mapped %>% nrow()} therapies do not have class or subclass mapped.\nUpdate "Therapy Mapping" tab in latest data dictionary googlesheet and recreate ca.map_therapy_class_subclass'))
    class_not_mapped %>% print(n=Inf)
  }

  systemic_therapy_oc %>%
    select(patientid,
           therapy_name = therapyname,
           therapy_start_date = start_date,
           therapy_start_date_granularity = start_date_granularity,
           therapy_start_date_year = start_date_year,
           therapy_start_reason_progression = startreason,
           therapy_start_reason_progression_source = startsource,
           therapy_end_date = end_date,
           therapy_end_date_granularity = end_date_granularity,
           therapy_end_date_year = end_date_year,
           therapy_stop_reason = stopreason,
           therapy_status_source = stopsource,
           therapy_ongoing_flag = therapyongoing,
           therapy_clinical_trial = clinicaltrial,
           therapy_class = class,
           therapy_subclass = subclass,
           therapy_route_of_administration = route_of_administration
    ) %>%
    {if (cohort_name %in% c("lung", "breast", "bladder")) select(.,-therapy_route_of_administration) else .}

}

# derive_systemic_therapy ----
derive_systemic_therapy <- function(patients, spmd = spmd_con(), cohort_name, followup, morphology, nest = TRUE){

  get_systemic_therapy_oc(patients = patients, cohort_name = cohort_name, spmd = spmd) %>%
    clean_systemic_therapy_oc(., cohort_name = cohort_name, spmd = spmd) %>%
    wrangle_systemic_therapy(., cohort_name = cohort_name, followup = followup, morphology = morphology, spmd = spmd) %>%
    {if (nest) nest(., systemic_therapy = c(!patientid)) else .}

}
