# Retrieve OpenClinica's recurrence data
get_recurrence_oc <- function(cancer, patients = NULL, spmd = spmd_con()) {
  
  if (is.null(cancer)) rlang::abort("cancer is missing")
  
  # check to make sure we can filter on patients by at least one of patientid or participantid
  patient_cols <- check_patient_arg(patients, 
                                    oc_warning = TRUE)
  
  run_query(glue::glue("SELECT
                          patientid,
                          participantid,
                          formdata::json -> 'recurrence' -> 'group' ->> 'ma_recurrence' AS ma_recurrence, 
                          formdata::json -> 'recurrence' -> 'recurrence' AS recurrence
                        FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}' AND
                        formdata::json -> 'recurrence' IS NOT NULL"),
            connection = spmd) %>%
    {if (!is.null(patients)) semi_join(., patients, by = patient_cols) else .} %>% 
    mutate(recurrence = purrr::map(recurrence, fromJSON_na)) %>%
    unnest(recurrence, keep_empty = TRUE) %>%
    rename_oc() %>%
    map_oc(recurrence, "recurrence") %>% 
    map_oc(recurrenceconfirmtype, "recurrenceConfirmType (r)") %>%
    map_oc(recurrencesite, "recurrenceSite (r)") %>%  
    convert_empty_string(recurrencedate) %>%
    convert_empty_string(recurrencedateym) %>%
    convert_empty_string(recurrencedatey) %>%
    select(-recurdate) %>%
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>%
    extract_year_oc(recurrencedate) %>% 
    select(-recurrencedateym, -recurrencedatey)
}

# Clean-up OC recurrence data
clean_recurrence_oc <- function(.data, patients, diagnosis = NULL, diagnosis_date = diagnosis_date, silent = FALSE) {
  
  if (is.null(patients)) rlang::abort("patients is missing")
  patients %>% 
    check_for_column(patientid, 
                     return_df = FALSE) %>% 
    invisible()
  
  diagnosis_date <- enquo(diagnosis_date)
  if (!is.null(diagnosis)) {
    diagnosis %>% 
      check_for_column(patientid, 
                       return_df = FALSE)
    diagnosis %>% 
      check_for_column(!!diagnosis_date, 
                       .class = "Date", 
                       return_df = FALSE)
  } else if (!silent) {
    rlang::warn("diagnosis data frame not supplied, including all recurrence data")
  }
  
  .data <- .data %>% 
    group_by(patientid) %>% 
    mutate(recurrence_flag = case_when(any(recurrence == "Yes") ~ "Yes",
                                       any(!is.na(recurrencesite)) ~ "Yes",
                                       any(recurrence == "No") ~ "No",
                                       all(is.na(recurrence)) ~ "Unknown",
                                       TRUE ~ "Unknown")) %>% 
    ungroup() %>% 
    # note: set reccurencedate_granularity to NA as well? extend to recurrence_flag == Unknown?
    # argument against doing it for Unknowns is that we can differentiate known Unknowns (some stuff non-NA)
    # from unknown Unknowns (all NA)
    # mutate(across(c(recurrencedate, recurrencedate_year, recurrencesite, recurrenceconfirmtype), ~ if_else(recurrence_flag == "No", NA, .))) %>% # doesn't work in airflow, not sure why
    mutate(recurrencedate = if_else(recurrence_flag == "No", as.Date(NA), recurrencedate),
           recurrencedate_year = if_else(recurrence_flag == "No", NA_real_, recurrencedate_year),
           across(c(recurrencesite, recurrenceconfirmtype), ~ if_else(recurrence_flag == "No", NA_character_, .))) %>% 
    mutate(across(c(recurrencesite, recurrenceconfirmtype), ~ if_else(recurrence_flag == "Yes", coalesce(., "Unknown"), .))) %>% 
    mutate(recurrencesite = ifelse(recurrencesite == "Unknown Depth of Invasion", "Unknown", recurrencesite))
  
  if (!is.null(diagnosis)) {
    .data <- .data %>% 
      left_join(diagnosis %>% 
                  select(patientid, !!diagnosis_date), 
                by = "patientid") %>% 
      filter(recurrence_flag %in% c("No", "Unknown") | 
               recurrencedate > !!diagnosis_date | 
               recurrencedate_granularity == "NONE" | 
               (recurrencedate_granularity == "YEAR" & recurrencedate_year >= year(!!diagnosis_date)))
  }
  
  .data <- .data %>% 
    join_replace_na(patients, recurrence_flag)
  
  if (!silent) {
    .data %>% 
      filter(recurrence_flag == "Unknown") %>% 
      count() %>% 
      pull() %>% 
      {if (. > 0) rlang::warn(glue::glue("recurrence_flag = 'Unknown' detected for {.} patients")) else invisible(.)}
  }
  
  .data %>% 
    select(patientid, 
           recurrence_flag,
           recurrence_site = recurrencesite, 
           recurrence_date = recurrencedate, 
           recurrence_date_year = recurrencedate_year,
           recurrence_date_granularity = recurrencedate_granularity, 
           recurrence_confirm_type = recurrenceconfirmtype) %>% 
    distinct()
}

# Retrieve and clean-up OC's recurrence data
get_recurrence <- function(cancer, patients, diagnosis = NULL, diagnosis_date = diagnosis_date, spmd = spmd_con(), silent = FALSE) {
  diagnosis_date <- enquo(diagnosis_date)
  
  if (is.null(cancer)) rlang::abort("cancer is missing")
  if (is.null(patients)) rlang::abort("patients is missing")
  
  get_recurrence_oc(cancer = cancer, 
                    patients = patients, 
                    spmd = spmd) %>% 
    clean_recurrence_oc(patients = patients, 
                        diagnosis = diagnosis,
                        diagnosis_date = !!diagnosis_date,
                        silent = silent) %>% 
    nest(recurrence = -c(patientid, recurrence_flag))
}

# Get OC's metastasis data
get_metastasis_oc <- function(cancer, patients = NULL, spmd = spmd_con()) {
  
  if (is.null(cancer)) rlang::abort("cancer is missing")
  
  # check to make sure we can filter on patients by at least one of patientid or participantid
  patient_cols <- check_patient_arg(patients, 
                                    oc_warning = TRUE)
  
  run_query(glue::glue("SELECT
                         patientid,
                         participantid,
                         formdata::json -> 'staging_performance_status_metastatic_status' ->> 'mets' AS mets
                        FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}' AND
                        formdata::json -> 'staging_performance_status_metastatic_status' ->> 'mets' IS NOT NULL"),
            connection = spmd) %>%
    {if (!is.null(patients)) semi_join(., patients, by = patient_cols) else .} %>% 
    mutate(mets = purrr::map(mets, jsonlite::fromJSON)) %>%
    unnest_wider(mets) %>%
    rename_oc() %>% 
    map_oc(hasmetastaticdisease, "hasMetastaticDisease") %>% 
    convert_oc_yes_no(metsbrain) %>% 
    convert_oc_yes_no(metsbone) %>% 
    convert_oc_yes_no(metsliver) %>% 
    convert_oc_yes_no(metslung) %>% 
    convert_oc_yes_no(metsdistantln) %>% 
    convert_oc_yes_no(metsother) %>% 
    convert_empty_string(dateoffirstmets) %>%
    convert_empty_string(dateoffirstmetsy) %>%
    convert_empty_string(dateoffirstmetsym) %>%
    select(-metsdate) %>%
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>%
    extract_year_oc(dateoffirstmets) %>% 
    select(-dateoffirstmetsy, -dateoffirstmetsym)
}

# Clean up OC's metastasis data
clean_metastasis_oc <- function(.data, patients, 
                                diagnosis, diagnosis_date = diagnosis_date, 
                                staging, stage_group = prioritized_stage_group_dx,
                                recurrence, 
                                demos, deceased_date = deceased_date, 
                                birth_date = birth_date, birth_date_deid = birth_date_deid,
                                silent = FALSE) {
  
  if (is.null(patients)) rlang::abort("patients is missing")
  if (is.null(diagnosis)) rlang::abort("diagnosis is missing")
  if (is.null(staging)) rlang::abort("staging is missing")
  if (is.null(recurrence)) rlang::abort("recurrence is missing")
  if (is.null(demos)) rlang::abort("demos is missing")
  
  # column defs
  diagnosis_date <- enquo(diagnosis_date)
  stage_group <- enquo(stage_group)
  deceased_date <- enquo(deceased_date)
  birth_date <- enquo(birth_date)
  birth_date_deid <- enquo(birth_date_deid)
  
  # check inputs
  diagnosis %>% 
    check_for_column(!!diagnosis_date, "Date") %>% 
    check_for_column(patientid,
                     return_df = FALSE) %>% 
    invisible()
  staging %>% 
    check_for_column(!!stage_group, "character") %>% 
    check_for_column(patientid,
                     return_df = FALSE) %>% 
    invisible()
  recurrence %>% 
    check_for_column(recurrence, "list") %>% 
    check_for_column(patientid,
                     return_df = FALSE) %>% 
    invisible()
  demos %>% 
    check_for_column(!!deceased_date, "Date") %>% 
    check_for_column(!!birth_date, "Date") %>% 
    check_for_column(!!birth_date_deid, "Date") %>% 
    check_for_column(patientid,
                     return_df = FALSE) %>% 
    invisible()
  
  # identify first distant recurrence
  first_distant_recurrence <- recurrence %>% 
    unnest(recurrence) %>% 
    filter(recurrence_site == "Distant") %>% 
    group_by(patientid) %>% 
    arrange(recurrence_date, 
            recurrence_date_year) %>% 
    slice_head() %>% 
    ungroup() %>% 
    distinct(patientid, recurrence_site, recurrence_date, recurrence_date_year, recurrence_date_granularity)
  
  # attach all necessary info
  .data <- .data %>% 
    left_join(diagnosis %>% 
                select(patientid, !!diagnosis_date),
              by = "patientid") %>% 
    left_join(staging %>% 
                select(patientid, !!stage_group), 
              by = "patientid") %>% 
    left_join(demos %>% 
                select(patientid, !!deceased_date, !!birth_date, !!birth_date_deid),
              by = "patientid") %>% 
    left_join(first_distant_recurrence,
              by = "patientid")
  
  # If patient is diagnosed stage IV, then metastasis occurred at diagnosis
  # Else if patient abstracted with mets, then metastasis occurs then
  # Else if patient abstracted with distant recurrence, then metastasis occurs then
  # Else, not metastatic
  .data <- .data %>% 
    mutate(metastasis_date = case_when(grepl("IV", !!stage_group) ~ !!diagnosis_date,
                                       hasmetastaticdisease == "Yes" ~ dateoffirstmets,
                                       recurrence_site == "Distant" ~ recurrence_date,
                                       TRUE ~ NA_Date_),
           metastasis_date_year = case_when(grepl("IV", !!stage_group) ~ year(!!diagnosis_date),
                                            hasmetastaticdisease == "Yes" ~ dateoffirstmets_year,
                                            recurrence_site == "Distant" ~ recurrence_date_year,
                                            TRUE ~ NA_real_),
           metastasis_date_granularity = case_when(grepl("IV", !!stage_group) ~ "DAY",
                                                   hasmetastaticdisease == "Yes" ~ dateoffirstmets_granularity,
                                                   recurrence_site == "Distant" ~ recurrence_date_granularity,
                                                   TRUE ~ NA_character_),
           metastasis_presentation = case_when(grepl("IV", !!stage_group) ~ "De novo",
                                               hasmetastaticdisease == "Yes" ~ "Recurrent",
                                               recurrence_site == "Distant" ~ "Recurrent",
                                               TRUE ~ NA_character_),
           metastasis_flag = if_else(metastasis_presentation %in% c("De novo", "Recurrent"), 
                                     "Yes", 
                                     "No", 
                                     "No"))
  
  # Warn about edge cases (send for CTR review)
  if (!silent) {
    .data %>% 
      filter(hasmetastaticdisease == "No" & recurrence_site == "Distant") %>% 
      count() %>% 
      pull() %>% 
      {if (. > 0) rlang::warn(glue::glue("{.} patients detected with OC metastasis flag set to 'No', but with a distant recurrence (distant recurrence implies metastasis)")) else invisible(.)}
    
    .data %>% 
      filter(hasmetastaticdisease == "Yes" & is.na(!!deceased_date) & (is.na(metsbone) | is.na(metslung) | is.na(metsbrain) | is.na(metsother) | is.na(metsdistantln))) %>% 
      count() %>% 
      pull() %>% 
      {if (. > 0) rlang::warn(glue::glue("{.} alive patients detected with OC metastasis flag set to 'Yes' and one or more metastatic sites are NA (metastasis sites should be Yes/No for patients that are alive)")) else invisible(.)}
    
  }
  
  # Backfill the metastatic site variables when appropriate
  .data <- .data %>% 
                                                   # if not metastatic, leave as NA
    mutate(across(starts_with("mets"), ~ case_when(metastasis_flag != "Yes" ~ NA_character_,
                                                   # if abstraction didn't find metastasis then they are a Stage IV or distant 
                                                   # recurrence patient and so we will leave specific sites as unknown
                                                   hasmetastaticdisease != "Yes" ~ "Unknown",
                                                   # if abstraction found a yes/no for a site, keep that
                                                   !is.na(.) ~ .,
                                                   # if alive, leave NAs as-is (warning above)
                                                   is.na(!!deceased_date) ~ .,
                                                   # if deceased and all sites are NA, set to Unknown
                                                   is.na(coalesce(metsbone, metslung, metsbrain, metsliver, metsother, metsdistantln)) ~ "Unknown",
                                                   # if deceased and at least one site not NA, set to No
                                                   TRUE ~ "No")))
  
  # Age at metastasis
  .data <- .data %>% 
    mutate(age_metastasis = age_at(metastasis_date, !!birth_date, truncate = TRUE),
           age_metastasis_deid = age_at(metastasis_date, !!birth_date_deid, truncate = TRUE))
  
  # Subset and return relevant
  .data %>% 
    rename_with(.cols = starts_with("mets"), .fn = ~str_replace(., "mets", "mets_")) %>% # renames mets
    select(patientid, metastasis_flag, metastasis_presentation, starts_with("metastasis_date"), starts_with("mets"), age_metastasis, age_metastasis_deid) %>% 
    join_replace_na(patients, metastasis_flag)
  
}

get_metastasis <- function(cancer, 
                           patients, 
                           diagnosis, diagnosis_date = diagnosis_date, 
                           staging, stage_group = prioritized_stage_group_dx,
                           recurrence, 
                           demos, deceased_date = deceased_date,
                           birth_date = birth_date, birth_date_deid = birth_date_deid,
                           spmd = spmd_con(),
                           silent = FALSE) {
  
  if (is.null(cancer)) rlang::abort("cancer is missing")
  if (is.null(patients)) rlang::abort("patients is missing")
  if (is.null(diagnosis)) rlang::abort("diagnosis is missing")
  if (is.null(staging)) rlang::abort("staging is missing")
  if (is.null(recurrence)) rlang::abort("recurrence is missing")
  if (is.null(demos)) rlang::abort("demos is missing")
  
  diagnosis_date <- enquo(diagnosis_date)
  stage_group <- enquo(stage_group)
  deceased_date <- enquo(deceased_date)
  birth_date <- enquo(birth_date)
  birth_date_deid <- enquo(birth_date_deid)
  
  get_metastasis_oc(cancer = cancer, 
                    patients = patients, 
                    spmd = spmd) %>% 
    clean_metastasis_oc(patients = patients,
                        diagnosis = diagnosis, 
                        diagnosis_date = !!diagnosis_date,
                        staging = staging,
                        stage_group = !!stage_group,
                        recurrence = recurrence,
                        demos = demos,
                        deceased_date = !!deceased_date,
                        birth_date = !!birth_date, 
                        birth_date_deid = !!birth_date_deid,
                        silent = silent)
}
