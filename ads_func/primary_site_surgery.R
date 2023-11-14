# Pipeline functions ----

# Function retrieves primary site surgery data from openclinica.formdata for a given cancer type and patient
# population (optional).
get_primary_site_surgery_oc <- function(cancer, patients = NULL, spmd = spmd_con()) {

  if (is.null(cancer)) rlang::abort("cancer is missing")

  # check to make sure we can filter on patients by at least one of patientid or participantid
  patient_cols <- check_patient_arg(patients,
                                    oc_warning = TRUE)

  # query openclinica.formdata to get data
  run_query(glue::glue("SELECT
                          patientid,
                          participantid,
                          formdata::json -> 'primary_site_surgery' -> 'group' ->> 'ma_reasonfornosurgery' AS ma_reasonfornosurgery,
                          formdata::json -> 'primary_site_surgery' -> 'surgery' AS surgery
                          FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}'
                        AND formdata::json -> 'primary_site_surgery' IS NOT NULL"),
            connection = spmd) %>%
    {if (!is.null(patients)) semi_join(., patients, by = patient_cols) else .} %>%
    mutate(surgery = purrr::map(surgery, fromJSON_na)) %>%
    unnest(surgery, keep_empty = TRUE) %>%
    rename_oc() %>% 
    rename(rxsummsurgicalmargins = summsurgicalmargins)
}

# Function retrieves the primary site surgery data from ca.registry_naaccr for a given cohort query and
# cancer type.
get_primary_site_surgery_registry <- function(cohort_query, cancer_type, cols = c()) {
  get_registry(cohort_query,
               cancer_type,
               c("reasonfornosurgery",
                 "rxdatedxstgproc",
                 # "rxdatesurgery", # does not actually correspond to rxsummsurgprimsite
                 "rxdatemostdefinsurg", # corresponds to rxsummsurgprimsite and used for definitive surgery
                 "rxsummsurgprimsite",
                 "rxsummsurgicalmargins",
                 "lymphvascularinvasion",
                 "rxsummreglnexamined",
                 "rxsummscopereglnsur",
                 "residualtumvolpostcytoreduction", 
                 "datecasecompleted",
                 "datecasecompletedcoc",
                 "dateoflastcontact",
                 "classofcase",
                 "registryid",
                 "npiphysicianprimarysurg",
                 cols))
}

# Cleans up primary site surgery data from OpenClinica using general and sites-specific rules.
clean_primary_site_surgery_oc <- function(.data, cancer = NULL, spmd = spmd_con()) {
  # pan-ads initial cleaning
  .data <- .data %>%
    mutate(across(where(is.character), ~ na_if(., ""))) %>%
    # filter to rows containing data (row 1 is not enumerated per patient)
    # all other rows with NA surgrow are assumed to be errant entries
    group_by(patientid) %>%
    filter(row_number() == 1 |
             !is.na(surgrow)) %>%
    ungroup() %>%
    impute_dates(impute_year = FALSE,
                 keep_partial_date_field = TRUE) %>%
    extract_year_oc(rxdatesurgery) %>%
    select(-rxdatesurgeryym, -rxdatesurgeryy) %>%
    # pre-mapping code corrections
    mutate(
      rxsummsurgprimsite = ifelse(rxsummsurgprimsite == "00",
                                  "0",
                                  rxsummsurgprimsite),
      rxsummsurgprimsite = na_if(rxsummsurgprimsite, "NULL"),
      # retain coded version
      rxsummsurgprimsite_code = rxsummsurgprimsite,
      rxsummreglnexamined = ifelse(rxsummreglnexamined == "unknown",
                                   "99",
                                   rxsummreglnexamined),
      # using the NAACCR map so use the NAACCR unknown code
      lymphvascularinvasion = ifelse(lymphvascularinvasion == "unk",
                                     "9",
                                     lymphvascularinvasion),
      rxsummsurgicalmargins = ifelse(rxsummsurgicalmargins == "unk",
                                     "9",
                                     rxsummsurgicalmargins),
      # calculate before mapping for simplicity of using the codes
      gross_residual_disease = case_when(rxsummsurgicalmargins == "3" ~ "Yes",
                                         rxsummsurgicalmargins %in% c("0", "1", "2") ~ "No",
                                         rxsummsurgicalmargins %in% c("7", "8", "9", "unk") ~ "Unknown",
                                         is.na(rxsummsurgicalmargins) ~ "Unknown",
                                         TRUE ~ NA_character_))
  # cancer-specific cleaning
  if (!is.null(cancer)) {
    if (cancer == "bladder") {
      .data <- .data %>%
        clean_primary_site_surgery_oc_bladder()
    } else if (cancer == "breast") {
      .data <- .data %>%
        mutate(
          # these are likely an openclinica ghost - in all cases, surgeries lacking a 'br' were removed by CTRs from reg prepop
          # but we are still seeing them in the openclinica schema
          # these are set to NA and will be back-filled post-alignment based on reason for no surgery
          rxsummsurgprimsite = if_else(!grepl("br", rxsummsurgprimsite),
                                       NA_character_,
                                       rxsummsurgprimsite,
                                       rxsummsurgprimsite),
        )
    } else if (cancer == "lung") {
      # no lung rules at the moment
    } else if (cancer == "ovarian") {
      .data <- .data %>%
        # Ronda looked into these patients and believes these codes are in error. She and Betsy recommend
        # filtering them out.
        filter(!grepl("20ov|40ov", rxsummsurgprimsite)) %>% 
        # map residual tumor volume for ovarian
        map_residual_tumor_volume()
    } else {
      rlang::warn(glue::glue("cancer = '{cancer}' not a recognized option, no further cleaning performed"))
    }
  } else {
    rlang::warn(glue::glue("cancer argument not provided, no cancer-specific cleaning performed"))
  }

  # pan-ads mappings
  # do not need to map residualtumvolpostcytoreduction b/c always empty
  .data %>%
    map_oc(reasonfornosurgery, "reasonForNoSurgery") %>%
    map_oc(rxsummsurgprimsite, "rxSummSurgPrimSite (r)") %>%
    coalesce_map_naaccr("lymphvascularinvasion", spmd) %>% # prefer NAACR map
    coalesce_map_naaccr("rxsummsurgicalmargins", spmd) %>% # prefer NAACR map
    map_oc(rxsummscopereglnsur, "rxSummScopeRegLnSur (r)") %>%
    # consolidate unknown categories of surgical margins and scope of regional LN surgery variables
    # there is a minor casing difference between the OC and NAACCR mappings for the unknow category in these
    # vars
    mutate(rxsummscopereglnsur = ifelse(rxsummscopereglnsur == "Unknown or Not applicable", 
                                        "Unknown or not applicable", 
                                        rxsummscopereglnsur)) %>% 
    # handling no medrio and NA cases of rxsummsurgprimsite
    mutate(
      rxsummsurgprimsite = ifelse(reasonfornosurgery == "Marked No in Medrio",
                                  "00 - No surgery",
                                  rxsummsurgprimsite),
      ## we also account for occasions of lymph node procedures which are not tied to a primary site surgery
      ## we don't want to remove these rows so assign them as "99 - Unknown if surgery performed"
      rxsummsurgprimsite = case_when(is.na(rxsummsurgprimsite) &
                                       rxsummscopereglnsur != "Unknown or Not applicable" ~ "99 - Unknown if surgery performed",
                                     is.na(rxsummsurgprimsite) &
                                       !is.na(rxsummreglnexamined) ~ "99 - Unknown if surgery performed",
                                     TRUE ~ rxsummsurgprimsite))

}

# Cleans up primary site surgery data from registry according to general and site-specific rules.
clean_primary_site_surgery_registry <- function(.data, cancer, spmd = spmd_con()) {
  
  cancer_code_lookup <- list("bladder" = "bl",
                          "breast" = "br",
                          "lung" = "lu",
                          "ovarian" = "ov")
  if (is.null(cancer)) {
    rlang::abort("cancer argument is missing and is required")
  } else if (!cancer %in% names(cancer_code_lookup)) {
    rlang::abort(glue::glue("cancer = '{cancer}' is not a valid input. Please use one of 'bladder', 'breast', 'lung', or 'ovarian'."))
  }
  cancer_code <- cancer_code_lookup[[cancer]]
  
  # NA empty strings and handle dates
  .data <- .data %>% 
    mutate(across(where(is.character), ~ na_if(., ""))) %>%
    create_partial_date_registry(rxdatemostdefinsurg) %>% 
    create_partial_date_registry(rxdatedxstgproc) %>% 
    mutate(datecasecompleted= ymd(datecasecompleted, quiet = TRUE),
           datecasecompletedcoc = ymd(datecasecompletedcoc, quiet = TRUE),
           dateoflastcontact = ymd(dateoflastcontact, quiet = TRUE),
           rxdatemostdefinsurg = ymd(rxdatemostdefinsurg, quiet = TRUE),
           rxdatedxstgproc = ymd(rxdatedxstgproc, quiet = TRUE)) %>% 
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>% 
    mutate(rxdatemostdefinsurg_year = ifelse(!is.na(rxdatemostdefinsurg),
                                             year(rxdatemostdefinsurg),
                                             as.numeric(str_extract(rxdatemostdefinsurg_partial, "^\\d{4}"))),
           rxdatedxstgproc_year = ifelse(!is.na(rxdatedxstgproc),
                                         year(rxdatedxstgproc),
                                         as.numeric(str_extract(rxdatedxstgproc_partial, "^\\d{4}")))) %>% 
    select(-stagingalgorithmschema, -dateofdiagnosis, -primarysite, -dateofdiagnosis_year) %>%
    distinct()
  # Mapping registry rxsummsurgprimsite to OC version
  # Registry data contains site specific codes we can use the OC mappings with. We need to differentiate 
  # between codes that need a cancer-specific (e.g., "bl" or "br") and an "ot" suffix though.
  oc_rxsummsurgprimsite_other_codes <- map_oc_value %>%
    filter(valueset_name == "rxSummSurgPrimSite (r)") %>%
    filter(grepl("ot", code)) %>%
    pull(code) %>%
    str_remove("ot")
  oc_rxsummsurgprimsite_cohort_codes <- map_oc_value %>%
    filter(valueset_name == "rxSummSurgPrimSite (r)") %>%
    filter(grepl(cancer_code, code)) %>% 
    pull(code) %>%
    str_remove("bln") %>% # these occur in bladder sometimes, removing to use bladder codes only
    str_remove(cancer_code)
  # append 'cancer_code' or "ot" to unmapped codes and apply OC mapping
  .data <- .data %>%
    mutate(rxsummsurgprimsite = case_when(cancer == "bladder" & rxsummsurgprimsite == "40" ~ "40bln", # special bladder code
                                          rxsummsurgprimsite %in% oc_rxsummsurgprimsite_cohort_codes ~ paste0(rxsummsurgprimsite, cancer_code),
                                          rxsummsurgprimsite %in% oc_rxsummsurgprimsite_other_codes ~ paste0(rxsummsurgprimsite, "ot"),
                                          TRUE ~ rxsummsurgprimsite)) %>% 
    map_oc(rxsummsurgprimsite, "rxSummSurgPrimSite (r)") %>%
    # in case any registry codes escape the OC mapping, in theory this could only be 98 = "Site specific codes; special"
    coalesce_map_naaccr("rxsummsurgprimsite", spmd) %>% 
    ## we also account for occasions of lymph node procedures which are not tied to a primary site surgery
    ## we don't want to remove these rows so assign them as "99 - Unknown if surgery performed"
    mutate(rxsummsurgprimsite = case_when(is.na(rxsummsurgprimsite) & 
                                            rxsummscopereglnsur != "Unknown or Not applicable" ~ "99 - Unknown if surgery performed",
                                          is.na(rxsummsurgprimsite) &
                                            !is.na(rxsummreglnexamined) ~ "99 - Unknown if surgery performed",
                                          TRUE ~ rxsummsurgprimsite))
  
  # Map remaining fields
  .data %>% 
    map_oc(reasonfornosurgery, "reasonForNoSurgery") %>% 
    coalesce_map_naaccr("lymphvascularinvasion", spmd) %>% # prefer NAACR map
    coalesce_map_naaccr("rxsummsurgicalmargins", spmd) %>% # prefer NAACR map
    map_oc(rxsummscopereglnsur, "rxSummScopeRegLnSur (r)") %>% 
    # consolidate unknown categories of surgical margins and scope of regional LN surgery variables
    # there is a minor casing difference between the OC and NAACCR mappings for the unknow category in these
    # vars
    mutate(rxsummscopereglnsur = ifelse(rxsummscopereglnsur == "Unknown or Not applicable", 
                                        "Unknown or not applicable", 
                                        rxsummscopereglnsur)) %>% 
    # map residual tumor volume for ovarian
    {if (cancer == "ovarian") map_residual_tumor_volume(.) else .}

}

# Wrangles primary site surgery data from OC and registry into a single, logically aligned data frame. Expects
# outputs of the clean_primary_site_surgery_oc and clean_primary_site_surgery_registry functions.
# Principally, this prioritizes MA over registry data and then applies logic to ensure that the reason for 
# no surgery (and eventual derived flag) and the discrete surgery event data (e.g., site of surgery) align. 
# - we need to account for patients with multiple, unique reasonForNoSurgery, and 
# - with disagreements between reasonForNoSurgery and rxSummSurgPrimSite
# Approach was discussed and aligned on in https://syapse.atlassian.net/browse/CA-3723
# See document outlining steps here: https://docs.google.com/document/d/1zo-Pe9EgvJBcNWl3zIoOouXpzcWTjcsOYyCHNNPK358/edit
wrangle_primary_site_surgery <- function(.data_oc, .data_registry, cancer, patients, spmd = spmd_con()) {
  if (is.null(cancer)) {
    rlang::abort("cancer argument is missing and is required")
  } else if (!cancer %in% c('bladder', 'breast', 'lung', 'ovarian')) {
    rlang::abort(glue::glue("cancer = '{cancer}' is not a valid input. Please use one of 'bladder', 'breast', 'lung', or 'ovarian'."))
  }
  # 1. Combine data and prioritize MA > registry
  .data <- .data_registry %>%
    mutate(sourceschema = "registry") %>% 
    # per https://syapse.atlassian.net/browse/MDRLM-218 rxdatemostdefinsurg corresponds to rxsummsurgprimsite, not rxdatesurgery
    rename_with(~ str_replace(., "rxdatemostdefinsurg", "rxdatesurgery")) %>%
    bind_rows(.data_oc %>%
                mutate(sourceschema = "ma")) %>%
    # MA > registry
    prioritize_source() %>% 
    filter(!(is.na(reasonfornosurgery) & is.na(rxsummsurgprimsite))) # these are empty - move this up or to clean functions?
  # 2. Prioritize amongst multiple, unique reasonForNoSurgery
  .data <- .data %>% 
    # backfill cases with missing reason for no surgery but primary site information available - move to clean functions?
    mutate(reasonfornosurgery = case_when(grepl("[1-8][0-9]|90|02", rxsummsurgprimsite) & is.na(reasonfornosurgery) ~ "0 - Surgery of the primary site was performed",
                                          grepl("00|99|01", rxsummsurgprimsite) & is.na(reasonfornosurgery) ~ "9 - It is unknown if surgery of the primary site was recommended or performed.",
                                          TRUE ~ reasonfornosurgery)) %>% 
    # prioritize known reasons (0-7) over unknown reasons (8-9)
    mutate(reasonfornosurgery_rank = ifelse(grepl("^[0-7]", reasonfornosurgery),
                                            1,
                                            0)) %>% 
    group_by(patientid) %>% 
    filter(reasonfornosurgery_rank == max(reasonfornosurgery_rank, na.rm = TRUE)) %>% 
    # then take the most recent case (registry only)
    # this is defined by using datecasecompletedcoc first (registry only) as this is filled in for "analytic 
    # cases" which area priority (priority), then the case with the most recent information (dateoflastcontact), 
    # and then datecasecompleted (max of date complete and date of last abstraction from MA)
    filter(datecasecompletedcoc == suppressWarnings(max(datecasecompletedcoc, na.rm = TRUE)) | 
             all(is.na(datecasecompletedcoc))) %>% 
    filter(dateoflastcontact == suppressWarnings(max(dateoflastcontact, na.rm = TRUE)) | 
             all(is.na(dateoflastcontact))) %>% 
    filter(datecasecompleted == suppressWarnings(max(datecasecompleted, na.rm = TRUE)) | 
             all(is.na(datecasecompleted))) %>% 
    ungroup()
  # note this isn't perfect and there are a handful of registry patients still duplicated, they will be handled 
  # in the next step
  # 3. Align mismatched reasonForNoSurgery and rxSummSurgPrimSite
  # There are several cases we need to account for (see above document for detail):
  # 1) If a patient was coded as 1-9, but there was evidence of surgery to the primary site, they were recoded 
  #    to 0. 
  # 2) If a patient was coded as 0, but there was no evidence of surgery to the primary site (either 00-none, 
  #    01-biopsy non-primary, 99-unknown, or missing/Null/None), the primary site surgery flag was set to 'unknown' and 
  #    the reason for no surgery was recoded to 9. 
  # 3) If a patient was coded as 0-7, and all surgery primary sites were coded as 99-unknown, the primary 
  #    site surgery flag was set to 'unknown' and the reason for no surgery was recoded to 9. 
  # 4) If a patient was coded as 8-9, and all surgery primary sites were coded as 00-none or 01-biopsy 
  #    non-primary, the primary site surgery flag was set to 'unknown' and the reason for no surgery was kept as 8 or 9 
  #    (i.e., no change to reasonForNoSurgery).
  .data <- .data %>% 
    group_by(patientid) %>% 
    # note case #4 isn't implemented here because reasonForNoSurgery doesn't change
    mutate(reasonfornosurgery = case_when(
      # If a patient was coded as 1-9, but there was evidence of surgery to the primary site, they were recoded 
      # to 0. Evidence of surgery = code of 10-90.
      any(grepl("[1-9]|Marked No in Medrio", reasonfornosurgery)) & 
        any(grepl("[1-8][0-9]|90|02", rxsummsurgprimsite)) ~ "0 - Surgery of the primary site was performed",
      # If a patient was coded as 0, but there was no evidence of surgery to the primary site (either 00-none, 
      # 99-unknown, or missing), the primary site surgery flag was set to 'unknown' and the reason for no surgery 
      # was recoded to 9. 
      any(grepl("0", reasonfornosurgery)) & 
        all(grepl("00|99|01", rxsummsurgprimsite) | 
              is.na(rxsummsurgprimsite) |
              rxsummsurgprimsite == "NULL" |
              rxsummsurgprimsite == "None"
        ) ~ "9 - It is unknown if surgery of the primary site was recommended or performed.",
      # If a patient was coded as 0-7, and all surgery primary sites were coded as 99-unknown, the primary 
      # site surgery flag was set to 'unknown' and the reason for no surgery was recoded to 9. 
      all(grepl("[0-7]", reasonfornosurgery)) & 
        all(grepl("99", rxsummsurgprimsite)) ~ "9 - It is unknown if surgery of the primary site was recommended or performed.",
      # If a patient was Marked No in Medrio and no surgery was performed at primary site, set reason for no surgery to 
      # Surgery of the primary site was not performed; reason unknown"
      all(reasonfornosurgery == "Marked No in Medrio") &
        all(rxsummsurgprimsite == "00 - No surgery") ~ "Surgery of the primary site was not performed; reason unknown",
      ## otherwise stays the same
      TRUE ~ reasonfornosurgery)) %>% 
    # catch all in case patients still manage to have multiple reasonfornosurgery
    mutate(reasonfornosurgery = ifelse(n_distinct(reasonfornosurgery) > 1,
                                       "9 - It is unknown if surgery of the primary site was recommended or performed.",
                                       reasonfornosurgery),
    ) %>% 
    ungroup() %>% 
    mutate(primary_site_surgery_flag = case_when(grepl("0", reasonfornosurgery) ~ "Yes",
                                                 grepl("[1-7]", reasonfornosurgery) ~ "No",
                                                 reasonfornosurgery == "Surgery of the primary site was not performed; reason unknown" ~ "No",
                                                 grepl("8|9", reasonfornosurgery) ~ "Unknown",
                                                 TRUE ~ "Unknown")) %>% 
    rename(reason_for_no_surgery = reasonfornosurgery)
  .data %>% 
    select(patientid, primary_site_surgery_flag, reason_for_no_surgery) %>% 
    distinct() %>% 
    confirm_unique(patientid, return_df = FALSE, error_msg = "Surgery - patients have multiple, unique reasonForNoSurgery still!") %>% 
    invisible()
  # We also need to to the reverse and align surgery data (rxSummSurgPrimSite, rxDateSurgery, and residualTumorVolumePostCytoreduction)
  # with the newly aligned reason_for_no_surgery and primary_site_surgery_flag
  .data <- .data %>% 
    # If a patient was coded as 1-9, but there was evidence of surgery to the primary site, they were recoded 
    # to 0. We want to remove the 00, 99, and 01 records for bladder, lung, and ovarian, and we want to remove 
    # the 00 and 01 records only - **leave 99 due to lymph node procedures**.
    mutate(rxsummsurgprimsite = case_when(cancer != "breast" & 
                                            grepl("0", reason_for_no_surgery) & 
                                            grepl("00|99|01", rxsummsurgprimsite) ~ "Remove",
                                          cancer == "breast" & 
                                            grepl("0", reason_for_no_surgery) & 
                                            grepl("00|01", rxsummsurgprimsite) ~ "Remove",
                                          # If a patient was coded as 0, but there was no evidence of surgery to the primary site (either 00-none, 01-biopsy non-primary
                                          # 99-unknown, or missing), the primary site surgery flag was set to 'unknown' and the reason for no surgery 
                                          # was recoded to 9. 
                                          # If a patient was coded as 0-7, and all surgery primary sites were coded as 99-unknown, the primary 
                                          # site surgery flag was set to 'unknown' and the reason for no surgery was recoded to 9. 
                                          # If a patient was coded as 8-9, and all surgery primary sites were coded as 00-none or 01-biopsy 
                                          # non-primary, the primary site surgery flag was set to 'unknown' and the reason for no surgery was kept as 8 or 9.
                                          # In all three cases, we want to set rxsummsurgprimsite to 99 and NA the dates
                                          grepl("8|9", reason_for_no_surgery) ~ "99 - Unknown if surgery performed",
                                          # NA cases
                                          # If coded as 0 or recoded to 0, remove NA rxsummsurgprimsite records
                                          grepl("0", reason_for_no_surgery) & 
                                            (is.na(rxsummsurgprimsite) | 
                                               rxsummsurgprimsite %in% c("NULL", "None")) ~ "Remove",
                                          # If coded as 1-7, recode NA rxsummsurgprimsite records to 00
                                          grepl("[1-7]", reason_for_no_surgery) & 
                                            (is.na(rxsummsurgprimsite) |
                                               rxsummsurgprimsite %in% c("NULL", "None")) ~ "00 - No surgery",
                                          # If coded as 8/9, recode NA rxsummsurgprimsite records to 99
                                          grepl("[8-9]", reason_for_no_surgery) & 
                                            (is.na(rxsummsurgprimsite) |
                                               rxsummsurgprimsite %in% c("NULL", "None")) ~ "99 - Unknown if surgery performed",
                                          # All others leave as is
                                          TRUE ~ rxsummsurgprimsite)) %>% 
    # If no evidence of surgery, NA the date fields and recode residualtumvolpostcytoreduction 
    # in theory, rxsummsurgprimsite == 01-biopsy non-primary or 99 - lymphnode procedure (breast) could be 
    # associated with a reasonfornosurgery == 1-7 at this point, and we are leaving these records (including the dates etc.)
    mutate(rxdatesurgery = case_when(cancer != "breast" & grepl("00|99", rxsummsurgprimsite) ~ NA_Date_,
                                     cancer == "breast" & grepl("00", rxsummsurgprimsite) ~ NA_Date_,
                                     TRUE ~ rxdatesurgery),
           rxdatesurgery_year = case_when(cancer != "breast" & grepl("00|99", rxsummsurgprimsite) ~ NA_real_,
                                          cancer == "breast" & grepl("00", rxsummsurgprimsite) ~ NA_real_,
                                          TRUE ~ rxdatesurgery_year),
           rxdatesurgery_granularity = case_when(cancer != "breast" & grepl("00|99", rxsummsurgprimsite) ~ "NONE",
                                                 cancer == "breast" & grepl("00", rxsummsurgprimsite) ~ "NONE",
                                                 TRUE ~ rxdatesurgery_granularity),
           # only for ovarian right now, but including same breast logic
           residualtumvolpostcytoreduction = case_when(grepl("00", rxsummsurgprimsite) ~ "97 - No cytoreductive surgery performed",
                                                       cancer != "breast" & grepl("99", rxsummsurgprimsite) ~ NA_character_,
                                                       TRUE ~ residualtumvolpostcytoreduction)
    )
  # Perform the removal flagged above
  if (.data %>% count_patients() %>% pull(patients) != .data %>% filter(rxsummsurgprimsite != "Remove") %>% count_patients() %>% pull(patients)) {
    rlang::warn("Patients lost from surgery df when removing mis-algined reasonForNoSurgery and rxSummSurgPrimSite! These will become unknowns.")
  }
  .data <- .data %>% 
    filter(rxsummsurgprimsite != "Remove")
  # Check for any lingering mismatches
  mismatched_surgery <- .data %>% 
    mutate(mismatch = case_when(grepl("[1-7]", reason_for_no_surgery) & 
                                  !grepl("00|01", rxsummsurgprimsite) ~ "reasonForNoSurgery in 1-7, known/unknown surgery in rxSummSurgPrimSite",
                                grepl("8|9", reason_for_no_surgery) & 
                                  !grepl("99", rxsummsurgprimsite) ~ "reasonForNoSurgery in 8 or 9, known/no surgery in rxSummSurgPrimSite",
                                grepl("[1-9]", reason_for_no_surgery) &
                                  !is.na(rxdatesurgery) ~ "reasonForNoSurgery in 1-9, non-missing date in rxDateSurgery",
                                grepl("0", reason_for_no_surgery) & 
                                  (grepl("00|99", rxsummsurgprimsite) | 
                                     is.na(rxsummsurgprimsite)) ~ "reasonForNoSurgery == 0, evidence of no surgery/unknown surgery in rxSummSurgPrimSite")) %>% 
    filter(rxsummsurgprimsite != "01 - Biopsy non-primary site" & grepl("[1-7]", reason_for_no_surgery)) %>%
    filter(!is.na(mismatch))
  if (dim(mismatched_surgery)[1] > 0) {
    rlang::warn("Patients with mismatched reasonForNoSurgery and surgery fields (rxSummSurgPrimSite, or rxDateSurgery)")
    mismatched_surgery %>% 
      count(sourceschema, mismatch) %>% 
      print(n = 40)
  }
  # 4. Create final variables
  # number of lymph nodes removed is a mixed string- numeric variable 
  # Split into a categorical and continuous variables as per
  # https://docs.google.com/document/d/11n4WHKUgruOd4YPP0iZx2662IOfG1vB7dVWECdtIU7k/edit?usp=sharing\
  .data <- .data %>% 
    mutate(
      rxsummreglnexamined = as.numeric(rxsummreglnexamined),
      surgery_regional_ln_examined_number = case_when(rxsummreglnexamined <= 90 ~ rxsummreglnexamined, # numeric version of variable
                                                      (rxsummreglnexamined > 90 & rxsummreglnexamined <= 94) ~ 90, #cleaning in case
                                                      TRUE ~ NA_real_), 
      rxsummreglnexamined = case_when(rxsummreglnexamined == 00 ~ "No nodes were examined",
                                      rxsummreglnexamined >= 1 &
                                        rxsummreglnexamined < 90 ~ "1-89 nodes were examined",
                                      rxsummreglnexamined == 90 ~ "90 or more nodes were examined",
                                      (rxsummreglnexamined > 90 & rxsummreglnexamined <= 94) ~ "90 or more nodes were examined",
                                      TRUE ~ as.character(rxsummreglnexamined))
    ) %>% 
    # map codes not covered above
    coalesce_map_naaccr("rxsummreglnexamined", spmd) %>% 
    rename(surgery_regional_ln_examined = rxsummreglnexamined)
  if (cancer == "ovarian") {
    # Assign and prioritize a residual tumor volume per surgery
    # If multiple residual disease statuses associated with the same surgery, select the most specific code:
    # 1. 00, 50, 60 would be prioritized over the other codes.
    # 2. 00 > 50 > 60 (since 00 is the most favorable outcome and 60 is the least favorable)
    # 3. 80 > 70 (since 80 offers more detail)
    residual_disease_priority <- c("00 - No gross residual tumor nodules",
                                   "50 - Residual tumor nodule(s) 1 centimeter (cm) or less",
                                   "60 - Residual tumor nodule(s) greater than 1 cm",
                                   "80 - Procedure described as optimal debulking and size of residual tumor nodule(s) not given",
                                   "70 - Macroscopic residual tumor nodule(s), size not stated",
                                   # the exact order of these is not that big of a deal
                                   "97 - No cytoreductive surgery performed",
                                   "98 - Not applicable: Information not collected for this case",
                                   "99 - Not documented in medical record. Residual tumor status after cytoreductive surgery not assessed or unknown if assessed")
    .data <- .data %>%
      add_sorting(residualtumvolpostcytoreduction, "custom", residual_disease_priority) %>% 
      group_by(patientid, rxsummsurgprimsite, rxdatesurgery, rxdatesurgery_year) %>% # identifying same surgeries
      arrange(residualtumvolpostcytoreduction) %>%
      slice_head() %>%
      ungroup() %>% 
      mutate(residualtumvolpostcytoreduction = as.character(residualtumvolpostcytoreduction)) %>%
      mutate(residualtumvolpostcytoreduction_category = case_when(grepl("^(00|50|80)", residualtumvolpostcytoreduction) ~ "1 - Optimal debulking: <=1 cm residual disease",
                                                                  grepl("^60", residualtumvolpostcytoreduction) ~ "2 - Sub-optimal debulking: >1 cm residual disease",
                                                                  grepl("^(70|99)", residualtumvolpostcytoreduction) ~ "9 - Unknown",
                                                                  grepl("^(97|98)", residualtumvolpostcytoreduction) ~ "8 - Not applicable/no cytoreductive surgery"))
  }
  # 5. Rename and subset to final variables
  .data <- .data %>% 
    select_rename_primary_site_surgery(cancer) %>% 
    distinct()
  # 6. Set the flags to unknown for patients without surgery data
  .data %>% 
    full_join(patients %>% select(patientid), by = "patientid") %>% 
    replace_na(list(primary_site_surgery_flag = "Unknown",
                    reason_for_no_surgery = "9 - It is unknown if surgery of the primary site was recommended or performed."))

}

# Pipeline function for deriving primary site surgery variables and df-in-df
derive_primary_site_surgery <- function(cancer, patients, cohort_query, cancer_type, spmd = spmd_con()) {
  oc <- get_primary_site_surgery_oc(cancer = cancer,
                                    patients = patients,
                                    spmd = spmd) %>%
    clean_primary_site_surgery_oc(cancer = cancer,
                                  spmd = spmd)
  registry <- get_primary_site_surgery_registry(cohort_query = cohort_query,
                                                cancer_type = cancer_type) %>%
    clean_primary_site_surgery_registry(cancer = cancer,
                                        spmd = spmd)
  wrangle_primary_site_surgery(.data_oc = oc, 
                               .data_registry = registry,
                               cancer = cancer,
                               patients = patients,
                               spmd = spmd) %>% 
    nest(primary_site_surgery = c(-patientid, -primary_site_surgery_flag, -reason_for_no_surgery))

}

# Helper functions ----

# Select and renamed variables for the primary site surgery ADS data
select_rename_primary_site_surgery <- function(.data, cancer) {
  
  if (!cancer %in% c('bladder', 'breast', 'lung', 'ovarian')) {
    rlang::abort(glue::glue("cancer = '{cancer}' is not a valid input. Please use one of 'bladder', 'breast', 'lung', or 'ovarian'."))
  }
  
  select_rename_vars <- function(.data, cols = NULL) {
    .data %>% 
      select(patientid,
             primary_site_surgery_flag,
             reason_for_no_surgery,
             primary_site_surgery = rxsummsurgprimsite,
             surgery_date = rxdatesurgery,
             surgery_date_granularity = rxdatesurgery_granularity,
             surgery_date_year = rxdatesurgery_year,
             surgery_lymphovascular_invasion = lymphvascularinvasion,
             surgery_surgical_margins = rxsummsurgicalmargins,
             surgery_scope_regional_ln = rxsummscopereglnsur,
             all_of(cols),
             surgery_source = sourceschema)
  }
  
  if (cancer == "bladder") {
    .data %>% 
      select_rename_vars(c("surgery_tnm_t" = "tnmt",
                           "surgery_tnm_n" = "tnmn",
                           "surgery_tnm_m" = "tnmm"))
  } else if (cancer %in% c("breast", "lung")) {
    .data %>% 
      select_rename_vars(c("surgery_regional_ln_examined", 
                           "surgery_regional_ln_examined_number"))
  } else if (cancer == "ovarian") {
    .data %>% 
      select_rename_vars(c("surgery_regional_ln_examined", 
                           "surgery_regional_ln_examined_number",
                           "surgery_residual_disease_status" = "residualtumvolpostcytoreduction",
                           "surgery_residual_disease_status_category" = "residualtumvolpostcytoreduction_category"))
  }
  
}

# Maps the residualTumVolPostCytoreduction field from registry according to the OC mapping. 
map_residual_tumor_volume <- function(.data) {
  # Mapping registry residualtumvolpostcytoreduction to OC version
  # We have older codes in residualtumvolpostcytoreduction field from the registry. The MA team is using the
  # latest codes. We decided to map the old codes to the new labels (and then we can apply the OC map). The
  # codes changed after version 18.
  # For future development, it could be worthwhile to ingest all NAACCR codes across all versions for each data
  # element. Then append columns describing the dx years this applies to.
  residual_tumor_naaccr_map_old_to_new <- tribble(~code, ~value,
                                                  "10", "50",
                                                  "20", "50",
                                                  "30", "60",
                                                  "40", "60",
                                                  "90", "70",
                                                  "91", "70",
                                                  "92", "80",
                                                  "93", "80")
  .data %>%
    left_join(residual_tumor_naaccr_map_old_to_new, by = c("residualtumvolpostcytoreduction" = "code")) %>%
    mutate(residualtumvolpostcytoreduction = coalesce(value, residualtumvolpostcytoreduction)) %>%
    select(-value) %>% 
    map_oc(residualtumvolpostcytoreduction, "residualTumVolPostCytoreduction (r)")
}

# Cleans up primary site surgery data from OC for the bladder cohort. Specifically handles mapping of TNM and
# corrects 40/40bl to 40bln in rxSummSurgPrimSite (the former are not real codes).
clean_primary_site_surgery_oc_bladder <- function(.data) {
  .data <- .data %>%
    # Map TNMs
    map_oc(tnmt, "TnmT (r)") %>%
    map_oc(tnmn, "TnmN (r)") %>%
    map_oc(tnmm, "TnmM (r)") %>%
    mutate(tnmt = toupper(tnmt),
           tnmt = case_when((tnmt == "P0" | tnmt == "PT0" | tnmt == "T0") ~ "T0",
                            (tnmt == "P1" | tnmt == "PT1" | tnmt == "T1") ~ "T1",
                            (tnmt == "P1A" | tnmt == "PT1A" | tnmt == "T1A") ~ "T1",
                            (tnmt == "P1B" | tnmt == "PT1B" | tnmt == "T1B") ~ "T1",
                            (tnmt == "P1C" | tnmt == "PT1C" | tnmt == "T1C") ~ "T1",
                            (tnmt == "P2" | tnmt == "PT2" | tnmt == "T2") ~ "T2",
                            (tnmt == "P2A" | tnmt == "PT2A" | tnmt == "T2A") ~ "T2a",
                            (tnmt == "P2B" | tnmt == "PT2B" | tnmt == "T2B" | tnmt == "P2C" | tnmt == "PT2C" | tnmt == "T2C") ~ "T2b",
                            (tnmt == "P3"| tnmt == "PT3" | tnmt == "T3" | tnmt == "T3C") ~ "T3",
                            (tnmt == "P3A" | tnmt == "PT3A" | tnmt == "T3A") ~ "T3a",
                            (tnmt == "P3B" | tnmt == "PT3B" | tnmt == "T3B") ~ "T3b",
                            (tnmt == "P4" | tnmt == "PT4" | tnmt == "T4") ~ "T4",
                            (tnmt == "PA" | tnmt == "PTA" | tnmt == "TA") ~ "Ta",
                            (tnmt == "P4A" | tnmt == "PT4A" | tnmt == "T4A") ~ "T4a",
                            (tnmt == "PT4B" | tnmt == "PT4B" | tnmt == "T4B") ~ "T4b",
                            (tnmt == "PX" | tnmt == "PTX" | tnmt == "TX") ~ "TX",
                            (tnmt == "PA" | tnmt == "PTA" | tnmt == "TA") ~ "Ta",
                            (tnmt %in% c("PIS", "PTIS", "TIS", "PIDS", "PTIDS", "TIDS", "DCIS", "LCIS", "TIS (LCIS)", "PTIS(DCIS)", "TIS (DCIS)")) ~ "Tis",
                            is.na(tnmt) ~ NA_character_,
                            TRUE ~ "Other - FLAG"),
           tnmn = toupper(tnmn),
           tnmn = case_when((tnmn == "C0" | tnmn == "CN0" | tnmn == "P0" | tnmn == "PN0" | tnmn == "N0" | tnmn == "N0(I-)") ~ "N0",
                            (tnmn == "N1" | tnmn == "N1A" | tnmn == "N1B" | tnmn == "N1C" | tnmn == "PN1" | tnmn == "PN1A") ~ "N1",
                            (tnmn == "CN2" | tnmn == "N2" | tnmn == "N2A" | tnmn == "N2B" | tnmn == "P2" | tnmn == "PN2") ~ "N2",
                            (tnmn == "N3" | tnmn == "P3" | tnmn == "PN3") ~ "N3",
                            (tnmn == "CNX" | tnmn == "NX" | tnmn == "PNX" | tnmn == "PTX") ~ "NX",
                            is.na(tnmn) ~ NA_character_,
                            TRUE ~ "Other - FLAG"),
           tnmm = toupper(tnmm),
           tnmm = case_when((tnmm == "C0" | tnmm == "CM0" | tnmm == "P0" | tnmm == "PM0" | tnmm == "M0") ~ "M0",
                            (tnmm == "M1" | tnmm == "PM1" | tnmm == "CM1" | tnmm == "P1" | tnmm == "C1") ~ "M1",
                            (tnmm == "M1A" | tnmm == "PM1A" | tnmm == "CM1A") ~ "M1a",
                            (tnmm == "M1B" | tnmm == "PM1B" | tnmm == "CM1B" | tnmm == "M1C" | tnmm == "CM1C") ~ "M1b",
                            (tnmm == "MX" | tnmm == "PX" | tnmm == "PMX" | tnmm == "CX") ~ "MX",
                            is.na(tnmm) ~ NA_character_,
                            TRUE ~ "Other - FLAG")) %>%
    # specific-code correction
    mutate(rxsummsurgprimsite = if_else(rxsummsurgprimsite %in% c("40", "40bl"),
                                        "40bln",
                                        rxsummsurgprimsite,
                                        rxsummsurgprimsite))
  # Check for normalization escapes
  tnm_check <- .data %>%
    filter(tnmt == "Other - FLAG" | tnmn == "Other - FLAG" | tnmm == "Other - FLAG")
  if (dim(tnm_check)[1] > 0) {
    rlang::warn(glue::glue("{nrow(tnm_check)} t or n stage escape normalization, review warranted"))
    tnm_check %>% print(50)}
  
  .data
}
