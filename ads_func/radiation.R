
get_radiation_oc <- function(cancer, patients = NULL, spmd = spmd_con()) {
  
  # Check patients argument
  if (!is.null(patients)) {
    if (!is.data.frame(patients)) {
      rlang::abort("patients is not a data frame")
    } else {
      filter_columns <- patients %>% 
        select(matches("patientid|participantid")) %>% 
        colnames()
    }
  }
  
  run_query(glue::glue("SELECT
                        patientid,
                        participantid,
                        formdata::json -> 'radiation_therapy' -> 'group' ->> 'ma_reasonfornoradiation' AS ma_reasonfornoradiation,
                        formdata::json -> 'radiation_therapy' -> 'radiation' AS radiation
                        FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}' AND
                        formdata::json -> 'radiation_therapy' IS NOT NULL"),
            connection = spmd) %>%
    {if (!is.null(patients)) semi_join(., patients, by = filter_columns) else .} %>%
    mutate(radiation = purrr::map(radiation, fromJSON_na)) %>%
    unnest(radiation, keep_empty = TRUE) %>%
    rename_oc() %>%
    convert_empty_strings() %>%
    separate_rows(radtreatmentvolume, sep = ",") %>%
    # encountered missing data in QC for ovarian cohort on 9/29 and added this workaround 
    {if (!"rxdateradiationym" %in% colnames(.)) mutate(., rxdateradiationym = NA_character_) else .} %>% 
    {if (!"rxdateradiationendedy" %in% colnames(.)) mutate(., rxdateradiationendedy = NA_character_) else .} %>% 
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>%
    extract_year_oc(rxdateradiation) %>%
    extract_year_oc(rxdateradiationended) %>%
    map_oc(radiationintent, "radiationIntent (r)") %>%
    map_oc(ongoingrt, "ongoingRT (r)") %>%
    # silly error
    mutate(radtreatmentvolume = if_else(radtreatmentvolume == "9", "09", radtreatmentvolume, radtreatmentvolume)) %>%
    map_oc(radtreatmentvolume, "radTreatmentVolume (r)") %>%
    # this fills in inappropriately I think
    # replace_na(list(ongoingrt = "Unknown / Not stated",
    #                 radiationintent = "Unknown")) %>%
    # filter out these obs. as could be affected by registry error without CTR correction
    filter(!radiationintent %in% c("Palliative"))
}

get_radiation_registry <- function(cohort_query, cancer_type, cols = c()) {
  get_registry(cohort_query,
               cancer_type,
               c("reasonfornoradiation",
                 "rxdateradiation",
                 "rxdateradiationended",
                 "radtreatmentvolume",
                 "phase1radiationprimarytxvolume",
                 "datecasecompleted",
                 "datecasecompletedcoc",
                 "dateoflastcontact",
                 "classofcase",
                 "registryid", 
                 cols)) %>%
    create_partial_date_registry(rxdateradiation) %>% 
    create_partial_date_registry(rxdateradiationended) %>% 
    ## phase1radiationprimarytxvolume is the current variable
    ## radtreamtne volume was the variable pre-2018
    ## they have different maps/options - most irritiatingly the 00 mapping, 
    ### so hack that first using the post-2018 verbiage
    mutate(across(.cols = c(radtreatmentvolume, phase1radiationprimarytxvolume),
                  .fns = ~ifelse(. == "00", "No radiation treatment", .))) %>% 
    coalesce_map_naaccr("radtreatmentvolume", cohort_query$src$con) %>% 
    coalesce_map_naaccr("phase1radiationprimarytxvolume", cohort_query$src$con) %>% 
    ## phase1radiationprimarytxvolume is the more detailed map
    ## therefore use that if both variables exist per patient
    ## but naming the variable "radtreatmentvolume" for ease of use with MA data later
    mutate(radtreatmentvolume = coalesce(phase1radiationprimarytxvolume, radtreatmentvolume),
           radtreatmentvolume = ifelse(radtreatmentvolume == "43", "Unknown", radtreatmentvolume), 
           datecasecompleted= ymd(datecasecompleted, quiet = TRUE),
           datecasecompletedcoc = ymd(datecasecompletedcoc, quiet = TRUE),
           dateoflastcontact = ymd(dateoflastcontact, quiet = TRUE),
           rxdateradiation = ymd(rxdateradiation, quiet = TRUE),
           rxdateradiationended = ymd(rxdateradiationended, quiet = TRUE)) %>% 
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>% 
    mutate(rxdateradiation_year = ifelse(!is.na(rxdateradiation),
                                         year(rxdateradiation),
                                         as.numeric(str_extract(rxdateradiation_partial, "^\\d{4}"))),
           rxdateradiationended_year = ifelse(!is.na(rxdateradiationended),
                                              year(rxdateradiationended),
                                              as.numeric(str_extract(rxdateradiationended_partial, "^\\d{4}"))))  %>% 
    select(-stagingalgorithmschema, -dateofdiagnosis, -primarysite, -dateofdiagnosis_year, -phase1radiationprimarytxvolume) %>%
    distinct()
}

# We need to account for patients with multiple, unique reasonForNoRadiation and with disagreements between 
# reasonForNoRadiation and radiation sites. This step must be accomplished before subsequent steps in order 
# to keep the data internally consistent.
# Mirroring approach from surgery
resolve_discordant_radiation_reason <- function(.data) {
  # Step 1: If multiple unique reasonForNoSurgery 
  # Prioritize 0-7 over 8-9, then take the most recent case
  step1 <-  .data %>% 
    mutate(reasonfornoradiation = case_when(!radtreatmentvolume %in% c("No radiation treatment",
                                                                       "Unknown") 
                                            & is.na(reasonfornoradiation) ~ "0",
                                            radtreatmentvolume %in% c("No radiation treatment",
                                                                      "Unknown") 
                                            & is.na(reasonfornoradiation) ~ "9",
                                            TRUE ~ reasonfornoradiation),
           reasonfornoradiation_rank = ifelse(grepl("^[0-7]", reasonfornoradiation),
                                              1,
                                              0)) %>% 
    group_by(patientid) %>% 
    # prioritize 0-7 over 8-9
    filter(reasonfornoradiation_rank == max(reasonfornoradiation_rank, na.rm = TRUE)) %>% 
    # then take the most recent case 
    # this is defined by using datecasecompletedcoc first (registry only) as this is filled in for "analytic 
    # cases" which area priority (priority), then the case with the most recent information (dateoflastcontact), 
    # and then datecasecompleted (max of date complete and date of last abstraction from MA)
    {if ("datecasecompletedcoc" %in% colnames(.)) filter(., datecasecompletedcoc == suppressWarnings(max(datecasecompletedcoc, na.rm = TRUE)) | all(is.na(datecasecompletedcoc))) else .} %>% 
    {if ("dateoflastcontact" %in% colnames(.)) filter(., dateoflastcontact == suppressWarnings(max(dateoflastcontact, na.rm = TRUE)) | all(is.na(dateoflastcontact))) else .} %>% 
    {if ("datecasecompleted" %in% colnames(.)) filter(., datecasecompleted == suppressWarnings(max(datecasecompleted, na.rm = TRUE)) | all(is.na(datecasecompleted))) else .} %>% 
    ungroup()
  # Step 2: Aligning mismatched reasonForNoRadiation and site
  # There are several cases we need to account for (see above document for detail):
  # 1) If a patient was coded as 1-9, but there was evidence of radiation, they were recoded 
  #    to 0. 
  # 2) If a patient was coded as 0, but there was no evidence of radiation (either 00-none, 
  #    99-unknown, or Null/None), the primary site radiation flag was set to 'unknown' and 
  #    the reason for no radiation was recoded to 9.  (ALLOWING FOR MISSING DUE TO HIGH MISSINGNESS).
  # 3) If a patient was coded as 0-7, and all radiation were coded as 99-unknown, the radiation 
  #    flag was set to 'unknown' and the reason for no radiation was recoded to 9. 
  # 4) If a patient was coded as 8-9, and all radiation were coded as 00,
  #    the radiation flag was set to 'unknown' and the reason for no radiation was kept as 8 or 9 
  #    (i.e., no change to reasonForNoRadiation).
  step1 %>% 
    group_by(patientid) %>% 
    # note case #4 isn't implemented here because reasonForNoRadiation doesn't change
    mutate(reasonfornoradiation = case_when(
      # If a patient was coded as 1-9, but there was evidence of radiation, they are recoded 
      # to 0. Evidence of radiation = anything other than "No radiation treatment" or "Unknown".
      # But need to account for NA radtreatmentvolume, too
      any(grepl("[1-9]|nomedrio", reasonfornoradiation)) & 
      all(is.na(radtreatmentvolume)) ~ reasonfornoradiation,
      any(grepl("[1-9]|nomedrio", reasonfornoradiation)) & 
      any(!radtreatmentvolume %in% c("No radiation treatment",
      "Unknown") &
      !is.na(radtreatmentvolume)) ~ "0",
      # If a patient was coded as 0, but there was no evidence of radiation (No rad, unknown, all NA), 
      # the radiation flag was set to 'unknown' and the reason for no radiation 
      # is recoded to 9. 
      any(grepl("0", reasonfornoradiation)) & 
        all(radtreatmentvolume %in% c("No radiation treatment",
                                      "Unknown") |
              is.na(radtreatmentvolume)) ~ "9",
      # If a patient was coded as 0-7, and all radiation were coded as 99-unknown, the  
      # radiation was set to 'unknown' and the reason for no radiation was recoded to 9. 
      all(grepl("[0-7]", reasonfornoradiation)) & 
        all(radtreatmentvolume == "Unknown") ~ "9",
      # If a patient was Marked No in Medrio and no radiation, set reason for no radiation to 
      # no radiation; reason unknown"
      all(reasonfornoradiation == "nomedrio") &
        all(radtreatmentvolume == "No radiation treatment" |
              is.na(radtreatmentvolume)) ~ "Radiation therapy was not administered; reason unknown",
      ## otherwise stays the same
      TRUE ~ reasonfornoradiation)) %>% 
    # nomedrio from case 1 in above logic needs to be filled in or else it won't be mapped/ flag with be wrong
    mutate(reasonfornoradiation = ifelse(reasonfornoradiation == "nomedrio", 
                                         "Radiation therapy was not administered; reason unknown",
                                         reasonfornoradiation),
           # catch all in case patients still manage to have multiple reasonfornosurgery
           reasonfornoradiation = ifelse(n_distinct(reasonfornoradiation) > 1,
                                         "9",
                                         reasonfornoradiation),
    ) %>%
    ungroup() %>%
    mutate(radiation_therapy_flag = case_when(grepl("0", reasonfornoradiation) ~ "Yes",
                                              grepl("[1-7]", reasonfornoradiation) ~ "No",
                                              reasonfornoradiation == "Radiation therapy was not administered; reason unknown" ~ "No",
                                              grepl("8|9", reasonfornoradiation) ~ "Unknown",
                                              TRUE ~ "Unknown")) %>%
    rename(reason_for_no_radiation = reasonfornoradiation)
}

# We also need to to the reverse and align surgery data (reasonfornoradiation, rxdateradiation, and rxdateradiationended)
# with the newly aligned reason_for_no_radiation and radiation_therapy_flag
resolve_discordant_radiation_site_date <- function(.data) {
  .data %>% 
    # If a patient was coded as 1-9, but there was evidence of radiation, they were recoded 
    # to 0. We want to remove the 00, 99 records.
    mutate(radtreatmentvolume = case_when(grepl("0", reason_for_no_radiation) & 
                                            (radtreatmentvolume == "No radiation treatment") ~ "Remove",
                                          # If a patient was coded as 0, but there was no evidence radiation (either 00
                                          # 99-unknown, or missing), the radiation flag was set to 'unknown' and the reason for no radiation
                                          # was recoded to 9. 
                                          # If a patient was coded as 0-7, and all radiation were coded as 99-unknown, the  
                                          # radiation was set to 'unknown' and the reason for no radiation was recoded to 9. 
                                          # If a patient was coded as 8-9, and all radiation were coded as 00,
                                          # the radiation flag was set to 'unknown' and the reason for no radiation was kept as 8 or 9.
                                          # In all three cases, we want to set radtreatmentvolume to 99 and NA the dates
                                          grepl("8|9", reason_for_no_radiation) ~ "Unknown",
                                          # NA cases
                                          # Unfortunately missing site is likely for patients with otherwise known dates
                                          grepl("0", reason_for_no_radiation) & 
                                            rxdateradiation_granularity != "NONE" &
                                            is.na(radtreatmentvolume) ~ "Unknown",
                                          # If coded as 0 or recoded to 0, remove NULL radtreatmentvolume records
                                          grepl("0", reason_for_no_radiation) & 
                                            rxdateradiation_granularity == "NONE" &
                                            is.na(radtreatmentvolume) ~ "Remove",
                                          # If coded as 1-7, recode NA radtreatmentvolume records to 00
                                          grepl("[1-7]|Radiation therapy was not administered; reason unknown", reason_for_no_radiation) & 
                                            (is.na(radtreatmentvolume) | 
                                               !radtreatmentvolume %in% c("No radiation treatment"))~ "No radiation treatment",
                                          # If coded as 8/9, recode NA radtreatmentvolume records to 99
                                          grepl("[8-9]", reason_for_no_radiation) & is.na(radtreatmentvolume) ~ "Unknown",
                                          # All others leave as is
                                          TRUE ~ radtreatmentvolume)) %>% 
    # If no evidence of radiation, NA the date fields - using the cleaned up reason_for_no_radiation since 
    # we can have "Unknown" radtreatmentvolumes for patients with radiation therapy
    mutate(across(.cols = c(rxdateradiation, rxdateradiationended),
                  .fns = ~if_else(grepl("[1-9]", reason_for_no_radiation),
                                  NA_Date_,
                                  .,
                                  .)),
           across(.cols = c(rxdateradiation_year, rxdateradiationended_year),
                  .fns = ~if_else(grepl("[1-9]", reason_for_no_radiation),
                                  NA_real_,
                                  .,
                                  .)),
           across(.cols = c(rxdateradiation_granularity, rxdateradiation_granularity),
                  .fns = ~if_else(grepl("[1-9]", reason_for_no_radiation),
                                  "NONE",
                                  .,
                                  .))
    )
}

resolve_discordant_radiation <- function(.data, check_dedup = TRUE, verbose = TRUE) {
  # resolve discordance between radiation reason and radiation treatment volume
  .data <- .data %>% 
    resolve_discordant_radiation_reason()
  # check to see if this resolved patients with multiple, unique reasons even after discordance resolution
  if (check_dedup) {
    .data %>% 
      select(patientid, radiation_therapy_flag, reason_for_no_radiation) %>% 
      distinct() %>% 
      confirm_unique(patientid, return_df = FALSE, error_msg = "Radiation - patients have multiple, unique reasonForNoRadiation still!") %>% 
      invisible()
  }
  # fix radiation sites and dates so they correspond to radiation reason
  .data <- .data %>% 
    resolve_discordant_radiation_site_date()
  if (verbose & .data %>% count_patients() %>% pull(patients) != .data %>% filter(radtreatmentvolume != "Remove") %>% count_patients() %>% pull(patients)) {
    rlang::warn("Patients lost from radiation df when removing mis-algined reasonForNoRadiation and radtreatmentvolume!")
  }
  .data <- .data %>% 
    filter(radtreatmentvolume != "Remove")
  # check for lingering mismatches
  radiation_mismatch <- .data %>% 
    mutate(mismatch = case_when(grepl("[1-7]", reason_for_no_radiation) & 
                                  !radtreatmentvolume %in% c("No radiation treatment") ~ "reason_for_no_radiation in 1-7, known/unknown radiation in radtreatmentvolume",
                                grepl("8|9", reason_for_no_radiation) & 
                                  !grepl("Unknown", radtreatmentvolume) ~ "reason_for_no_radiation in 8 or 9, known/no radiation in radtreatmentvolume",
                                grepl("[1-9]", reason_for_no_radiation) &
                                  !is.na(rxdateradiation) ~ "reason_for_no_radiation in 1-9, non-missing date in rxDateRadiation",
                                grepl("0", reason_for_no_radiation) & 
                                  (grepl("00|99", radtreatmentvolume) | 
                                     is.na(radtreatmentvolume)) ~ "reason_for_no_radiation == 0, evidence of no surgery/unknown surgery in radtreatmentvolume")
    ) %>% 
    filter(!is.na(mismatch))
  if (verbose & dim(radiation_mismatch)[1] > 0) {
    rlang::warn("Patients with mismatched reason_for_no_radiation and radiation fields (radtreatmentvolume, or Date)")
    radiation_mismatch %>% 
      count(sourceschema, mismatch) %>% 
      print(n = 40)
  }
  
  .data
}

clean_radiation_df <- function(.data, patients = NULL, spmd = spmd_con()) {
  # Check patients argument
  if (!is.null(patients)) {
    if (!is.data.frame(patients)) {
      rlang::abort("patients is not a data frame")
    } else {
      patients %>% check_for_column(patientid, return_df = FALSE) %>% invisible()
    }
  }
  
  .data %>% 
    map_oc(reason_for_no_radiation, "reasonForNoRadiation") %>% 
    # in case any registry codes escape the OC mapping, in theory this could only be 98 = "Site specific codes; special"
    coalesce_map_naaccr("radtreatmentvolume", spmd)  %>% 
    rename_with(~ str_replace(., "rxdateradiationended", "end_date")) %>%
    rename_with(~ str_replace(., "rxdateradiation", "start_date")) %>% 
    rename(site = radtreatmentvolume,
           intent = radiationintent,
           radiation_ongoing = ongoingrt) %>% 
    select(patientid, radiation_therapy_flag, reason_for_no_radiation, 
           site, intent, radiation_ongoing,
           start_date, start_date_granularity, start_date_year,
           end_date, end_date_granularity, end_date_year) %>% 
    # move these up into the "resolve" functions? start/end date are handled there
    mutate(site = case_when(!is.na(site) ~ site,
                            radiation_therapy_flag == "Yes" ~ "Unknown", 
                            radiation_therapy_flag == "No" ~ "No radiation treatment",
                            radiation_therapy_flag == "Unknown" ~ "Unknown"),
           intent = ifelse(radiation_therapy_flag == "Yes", 
                           coalesce(intent, "Unknown"),
                           NA_character_),
           radiation_ongoing = ifelse(radiation_therapy_flag == "Yes", 
                                      coalesce(radiation_ongoing, "Unknown / Not stated"),
                                      NA_character_)) %>% 
    distinct() %>% 
    # set NAs to unknown - easiest to do this now upstream of the nest() and most definitive surgery 
    # definitions
    {if (!is.null(patients)) full_join(., patients %>% select(patientid), by = "patientid") else .} %>% 
    # for primary site surgery, we back-filled the flag and the reason, and then left other fields NA to 
    # differentiate between unknown unknowns and known unknowns
    replace_na(list(radiation_therapy_flag = "Unknown",
                    reason_for_no_radiation = "It is unknown if radiation therapy was recommended or administered"))
    #                 site = "Unknown",
    #                 start_date_granularity = "NONE",
    #                 end_date_granularity = "NONE"
    # )) 
}

nest_radiation_df <- function(.data) {
  .data %>% 
    select(patientid,
           radiation_therapy_flag,
           reason_for_no_radiation,
           radiation_site = site, 
           radiation_start_date= start_date,
           radiation_start_date_granularity= start_date_granularity,
           radiation_start_date_year= start_date_year,
           radiation_end_date= end_date,
           radiation_end_date_granularity= end_date_granularity,
           radiation_end_date_year= end_date_year,
           radiation_intent= intent, 
           radiation_ongoing) %>% 
    nest(radiation_therapy = c(-patientid, -radiation_therapy_flag, -reason_for_no_radiation))
}
