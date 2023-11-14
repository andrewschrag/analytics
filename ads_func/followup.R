### Date of last abstraction lives in the dateoflastabstraction and datecompleted across two forms.
### - dateoflastabstraction is populated once a patient is fully refreshes
### - datecomplete is static and is only populated the first time a patient is abstracted
### patients can sometimes be abstracted under two different ids, to protect against this
### we deduplicate alongside ensuring we are using only those patientids with OC cases marked 
### Okay to Release (0) to limit our dataset to complete data only among the enriched cohort.

get_followup_oc <- function(cancer, spmd = spmd_con()) {
  run_query(glue::glue("SELECT
                          patientid,
                          participantid,
                          formdata::json -> 'follow_up' ->> 'group' AS follow_up
                          FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}'
                        AND formdata::json -> 'follow_up' IS NOT NULL"),
            connection = spmd) %>%
    mutate(follow_up = purrr::map(follow_up, jsonlite::fromJSON)) %>%
    unnest_wider(follow_up) %>%
    rename_oc() %>%
    mutate(last_abstraction_date = ymd(dateoflastabstraction),
           dateoflastcontact = ymd(dateoflastcontact)) %>%
    select(patientid, participantid, last_abstraction_date, vitalstatus, dateoflastcontact) %>%
    # prioritize most recent date
    group_by(patientid) %>%
    arrange(is.na(last_abstraction_date),
            desc(last_abstraction_date)) %>%
    slice_head() %>%
    ungroup() %>%
    distinct()
}

get_case_release <- function(cancer, spmd = spmd_con()) {
  run_query(glue::glue("SELECT
                          patientid,
                          participantid,
                          formdata::json -> 'case_release' ->> 'group1' AS release
                          FROM openclinica.formdata
                        WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}'
                        AND formdata::json -> 'case_release' IS NOT NULL"),
                            connection = spmd) %>%
    mutate(release = purrr::map(release, jsonlite::fromJSON)) %>%
    unnest_wider(release) %>%
    rename_oc() %>%
    mutate(datecompleted = ymd(datecompleted)) %>%
    group_by(patientid) %>%
    filter(abstractfinal == 0) %>% 
    arrange(is.na(abstractfinal),
            is.na(datecompleted),
            desc(datecompleted)) %>%
    slice_head() %>%
    ungroup() %>%
    select(patientid, participantid, datecompleted) %>%
    distinct()
}

get_followup <- function(cancer, spmd = spmd_con()) {
  
  followup_ma <- get_followup_oc(cancer = cancer, 
                                 spmd = spmd)
  case_release <- get_case_release(cancer = cancer,
                                   spmd = spmd)
  
  case_release %>%
    # patients not marked for release will still show up
    full_join(followup_ma, by = c("patientid", "participantid")) %>%
    mutate(last_abstraction_date = pmax(last_abstraction_date, datecompleted, na.rm = TRUE)) %>% 
    group_by(patientid) %>% 
    arrange(last_abstraction_date) %>% 
    slice_head() %>% 
    ungroup() %>% 
    select(-vitalstatus)
}
