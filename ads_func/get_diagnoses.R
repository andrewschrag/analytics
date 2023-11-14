#' Get diagnosis codes
#' 
#' Functions returns all diagnosis records for a patient after attempting to map them to ICD10 codes
#' whether or not a patient had 
#'
#' @param .data data frame containing `patientid` column of SPMD patientids
#' @param spmd connection to spmd, defaults to spmd_con()
#' @param schema name of the schema in SPMD to pull the diagnosis table from
#' @param cols vector of additional columns to pull from diagnosis table
#'
#' @return a data frame containing all diagnosis records mapped to ICD10 codes for the patients, columns are `patientid`, `problem_date`, and `icd10`
#' @export
#'
#' @examples
#' tbl(spmd_con(), in_schema("mdr", "tumor")) %>% head() %>% collect() %>% get_diagnoses()
#' 
get_diagnoses <- function(.data, spmd = spmd_con(), schema = "mdr", cols = c()) {
  patients <- .data %>%
    check_for_column(patientid) %>%
    pull(patientid) %>% 
    unique()
  # Mapping from name to code
  icd10_map <- tbl(spmd, in_schema("maps", "concept")) %>%
    filter(domain == "diagnosis" & vocabulary == "ICD10CM") %>%
    select(name, code) 
  # get patient diagnoses, apply ICD10 mapping, filter to ICD10 codes and return
  tbl(spmd, in_schema(schema, "diagnosis")) %>%
    filter(patientid %in% patients) %>% 
    left_join(icd10_map, by = c("diagnosis" = "name")) %>%
    left_join(icd10_map, by = c("diagnosis_rawvalue" = "name")) %>%
    mutate(problem_date = as.Date(coalesce(diagnosisdate, startdts, enddts, noteddate)),
           icd10 = coalesce(code.x, code.y, diagnosis_code, diagnosis_rawvalue)) %>%
    select(patientid, problem_date, icd10, all_of(cols)) %>% 
    collect()
  # distinct()? Unsure if we want to make this distinct (e.g., would we ever want to see if a patient had the
  # same ICD10 code multiple times with NA dates?)
}