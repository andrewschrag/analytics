# TEMPORARY UNTIL https://syapse.atlassian.net/browse/CA-3892 is completed

#' Get Charlson or Elixhauser comorbidity score
#'
#' Generates a weighted Exlixhauser score, AHRQ, or the Charlson comorbidity score from a data frame.
#'
#' @param .data a data frame
#' @param patient_col the column in the associated data frame containing the patient identifier
#' @param index_date_col the column in associated data frame containing index date (e.g. diagnosis or metastatic date). Default is diagnosisdate.
#' @param remove_icd10 a character vector in regex syntax for icd10 codes to exclude, typically primary site malignancies
#' @param cci_baseline_period dbl; number of years prior to initial diagnosis date that we want CCI for. Default is 1 year.
#' @param infer_cci_outside_window flag to infer CCI to be 0 when the patient has ICD10 codes but none of them qualify for the time window. Default is TRUE.
#' @param algorithm a character vector, either 'charlson' or 'elixhauser'. 
#' See [the Insight Analytics comorbidity calculation wiki](https://syapse.atlassian.net/wiki/spaces/DI/pages/2304540716/Comorbidity+Calculation)
# for more details.
#' @param apply_hierarchy_logic Passed to the [comorbidity::comorbidity()] `assign0` argument. 
#' Apply a hierarchy of comorbidities: should a comorbidity be present in a patient with different 
#' degrees of severity, then the milder form will be assigned a value of 0. By doing this, a type of
#' comorbidity is not counted more than once in each patient.
#' @param weights Character vector passed to [comorbidity::score()]
#' @param spmd connection to spmd, defaults to spmd_con()
#' @param schema name of the schema in SPMD to pull the diagnosis table from
#' @param ... Additional variables passed to [comorbidity::comorbidity()] 
#' 
#' @examples
#' dx = tbl(spmd_con(),in_schema("mdr","diagnosis")) %>% head(10000) %>% collect
#' dx %>% get_comorbidity_score(index_date_col = diagnosisdate, remove_icd10 = c("C50"), algorithm='charlson')
#' dx %>% get_comorbidity_score(index_date_col = diagnosisdate, algorithm='elixhauser')
#' 
#' @export
#'
get_comorbidity_score <- function(.data, patient_col = patientid, index_date_col = diagnosisdate,  remove_icd10 = NULL,cci_baseline_period = 1,
                                  infer_cci_outside_window = TRUE, algorithm = "charlson", apply_hierarchy_logic = TRUE, weights = NULL, 
                                  spmd = spmd_con(), schema = "mdr", ...) {
  patient_col_quo = rlang::enquo(patient_col)
  index_date_col_quo = rlang::enquo(index_date_col)
  
  # This hardcoding is regrettable, but avoids a warning message within the comorbidity package.
  .data = .data %>%
    check_for_column(!!patient_col_quo) %>%
    check_for_column(!!index_date_col_quo) %>% 
    select(
      patientid = !!patient_col_quo,
      diagnosisdate = !!index_date_col_quo
    ) 
  pt_list <- .data %>% distinct(patientid) %>% pull
  
  ## Prepare the data
  # Mapping from name to code
  icd10_map <- tbl(spmd, dbplyr::in_schema("maps", "concept")) %>%
    filter(domain == "diagnosis" & vocabulary == "ICD10CM") %>%
    select(name, code)
  # Get patient diagnoses and apply ICD10 mapping
  diagnoses <- tbl(spmd, dbplyr::in_schema(schema,"diagnosis")) %>%
    filter(patientid %in% pt_list) %>%
    left_join(icd10_map, by = c("diagnosis" = "name")) %>%
    left_join(icd10_map, by = c("diagnosis_rawvalue" = "name")) %>%
    mutate(problem_date = as.Date(coalesce(diagnosisdate, startdts, enddts, noteddate)),
           icd10 = coalesce(code.x, code.y)) %>%
    select(patientid, diagnosis, diagnosis_code, diagnosis_rawvalue, icd10, problem_date) %>%
    collect()
  # To infer CCI = 0 when patient has diagnosis date, but none qualifying for CCI calculation
  patient_has_diagnoses <- diagnoses %>%
    select(patientid) %>%
    distinct() %>%
    mutate(has_diagnosis_data = TRUE)
  # Filter to valid ICD10 code format and exclude primary site malignancies
  diagnoses <- diagnoses %>%
    mutate(icd10 = coalesce(icd10, diagnosis_code, diagnosis_rawvalue)) %>%
    filter(grepl("[A-Z]",icd10) & grepl("\\d",icd10)) %>%
    {`if`(!is.null(remove_icd10),filter(.,!grepl(remove_icd10, icd10)),.)} %>%
    left_join(.data %>%
                select(patientid, diagnosisdate),
              by = "patientid")
  # Filter to ICD10 codes occurring in the specified time prior to initial diagnosis date ("baseline")
  diagnosis_baseline <- diagnoses %>%
    filter(problem_date <= diagnosisdate &
             problem_date >= diagnosisdate - lubridate::years(cci_baseline_period)) %>%
    select(patientid, icd10) %>%
    distinct()
  
  ## Apply the CCI algorithm
  map_val = NULL
  if(algorithm == "charlson") {
    map_val = "charlson_icd10_quan"
  } else if(algorithm == "elixhauser") {
    map_val = "elixhauser_icd10_quan"
  } else {
    rlang::abort(stringr::str_glue("Unrecognized algorithm: '{algorithm}'"))
  }
  
  c = suppressWarnings(comorbidity::comorbidity(
    diagnosis_baseline, 
    id= "patientid", 
    code = "icd10", 
    map = map_val, 
    assign0 = apply_hierarchy_logic
  ))
  
  s = as.numeric(comorbidity::score(
    c,
    weights = weights,
    assign0 = apply_hierarchy_logic
  ))
  
  cmb_by_patient = tibble(patientid = c$patientid, comorbidity_score = s)
  
  if (infer_cci_outside_window) {
    cci_diagnosis = .data %>% 
      select(patientid) %>% 
      left_join(patient_has_diagnoses, by = "patientid") %>% 
      tidyr::replace_na(list(has_diagnosis_data = FALSE)) %>% 
      left_join(cmb_by_patient, by = "patientid") %>% 
      mutate(comorbidity_score = if_else(is.na(comorbidity_score) & has_diagnosis_data, 0, comorbidity_score)) %>% # they have diagnosis data, but none in the window of interest and therefore we assume CCI = 0
      select(-has_diagnosis_data) %>% 
      rename(!!patient_col_quo := patientid)
    
  } else {
    cci_diagnosis = .data %>% 
      select(patientid) %>% 
      left_join(cmb_by_patient, by = "patientid") %>%
      rename(!!patient_col_quo := patientid)
  }
  
  return(cci_diagnosis)
}