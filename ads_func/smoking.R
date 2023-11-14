# Retrieve smoking data from ca.smokinghistory data ----
get_smoking_spmd <- function(patients = NULL, spmd = spmd_con()) {

  # Check patients argument
  if (is.data.frame(patients)) {
    patients <- patients %>%
      check_for_column(patientid) %>%
      pull(patientid)
  }
  
  patients <- enquo(patients)
  
  # smoking data
  tbl(spmd, in_schema("ca", "smokinghistory")) %>%
    filter(., patientid %in% !!patients) %>%
    select(-diagnosisdate) %>%
    collect()
}

# clean smoking status ----
## pack years are expected to be non-negative
## years since quit are reasonably expected be under 100, 
### and protects against entries of the _year_ the patients quit
clean_smoking_spmd <- function(.data) {

  if (is.null(.data)) rlang::abort("no data provided")

.data %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>% 
  mutate(smokingstatus = str_remove(smokingstatus, " smoker"),
         packyears = as.numeric(packyears),
         packyears = if_else(packyears <= 0, NA_real_, packyears),
         yearssincequit = as.numeric(yearssincequit),
         yearssincequit = if_else(yearssincequit > 100, NA_real_, yearssincequit),
         packsperday = as.numeric(packsperday),
         smokingyears = as.numeric(smokingyears)
  )
}

# wrangle/format the patients smoking status vars for set population ----
wrangle_smoking <- function(.data, patients){
  .data %>% 
    rename(smoking_status = smokingstatus,
           smoking_years = smokingyears,
           smoking_packs_per_day = packsperday,
           smoking_pack_years = packyears,
           smoking_years_since_quit = yearssincequit) %>%
    join_replace_na(patients, smoking_status)
    
}
 
# derive a patients smoking status and associated variables ----
## fetch and clean smoking variables
derive_smoking <- function(patients = NULL, spmd = spmd_con()){
  
  get_smoking_spmd(patients, spmd) %>%
    clean_smoking_spmd(.) %>% 
    wrangle_smoking(., patients)
}





