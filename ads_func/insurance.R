## Insurance at diagnosis function

# Retrieve insurance status data from registry ----
get_payeratdx_registry <- function(cohort_query, cancer_type) {
  get_registry(cohort_query,
               cancer_type,
               c("primarypayeratdx")) %>%
    select(-stagingalgorithmschema, -dateofdiagnosis, -primarysite, -dateofdiagnosis_year) %>%
    distinct()
}

# clean insurance status ----
clean_payeratdx_registry <- function(.data, spmd = spmd_con()) {
  
  .data %>%
    mutate(across(where(is.character), ~ na_if(., ""))) %>% 
    coalesce_map_naaccr("primarypayeratdx", spmd) %>% 
    mutate(rank = case_when(primarypayeratdx == "Not insured" ~ 5, 
                            primarypayeratdx == "Insurance status unknown" ~ 4,
                            primarypayeratdx == "Not insured, self-pay" ~ 3,
                            primarypayeratdx == "Insurance, NOS" ~ 2,
                            TRUE ~ 1)) %>% 
    group_by(patientid) %>% 
    slice_min(order_by = rank) %>% 
    mutate(n_payer = n_distinct(primarypayeratdx),
           primarypayeratdx = case_when(n_payer > 1~ "Multiple insurance types reported",
                                        TRUE ~ primarypayeratdx)) %>% 
    ungroup() %>% 
    distinct(patientid, primarypayeratdx)
}

# wrangle/format the patients insurance status ----
## roll up insurance to categories of private, medicare, medicaid, other, not insured, unknown
## if multiple, remove "unknown" or "not insured" for those with other, insured categories
## if still multiple, combine into a "multiple" category
wrangle_payeratdx <- function(.data, patients){
  .data %>% 
   rename(insurance_status_dx = primarypayeratdx)
}

# derive a patients insurance status at diagnosis ----
derive_insurance_status <- function(cohort_query, cancer_type, spmd = spmd_con()){
  
  get_payeratdx_registry(cohort_query, cancer_type) %>%
    clean_payeratdx_registry(., spmd) %>% 
    wrangle_payeratdx(.)
}

