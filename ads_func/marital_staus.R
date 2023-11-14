# Creates a de-identified copy of marital_status (marital_status_deid)
deid_marital_status <- function(.data) {
  .data %>% 
    check_for_column(marital_status, 
                     return_df = FALSE) %>% 
    invisible()
  
  .data %>% 
    mutate(marital_status_deid = case_when(marital_status == "Single (never married)" ~ "Single",
                                           marital_status == "Married (including common law)" ~ "Married",
                                           TRUE ~ "Other or Unknown"))
}