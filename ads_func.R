get_ads_dd <- function(cohort = NULL) {
  query = tbl(spmd_con('prod'), in_schema('ca', 'ads_data_dictionary'))
  if (!is.null(cohort)) {
    .output = query %>%
      filter(!!as.symbol(cohort)) %>%
      collect %>%
      mutate(data_element = trimws(gsub('\\\n<.*>', '', data_element)))
  } else {
    .output = query %>% filter(tumor_agnostic) %>% collect
  }
  return(.output)
}



.get_ads_dd<- function(cohort = NULL, ...){
  syapi_call(path = 'get_data_dict', list(cohort = cohort, ...))
}

get_dd_elements <- function(cohort = NULL){
  syapi_call(path = 'get_elements', list(cohort = cohort))
}

# .get_ads_data <- function(cohort = NULL, variables){
#   syapi_call(path = 'get_ads_data', list(cohort = cohort, variables = variables))
# }


.get_ads_data <- function(cohort, variables, spmd = spmd_con()){
  .cohort = tolower(cohort)
  .variables = c(strsplit(variables, split = ', {0,1}'))[[1]]
  
  ads_type_map = list('enriched'   = syhelpr::list_ads(type='enriched'),
                      'essentials' = syhelpr::list_ads(type='essentials'))
  
  if(.cohort %in% ads_type_map$enriched){
    ads_table_name = glue::glue('ads_{.cohort}_enriched')
  } else {
    ads_table_name = glue::glue('ads_{.cohort}_essentials')
  }
  
  tbl(spmd, dbplyr::in_schema('ca', ads_table_name)) %>%
    select(patientid, sourcename, suborg, any_of(c(.variables))) %>% 
    collect
}
