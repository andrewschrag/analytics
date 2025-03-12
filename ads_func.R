## get_next_element_id
get_next_element_id <- function(){
  last_id = tbl(spmd, dbplyr::in_schema('ca', 'ads_data_dictionary')) %>% 
    filter(element_id == max(element_id, na.rm = TRUE)) %>% 
    distinct(element_id) %>% 
    collect %>% 
    pull(element_id)
  
  return(last_id+1)
}

get_next_unique_element_id <- function(elementid){
  last_id = tbl(spmd, dbplyr::in_schema('ca', 'ads_data_dictionary')) %>% 
    filter(element_id == elementid) %>% 
    group_by(element_id) %>% 
    filter(unique_element_id == max(unique_element_id, na.rm = T)) %>% 
    collect %>% 
    distinct(unique_element_id) %>% 
    pull(unique_element_id)
  
  return(last_id+1)
}

get_unique_element_ids <- function(element_id){
  get_ads_dd(element_id = element_id) %>% 
    distinct(unique_element_id) %>% 
    collect %>% 
    pull(unique_element_id)
}

get_latest_record_ids <- function(element_id){
  get_ads_dd(element_id = element_id, version = 'recent') %>% 
    group_by(element_id, cohort) %>% 
    filter(record_id == max(record_id, na.rm = T)) %>% 
    collect %>% 
    pull(record_id)
}

# Get ADS Data Dictionary from CA schema
get_ads_dd <- function(..., version = NULL, collect = TRUE) {
  args = list(...)
  query = tbl(spmd, dbplyr::in_schema('ca', 'ads_data_dictionary'))
  
  if(!is.null(args)){
    for(arg in names(args)){
      variable = as.symbol(arg)
      value = args[[arg]]
      
      if(!is.null(value)){
        query = query %>% filter({{ variable }} %in% c(value))
      }
    }
  }
  
  if (!is.null(version)){
    if(tolower(version) == 'recent'){
      query = query %>% 
        collect %>% 
        group_by(element_id, cohort) %>% 
        filter(record_id == max(record_id, na.rm = T)) %>% 
        ungroup
    } else if (tolower(version) == 'latest'){
      query = query %>% 
      filter(tolower(status) == 'production') %>% 
      group_by(element_id, cohort) %>% 
      filter(record_id == max(record_id)) %>% 
      ungroup
    } else if (tolower(version) == 'planned'){
      query = query %>% 
        filter(!is.na(planned_version)) 
    } else {
      query = query %>% filter(version == version)
    }
  }
  
  if(collect){
    return(query %>% collect)
  } else{
    return(query)
  }
}


# Get ADS Planned Version
get_ads_dd_planned <- function(dd = NULL, ...){
  if(is.null(dd)){
    df = get_ads_dd(
      ...,
      version = 'planned')
  } else {
    df = dd %>% 
      {if(!is.null(status)) filter(., status == status) else .} %>%  
      {if(!is.null(cohort)) filter(., cohort == cohort) else .}
  }
  
  dd = df %>%  
    filter(!is.na(planned_version)) %>% 
    group_by(element_id, cohort) %>% 
    arrange(desc(record_id)) %>% 
    slice(1) %>% 
    ungroup 

  return(dd)
}


## Get Data Dict Comment ----
get_ads_dd_comment <- function(...){
  args = list(...)
  query = tbl(spmd, in_schema('ca', 'ads_data_dictionary_comment'))
  
  if(!is.null(args)){
    for(arg in names(args)){
      variable = as.symbol(arg)
      value = args[[arg]]
      
      query = query %>% filter({{ variable }} == value)
    }
  }
  return(query %>% 
           arrange(desc(comment_dts)) %>% 
           collect)
}


dd_unique_elements <- function(dd){
  .ads_options = ads_options[ads_options != 'tumor agnostic']
  all_cohorts = paste0(ifelse(nchar(.ads_options)>3, str_to_title(.ads_options), toupper(.ads_options)), collapse = ', ') 
  
  output = dd %>% 
    group_by(
      unique_element_id,
      cohort
    ) %>% 
    arrange(desc(record_id)) %>% 
    slice(1) %>% 
    ungroup %>% 
    group_by(unique_element_id) %>% 
    mutate(
      cohorts_label = paste0(ifelse(nchar(cohort)>3, str_to_title(cohort), toupper(cohort)), collapse = ', '),
      cohorts = list(cohort),
      record_ids = list(record_id),
      element_ids = list(element_id)) %>% 
    mutate(
      cohorts_label = ifelse(cohorts_label == all_cohorts, 'Tumor Agnostic', cohorts_label)
    ) %>% 
    group_by(cohorts) %>% 
    slice(1) %>% 
    ungroup %>% 
    distinct(unique_element_id, cohorts_label, record_ids, element_ids, cohorts) %>% 
    arrange(desc(nchar(cohorts)))
  
  output %>%
    rowwise %>%
    filter(length(record_ids)>1) %>%
    arrange(desc(length(record_ids))) %>%
    bind_rows(
      output %>%
        rowwise %>%
        filter(length(record_ids)==1) %>%
        arrange(cohorts_label)
    ) %>% 
    ungroup
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


get_ads_data <- function(cohort, variables, spmd = spmd_con()){
  .cohort = tolower(cohort)
  #.variables = c(strsplit(variables, split = ', {0,1}'))[[1]]
  
  ads_type_map = list('enriched'   = syhelpr::list_ads(type='enriched'),
                      'essentials' = syhelpr::list_ads(type='essentials'))
  
  if(.cohort %in% ads_type_map$enriched){
    ads_table_name = glue::glue('ads_{.cohort}_enriched')
  } else {
    ads_table_name = glue::glue('ads_{.cohort}_essentials')
  }
  
  tbl(spmd, dbplyr::in_schema('ca', ads_table_name)) %>%
    select(patientid, sourcename, suborg, any_of(c(variables))) %>% 
    collect
}
