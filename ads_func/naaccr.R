# Supported functions ----

#' Search for NAACCR Columns and Metadata - NOT IN SYHELPR
#' 
#' Searches for NAACCR columns containing the "pattern" and returns their rows from ca.map_naaccr_metadata.
#'
#' @param pattern a character vector containing the regex pattern to search for
#' @param spmd a DBI connection to the spmd database. Defaults to clone.

#' @return ca.map_naaccr_metadata filtered to rows where naaccr_item_id matches the pattern
#' @export
#'
#' @examples 
#' search_naaccr_cols("surg")
#' 
search_naaccr_metadata <- function(pattern, spmd = spmd_con()) {
  # this should maybe eventually replace search_naaccr_cols, or be added to it
  tbl(spmd, in_schema("ca", "map_naaccr_metadata")) %>% 
    filter(grepl(pattern, naaccr_item_id, ignore.case = TRUE))
}

# Wrapper for map_naaccr - NOT IN SYHELPR
coalesce_map_naaccr <- function(.data, column, spmd = spmd_con()) {
  column_description <- rlang::sym(paste0(column,"_description"))
  column_name <- rlang::sym(column)
  .data %>% 
    map_naaccr(column, spmd) %>%
    mutate(!!column_name := coalesce(!!column_description, !!column_name)) %>% 
    select(-!!column_description)
}

# Not yet supported functions ----

# Returns individual cancer types available in various registry.<customer>_naaccr schemas
# Convert to using ca.registry_naaccr? 
get_naaccr_cancer_types <-
  function(cancer_type = '[a-z]') {
    results <- list()
    orgs <- tbl(spmd, in_schema('mdr', 'patient')) %>%
      distinct(sourcename) %>%
      collect %>%
      pull(sourcename)
    
    for (i in 1:length(orgs)) {
      table <- paste0(orgs[[i]], '_naaccr')
      
      results[[orgs[[i]]]] <-
        tbl(spmd, in_schema('registry', table)) %>%
        filter(grepl(cancer_type, stagingalgorithmschema, ignore.case = T)) %>%
        distinct(cancer_type = stagingalgorithmschema) %>%
        collect
    }
    
    return(results %>% bind_rows %>% distinct)
  }
# Returns all columns and all rows from individual registry.<customer>_naaccr schemas matching a given 
# cancer_type
get_naaccr_data <-
  function(cancer_type,
           cols = '*',
           list_cancer_types = FALSE) {
    results <- list()
    orgs <- tbl(spmd, in_schema('mdr', 'patient')) %>%
      distinct(sourcename) %>%
      collect %>%
      pull(sourcename)
    
    for (i in 1:length(orgs)) {
      table <- paste0(orgs[[i]], '_naaccr')
      
      results[[orgs[[i]]]] <-
        tbl(spmd, in_schema('registry', table)) %>%
        filter(grepl(cancer_type, stagingalgorithmschema, ignore.case = T)) %>%
        select(mrn = clientid, matches(cols)) %>%
        collect %>%
        mutate(source = table, mrn = as.character(as.numeric(mrn)))
    }
    
    
    return(results)
  }
# Similar to the above but uses ca.registry_naaccr
get_naaccr_cancer_patients <- function(cancertype) {
  tbl(spmd_con(), in_schema("ca","registry_naaccr")) %>% 
    filter(grepl(cancertype, stagingalgorithmschema, ignore.case = TRUE))
}
