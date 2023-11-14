# Functions to test for PHI and other sensistive data in the De-ID ADS
# Tests are documented in plain text here: https://docs.google.com/document/d/1ZqIAmYcIQFbck7oxLYW8QbcnMZaqJ5zjDYPjodC_wgY/edit#heading=h.394a61injza1

# Pipeline function to perform all PHI/sensitive data checks
check_deid_explicitly <- function(.data, ignore_columns = c("ethnicity"), marital_status_allowed_values = c("Single", "Married", "Other or Unknown")) {
  check_restricted_columns(.data, 
                           ignore_columns = ignore_columns)
  check_deceased_date(.data)
  check_age_89(.data)
  check_birth_date_year(.data)
  check_marital_status(.data, 
                       marital_status_allowed_values = marital_status_allowed_values)
}

# Checks for restricted (PHI or sensitive data) columns. These are columns that match a regular expression: 
# city|zip|postal|street|sourcename|suborg|sourceschema|npi|^birth_date$|_deid$
check_restricted_columns <- function(.data, ignore_columns = c("ethnicity"), restricted_column_regex = "city|zip|postal|street|sourcename|suborg|sourceschema|npi|^birth_date$|_deid$") {
  restricted_columns <- .data %>% 
    select(matches(restricted_column_regex)) %>% 
    colnames()
  restricted_columns <- restricted_columns[!restricted_columns %in% ignore_columns]
  if (length(restricted_columns) > 0) {
    rlang::abort(glue::glue("Restricted columns {paste0(restricted_columns, collapse = ' & ')} are present!"))
  }
}

# Checks for date of death not ending in "-01" as proxy for testing if de-id rules were applied
check_deceased_date <- function(.data) {
  dod_not_first <- .data %>% 
    filter(day(deceased_date) != 1)
  if (dim(dod_not_first)[1] > 0) {
    rlang::abort("Date of death that does not end in '-01'! Per DE-ID rules, all dates of death should be set two months ahead of actual date and set to the first of the month.")
  }
}

# Checks that all ages <= 89
check_age_89 <- function(.data) {
  age_great_89 <- .data %>% 
    select(matches("^age_")) %>% 
    filter(if_any(everything(), ~ . > 89))
  if (dim(age_great_89)[1] > 0) {
    rlang::abort("Age >89! Per DE-ID rules, age should be <=89.")
  }
}

# Checks that year of birth is not more than 89 years prior to death/today's date
check_birth_date_year <- function(.data) {
  # release_date <- curated$list_versions() %>% tail(1) %>% pull(version) %>% ymd()
  dob_greater_89 <- .data %>% 
    filter(year(coalesce(deceased_date, today())) - birth_date_year > 89)
  if (dim(dob_greater_89)[1] > 0) {
    rlang::abort("Year of birth >89 years prior to death/today's date! Per DE-ID rules, Year of birth should be <=89 prior to death/today\n")
  }
}

# Checks for marital status values that break de-id rules
check_marital_status <- function(.data, marital_status_allowed_values = c("Single", "Married", "Other or Unknown")) {
  if ("marital_status" %in% colnames(.data)) {
    marial_status_not_deid <- .data %>% 
      filter(!marital_status %in% marital_status_allowed_values)
    if (dim(marial_status_not_deid)[1] > 0) {
      rlang::abort("Marital status values not de-identified! Per DE-ID rules, marital status should be one of 'Married', 'Single', or 'Other or Unknown'.\n")
    }
  }
}
