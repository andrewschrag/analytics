# State/ZIP/Region Map ----
# region breakdown according to US Census
# https://www2.census.gov/geo/docs/maps-data/maps/reg_div.txt
# library(zipcode)
  state_map  <- list(state.name %>% tibble::enframe(name = NULL) %>% rename(state=value),
                     state.abb %>% tibble::enframe(name = NULL) %>% rename(abbr=value),
                     state.region %>% tibble::enframe(name = NULL) %>% rename(region=value)) %>%
    reduce(cbind) %>%
    as_tibble() %>%
    mutate(region = as.character(region)) %>%
    mutate(state = tolower(state),
           abbr = tolower(abbr),
           region = ifelse(region=="North Central","Midwest",region)) %>%
    left_join(zipcodeR::zip_code_db %>%
                select(zipcode, state) %>%
                as_tibble() %>%
                mutate(state = tolower(state)),
              by = c("abbr" = "state"))

  # Geographic region (at diagnosis)
# Geographic region as recorded on file at diagnosis. Categorized as Northeast, Midwest, South, West, unknown.
# Derive region for all addresses on file for each patient and then resolve conflicting addresses, using
# zipcode or state.

derive_region <- function(data){
  region <- tbl(spmd, in_schema("mdr", "address")) %>%
    select(-provenance) %>%
    inner_join(data, by = "patientid") %>%
    select(patientid, state, postalcode) %>%
    collect() %>%
    mutate(postalcode = zipcodeR::normalize_zip(postalcode),
           state = tolower(state)) %>%
    distinct() %>%
    # first join by state full name
    left_join(state_map %>%
                distinct(state, region),
              by = "state") %>%
    rename(region_state = region) %>%
    # then try by abbreviated name
    left_join(state_map %>%
                distinct(abbr,
                         region),
              by = c("state" = "abbr")) %>%
    rename(region_abbr = region) %>%
    # then try y zipcode
    left_join(state_map %>%
                distinct(zipcode, region), by = c("postalcode" = "zipcode")) %>%
    rename(region_zip = region) %>%
    mutate(region = coalesce(region_zip, region_state, region_abbr)) %>%
    filter(!is.na(region)) %>%
    group_by(patientid) %>%
    summarize(region = ifelse(n_distinct(region) == 1,
                               unique(region),
                               "Conflicting")) %>%
    ungroup()

  region <- data %>%
    select(patientid, sourcename, suborg) %>%
    distinct() %>%
    collect() %>%
    full_join(region, by = "patientid") %>%
    mutate(region = case_when(region != "Conflicting" ~ region,
                              grepl("aah\\-(il|wi)|great lakes|henry ford|illinois|indiana|kansas|michigan|wisconsin", suborg, ignore.case = TRUE) ~ "Midwest",
                              grepl("alabama|bayhealth|florida|maryland|oklahoma|tennessee|texas", suborg, ignore.case = TRUE) ~ "South",
                              grepl("main line|new york", suborg, ignore.case = TRUE) ~ "Northeast",
                              TRUE ~ region)) %>%
    replace_na(list(region = "Unknown")) %>%
    mutate(region = ifelse(region != "Conflicting", region, "Unknown")) %>%
    select(patientid, region) %>%
    distinct()

    return (region)
}
