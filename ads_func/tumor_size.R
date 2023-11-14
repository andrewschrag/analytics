## Tumor size function

# Retrieve tumor size data from registry ----
### Clinical, pathologic, and summary tumor size (NAACCR items )

get_tumorsize_registry <- function(cohort_query, cancer_type) {
  get_registry(cohort_query,
               cancer_type,
               c("tumorsizeclinical", "tumorsizepathologic", "tumorsizesummary")) %>%
    select(-stagingalgorithmschema, -dateofdiagnosis, -primarysite, -dateofdiagnosis_year) %>%
    distinct()
}

# clean tumor size ----
### Split into a categorical and continuous variables as per
### https://docs.google.com/document/d/11n4WHKUgruOd4YPP0iZx2662IOfG1vB7dVWECdtIU7k/edit?usp=sharing

clean_tumorsize_registry <- function(.data, spmd = spmd_con()) {
  
  .data %>%
    mutate(
      across(.cols = c(tumorsizeclinical, tumorsizepathologic, tumorsizesummary),  # converts to numeric so we can handle measurements
             .fns = ~as.numeric(.)),
      across(.cols = c(tumorsizeclinical, tumorsizepathologic, tumorsizesummary),  # retain exact measurements as their own variables (_number)
             .fns = ~case_when(. <= 989 ~ .,
                               (. > 990 & . < 998 ) ~ 989, # fix some cleaning issues - these mistakes would be measurements larger than 988
                               TRUE ~ NA_real_),
             .names = "{col}_number"),
      across(.cols = c(tumorsizeclinical, tumorsizepathologic, tumorsizesummary),  # assign slightly tidied versions of the NAACCR map 
             .fns = ~case_when(is.na(.) ~ "Unknown; size not stated; Not documented in patient record; Not applicable",
                               . == 999 ~ "Unknown; size not stated; Not documented in patient record; Not applicable",
                               . == 0 ~ "No mass/tumor found",
                               . > 0 & # this cleans messy data
                                 . <= 001 ~ "1 mm or described as less than 1 mm",
                               . > 001 &
                                 . <= 988 ~ "2 mm to 988 mm",
                               . == 989 ~ "989 millimeters or larger",
                               (. > 990 & . < 998 ) ~"989 millimeters or larger",
                               . == 998 ~ "Diffuse",
                               TRUE ~ as.character(.)) # all other cases will map to the available NAACCR map
      )
    ) %>% 
    coalesce_map_naaccr("tumorsizeclinical", spmd) %>%  
    coalesce_map_naaccr("tumorsizepathologic",spmd) %>% 
    coalesce_map_naaccr("tumorsizesummary", spmd)
}

# wrangle/format the patients tumor size ----
### For small sub set of multiple records, prioritization of (1) most complete and (2) largest summary > path > clinical
### and Sorting of categorical consideration in order of largest to smallest
### per discussion https://syapse.atlassian.net/browse/CA-3971
wrangle_tumorsize <- function(.data, cohort_query){
  patients <- cohort_query %>% 
    select(patientid) %>% 
    collect()
  
  .data %>% 
    rename_with(.cols = starts_with("tumorsize"),     # rename to fit ADS/ best practices naming style
                .fn = ~str_replace(., "tumorsize", 
                                   "tumor_size_")) %>% 
    select(patientid, starts_with("tumor_size")) %>% 
    distinct() %>%
    add_sorting(tumor_size_summary, by = "custom", sort_vector = sort_tumor_size) %>%   # very small tumors won't have numeric value in the _number col
    add_sorting(tumor_size_pathologic, by = "custom", sort_vector = sort_tumor_size) %>%# this sorting allows us to respect that when choosing
    add_sorting(tumor_size_clinical, by = "custom", sort_vector = sort_tumor_size) %>%  # the larger of two records
    group_by(patientid) %>%
    arrange(tumor_size_summary,             # want the largest summary size, first - categorical so catch those without a numeric value
            desc(tumor_size_summary_number), # take the largest exact summary size, if available
            tumor_size_pathologic,          # and so on for pathologic and clinical
            desc(tumor_size_pathologic_number),
            tumor_size_clinical,
            desc(tumor_size_clinical_number)) %>%
    slice_head() %>%
    ungroup() %>% 
    join_replace_na(patients, tumor_size_summary, na_fill = "Unknown; size not stated; Not documented in patient record; Not applicable") %>% 
    join_replace_na(patients, tumor_size_pathologic, na_fill = "Unknown; size not stated; Not documented in patient record; Not applicable") %>% 
    join_replace_na(patients, tumor_size_clinical, na_fill = "Unknown; size not stated; Not documented in patient record; Not applicable") 
  
}

# derive tumor size variables ----
derive_tumor_size <- function(cohort_query, cancer_type, spmd = spmd_con()){
  
  get_tumorsize_registry(cohort_query, cancer_type) %>%
    clean_tumorsize_registry(., spmd) %>% 
    wrangle_tumorsize(., cohort_query)
}
