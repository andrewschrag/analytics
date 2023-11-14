## HER2, ER, PR, and BRCA Status for Breast



## ASCO/CAP HER2 Status
## written down stream of HER2 data pull and 
## date prep due to difference in mets versus diagnosis
## date window handling. Will think about how to make
## a single function with more time
## Probably a window argument plus argument to allow 
## take the most recent previous specimen/report date if
## no tests in window?
prioritize_her2 <- function(.data, specimen_date = specimen_date, report_date = report_date,
                            her2_status = her2_status, platform_technology = platform_technology,
                            copy_number, intensity_score) {
  specimen_date <- rlang::enquo(specimen_date)
  report_date <- rlang::enquo(report_date)
  her2_status <- rlang::enquo(her2_status)
  copy_number <- rlang::enquo(copy_number)
  intensity_score <- rlang::enquo(intensity_score)
  platform_technology <- rlang::enquo(platform_technology)
  
  .data %>% 
    check_for_column(!!specimen_date, .class = "Date") %>% 
    check_for_column(!!report_date, .class = "Date", return_df = FALSE)
  .data %>% 
    filter(!!platform_technology != "other") %>% 
    ## most recent specimen because if re-biopsied, that's important
    group_by(patientid) %>% 
    mutate(last_specimen = max(!!specimen_date, na.rm = TRUE)) %>% 
    filter(!!specimen_date == last_specimen) %>% 
    ## most recent test, too, because re-testing is also important
    mutate(most_recent_report = max(!!report_date, na.rm = TRUE)) %>% 
    filter(!!report_date == most_recent_report) %>% 
    add_sorting(!!her2_status, by = "custom", sort_vector = c("positive", "negative", "equivocal")) %>% 
    # When a patient has both an IHC and an ISH result, and both are non-equivocal, the ISH result is used.
    # When a patient has either IHC or ISH testing, the available result is used.
    # add_sorting(!!copy_number, by = "numeric") %>%
    # add_sorting(!!intensity_score, by = c("0", "1+", "2+", "3+")) %>%
    # When a patient has only the summary result or registry summary value available, the summary result is used.
    # The summary result does not overwrite an individual test result, if an individual test result is available. 
    # Non-equivocal prioritized over equivocal results
    # Positive results are prioritized over negative results 
    add_sorting(!!platform_technology, by = "custom", sort_vector = c("ish", "ish/ihc", "ihc", "overall summary")) %>% 
    arrange(!!platform_technology,
            !!her2_status
            # ,
            # desc(!!copy_number),
            # desc(!!intensity_score)
            ) %>% 
    slice_head() %>% 
    ungroup()
}
  