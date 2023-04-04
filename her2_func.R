#' Apply Source Factor
#' Applies factor sorting based on data source
#' mdx > ma > registry > health_system
#' For use in Sandbox
apply_source_factor <- function(df, cols = c(), .levels = c('mdx', 'ma', 'registry', 'health_system')) {
  for(col in cols){
    col = as.symbol(col)
    df <- df %>% 
      mutate({{ col }} := factor({{ col }}, levels = .levels))
  }
  return(df)
}


#' Closest to Index
#' Joins index date of user choice to a dataset and 
#' filters to row where the event_date that is closest to index_date
#' Useful for finding closest event to diagnosis date, for example
closest_to_index <- function(df, index_date, event_date, suffix = NA, prefix = NA) {
  .event_date = enquo(event_date)
  .index_date = enquo(index_date)
  df <- df  #%>% rename_with(~ tolower(.), -patientId) %>%
  mutate(days_from_index = abs(day_difference(!!.index_date, !!.event_date))) %>%
    arrange(patientId, desc(days_from_index)) #%>%
  #filter_first_patient_row()
  
  if(!is.na(suffix)) {
    df <- df %>% dplyr::rename_with(~ paste0(replace_na(prefix, ''), str_remove_all(., paste0(prefix, '|', suffix)), suffix), -patientId)
  }
  return(df)
}



#' Get HER2 IHC
#' Function for returning all HER2 IHC results
get_her2_ihc <-
  function(cohort,
           cols = c(),
           con = spmd_con(),
           schema = 'mdr') {
    dxdate_col_nm = intersect(
      c(
        "diagnosisdate",
        "diagnosisDate",
        "dx_index_date",
        "diagnosis_date"
      ),
      names(cohort)
    )[1]
    cols = c(cols, dxdate_col_nm)
    if (!is.na(dxdate_col_nm))
      dxdate_col = sym(dxdate_col_nm)
    
    data = dplyr::tbl(con, dbplyr::in_schema(schema, "molecularbiomarker"))  %>%
      filter(
        biomarkertype %in% c('Protein Expression'),
        grepl('ERBB2', as.character(genes)),
        patientid %in% !!cohort$patientid
      ) %>%
      left_join(
        dplyr::tbl(.$src$con, dbplyr::in_schema(schema, 'molecularreport')) %>% select(
          patientid,
          molecularreportid = id,
          orderdate,
          reportdate,
          reportsignoffdts
        ),
        by = c('patientid', 'molecularreportid'),
        suffix = c('_mb', '_mr')
      ) %>%
      mutate(report_date = pmax(
        orderdate,
        reportdate_mr,
        reportdate_mb,
        reportsignoffdts,
        na.rm = T
      )) %>%
      collect
    
    
    output = data %>%
      inner_join(cohort %>% distinct(patientid, {
        {
          dxdate_col
        }
      })) %>%
      mutate(
        ihc_result = case_when(
          grepl('1\\+',
                paste(biomarkername, call),
                ignore.case = TRUE) ~ '1+',
          grepl('equivocal|2\\+',
                paste(biomarkername, call),
                ignore.case = TRUE) ~ '2+',
          grepl('positive|3\\+',
                paste(biomarkername, call),
                ignore.case = TRUE) ~ '3+',
          grepl('0', biomarkername, ignore.case = TRUE) |
            grepl('0\\+', intensityscore, ignore.case = TRUE) |
            tolower(call) == '0 negative' |
            tolower(call) == 'no_gain_loss' ~ '0',
          !is.na(intensityscore) ~ intensityscore,
          biomarkername %in% c(
            'HER2 Expression Negative',
            '17: Stated as negative, but score not stated',
            '17: Stated as negative, but result not stated',
            '22: Negative/normal; within normal limits'
          ) | tolower(call) == 'negative' ~ 'Negative, Unknown',
          TRUE ~ NA_character_
        ),
        ihc_result = tolower(ihc_result),
        ihc_call = tolower(call),
        her2_ihc_source = gsub('^[^_]*_', '', sourceschema)
      ) %>%
      filter(!is.na(ihc_result)) %>%
      convert_dates %>%
      select(
        patientid,
        her2_ihc_report_date = report_date,
        ihc_result,
        her2_ihc_source,
        tidyselect::any_of(cols)
      ) %>%
      arrange(patientid, her2_ihc_report_date) %>%
      distinct
    
    if (is.symbol(dxdate_col)) {
      output = output %>%
        mutate(her2_ihc_report_date = if_else(
          is.na(her2_ihc_report_date) &
            grepl('registry', her2_ihc_source),
          {
            {
              dxdate_col
            }
          },
          her2_ihc_report_date
        ))
    }
    
    return(output)
  }



#' Get HER2 FISH Results
#'
#'
get_her2_fish <-
  function(cohort,
           cols = c(),
           con = spmd_con(),
           schema = 'mdr') {
    dxdate_col_nm =
      intersect(
        c(
          "diagnosisdate",
          "diagnosisDate",
          "dx_index_date",
          "diagnosis_date"
        ),
        names(cohort)
      )[1]
    cols = c(cols, dxdate_col_nm)
    if (!is.na(dxdate_col_nm))
      dxdate_col = sym(dxdate_col_nm)
    
    
    data = dplyr::tbl(con, dbplyr::in_schema(schema, "molecularbiomarker"))  %>%
      filter(
        platformtechnology %in% c('ISH', 'FISH', 'CISH'),
        grepl('ERBB2', as.character(genes)),
        patientid %in% !!cohort$patientid
      ) %>%
      left_join(
        dplyr::tbl(.$src$con, dbplyr::in_schema('mdr', 'molecularreport')) %>%
          select(
            patientid,
            molecularreportid = id,
            orderdate,
            reportdate,
            reportsignoffdts
          ),
        by = c('patientid', 'molecularreportid'),
        suffix = c('_mb', '_mr')
      ) %>%
      mutate(report_date = pmax(
        orderdate,
        reportdate_mr,
        reportdate_mb,
        reportsignoffdts,
        na.rm = T
      )) %>%
      collect
    
    
    output = data %>%
      inner_join(cohort %>% distinct(patientid, {
        {
          dxdate_col
        }
      })) %>%
      process_call %>%
      rename(fish_result = call) %>%
      mutate(
        fish_result = tolower(fish_result),
        her2_fish_source = gsub('^[^_]*_', '', sourceschema)
      ) %>%
      filter(!is.na(fish_result)) %>%
      convert_dates %>%
      select(
        patientid,
        her2_fish_report_date = report_date,
        fish_result,
        her2_fish_source,
        tidyselect::any_of(cols)
      ) %>%
      distinct
    
    if (is.symbol(dxdate_col)) {
      output = output %>%
        mutate(her2_fish_report_date = if_else(
          is.na(her2_fish_report_date) &
            grepl('registry', her2_fish_source),
          {
            {
              dxdate_col
            }
          },
          her2_fish_report_date
        ))
    }
    
    return(output)
  }


#' Compute HER2 Status
#'
#'
compute_her2_status <-
  function(cohort,
           her2_logic = .her2_logic,
           cols = c(),
           con = spmd_con(),
           schema = 'mdr',
           delta_limit = 45) {
    her2_levels <- c(
      'HER2-low',
      'HER2-',
      'HER2-, unspecified',
      'HER2+',
      'IHC 2+, FISH equivocal',
      'IHC 2+, FISH unknown',
      'IHC- score unknown, FISH-',
      'IHC- score unknown, FISH equivocal',
      'IHC- score unknown, FISH unknown',
      'IHC unknown, FISH-',
      'IHC unknown, FISH equivocal',
      'IHC Equivocal, No FISH Testing',
      'HER2 Uncertain'
    )
    
    her2_ihc <- cohort %>% get_her2_ihc(cols, con, schema)
    her2_fish <- cohort %>% get_her2_fish(cols, con, schema)
    her2_overall <-
      dplyr::tbl(con, dbplyr::in_schema(schema, "molecularbiomarker"))  %>%
      filter(
        !patientid %in% !!c(her2_ihc$patientid, her2_fish$patientid),
        grepl('ERBB2', as.character(genes)),
        biomarkertype == 'HER2 Overall Status'
      )  %>%
      left_join(
        dplyr::tbl(.$src$con, dbplyr::in_schema(schema, 'molecularreport')) %>% select(
          patientid,
          molecularreportid = id,
          orderdate,
          reportdate,
          reportsignoffdts
        ),
        by = c('patientid', 'molecularreportid'),
        suffix = c('_mb', '_mr')
      ) %>%
      mutate(report_date = pmax(
        orderdate,
        reportdate_mr,
        reportdate_mb,
        reportsignoffdts,
        na.rm = T
      )) %>%
      collect %>%
      process_call() %>%
      filter(call != 'Indeterminate') %>%
      mutate(
        her2_overall = case_when(
          call == 'Positive' ~ 'HER2+',
          call == 'Negative' ~ 'HER2-',
          TRUE ~ 'Unknown'
        )
      )
    
    
    output = cohort  %>%
      distinct(patientid) %>%
      full_join(her2_ihc, multiple = "all") %>%
      full_join(her2_fish, multiple = "all") %>%
      mutate(
        delta = abs(
          day_difference(her2_ihc_report_date, her2_fish_report_date)
        ),
        ihc_result = tolower(ihc_result),
        fish_result = tolower(fish_result)
      ) %>%
      filter(delta <= delta_limit | is.na(delta)) %>%
      bind_rows(her2_ihc) %>%
      group_by(patientid, her2_ihc_report_date) %>%
      arrange(patientid, her2_ihc_report_date, delta) %>%
      slice(1) %>%
      ungroup %>%
      bind_rows(her2_fish) %>%
      group_by(patientid, her2_fish_report_date) %>%
      arrange(patientid, her2_fish_report_date, delta) %>%
      slice(1) %>%
      ungroup  %>%
      filter(!is.na(her2_ihc_report_date) |
               !is.na(her2_fish_report_date)) %>%
      her2_logic %>%
      mutate(her2_status = factor(her2_status, levels = her2_levels)) %>%
      apply_source_factor(cols = c('her2_ihc_source', 'her2_fish_source')) %>%
      convert_dates() %>%
      group_by(patientid) %>%
      mutate(
        delta = coalesce(abs(
          day_difference(her2_ihc_report_date, lag(her2_ihc_report_date, 1))
        ),
        abs(
          day_difference(her2_ihc_report_date, lead(her2_ihc_report_date, 1))
        )),
        group_var = ifelse(delta <= delta_limit, delta, row_number()),
        her2_report_date = pmax(her2_ihc_report_date, her2_fish_report_date, na.rm = T),
        her2_days_from_diagnosis = day_difference(diagnosisdate, her2_report_date)
      )  %>%
      group_by(patientid, group_var) %>%
      arrange(patientid,
              her2_status,
              her2_ihc_source,
              her2_days_from_diagnosis) %>%
      slice(1) %>%
      ungroup %>%
      select(patientid, diagnosisdate, everything(), -c(group_var, delta))
    
    return(output)
  }




#' Compute HR Status
#' 
#' 
compute_hr_status <-
  function(cohort,
           hr_logic = .hr_logic,
           cols = c(),
           con = spmd_con(),
           schema = 'mdr',
           delta_limit = 45) {
    
    dxdate_col_nm = intersect(
      c(
        "diagnosisdate",
        "diagnosisDate",
        "dx_index_date",
        "diagnosis_date"
      ),
      names(cohort)
    )[1]
    cols = c(cols, dxdate_col_nm)
    if (!is.na(dxdate_col_nm))
      dxdate_col = sym(dxdate_col_nm)
    
    
    hr.raw =
      dplyr::tbl(con, dbplyr::in_schema(schema, "molecularbiomarker"))  %>%
      filter(grepl('PGR|ESR1', as.character(genes)),
             patientid %in% !!cohort$patientid) %>%
      left_join(
        dplyr::tbl(.$src$con, dbplyr::in_schema(schema, 'molecularreport')) %>% select(
          patientid,
          molecularreportid = id,
          orderdate,
          reportdate,
          reportsignoffdts
        ),
        by = c('patientid', 'molecularreportid'),
        suffix = c('_mb', '_mr')
      ) %>%
      mutate(report_date = pmax(
        orderdate,
        reportdate_mr,
        reportdate_mb,
        reportsignoffdts,
        na.rm = T
      )) %>%
      collect %>%
      process_call %>%
      mutate(status = tolower(call),
             source = gsub('^[^_]*_', '', sourceschema)) %>%
      select(
        patientid,
        biomarkername,
        genes,
        platformtechnology,
        status,
        status_flag = call_flag,
        report_date,
        source,
        tidyselect::any_of(cols)
      ) %>% 
      inner_join(cohort %>% select(patientid, {{ dxdate_col }}), multiple = 'all') %>% 
      convert_dates
    
    
    pr <- hr.raw %>%
      filter(grepl('PGR', as.character(genes)), )  %>%
      mutate(status = factor(status, levels = c('positive', 'negative'))) %>%
      arrange(patientid, status, desc(report_date)) %>%
      select(
        patientid,
        pr_positive = status_flag,
        pr_report_date = report_date,
        pr_source = source,
        tidyselect::any_of(cols)
      ) %>%
      distinct
    er <- hr.raw %>%
      filter(grepl('ESR1', as.character(genes)), ) %>%
      mutate(status = factor(status, levels = c('positive', 'negative'))) %>%
      arrange(patientid, status, desc(report_date)) %>%
      select(
        patientid,
        er_positive = status_flag,
        er_report_date = report_date,
        er_source = source,
        tidyselect::any_of(cols)
      ) %>%
      distinct
    
    if (is.symbol(dxdate_col)) {
      er = er %>% mutate(er_report_date = if_else(er_source == 'registry', {
        {
          dxdate_col
        }
      }, er_report_date))
      
      pr = pr %>% mutate(pr_report_date = if_else(pr_source == 'registry', {
        {
          dxdate_col
        }
      }, pr_report_date))
    }
    
    
    
    output = cohort %>%
      distinct(patientid, {{ dxdate_col }}) %>%
      left_join(pr, multiple = "all") %>%
      left_join(er, multiple = "all") %>%
      hr_logic %>%
      mutate(
        hr_report_date = pmax(pr_report_date, er_report_date, na.rm = T),
        hr_source = coalesce(er_source, pr_source),
        delta = abs(day_difference(pr_report_date, er_report_date)),
        hr_days_from_diagnosis = day_difference(diagnosisdate, hr_report_date)
      ) %>%
      filter(
        delta <= delta_limit, 
        !is.na(hr_status)) %>%
      select(
        patientid,
        hr_report_date,
        hr_status,
        hr_positive,
        hr_days_from_diagnosis,
        pr_positive,
        pr_report_date,
        er_positive,
        er_report_date,
        hr_source,
        tidyselect::any_of(cols)
      ) %>%
      distinct
    

    return(output)
  }



#' Computer Breast Subtype
#'
#'
compute_breast_subtype <- function(cohort,
                                   her2_logic = .her2_logic,
                                   hr_logic = .hr_logic,
                                   cols = c(),
                                   con = spmd_con(),
                                   schema = 'mdr',
                                   delta_limit = 45) {
  
  her2 = compute_her2_status(cohort, her2_logic, cols, con, schema, delta_limit)
  hr   = compute_hr_status(cohort, hr_logic, cols, con, schema, delta_limit)
  
  
  output = cohort  %>%
    left_join(her2, multiple = "all") %>%
    left_join(hr, multiple = "all") %>%
    mutate(
      subtype_report_date = pmax(her2_report_date, hr_report_date, na.rm = T),
      delta = abs(day_difference(her2_report_date, hr_report_date))
    ) %>%
    filter(delta <= delta_limit | is.na(delta),
           !is.na(hr_status) | !is.na(her2_status)) %>%
    mutate(
      subtype = case_when(
        hr_positive & her2_status == 'HER2+'  ~ 'HR+ / HER2+',
        hr_positive & her2_status == 'HER2-'  ~ 'HR+ / HER2-',
        hr_positive & her2_status == 'HER2-low' ~ 'HR+ / HER2-low',
        hr_positive & her2_status %in% c('HER2+', 'HER2-', 'HER2-low') ~ 'HR+ / HER2 Uncertain',
        !hr_positive & her2_status == 'HER2+' ~ 'HR- / HER2+',
        !hr_positive & her2_status == 'HER2-' ~ 'TNBC',
        !hr_positive & her2_status == 'HER2-low' ~ 'HR- / HER2-low',
        !hr_positive & !her2_status %in% c('HER2+', 'HER2-', 'HER2-low') ~ 'HR- / HER2 Uncertain',
        TRUE ~ paste0(replace_na(as.character(hr_status), 'HR Unknown'), ' / ', replace_na(as.character(her2_status), 'HER2 Unknown')) #'Unknown'
      )
    ) %>%
    group_by(patientid) %>%
    mutate(
      subtype_days_from_diagnosis = day_difference(diagnosisdate, subtype_report_date),
      subtype_at_diagnosis_flag = subtype_days_from_diagnosis == min(subtype_days_from_diagnosis, na.rm = T)) %>%
    ungroup %>%
    select(
      patientid,
      subtype,
      subtype_report_date,
      subtype_days_from_diagnosis,
      subtype_at_diagnosis_flag,
      her2_status,
      her2_report_date,
      hr_status,
      hr_report_date,
      any_of(cols),
      matches('source'),
      -delta
    ) %>% 
    convert_dates %>%
    suppressWarnings
  
  return(output)
}



#' HER2 Logic
#'
.her2_logic <- function(data) {
  data %>%
    mutate(
      ihc_result = tolower(ihc_result),
      fish_result = tolower(fish_result),
      her2_results = case_when(
        ihc_result == '2+' & fish_result == 'negative' ~ 'IHC 2+, FISH-',
        ihc_result == '2+' &
          fish_result == 'positive' ~ 'IHC 2+, FISH+',
        ihc_result == '2+' &
          is.na(fish_result) ~ 'IHC 2+, FISH equivocal or unknown',
        ihc_result == '2+' &
          fish_result == 'positive' ~ 'IHC 2+, FISH+',
        ihc_result == '2+' &
          fish_result == 'equivocal' ~ 'IHC 2+, FISH equivocal or unknown',
        ihc_result == '1+' ~ 'IHC 1+',
        ihc_result == '0' ~ 'IHC 0',
        ihc_result == '3+' | fish_result == 'positive' ~ 'IHC 3+',
        ihc_result == 'negative, unknown' &
          fish_result == 'negative' ~ 'IHC- score unknown, FISH-',
        ihc_result == 'negative, unknown' &
          is.na(fish_result) ~ 'IHC- score unknown, FISH unknown',
        ihc_result == 'negative, unknown' &
          fish_result == 'equivocal' ~ 'IHC- score unknown, FISH equivocal',
        is.na(ihc_result) &
          fish_result == 'negative' ~ 'HER2-',
        is.na(ihc_result) &
          fish_result == 'equivocal' ~ 'IHC unknown, FISH equivocal',
        is.na(fish_result) &
          ihc_result == '3+' ~ 'IHC 3+, No FISH Testing',
        is.na(fish_result) &
          ihc_result == '2+' ~ 'IHC 2+, No FISH Testing',
        is.na(fish_result) &
          ihc_result == '1+' ~ 'IHC 1+, No FISH Testing',
        is.na(fish_result) &
          ihc_result == 'negative, unknown' ~ 'IHC- score unknown, No FISH Testing'
      ),
      her2_status = case_when(
        ihc_result == '2+' & fish_result == 'negative' ~ 'HER2-low',
        ihc_result == '2+' & fish_result == 'positive' ~ 'HER2+',
        ihc_result == '2+' & is.na(fish_result) ~ 'HER2-low',
        ihc_result == '2+' & fish_result == 'positive' ~ 'HER2+',
        ihc_result == '2+' &
          fish_result == 'equivocal' ~ 'HER2-low',
        ihc_result == '1+' ~ 'HER2-low',
        ihc_result == '0' ~ 'HER2-',
        ihc_result == '3+' | fish_result == 'positive' ~ 'HER2+',
        ihc_result == 'negative, unknown' &
          fish_result == 'negative' ~ 'HER2-',
        ihc_result == 'negative, unknown' &
          is.na(fish_result) ~ 'HER2-',
        ihc_result == 'negative, unknown' &
          fish_result == 'equivocal' ~ 'HER2-',
        is.na(ihc_result) &
          fish_result == 'negative' ~ 'HER2-low',
        is.na(ihc_result) &
          fish_result == 'equivocal' ~ 'HER2 Uncertain',
        is.na(fish_result) & ihc_result == '3+' ~ 'HER2+',
        is.na(fish_result) &
          ihc_result == '2+' ~ 'IHC Equivocal, No FISH Testing',
        is.na(fish_result) & ihc_result == '1+' ~ 'HER2-low',
        is.na(fish_result) &
          ihc_result == 'negative, unknown' ~ 'HER2-, unspecified',
        TRUE ~ her2_results
      )
    )
}



#' HR Logic
#'
.hr_logic <- function(data) {
  data %>%
    mutate(
      hr_positive = pr_positive | er_positive,
      hr_status = ifelse(hr_positive, 'HR+', 'HR-')
    )
}
