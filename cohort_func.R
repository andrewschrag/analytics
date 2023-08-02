## Custom Functions used in Cohort datasets

## Get Cohort Base Data ====
get_cohort <-
  function(cohort,
           cohort_type = 'Curated',
           con = spmd_con(), 
           schema = 'mdr') {
    tbl(con, in_schema("cohorts", "cohort")) %>%
      filter(tolower(name) == tolower(cohort),
             tolower(cohorttype) == tolower(cohort_type)) %>%
      select(cohortid = id, name, cohorttype) %>%
      inner_join(tbl(.$src$con, in_schema("cohorts", "cohortmember"))) %>%
      select(patientid, tumorid, cohorttype, entrancedate, exitdate, addeddts, updateddts) %>%
      inner_join(
        tbl(.$src$con, in_schema(schema, 'tumor')) %>%
          select(
            patientid,
            tumorid = id,
            diagnosisdate = diagnosisdate,
            tumortype,
            primarysite,
            histology,
            sourcename
          )
      ) %>%
      inner_join(tbl(.$src$con, in_schema(schema, 'patient')) %>% rename(patientid = id) %>% select(patientid, suborg)) %>%
      distinct
  }


get_stage <-
  function(cohort_query,
           cols = c(),
           schema = 'mdr',
           one_row = TRUE,
           con = spmd_con()) {
    patients = cohort_query %>% select(patientid) %>% collect %>% distinct %>% pull(patientid)
    
    df <- tbl(con, in_schema(schema, "stage")) %>%
      filter(patientid %in% !!patients) %>%
      #inner_join(cohort_query %>% select(patientid)) %>%
      filter(
        !is.na(stagegroup),!grepl('[0-9]{2}', stagegroup),!stagegroup %in% c('Unknown', 'Not Applicable')
      ) %>%
      select(patientid, stagegroup, designator, stagingdate, any_of(cols)) %>%
      rename(stage = stagegroup) %>%
      collect %>%
      mutate(
        designator = factor(
          designator,
          labels = c('p', 'c'),
          levels = c('p', 'c')
        ),
        stage = ifelse(
          grepl('.*Occult.*|OC|in situ', stage, ignore.case = T),
          '0',
          stage
        ),
        stage_group = str_extract(stage, '(IV|I{1,3}|0)'),
        stage_subgroup = toupper(str_extract(stage, '[ABCa]|is')),
        stage = na_if(paste0(
          replace_na(stage_group, ''), replace_na(stage_subgroup, '')
        ), '')
      ) %>%
      filter(!is.na(stage_group)) %>%
      distinct %>%
      arrange(patientid, designator, desc(stage)) %>%
      {
        if (one_row)
          filter_first_patient_row(.)
        else
          .
      }
    
    df %>%
      # mutate_at(.vars = c('stage', 'stage_group'), ~ replace_na(., 'Unknown')) %>%
      select(patientid,
             stagingdate,
             stage,
             stage_group,
             stage_subgroup,
             any_of(cols))
  }


## Join Index Dates ====
join_index_dates <- function(df,
                             cohort,
                             cols = c('diagnosisdate',
                                      'metastasisdate',
                                      'advanceddate')) {
  cols <- c('patientid', cols)
  df %>% inner_join(cohort %>% select(one_of(cols)) %>% distinct)
}


## Get Antineoplastics using Cohort Query ====
get_antineoplastics <-
  function (cohort = NULL,
            con = spmd_con(),
            silent = TRUE) {
    # if (!is.null(cohort)) {
    #  con = cohort$src$con
    # }
    antineoplastics_regex <-
      dplyr::tbl(con, dbplyr::in_schema("ca",
                                        "map_atc")) %>%
      dplyr::filter(class.level2 %in% c("L01: antineoplastic agents",
                                        "L02: endocrine therapy")) %>% dplyr::pull(all_names_regex) %>%
      unique() %>% paste0(collapse = "|")
    antineoplastics_regex <-
      build_med_regex(antineoplastics_regex, spmd = spmd_con())
    if (!silent) {
      rlang::inform(glue::glue("Using med_regex = '{antineoplastics_regex}'"))
    }
    antineoplastic_ts = str_split(antineoplastics_regex, "\\|")
    antineoplastic_ts = str_c(paste0("(", unlist(antineoplastic_ts),
                                     ")"), collapse = "|")
    antineoplastic_ts = str_replace_all(antineoplastic_ts, " ",
                                        "&")
    antineoplastics = dplyr::tbl(con, dbplyr::in_schema("ca", "medications")) %>%
      dplyr::filter(dbplyr::sql(
        glue::glue(
          "drugsearch_index_col @@ to_tsquery('simple','{antineoplastic_ts}')"
        )
      ))
    if (is.null(cohort)) {
      return(antineoplastics)
    } else if (typeof(cohort) == 'character') {
      antineoplastics = antineoplastics %>% dplyr::filter(patientid %in% !!cohort)
    } else if (typeof(cohort) == 'list') {
      cohort$src$con = antineoplastics$src$con
      antineoplastics = antineoplastics %>%
        inner_join(cohort %>% select(patientid))
    }
    return(antineoplastics)
  }



## Get Metastasis Date ====
get_metastasis_date <- function(cohort_query,
                                return = 'first',
                                con = cohort_query$src$con,
                                schema = 'mdr') {
  patients = cohort_query %>% select(patientid) %>% collect %>% distinct %>% pull(patientid)
  secondary_cancer =
    tbl(cohort_query$src$con, in_schema(schema, 'diagnosis')) %>%
    filter(patientid %in% !!patients) %>%
    # inner_join(cohort_query %>% select(patientid, diagnosisdate)) %>%
    filter(grepl('C7[7-9]', diagnosis),
           diagnosisdate >= diagnosisdate,
           diagnosis != 'C79.81') %>%
    select(patientid, diagnosis, diagnosisdate, diagnosisdate_partial) %>%
    collect %>%
    impute_dates %>%
    rename(metastasisdate = diagnosisdate) %>%
    mutate(icd10 = gsub('\\.', '', diagnosis)) %>%
    left_join(icd10() %>%
                distinct(
                  icd10,
                  mets_site = gsub('Secondary malignant neoplasm of ', '', icd_description)
                )) %>%
    distinct %>%
    group_by(patientid, metastasisdate) %>%
    summarise(mets_site = paste0(mets_site, collapse = ', ')) %>%
    ungroup
  
  metastasis =
    tbl(cohort_query$src$con, in_schema(schema, 'metastasis')) %>%
    distinct(patientid,
             metastasisdate,
             metastasisdate_partial,
             mets_site = bodysite) %>%
    collect %>%
    impute_dates
  
  output = bind_rows(secondary_cancer, metastasis) %>%
    select(-metastasisdate_granularity) %>%
    arrange(patientid, metastasisdate)
  
  if (return == 'first') {
    output = output %>% filter_first_patient_row()
  }
  if (output %>% check_for_dupes() %>% nrow > 0) {
    warning('!!!! Duplicate Records Present - Please Review')
  }
  
  return(output)
}


## Get Recurrence Date ====
get_recurrence_date <-
  function(cohort_query,
           recurrence_type = 'distant',
           return = 'first',
           con = cohort_query$src$con,
           schema = 'mdr') {
    patients = cohort_query %>% select(patientid, diagnosisdate) %>% collect %>% distinct
    recurrence_type = tolower(recurrence_type)
    map_bodysite <- tribble(
      ~ bodysite,
      ~ recurrencetype,
      'Nervous system, NOS',
      'Distant',
      'Lung, NOS',
      'Distant',
      'Bone, NOS',
      'Distant',
      'Pleura, NOS',
      'Distant',
      'Distant',
      'Distant',
      'Skin, NOS',
      'Distant',
      'Local',
      'Local',
      'Lymph node, NOS',
      'Regional',
      'Liver',
      'Distant',
      'Peritoneum, NOS',
      'Distant',
      'Regional',
      'Regional'
    )
    
    recurrence <-
      tbl(cohort_query$src$con, in_schema(schema, 'recurrence')) %>%
      filter(patientid %in% !!patients$patientid) %>%
      # inner_join(cohort_query %>% select(patientid, diagnosisdate)) %>%
      select(patientid, recurrencedate, bodysite) %>%
      collect %>%
      distinct %>%
      inner_join(patients) %>%
      filter(recurrencedate > diagnosisdate) %>%
      left_join(map_bodysite) %>%
      mutate(bodysite = tolower(gsub(', NOS', '',  bodysite))) %>%
      arrange(patientid, recurrencedate)
    
    if (recurrence_type %in% c('distant', 'regional')) {
      output = recurrence %>%
        filter(tolower(recurrencetype) == recurrence_type) %>%
        select(-recurrencetype) %>%
        dplyr::rename_with( ~ paste0(
          "recurrence_",
          recurrence_type,
          '_',
          gsub('recurrence', '', .)
        ),
        -patientid) %>%
        distinct
    } else if (recurrence_type == 'first') {
      output = recurrence %>%
        select(-recurrencetype) %>%
        dplyr::rename_with( ~ paste0(
          "recurrence_",
          recurrence_type,
          '_',
          gsub('recurrence', '', .)
        ),
        -patientid) %>%
        distinct %>%
        filter_first_patient_row()
    }
    
    if (return == 'first') {
      output = output %>% filter_first_patient_row()
    }
    
    if (output %>% check_for_dupes() %>% nrow > 0) {
      warning('!!!! Duplicate Records Present - Please Review')
    }
    
    return(output)
  }




get_demographics <- function (cohort,
                              con = spmd_con('clone'),
                              cols = c(),
                              schema = "mdr")
{
  patients = cohort %>% select(patientid) %>% collect %>% distinct %>% pull(patientid)
  demo_query <-
    dplyr::tbl(con, in_schema(schema, "patient")) %>%
    filter(id %in% !!patients) %>%
    dplyr::select(
      patientid = id,
      birthdate,
      deceaseddate,
      isdeceased,
      sex,
      race,
      ethnicity,
      tidyselect::all_of(cols)
    )
  
  
  
  # if (typeof(cohort) == 'character') {
  #   demo_query = demo_query %>% dplyr::filter(patientid %in% !!cohort)
  # } else if (typeof(cohort) == 'list') {
  #   demo_query$src$con <- cohort$src$con
  #   demo_query = demo_query %>% inner_join(cohort %>% select(patientid))
  # }
  
  demo_query %>%
    collect %>%
    distinct %>%
    dplyr::group_by(patientid) %>%
    dplyr::mutate(
      sex = replace_na(sex, "Unknown"),
      race = replace_na(race, "Not Provided"),
      ethnicity = replace_na(ethnicity, "Unknown"),
      race_ethnicity = case_when(
        ethnicity ==
          "Hispanic/Latino" ~ "Hispanic/Latino",
        race == "Black or African American" ~
          "Black",
        race %in% c("Asian", "Native Hawaiian or Other Pacific Islander") ~
          "Asian/Hawaiian/Pacific Islander",
        race == "White" ~
          "Non-Hispanic White",
        race == "American Indian or Alaska Native" ~
          "American Indian or Alaska Native",
        TRUE ~ "Other or Unknown"
      )
    ) %>%
    ungroup
  
}



get_biomarker <- function (gene = NULL,
                           type = NULL,
                           cols = c(),
                           schema = "mdr")
{
  if (!is.null(gene)) {
    gene_list <- paste("'", paste0(gene, collapse = "','"),
                       "'", sep = "")
  }
  dplyr::tbl(spmd_con(),
             dbplyr::in_schema(schema, "molecularbiomarker")) %>%
    dplyr::select(
      patientid,
      molecularreportid,
      reportdate,
      reportdate_partial,
      genes,
      biomarkertype,
      call,
      biomarkername,
      tidyselect::all_of(cols)
    ) %>% {
      if (!is.null(type))
        dplyr::filter(., biomarkertype %in% type)
      else
        .
    } %>% {
      if (!is.null(gene))
        dplyr::filter(., dbplyr::sql(stringr::str_glue(
          "ARRAY[{gene_list}] && genes::text[]"
        )))
      else
        .
    } %>% dplyr::collect() %>% dplyr::mutate(call = toupper(call)) %>%
    dplyr::mutate(
      reportdate_partial = ifelse(
        reportdate_partial ==
          "" |
          reportdate_partial == "UNKN-UN-UN",
        as.character(NA),
        reportdate_partial
      )
    ) %>% dplyr::mutate(reportdate_year = ifelse(
      !is.na(reportdate),
      suppressWarnings(lubridate::year(reportdate)),
      as.numeric(stringr::str_extract(reportdate_partial,
                                      "^\\d{4}"))
    )) %>% dplyr::select(
      patientid,
      molecularreportid,
      reportdate,
      reportdate_partial,
      reportdate_year,
      genes,
      biomarkertype,
      call,
      biomarkername,
      tidyselect::all_of(cols)
    )
}




build_cohort <-
  function(cohort_name,
           cohort_type = c('Structured', 'Curated', 'Enriched'),
           followup = FALSE,
           age_breaks = c(0, 50, 64, 74, Inf),
           age_labels = c('<50', '50-64', '65-74', '75+'),
           con = spmd_con('prod'),
           write_table = F,
           ...) {
    tictoc::tic('====>> build_cohort() run time')
    message(glue::glue('{syhelpr::timestamp()} - building {cohort_name} cohort'))
    
    output = list()
    .cohort_type = sort(factor(
      cohort_type,
      levels = c('Structured', 'Curated', 'Manual Abstraction', 'Enriched')
    ))
    
    
    index_cols =
      c('patientid',
        'diagnosisdate',
        'metastasisdate',
        'advanceddate')
    
    cohorttypes = list_cohorts(cohort_name, con) %>% pull(cohorttype)
    if ('Structured' %in% cohorttypes &
        'Structured' %in% .cohort_type) {
      cohort_base = 'Structured'
    } else {
      cohort_base = 'Curated'
    }
    
    #### Cancer Cohort ====
    cohorts = list()
    output$queries  = list()
    for (cohort_type in .cohort_type) {
      cohort.query = get_cohort(cohort_name, cohort_type, 'mdr', con = con)
      output$queries[[trimws(cohort_type)]] = cohort.query
      .cohort = cohort.query %>%
        select(-c(cohorttype)) %>%
        collect %>%
        distinct %>%
        mutate({
          {
            cohort_type
          }
        } := TRUE) %>%
        setNames(tolower(names(.))) %>%
        arrange(patientid, diagnosisdate) %>%
        filter_first_patient_row()
      if (nrow(.cohort) > 0)
        cohorts[[trimws(cohort_type)]] = .cohort
    }
    
    cohort_query = output$queries[[cohort_base]] %>% select(patientid, diagnosisdate)
    
    output$data$cancer = cohorts %>% reduce(full_join) %>% mutate_if(is.logical, ~ replace_na(., FALSE))
    output$data$demographics =
      get_demographics(cohort_query,
                       cols = c('birthdate')) %>%
      select(-deceaseddate)
    output$data$stage =
      get_stage(cohort_query, cols = c('t', 'n', 'm'))
    output$data$metastasis =
      get_metastasis_date(cohort_query, 'first')
    output$data$recurrence_regional =
      get_recurrence_date(cohort_query, 'regional', 'first')
    output$data$recurrence_distant =
      get_recurrence_date(cohort_query, 'distant', 'first')
    
    if (followup) {
      output$data$last_contact = get_last_contact(cohort_query)
    }
    
    # output %>% reduce(left_join) %>% check_for_dupes()
    output$cohort = output$data %>%
      reduce(left_join) %>%
      mutate(
        diagnosisyear = year(diagnosisdate),
        age_at_dx = age_at(diagnosisdate, birthdate),
        age_at_dx_group := cut(
          age_at_dx,
          breaks = age_breaks,
          right = F,
          labels = age_labels
        ),
        adult = age_at_dx >= 18,
        stage_category = case_when(
          grepl('IV', stage) ~ 'Metastatic',
          grepl('III[B-C]', stage) ~ 'Locally Advanced',
          is.na(stage) ~ 'Unknown',
          stage_group == '0' ~ 'In Situ',
          TRUE ~ 'Early Stage'
        ),
        metastatic_at_dx = stage_category == 'Metastatic',
        metastatic = (
          metastatic_at_dx |
            patientid %in% c(
              output$data$metastasis$patientid,
              output$data$recurrence_distant$patientid
            )
        ),
        metastasisdate = if_else(
          grepl('IV', stage_group),
          diagnosisdate,
          pmin(metastasisdate, recurrence_distant_date, na.rm = TRUE)
        ),
        index_mets_date = pmin(metastasisdate, recurrence_distant_date, na.rm = TRUE),
        advanced = metastatic |
          stage_category %in% c('Metastatic', 'Locally Advanced'),
        advanceddate = case_when(
          grepl('III[B-C]|IV', stage_group) ~ diagnosisdate,
          metastatic ~ metastasisdate,
          TRUE ~ NA_Date_
        )
      ) %>%
      rename(first_mets_site = mets_site) %>%
      mutate_if(is.logical, ~ replace_na(., FALSE)) %>%
      mutate_at(c('t', 'n', 'm'), ~ gsub(' NOS', '', na_if(., 'Not Applicable')))
    
    if (followup) {
      output$cohort = output$cohort %>%
        mutate(
          last_contact_date = coalesce(deceaseddate, last_contact_date),
          followup = round(interval(diagnosisdate, last_contact_date) / months(1), 1),
          followup = ifelse(followup < 0, NA_real_, followup)
        )
    }
    
    if (write_table) {
      message(glue::glue("Writing table to ca.cohort_{tolower(cohort_name)}..."))
      dplyr::copy_to(
        spmd_con('prod', write = T),
        output$cohort,
        in_schema('ca', glue::glue('cohort {cohort_name}')),
        indexes = c('patientid'),
        temporary = FALSE,
        overwrite = TRUE,
        nrows = 1000000
      )
      # DBI::dbWriteTable(
      #   conn = spmd_con('prod', write = T),
      #   name = in_schema('ca', glue::glue('cohort {cohort_name}')),
      #   value = as.data.frame(output$cohort),
      #   overwrite = overwrite
      # )
      message("Writing table complete.")
    }
    
    message(glue::glue('{syhelpr::timestamp()} - build_cohort() complete'))
    tictoc::toc()
    return(output[c('cohort', 'queries')])
  }



build_structured_cohort <-
  function(cohort_query,
           followup = FALSE,
           age_breaks = c(0, 50, 64, 74, Inf),
           age_labels = c('<50', '50-64', '65-74', '75+'),
           con = spmd_con(),
           schema = 'mdr') {
    output = list()
    
    index_cols =
      c('patientid',
        'diagnosisdate',
        'metastasisdate',
        'advanceddate')
    
    #### Cancer Cohort ====
    tic('==> cohort query')
    cohorts = list()
    output$queries$Structured = cohort_query
    .cohort = cohort_query %>%
      collect %>%
      setNames(tolower(names(.))) %>%
      arrange(patientid, diagnosisdate) %>%
      filter_first_patient_row()
    
    if (nrow(.cohort) > 0)
      cohorts$Structured = .cohort
    
    cohort_query =
      output$queries$Structured %>% select(patientid, diagnosisdate)
    
    output$data$cancer =
      cohorts %>% reduce(left_join) %>% mutate_if(is.logical, ~ replace_na(., FALSE))
    toc()
    
    tic('==> demographics')
    output$data$demographics =
      get_demographics(cohort_query,
                       schema = schema,
                       cols = c('birthdate')) %>%
      select(-deceaseddate)
    toc()
    
    tic('==> stage')
    output$data$stage =
      get_stage(cohort_query,
                schema = schema,
                cols = c('t', 'n', 'm'))
    toc()
    
    tic('==> metastasis')
    output$data$metastasis =
      get_metastasis_date(cohort_query, schema = schema, 'first')
    toc()
    
    tic('==> recurrence regional')
    output$data$recurrence_regional =
      get_recurrence_date(cohort_query, 'regional', 'first', schema = schema)
    toc()
    
    tic('==> recurrence distant')
    output$data$recurrence_distant =
      get_recurrence_date(cohort_query, 'distant', 'first', schema = schema)
    toc()
    
    if (followup) {
      output$data$last_contact = get_last_contact(cohort_query)
    }
    
    # output %>% reduce(left_join) %>% check_for_dupes()
    output$cohort = output$data %>%
      reduce(left_join) %>%
      convert_dates %>%
      mutate(
        age_at_dx = age_at(diagnosisdate, birthdate),
        age_at_dx_group := cut(
          age_at_dx,
          breaks = age_breaks,
          right = F,
          labels = age_labels
        ),
        adult = age_at_dx >= 18,
        stage_category = case_when(
          grepl('IV', stage) ~ 'Metastatic',
          grepl('III[B-C]', stage) ~ 'Locally Advanced',
          is.na(stage) ~ 'Unknown',
          stage_group == '0' ~ 'In Situ',
          TRUE ~ 'Early Stage'
        ),
        metastatic_at_dx = stage_category == 'Metastatic',
        metastatic = (
          metastatic_at_dx |
            patientid %in% c(
              output$data$metastasis$patientid,
              output$data$recurrence_distant$patientid
            )
        ),
        metastasisdate = if_else(
          grepl('IV', stage_group),
          diagnosisdate,
          pmin(metastasisdate, recurrence_distant_date, na.rm = TRUE)
        ),
        index_mets_date = pmin(metastasisdate, recurrence_distant_date, na.rm = TRUE),
        advanced = metastatic |
          stage_category %in% c('Metastatic', 'Locally Advanced'),
        advanceddate = case_when(
          grepl('III[B-C]|IV', stage_group) ~ diagnosisdate,
          metastatic ~ metastasisdate,
          TRUE ~ NA_Date_
        )
      ) %>%
      rename(first_mets_site = mets_site) %>%
      mutate_if(is.logical, ~ replace_na(., FALSE)) %>%
      mutate_at(c('t', 'n', 'm'), ~ gsub(' NOS', '', na_if(., 'Not Applicable')))
    
    if (followup) {
      output$cohort = output$cohort %>%
        mutate(
          last_contact_date = coalesce(deceaseddate, last_contact_date),
          followup = round(interval(diagnosisdate, last_contact_date) / months(1), 1)
        )
    }
    
    return(output[c('cohort', 'queries')])
  }



calc_tx_setting <-
  function(meds_df,
           cohort_df,
           surgery_df,
           radiation_df,
           reasonfornosurgery_df,
           reasonfornorad_df) {
    meds_df %>%
      join_index_dates(cohort_df,
                       cols = c('patientid',
                                'diagnosisdate',
                                'advanceddate')) %>%
      filter(startdate >= diagnosisdate | is.na.Date(startdate)) %>%
      left_join(surgery_df) %>%
      left_join(radiation_df) %>%
      mutate_if(is.logical, ~ replace_na(., FALSE)) %>%
      left_join(reasonfornosurgery_df) %>%
      left_join(reasonfornorad_df) %>%
      mutate(
        startdate_ = pmin(startdate, enddate - days(7), na.rm = T),
        tx_setting = case_when(
          is.na.Date(startdate_) ~ 'Unknown Sys. Therapy Date',!is.na.Date(advanceddate) &
            startdate_ >= advanceddate ~ "Advanced",
          (!has_surgery & !has_rad_therapy) &
            (!reasonfornosurgeryflag & !reasonfornoradiationflag) &
            startdate_ < coalesce(advanceddate, today()) &
            startdate_ >= diagnosisdate ~ "Definitive Systemic Therapy",
          (!has_surgery & !has_rad_therapy) &
            !is.na(advanceddate) &
            startdate_ < advanceddate &
            startdate_ >= diagnosisdate ~ "Definitive Systemic Therapy",
          startdate_ >= rad_startdate & startdate_ < rad_
          (!has_surgery & !has_rad_therapy) &
            (reasonfornosurgeryflag | reasonfornoradiationflag) &
            startdate_ < coalesce(advanceddate, today()) &
            startdate_ >= diagnosisdate ~ "Neoadjuvant",
          startdate_ < pmin(surgery_date, rad_startdate, na.rm = T) ~ "Neoadjuvant",
          startdate_ >= pmin(surgery_date, rad_startdate, na.rm = T) &
            startdate_ < coalesce(advanceddate, today()) ~ "Adjuvant",
          TRUE ~ 'Uncategorized'
        )
      )
  }
