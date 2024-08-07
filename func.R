## Functions for building datasets and supporting quick reports.
not_all_na <- function(x) any(!is.na(x))
not_any_na <- function(x) all(!is.na(x))

is.na.Date <- function(x)
  is.na(as.POSIXlt(x))


agebreaks <<-
  c(0,
    18,
    30,
    50,
    65,
    79,
    Inf)


agelabels <<-
  c('< 18',
    '18-29',
    '30-49',
    '50-64',
    '65-79',
    '80+')

# Define our Table output for easy coding
make_table <- function (df,
                        ...,
                        add_overall = F,
                        statistic = all_continuous() ~ c(
                          "{N_nonmiss} ({p_nonmiss}%)",
                          "{median} ({p25}, {p75})",
                          "{mean} [{min}, {max}]"),
                        type = all_continuous() ~ "continuous2",
                        sort = c(all_categorical() ~ "frequency"),
                        missing = 'ifany'){
  args = enquos(...)
  .label = ifelse(length(args)>0, as_label(args[[1]]), '')
  gttable = df %>%
    tbl_summary(
      missing = missing,
      sort = sort,
      type = type,
      statistic = statistic,
      ...
    ) %>%
    suppressMessages() %>%
    modify_header(all_stat_cols() ~ "**{level}**<br>N = {prettyNum(n, big.mark = ',')} ({style_percent(p)}%)",
                  label = .label) %>%
    bold_labels()
  if(add_overall) gttable = gttable %>% add_overall()
  gttable %>%
    # as_gt() %>%
    # gt:::as.tags.gt_tbl() %>%
    suppressWarnings()%>%
    suppressMessages()
}

# Make table1 **Deprecated**
make_table1 <-
  function(df,
           cols,
           labels,
           groups = list(""),
           groupspan = c(0),
           digits = 2,
           missing_name_map = 'Unknown',
           ...) {
    .variables = as.list(setNames(labels, cols))
    .df = df %>% apply_factors(cols, ...)

    if ('cohort' %in% names(.df)) {
      strata = c(split(.df, .df$cohort), list(Overall = .df))
    } else {
      strata = c(list(Overall = .df))
    }

    labels = list(variables = .variables,
                  groups = groups)

    table1(
      strata,
      labels,
      groupspan = groupspan,
      render.continuous = render.continuous.syapse,
      render.logical = render.logical.syapse,
      render.missing = render.missing.syapse,
      missing_name_map = missing_name_map,
      digits = digits
    )
  }



## Process Biomarker Calls
process_call <- function(df, call_col = call) {
  .call_col <- enquo(call_col)
  call.map <- tbl(spmd_con('prod'), in_schema('ca', 'map_biomarker_call')) %>%
      collect %>%
      mutate(call = tolower(call))
  df %>%
    mutate(call = ifelse(
      tolower(biomarkertype) != "wild type",
      tolower(call),
      "negative"
    )) %>%
    left_join(call.map, by = "call") %>%
    mutate(call = call_simple
    #factor(
    #  call_simple,
    #  levels = c(
    #    "Positive",
    #    "Equivocal",
    #    "Negative",
    #    "Low",
    #    "Indeterminate",
    #    "QNS",
    #    "Unknown"
    #  )
    #)
    ) %>%
    select(-call_simple)
}





decompress_file <- function(directory, file, .file_cache = FALSE) {
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    # Set working directory for decompression
    # simplifies unzip directory location behavior
    wd <- here::here('data')
    setwd(directory)

    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)

    # uncomment to delete archive once decompressed
    # file.remove(file)

    # Reset working directory
    setwd(wd)
    rm(wd)

    # Test for success criteria
    # change the search depending on
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
}



elapsed_months <- function(start_date, end_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}


## NAACCR Functions ====
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
        select(patientid = clientid,
               mrn = patientidnumber,
               matches(cols)) %>%
        collect %>%
        mutate(source = table, mrn = as.character(as.numeric(mrn)))
    }


    return(results)
  }

# Clean Registry Stage
clean_registry_stage <- function(df, stage_col, prefix = NULL) {
  .stage_col <- enquo(stage_col)
  stage_group_col <- sym(paste0(prefix, '_stage_group'))
  stage_group_num_col <- sym(paste0(prefix, '_stage_group_num'))
  stage_subgroup_col <- sym(paste0(prefix, '_stage_subgroup'))
  stage_rollup_col <- sym(paste0(prefix, '_stage_rollup'))

  .stage_group_col <- enquo(stage_group_col)
  .stage_group_num_col <- enquo(stage_group_num_col)
  .stage_subgroup_col <- enquo(stage_subgroup_col)
  .stage_rollup_col <- enquo(stage_rollup_col)

  df <- df %>%
    mutate(
      !!.stage_col := case_when(
        !!.stage_col == '88' ~ 'Not Applicable',
        !!.stage_col == '99' ~ 'Unknown',
        TRUE ~ !!.stage_col
      )
    ) %>%
    mutate(
      !!.stage_group_num_col := str_extract(!!.stage_col, '(^[0-4]|I{0,3})'),!!.stage_subgroup_col := str_extract(!!.stage_col, '([A-Da-dSs][1-2]{0,1}|is)'),
      stage_group_r = as.character(as.roman(!!.stage_group_num_col)),!!.stage_group_col := coalesce(stage_group_r,!!.stage_group_num_col),!!.stage_rollup_col := paste0(!!.stage_group_col, replace_na(!!.stage_subgroup_col, '')),!!.stage_group_col := ifelse(
        !!.stage_col %in% c('Unknown', 'Not Applicable'),
        !!.stage_col,!!.stage_group_col
      ),!!.stage_rollup_col := ifelse(
        !!.stage_col %in% c('Unknown', 'Not Applicable'),
        !!.stage_col,!!.stage_rollup_col
      )
    )
  return(df)
}


naaccr_med_search <- function(pattern) {
  med.generic_regex <-
    build_med_regex(pattern)
  med.atc <- tbl(spmd, in_schema('ca', 'map_atc')) %>%
    filter(grepl(med.generic_regex, all_names_regex, ignore.case = T)) %>%
    select(generic_name, brand_names_regex, all_names_regex) %>%
    collect
  med_table <- bind_rows(
    med.atc %>% distinct(match_name = generic_name, generic_name),
    med.atc %>% distinct(match_name = brand_names_regex, generic_name)
  )

  naaccr.rx <-
    bind_rows(
      search_naaccr('rxtextchemo', cols = c('stagingalgorithmschema')) %>% collect,
      search_naaccr('rxtexthormone', cols = c('stagingalgorithmschema')) %>% collect,
      search_naaccr('rxtextbrm', cols = c('stagingalgorithmschema')) %>% collect,
      search_naaccr('rxtextother', cols = c('stagingalgorithmschema')) %>% collect
    ) %>%
    mutate(rxtext = tolower(
      coalesce(rxtextchemo, rxtexthormone, rxtextbrm, rxtextother)
    )) %>%
    filter(!is.na(rxtext)) %>%
    collect %>%
    tidytext::unnest_tokens(word, rxtext) %>%
    filter(word != 'none') %>%
    anti_join(tidytext::stop_words) %>%
    stringdist_inner_join(
      med_table,
      by = c('word' = 'match_name'),
      max_dist = 3,
      distance_col = "dist"
    ) %>%
    group_by(patientid,
             stagingalgorithmschema,
             sourcename,
             word) %>%
    filter(dist == min(dist, na.rm = T)) %>%
    ungroup %>%
    distinct(patientid,
             cancer_type = stagingalgorithmschema,
             sourcename,
             generic_name) %>%
    group_by(patientid, cancer_type, sourcename) %>%
    arrange(patientid, generic_name) %>%
    summarise(regimen = paste0(generic_name, collapse = ', ')) %>%
    ungroup

  return(naaccr.rx)
}




attrition_table <-
  function(df,
           levels = NULL,
           labels = NULL,
           html_output = T) {
    attr.table <- tibble()

    if (is.null(levels)) {
      levels <- names(df %>% select_if(is.logical))
    } else if (is.list(levels)) {
      if (names(levels)[1] != '')
        levels <- c(list('Top of Funnel'), levels)
      labels <- paste(unlist(levels))
      levels <- names(levels)
    }
    if (is.null(labels) | length(labels) < length(levels)) {
      labels <- c('Top of Funnel', levels)
      levels <- c('', levels)
    }

    for (i in 1:length(labels)) {
      if (i == 1) {
        attr.table <- df %>%
          summarise(cohort = labels[[i]],
                    n = n_distinct(patientid))
      } else {
        df <- df %>% filter(!!sym(levels[[i]]))
        attr.table <-
          bind_rows(attr.table,
                    df %>%
                      summarise(cohort = labels[[i]],
                                n = n_distinct(patientid)))
      }
    }

    attr.table <- attr.table %>%
      mutate(
        n_excluded = lag(n, default = n[1]) - n,
        percent_excluded =  ifelse(row_number() == 1, '-', as_percent(n_excluded / lag(n)))
      )

    if (html_output)
      attr.table <- attr.table %>% sykable
    return(attr.table)
  }


clean_string_dates <- function(df, date_col) {
  .date_col <- enquo(date_col)
  .date_col_p <-
    as.symbol(paste0(deparse(substitute(date_col)), '_p'))

  df %>%
    mutate(!!.date_col_p := ifelse(nchar(!!.date_col) < 8,!!.date_col, NA),!!.date_col := ymd(trimws(ifelse(
      nchar(!!.date_col) == 8,!!.date_col, NA
    )))) %>%
    convert_dates %>%
    impute_dates(keep_partial_date_field = TRUE)

}


filter_first_patient_row <- function(df, .group = NULL) {
  df %>% group_by(patientid) %>% filter(row_number() == 1) %>% ungroup
}


## Label a dataset for Table1 usage
apply_table1_labels <- function(df, cols, labels, factor = T, ...) {
  if (length(cols) != length(labels))
    stop("COLS and LABELS not of equal length")
  var_labels <- setNames(as.list(labels), cols)
  df_labelled <- df %>%
    {
      if (factor)
        apply_factors(., cols, ...)
      else
        .
    } %>%
    labelled::set_variable_labels(.labels = var_labels, .strict = FALSE)
}


apply_factors <-
  function (df,
            cols,
            pat_col = patientid,
            label_logical = T)
  {
    # df = dataset.table5
    # cols = cols.table5b
    #cols = c(cols)
    cols.numeric <<-
      paste0(df %>% select_if(is.numeric) %>% names(.), collapse = '|')
    exclude_regex <<-
      paste0(cols.numeric,
             ifelse(cols.numeric == '', '', '|'),
             'date|age|freq',
             collapse = '|')
    .cols <<- cols[!grepl(exclude_regex, cols)]

    if (length(.cols) == 0) {
      return(df)
    } else {
      for (col in 1:length(.cols)) {
        col_nm <- .cols[col]
        if (grepl("*stage*", col_nm, ignore.case = TRUE)) {
          stage.levels <- c(
            "Unknown",
            "In Situ",
            "0",
            "0is",
            "0a",
            "I",
            "IA",
            "IB",
            "IC",
            "II",
            "IIA",
            "IIB",
            "IIC",
            "III",
            "IIIA",
            "IIIB",
            "IIIC",
            "IV",
            "IVA",
            "IVB",
            "IVC"
          )
          stage.levels <-
            stage.levels[(stage.levels %in% (df %>% distinct(!!rlang::sym(col_nm)) %>% pull(!!rlang::sym(col_nm))))]
          df <-
            df %>% mutate(`:=`(!!col_nm, factor(!!rlang::sym(col_nm),
                                                levels = stage.levels))) %>% arrange(!!rlang::sym(col_nm))
        } else if (grepl("ipss", col_nm, ignore.case = TRUE)) {
          if (grepl("ipssr", col_nm, ignore.case = TRUE)) {
            .levels <- ipssr_levels
          } else {
            .levels <- ipss_levels
          }
          df <- df %>%
            mutate(`:=`(!!col_nm, factor(!!rlang::sym(col_nm), levels = .levels))) %>%
            arrange(!!rlang::sym(col_nm))
        } else {
          df_c <- df %>% count_patients(
            patients = !!enquo(pat_col),!!rlang::sym(col_nm),
            row_total = FALSE
          )
          assign(paste0(col_nm, ".levels"), pull(df[, 1]))
          if (df %>% pull(!!rlang::sym(col_nm)) %>% is.logical) {
            if (label_logical) {
              df <-
                df %>% mutate(`:=`(
                  !!col_nm,
                  factor(
                    !!rlang::sym(col_nm),
                    levels = c(TRUE, FALSE),
                    labels = c('Yes', 'No')
                  )
                )) %>% arrange(!!rlang::sym(col_nm))
            }
          } else {
            df <-
              df %>% mutate(`:=`(!!col_nm, factor(
                !!rlang::sym(col_nm),
                levels = pull(df_c[, 1])
              ))) %>% arrange(!!rlang::sym(col_nm))
          }
        }
      }
    }
    return(df)
  }



closest_to_index <-
  function(df,
           index_date,
           event_date,
           suffix = NA,
           prefix = NA) {
    .event_date = enquo(event_date)
    .index_date = enquo(index_date)
    df <- df %>%
      mutate(days_from_index = abs(day_difference(!!.index_date,!!.event_date))) %>%
      arrange(patientid, days_from_index) %>%
      filter_first_patient_row()

    if (!is.na(suffix)) {
      df <-
        df %>% dplyr::rename_with( ~ paste0(replace_na(prefix, ''), str_remove_all(., paste0(
          prefix, '|', suffix
        )), suffix),-patientid)
    }
    return(df)
  }


list_cohorts <-
  function(pattern = '[a-z]',
           con = spmd_con('prod'),
           realtime = FALSE,
           schema = 'cohorts') {
    if (realtime) {
      tbl(con, in_schema(schema, "cohort")) %>%
        select(cohortid = id, name, cohorttype, updateddts) %>%
        filter(grepl(pattern, name, ignore.case = TRUE)) %>%
        inner_join(tbl(.$src$con, in_schema(schema, "currentcohortmember")) %>% select(patientid, cohortid)) %>%
        select(cohortid, name, cohorttype, updateddts, patientid) %>%
        distinct %>%
        group_by(name, cohorttype, cohortid, updateddts) %>%
        summarise(patients = n_distinct(patientid)) %>%
        ungroup %>%
        mutate(updateddts = as.Date(updateddts)) %>%
        arrange(name, cohorttype) %>%
        collect
    } else {
      tbl(spmd, in_schema('ca', 'cohort_counts')) %>% filter(grepl(pattern, name, ignore.case = TRUE)) %>% collect
    }
  }



find_dupe_cols <- function(df, summary = T) {
  output_ <- df %>%
    check_for_dupes() %>%
    group_by(patientid) %>%
    summarise_all( ~ n_distinct(.)) %>%
    pivot_longer(!patientid, names_to = "field", values_to = "count") %>%
    filter(count > 1)

  if(summary) {
    output <- output_ %>%
      distinct(dupe_records=field) %>%
      arrange(dupe_records) #%>%
      #pull(dupe_records)
  } else {
    output <- output_ %>%
      group_by(patientid) %>%
      summarise(dupe_records = paste0(unique(field), collapse = ', '))
  }

  return(output)
}


refresh_cohort_counts <- function(recreate_mat_view = FALSE) {
  spmd_write = syhelpr::spmd_con("prod", write = TRUE)
  tictoc::tic('====>> run time')
  rlang::inform(as.character(syhelpr::timestamp()))

  if (!recreate_mat_view) {
    syhelpr::run_query("REFRESH MATERIALIZED VIEW ca.cohort_counts;", spmd_write)
  } else {
    message('Dropping Existing Materialized View - ca.cohort_counts')
    run_query("DROP MATERIALIZED VIEW ca.cohort_counts", spmd_write)
    query = "
      CREATE MATERIALIZED VIEW ca.cohort_counts as
      SELECT
        c.name,
        c.cohorttype,
        cm.cohortid,
        (c.updateddts)::date AS updateddts,
        count(DISTINCT cm.patientid) AS patients
      FROM (cohorts.cohort c
            JOIN cohorts.currentcohortmember cm ON ((c.id = cm.cohortid)))
      GROUP BY cm.cohortid, c.name, c.cohorttype, c.updateddts
      ORDER BY c.name, (count(DISTINCT cm.patientid)) DESC"
    syhelpr::run_query(query, spmd_write)

    # indexes
    DBI::dbExecute(spmd_write,
                   glue::glue('GRANT SELECT ON ca.cohort_counts TO PUBLIC'))

    tictoc::toc()
  }
}





build_all_encounters <-
  function(con = spmd_con('replica', write = T),
           OVERWRITE = FALSE) {
    tictoc::tic('====>> run time')
    rlang::inform(as.character(syhelpr::timestamp()))
    if (DBI::dbExistsTable(con, 'all_encounters') & !OVERWRITE) {
      message('Refreshsing Materialized View - ca.all_encounters')
      run_query("REFRESH MATERIALIZED VIEW ca.all_encounters", con)
    }
    else {
      if (OVERWRITE) {
        message('Dropping Existing Materialized View - ca.all_encounters')
        run_query("DROP MATERIALIZED VIEW ca.all_encounters", con)
      }
      message('Creating Materialized View - ca.all_encounters')
      run_query(
        "CREATE MATERIALIZED VIEW ca.all_encounters AS
            SELECT patientid,
            orderdts as encounter_date,
            'mdr_clinicallab_orderdts' as source
            FROM msmdro.clinicallab
            UNION
            SELECT patientid,
            diagnosisdate as encounter_date,
            'mdr_diagnosis_diagnosisdate' as source
            FROM msmdro.diagnosis
            UNION
            SELECT patientid,
            startdts as encounter_date,
            'mdr_encounter_startdts' as source
            FROM msmdro.encounter
            UNION
            SELECT patientid,
            startdts as encounter_date,
            'mdr_imaging_startdts' as source
            FROM msmdro.imaging
            UNION
            SELECT patientid,
            startdts as encounter_date,
            'ca_medications_startdts' as source
            FROM ca.medications
            UNION
            SELECT patientid,
            metastasisdate as encounter_date,
            'mdr_metastasis_metastasisdate' as source
            FROM msmdro.metastasis
            UNION
            SELECT patientid,
            orderdate as encounter_date,
            'mdr_molecularreport_orderdate' as source
            FROM msmdro.molecularreport
            UNION
            SELECT patientid,
            performancedate as encounter_date,
            'mdr_performancestatus_performancedate' as source
            FROM msmdro.performancestatus
            UNION
            SELECT patientid,
            startdts as encounter_date,
            'mdr_procedure_startdts' as source
            FROM msmdro.procedure
            UNION
            SELECT patientid,
            orderdate as encounter_date,
            'mdr_radiation_orderdate' as source
            FROM msmdro.radiation
            UNION
            SELECT patientid,
            recurrencedate as encounter_date,
            'mdr_recurrence_recurrencedate' as source
            FROM msmdro.recurrence
            UNION
            SELECT patientid,
            stagingdate as encounter_date,
            'mdr_stage_stagingdate' as source
            FROM msmdro.stage
            UNION
            SELECT patientid,
            diagnosisdate as encounter_date,
            'mdr_tumor_diagnosisdate' as source
            FROM msmdro.tumor",
        connection = con
      )
    }

    DBI::dbExecute(spmd,
                   glue::glue('GRANT SELECT ON ca.all_encounters TO PUBLIC'))

    tictoc::toc()


    run_query('CREATE INDEX patientid_idx ON ca.all_encounters (patientid);',
              con)
  }






get_ads_deid <- function(cohort='lung', env = 'prod'){
  s3_root = 'syapse-deidentify-emr-data'
  if(!env %in% c('prod')) s3_root = glue::glue('syapse-deidentify-{ env }-emr-data')
  .cohort = tolower(cohort)
  bucket = arrow::s3_bucket(glue::glue('{ s3_root }/deidentify/ads'))
  available_cohort = gsub('\\.parquet', '', bucket$ls())

  if(!.cohort %in% available_cohort){
    stop(glue::glue('{ cohort } ADS not available in { env } environment. ADS currently available for { paste(available_cohort, collapse = ", ")}'))
  }

  ads_path = bucket$path(glue::glue("{ .cohort }.parquet"))
  message(glue::glue('Getting ADS from { s3_root }/{ .cohort }.parquet'))
  ads <- arrow::read_parquet(ads_path) %>% as_tibble
  return(ads)
}


get_health_system_ads <- function(ads){
  hs_data = tbl(spmd_con('prod'), in_schema('mdr', 'patient')) %>%
    filter(id %in% !!ads$patientid) %>%
    select(id, sourcename, suborg) %>%
    mutate(id = tolower(as.character(id))) %>%
    rename(patientid = id, health_system = sourcename) %>%
    collect %>%
    distinct
  output = ads %>% left_join(hs_data, by = 'patientid')
  return(output)
}

#' Unnest Embedded Dataframe in ADS
#' Provides a list of columns available in each table, by cohort
#'
#' @param df dataframe/tibble containing Analytical Dataset
#' @param col name of column where dataframe is embedded
#'
#' @examples
#' unnest_ads(systemic_therapy)
#'
unnest_ads <- function(df, col){
  .col = enquo(col)
  patient_col = intersect(c("patientId", "patientid", "patient_id", "patientID", "patients", "patient", "subjectid"), names(df))[1]
  patient_col = sym(patient_col)

  df %>%
    select({{ patient_col }}, {{ .col }}) %>%
    filter(!is.na({{ .col }})) %>%
    unnest({{ .col }})
}


list_nested_cols <- function(df, include_node = FALSE){
  patient_col = intersect(c("patientId", "patientid", "patient_id", "patientID", "patients", "patient", "subjectid"), names(df))[1]
  patient_col = sym(patient_col)
  list_cols = df %>% select(where(is.list)) %>% colnames
  unnest_list = list()
  for(col in list_cols){
      .col = sym(col)
      unnest_list[[col]] = df %>% 
      select({{patient_col}}, {{ .col }})  %>% 
      filter(!is.na({{ .col }})) %>% 
      head(0) %>% 
      unnest_ads({{ .col }}) %>%
      select(-{{patient_col}}) %>% 
      colnames
    }
    
    output = unnest_list %>% unlist(use.names = FALSE)
    if(include_node) output = c(output, list_cols)
    return(output)
}
