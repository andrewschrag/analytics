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


format_year <- function(x) style_number(as.numeric(x), digits = 0, big.mark = '') %>% suppressWarnings
convert_dates <- function (df, regex = "(^date_of_|_dt$|dts$|^dob$|date$)", ...) 
{
  if (length(list(...)$key) > 0) 
    regex = list(...)$key
  df %>% mutate_at(vars(matches(regex)), list(~lubridate::date(.)))
}

# Define our Table output for easy coding
# Conveneince wrapper for making table with gtsummary
make_table <-  function (df,
                         ...,
                         add_overall = F,                         
                         statistic = list(
                           all_continuous() ~ c(
                             #"{N_nonmiss} ({p_nonmiss}%)",
                             "{median} ({p25}, {p75})",
                             "{mean} [{min}, {max}]"
                           ),
                           all_categorical() ~ c("{n} ({p}%)")
                         ),
                         type = list(
                           where(is.logical) ~ "dichotomous",
                           all_continuous() ~ "continuous2"
                         ),
                         digits =  list(
                           all_categorical() ~ 0,
                           ends_with('year') ~ list(format_year)
                         ),
                         sort = 'freq',  
                         missing = 'ifany',
                         last = TRUE) {
  args = enquos(...)
  .label = ifelse(length(args) > 1, glue::glue(as_label(args[[1]])), ' ')

  if(sort %in% c('freq', 'frequency') | is.null(sort)){
    .sort = c(all_categorical() ~ "frequency")
  } else if (sort %in% c('alpha', 'alphanumeric')) {
    .sort = c(all_categorical() ~ "alphanumeric")
  } else {
    .sort = sort
  } 
  
  gttable = df %>%
    tbl_summary(
      missing = missing,
      sort = .sort,
      type = type,
      statistic = statistic,
      digits = digits,
      ...
    ) %>%
    suppressMessages() %>%
    modify_header(
      all_stat_cols() ~ "**{level}**<br>N = {prettyNum(n, big.mark = ',')} ({style_percent(p)}%)",
      label = ' '
    ) %>%
    modify_footnote(c(all_stat_cols()) ~ NA) %>%
    bold_labels() 
  
  if (add_overall)
    gttable = gttable %>% add_overall(last=last) %>% modify_footnote(c(all_stat_cols()) ~ NA) 
  
  gttable %>%
    # as_gt() %>%
    # gt:::as.tags.gt_tbl() %>%
    suppressWarnings() %>%
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
process_call <- function(df, call_col = call, spmd = spmd_con('prod', write = T)) {
  .call_col <- enquo(call_col)
  call.map <- tbl(spmd, in_schema('ca', 'map_biomarker_call')) %>%
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



attrition_table <- function(data, filters, strat = NULL) {

  if(!is.list(filters)){
    .filt = filters
    filters <- data %>% colnames %>% tail(length(.filt)) 
    names(filters) <- .filt
  }
  
  patients <- list() 
  patients$all <- data$patientid
  patients$included <- data$patientid
  table <- tibble()
  
  if (strat %>% length < 1) {
    for (filt in 1:length(filters)) {
      filt_pats = data %>%
        filter(patientid %in% patients$included, !!as.symbol(filters[[filt]])) %>%
        select(patientid, sourcename)
      
      table <- bind_rows(
        table,
        tibble(
          'Criteria' = names(filters[filt]),
          'Total' = filt_pats %>%
            count_patients() %>%
            .$patients
        )
      )
      
      patients$included <- filt_pats$patientid
    }
  } else {
    strat_tables <- list()
    strat_var = sym(strat)
    
    #unique(data[[strat]]) %>% print
    
    for (hs in c(sort(unique(data[[strat]])), 'Total')) {
      strat_tables[[hs]] <- table
      for (filt in 1:length(filters)) {
        filt_pats <- data %>%
          filter(patientid %in% patients$included,!!as.symbol(filters[[filt]])) %>%
          select(patientid, sourcename, {{strat_var}})
        
        patients$included <- filt_pats$patientid
        
        .label =  ifelse(hs == 'Total', hs, str_to_upper(hs))
        strat_tables[[hs]] <-  bind_rows(
          strat_tables[[hs]],
          tibble(
            'Criteria' = names(filters[filt]),
            !!as.symbol(.label) := filt_pats %>%
              { `if`(hs != 'Total', filter(., {{strat_var}} == hs), .) } %>%
              count_patients() %>%
              .$patients
          )
        )
      }
      patients$included <- patients$all
    }
    
    table <- strat_tables %>% reduce(left_join)
  }
  
  table %>% 
    mutate(
      `Included(%)` =  ifelse(row_number() == 1, '-', as_percent(`Total` / lag(Total), 1))
      ,`Excluded(n)` = lag(Total, default = Total[1]) - Total
    ) %>%
    gt(rowname_col = 'Criteria') %>% 
    fmt_number(
      decimals = 0,
      sep_mark = ","
    ) %>% 
    tab_options(
      table.width = '85%',
      table_body.border.top.color = '#000',
      table_body.border.top.style = "solid",
      table_body.border.top.width = "3px",
      table.border.top.style = 'hidden',
      table.border.bottom.color = '#fff',
      table.font.size = px(16),
      column_labels.border.bottom.color = '#f5f5f5',
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = "1px",
      column_labels.font.weight = '600',
      data_row.padding = px(35),
      column_labels.padding = px(35),
      heading.padding = px(35)
    ) %>% 
    tab_style(
      style = list(cell_text(weight = "500"), cell_borders(sides = c("right"), style = 'hidden')),
      locations = cells_stub()
    ) %>% 
    tab_style(
      style = cell_text(size = '1.2rem'),
      locations = cells_body()
    ) %>% 
    tab_style(
      style = cell_text(size = '1.2rem', weight = '600'),
      locations = cells_body(columns = Total)
    ) %>% 
    opt_horizontal_padding(scale = 3) %>%
    tab_stubhead(label = label) %>%
    return()
}


## filter_study_pop with attrrition columns
filter_study_pop <- function(.data, cols){
  .data %>% filter_at(vars(all_of(cols)), all_vars(.==TRUE))
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

## ADS + DD functions ----

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



# pivot_cohorts <- function(df){
#   df %>%
#     pivot_longer(all_of(unname(list_ads('enriched')))) %>%
#     filter(value) %>%
#     group_by(across(-c(name, value))) %>%
#     mutate(cohorts = paste0(name, collapse = '; ')) %>%
#     ungroup %>%
#     select(-c(name, value)) %>%
#     distinct
# }


# Get NAACCR item definition using SEER API
get_naaccr_item_def <- function(item){
  # SEER API key
  seer_api_key <- "39c41a32868d4570ae97c3112ad4459e"
  version = "latest"
  
  request <- paste0("https://api.seer.cancer.gov/rest/naaccr/xml/22/item/", item)
  response = httr::GET(request, httr::add_headers("X-SEERAPI-Key" = seer_api_key))
  httr::stop_for_status(response)
  output = fromJSON_na(httr::content(response, as="text", encoding = 'UTF-8')) %>% 
    map_df(.f = ~.x)
  
  return(output[1,]$documentation)
}

# Create a Histogram
syhistogram <- function (.data, col_name, color = 1, ...)
{
  col_quo = sym(col_name)
  check_for_column(.data,!!col_quo)
  plot_title = col_name
  na_cnt = .data %>% filter(is.na(!!col_quo)) %>% summarise(cnt = n()) %>%
    pull(cnt)
  if (na_cnt > 0)
    plot_title = paste0(plot_title, " (", na_cnt, " NA values)")
  x = plotly::plot_ly(
    x = .data %>% pull(!!col_quo),
    type = "histogram",
    marker = list(color = sypalette$hex[color]),
    ...
  ) %>%
    plotly::layout(title = plot_title)
  return(x)
}


# Compute Missing
compute_missing <- function(df, cols) {
  df %>%
    mutate_at(vars(cols),
              ~ if_else(
                is.na(.) |
                  grepl('unknown|missing', ., ignore.case = T),
                FALSE,
                TRUE
              ))
}


# Compute variable completeness
compute_completeness <- function(df, cols, process = T) {
  if (process) {
    df %>%
      select(cohort, all_of(cols)) %>%
      compute_missing(names(.)[names(.) != 'cohort']) %>%
      make_table(cohort)
  } else {
    df %>%
      select(cohort, all_of(cols)) %>%
      make_table(cohort)
  }
}


# Create Groupings for SVI 
group_svi <- function(df, svi_col = 'svi'){
  svi_col = as.symbol(svi_col)
  df %>% 
    mutate(
      {{ svi_col }} := case_when(
        {{ svi_col }} < 0.299 ~ 'Quartile 1 (<0.299)',
        {{ svi_col }} >= 0.299 & {{ svi_col }} < 0.448 ~ 'Quartile 2 (>=0.299 and <0.448)',
        {{ svi_col }} >= 0.448 & {{ svi_col }} < 0.604 ~ 'Quartile 3 (>=0.448 and <0.604)',
        {{ svi_col }} >= 0.604 ~ 'Quartile 4 (>=0.604)',
        TRUE ~ NA)
    )
}


# Conveneince wrapper for making table with gtsummary



ads_unique_values <- function(ads){
  output = list()
  vars = ads %>% colnames
  
  for(var in vars){
    .var = as.symbol(var)
    
    var_data <- ads %>% 
      distinct(values = {{ .var }})  %>% 
      filter(!is.na(.[[1]])) %>% 
      suppressMessages(readr::type_convert())
    
    var_type = var_data[['values']] %>% class
    
    if(var_type == 'Date'){
      output[[var]] <- tibble(values = 'YYYY-MM-DD')
    #} #else if(grepl('_year$', var)){
      #output[[var]] <- tibble(values = 'YYYY')
    } else if(var_type == 'numeric') {
      output[[var]] <- var_data %>% head(5)
    } else {
      output[[var]] <- var_data %>% arrange(.[[1]])
    }
  }
  
  return(output)
}


kms_api_call <- function(api_root, params = list()){
  output = NULL
  tryCatch({
    response = httr::GET(api_root, query = params)
    httr::stop_for_status(response)
    output = RJSONIO::fromJSON(httr::content(response, as="text", encoding = 'UTF-8'), warn = F)$data %>% map_df(.f = ~.x)
    ## ...
  }, http_error=function(e) {
    output = NULL
  })
  
  return(output)
}


syapi_call <- function(path, params = list(), api_key = Sys.getenv("CONNECT_API_KEY")){
  api_root = 'https://analytics.syapse.com/'
  
  output = NULL
  tryCatch({
    response = httr::GET(api_root,
                         path = paste0('syapi/', path),
                         query = params,
                         httr::add_headers(Authorization = paste0("Key ", api_key)))
    httr::stop_for_status(response)
    output = RJSONIO::fromJSON(httr::content(response, as="text", encoding = 'UTF-8'), warn = F) %>% 
      map_df(.f = ~.x)
    ## ...
  }, http_error=function(e) {
    output = NULL
  })
  
  return(output)
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


.get_ads_data <- function(cohort, variables){
  .cohort = tolower(cohort)
  .variables = variables #c(strsplit(variables, split = ', {0,1}'))[[1]]
  
  ads_type_map = list('enriched'   = syhelpr::list_ads(type='enriched'),
                      'essentials' = syhelpr::list_ads(type='essentials'))
  
  if(.cohort %in% ads_type_map$enriched){
    ads_table_name = glue::glue('ads_{.cohort}_enriched')
  } else {
    ads_table_name = glue::glue('ads_{.cohort}_essentials')
  }
  
  tbl(syhelpr::spmd_con('prod', write = T), dbplyr::in_schema('ca', ads_table_name)) %>%
    select(patientid, sourcename, suborg, any_of(c(.variables))) %>% 
    collect
}

search_ads_col <- function(regex){
  ads %>% head(0) %>% 
    select(matches(regex)) %>% 
    colnames
}


# EMR access ----
get_emr_access <- function(.data, cancer_type_regex, emr_access_vector = emr_access_site){
  
  # Class of case definition for prioritization
  # https://apps.naaccr.org/data-dictionary/data-dictionary/version=24/data-item-view/item-number=610/
  class_of_case_priority <- c(14, 12, 10, 22, 20, 13, 11, 21, 0, 34:37, 40:42, 30:33, 43, 38, 49, 99)
  
  # Ascension site facility map
  # Most disagreement was a result of facilities changing emr facilities and resulting in an updated site name. Some may have been due to move from suborgs
  # We will take the most recent event crf to inform the accurate site_name.
  # This is only necessary for Ascension Michigan, Illinois, and Wisconsin since their are multiple site-names tied to a single suborg
  ascension_site_facility_map <- tbl(spmd_con("clone"), in_schema("ca", "ascension_oc_site_facility_map")) %>% collect()
  
  .data %>%
    distinct(patientid, sourcename, suborg) %>%
    left_join(get_registry(.data,
                           cancer_type = cancer_type_regex,
                           spmd = spmd,
                           c("reportingfacility",
                             "classofcase")) %>%
                left_join(ascension_site_facility_map %>%
                            filter(suborg != "ILLINOIS" |
                                     (suborg == "ILLINOIS" &
                                        reportingfacility %in% c("0006431613", #Alexian Brothers
                                                                 "0000000145", #Alexian Brothers
                                                                 "0006431926", #Saint Alexius
                                                                 "0000000155", #Saint Alexius
                                                                 "0006431070", #Saint Joseph Chicago
                                                                 "0000006065"))) %>% #Saint Joseph Chicago
                            select(reportingfacility, site_name),
                          by = "reportingfacility") %>%
                distinct(patientid, site_name, classofcase),
              by = "patientid") %>%
    mutate(site_name = case_when(sourcename == "advent" & suborg == "Advent Central Florida Division - North" ~ "advent-advent central florida division-north",
                                 sourcename == "advent" & suborg == "Advent Central Florida Division - South" ~ "advent-advent central florida division-south",
                                 sourcename == "advent" & suborg == "Advent Great Lakes" ~ "advent-advent great lakes-cerner",
                                 sourcename == "advent" & suborg == "Advent Kansas" ~ "advent-advent kansas-shawnee mission",
                                 suborg == "aah-il" ~ "aurora-aah-il",
                                 suborg == "aah-wi" ~ "aurora-aah-wi",
                                 sourcename == "bayhealth" ~ "bayhealth-bayhealth",
                                 sourcename == "mlh" ~ "mlh-main line health",
                                 TRUE ~ site_name),
           emr_access = if_else(site_name %in% emr_access_vector, TRUE, FALSE, FALSE),
           classofcase = as.numeric(classofcase)) %>%
    add_sorting(classofcase, by = "custom", sort_vector = class_of_case_priority) %>%
    arrange(
      desc(emr_access),
      # Cole made the below up. Work with the CTRs to determine the best ordering, and perhaps 
      # use a spreadsheet to communicate (with a googlesheets4::read_sheet in data-raw)
      desc(!is.na(site_name)),
      desc(site_name == 'ascension-illinois-aurora'),
      desc(site_name == 'ascension-illinois-epic'),
      desc(site_name == 'ascension-illinois-cerner'),
      desc(site_name == 'advent-advent central florida division-south'),
      desc(site_name == 'advent-advent central florida division-north'),
      desc(site_name == 'ascension-illinois-meditech'), # per Sheryl
      classofcase
    ) %>% 
    slice_head(by = patientid) %>%
    distinct(patientid, sourcename, suborg, site_name, emr_access)
}


# Map patientid to deid patientid ----
map_deid_patient <- function(.data){
  .data %>%
    mutate(patientid_deid = purrr::map_chr(patientid, syhelpr::deid_convert_spmd))
}



    #' Timestamp
#'
#' Return a timestamp in appropriate timezone. Useful for printing timstap statuses in function
#'
#' @param timezone string. Defaults to Pacific 'America/Los_Angeles'
#' @param time_only bool. Controls output format. TRUE gives only time, while FALSE reutnrs date and time. Defaults to TRUE
#'
#' @return a tibble with dates converted to YYYY-MM-DD format
#' @export
#'
timestamp <-
  function(timezone = 'America/Los_Angeles',
           time_only = TRUE) {
    ts <- lubridate::with_tz(lubridate::now(tzone = 'GMT'), tzone = timezone)
    if (time_only) {
      ts <- hms::round_hms(hms::as_hms(ts), 2)
    }
    return(ts)
  }

#' Syapse Shared Path
#'
#' File Path to Syapse CA Shared Directory
#' `/var/lib/rstudio-server/rstudio-users/syapse-shared`
#' @export
#'
syapse_shared_path <-
  "/var/lib/rstudio-server/rstudio-users/syapse-shared"

#' Not All NA
#'
#' Function for removing columns that contain all `NA` values
#'
#' @examples
#' df %>% select_if(not_all_na)
#'
#' @export
#'
not_all_na <- function(x) any(!is.na(x))

#' Check for Dupes
#'
#' Checks for instances where a patient is represented by multiple rows
#'
#' @param df a dataframe containing patient data that is intended to be 1 row per patient
#'
#' @return tibble of records with multiple rows per patient
#'
#' @export
#'
check_for_dupes <- function(df, pat_col = patientid) {
  pat_col = enquo(pat_col)
  df <- df %>%
    group_by(!!pat_col) %>%
    filter(n() > 1) %>%
    ungroup %>%
    arrange(!!pat_col)
  
  if (length(df) == 0) {
    print("==== No Duplicates Found ====")
  } else {
    return(df)
  }
}

#' Syapse Color Palette
#'
#' A Selection of colors consistent with Syapse branding
#'
#' @export
#'
sypalette <- tibble::tribble(
  ~ hex,
  ~ rgb,
  '#0050B3',
  paste0(col2rgb('#0050B3', alpha = F), collapse = ', '),
  '#36cfc9',
  paste0(col2rgb('#36cfc9', alpha = F), collapse = ', '),
  '#bae637',
  paste0(col2rgb('#bae637', alpha = F), collapse = ', '),
  '#acd533',
  paste0(col2rgb('#acd533', alpha = F), collapse = ', '),
  '#ffec3d',
  paste0(col2rgb('#ffec3d', alpha = F), collapse = ', '),
  '#ffa940',
  paste0(col2rgb('#ffa940', alpha = F), collapse = ', '),
  '#ff7a45',
  paste0(col2rgb('#ff7a45', alpha = F), collapse = ', '),
  '#9254de',
  paste0(col2rgb('#9254de', alpha = F), collapse = ', '),
  '#597ef7',
  paste0(col2rgb('#597ef7', alpha = F), collapse = ', '),
  '#40a9ff',
  paste0(col2rgb('#40a9ff', alpha = F), collapse = ', '),
  '#161cff',
  paste0(col2rgb('#161cff', alpha = F), collapse = ', ')
) %>% as.list


#' Syapse Theme Blue
#' @export
theme_blue <- '#0050b3'

#' Syapse Color Palette Ramp
#'
#' A function for imputing color palettes of N length using the sypalette colors
#'
#' @examples 
#' sypalette_ramp(n)
#' @export
#'
sypalette_ramp <- colorRampPalette(sypalette$hex)

#' Day Difference
#'
#' Difference in days between two times
#'
#' @examples
#' day_difference(as.Date("2020-05-30"),as.Date("2020-06-01"))
#' 
#' @return numeric
#'
#' @export
#'
day_difference <- function(time1,time2) {
  return(as.double(difftime(time2,time1,units="days")))
}

#' As percent
#'
#' Format numeric as a percent
#' 
#' @param prop_ a numeric vector
#' @param round integer indicating the number of decimal places to be used.
#'
#' @examples
#' as_percent(0.9)
#' 
#' @return character vector
#'
#' @export
#'
as_percent <- function(prop_,round=0) {
  paste0(round(prop_*100,round),'%')
}

    
#' tbl condition
#'
#' Output a message if a table condition is met
#' 
#' @param .data a data frame
#' @param condition a condition statement to be evaluated on .data (e.g. !is.na(subjectid))
#' @param message a message to be delivered. Note that you can access the number of rows where the 
#' condition is met using '{n}' (uses glue::glue())
#' @param level c("message","warning","error")
#' 
#' @examples
#' tibble(x=c(1,2,3,NA_integer_)) %>% tbl_condition(is.na(x),'{n} row(s) are missing column x')
#'
#' @export
#'
tbl_condition <- function(.data,condition,message,level='message') {
  condition = enquo(condition)
  
  n = nrow(.data %>% filter(!!condition))
  if(n>0) {
    if(level == 'message') {
      inform(paste('Info:',glue(message)))
    } else if(level == 'warning') {
      warn(glue(message))
    } else if(level == 'error') {
      abort(glue(message))
    } else {
      abort('tbl_condition() level must be one of "message", "warning", or "error".')
    }
  }
}

#' Check for column 
#'
#' Confirm a column exists, and if applicable, whether it is the correct type. Returns the unedited data frame.
#' 
#' @param .data a data frame
#' @param .col a variable
#' @param .class character vector of class of the variable. If a vector, make sure the class is one of the values
#' @param return_df a logical. If FALSE, return NULL (for performance)
#' 
#' @return a data frame, if return_df = TRUE, or NULL
#' 
#' @examples 
#' # Both checks pass so the original tibble is returned
#' tibble::tibble(x=1,y="a") %>% check_for_column(x,.class="numeric") %>% 
#'   check_for_column(y,.class=c("numeric","character"))
#'
#' @export
#'
check_for_column <- function(.data,.col,.class=NULL,return_df=TRUE) {
  .col = enquo(.col)
  .vars = NULL
  if(inherits(.data,'tbl_lazy')) {
    .vars = sapply(.data %>% head %>% collect,class)
  } else {
    .vars = sapply(.data,class)
  }
  
  if(as_label(.col) %in% names(.vars)) {
    if(!is.null(.class)) {
      if(!.vars[[as_label(.col)]] %in% .class) {
        abort(glue("{as_label(.col)} must be of class {paste(.class,collapse=' or ')}"))
      }
    }
  } else {
    abort(glue("Can't find column {as_label(.col)} in data frame"))
  }
  
  if(return_df) {
    return(.data)
  }
}

#' base_class
#'
#' Get the base class of an object
#'
#' @param x an object
#' 
#' @examples 
#' class(tibble::tibble())
#' base_class(tibble::tibble())
#' 
#' @return a character vector
#' 
#' @export
base_class <- function(x) {
  tail(class(x),n=1)
}

#' remove_na
#'
#' Remove na values from a vector
#'
#' @param x a vector
#' 
#' @examples 
#' remove_na(c(1,2,NA_integer_))
#' 
#' @return a character vector
#' 
#' @export
remove_na <- function(x) {
  return(x[!is.na(x)])
}

#' Confirm unique
#'
#' Confirm a column is unique
#' 
#' @param .data a data frame
#' @param .col a variable
#' @param allow_na a logical
#' @param return_df a logical. If FALSE, return NULL (for performance)
#' @param error_msg a glue-able character vector
#' 
#' @return a data frame, if return_df = TRUE, or NULL
#' 
#' @examples 
#' # Both checks pass so the original tibble is returned
#' tibble::tibble(x=c(1,2),y=c("a","a")) %>% confirm_unique(x) 
#' tibble::tibble(x=c(1,2),y=c("a","a")) %>% confirm_unique(y) # fails 
#' tibble::tibble(x=c(1,2),y=c("a","a")) %>% confirm_unique(z) # fails
#' tibble::tibble(x=c(1,2),y=c("a",NA_character_)) %>% confirm_unique(y) # passes 
#' tibble::tibble(x=c(1,2),y=c("a",NA_character_)) %>% confirm_unique(y,allow_na=F) # fails 
#' tibble::tibble(x=c(1,2),y=c("a","a")) %>% group_by(x) %>% confirm_unique(y) # fails 
#'
#' @export
#'
confirm_unique <- function(.data,.col,allow_na=TRUE,return_df=TRUE,
                           error_msg='{as_label(.col)} must be unique.') {
  .data = .data %>% ungroup()
  .col = enquo(.col)
  check_for_column(.data,!!.col,return_df=FALSE)
  
  any_dupes = .data %>% 
    count(!!.col) %>% 
    filter(!is.na(!!.col)) %>% 
    summarise(any(n>1,na.rm=T)) %>% 
    pull
  
  if(any_dupes) {
    abort(glue('confirm_unique: {as_label(.col)} must be unique.'))
  }
  
  if(!allow_na) {
    any_na = .data %>% summarise(any(is.na(!!.col),na.rm=T)) %>% pull
    if(any_na) {
      abort(glue('confirm_unique: {as_label(.col)} must not have NA values.'))
    }
  }
  
  if(return_df) {
    return(.data)
  } else {
    return(NULL)
  }
  
}

#' Cross join
#'
#' @param x a data frame
#' @param y a data frame
#' @param copy a logical. See [mutate-joins][dplyr::mutate-joins]
#' @param suffix a named character vector. See [mutate-joins][dplyr::mutate-joins]
#' @param ... See [mutate-joins][dplyr::mutate-joins]
#' @param keep a logical. See [mutate-joins][dplyr::mutate-joins]
#' 
#' @return a data frame
#' 
#' @examples tibble(x=1) %>% cross_join(tibble(y=1))
#'
#' @export
#'
cross_join <- function(x, y, copy = FALSE, suffix = c(".x", ".y"), ..., keep = FALSE) {
  .colbase = "tmp"
  i=0
  
  .col = paste0(.colbase,i)
  while(.col %in% colnames(x) | .col %in% colnames(y)) {
    i = i+1
    .col = paste0(.colbase,i)
  }
  
  full_join(
    x %>% mutate(!!.col := 1), 
    y %>% mutate(!!.col := 1), 
    by=.col, 
    copy = copy, 
    suffix = suffix, 
    #..., 
    keep = FALSE
  ) %>% select(-!!.col)
}
  
#' Expect No Error
#' 
#' Per documentation of expect_error(regexp), "If NA, asserts that there should be no errors."
#' 
expect_no_error <- function(object) {
  do.call(expect_error,list(object=object,regexp=NA))
}

#' Invert named character vector
#' 
#' @export
invert_named_vector <- function(y) setNames(names(y), y)


#' Convert Empty String
#' 
#' Convert empty strings ("") to NA
#'
#' @param .data data frame, tibble, or lazy tibble 
#' @param column character column in .data
#' @param to character string to convert instances of "" to, defaults to NA_character_
#'
#' @return data frame, tibble, or lazy tibble with the value "" in `column` mapped to the value `to`
#' @export
#'
#' @examples
#' example <- as_tibble(list(strings = c("hello", "world", "", NA_character_)))
#' example %>% convert_empty_string(strings)
#' example %>% convert_empty_string(strings, to = "I was empty")
convert_empty_string <- function(.data, column, to = NA_character_) {
  column <- enquo(column)
  .data %>% 
    check_for_column(!!column, .class = "character") %>%
    mutate(!!column := if_else(!!column == "", 
                               to, 
                               !!column, 
                               !!column))
}
  
#' Data source prioritization for combining data from OpenClinica and Registry
#'
#' @param .data Input data-frame containing patient information from multiple sources. The `sourceschema` variable must be present in the 
#' data-frame in order to apply the prioritization logic.
#'
#' @return
#' The input data-frame is returned prioritizing abstracted information from OpenClinica over registry data for each patient record. 
#' 
#' @examples
#' example_multi_source_data <- tribble(~patientid, ~mets_flag, ~mets_date, ~sourceschema,
#' "1", "Yes", "2022-01-16", "ma",
#' "1", "Yes", NA_character_, "registry",
#' "2", "Yes", NA_character_, "registry",
#' "2", "No", NA_character_, "ma",
#' "3", "Yes", "2023-05-05", "ma",
#' "3", "No", NA_character_, "other")
#' 
#' example_multi_source_data %>% prioritize_source()
prioritize_source <- function(.data) {
  .data %>% 
    # MA > registry > all others
    mutate(source_ranking = case_when(sourceschema == "ma" ~ 1,
                                      grepl("registry", sourceschema) ~ 2,
                                      TRUE ~ 3)) %>% 
    group_by(patientid) %>% 
    filter(source_ranking == min(source_ranking, na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(-source_ranking)
}


#' Join data-frames and replace missing values with NAs
#'
#' Full join one data frame to another and replace missing values in column with na_fill 
#' Intended for filling NAs "at the site of definition" of ADS variables as opposed to at the end when all 
#' data frames are joined together. Both `.data` and `ads` must have patientid.
#' 
#' @param .data Input data-frame, typically containing a subset of patientids from `ads`.
#' @param ads Input data-frame to be used for the full_join() that contains all the patientids for the cohort of interest.
#' @param column Column names from `.data` that have NA values that need to be replaced during the full_join().
#' @param na_fill String for replacing the NA values, defaults to "Unknown".
#'
#' @return Data-frame containing all the patientids from `ads`. Additionally, `na_fill` values are used to replace the NAs introduced   
#' in the `column` variables during the full_join() between `.data` and `ads`.
#' 
#' @examples
#' example_ads <- tribble(~patientid, ~diagnosis_date,
#' "id1", "2023-04-21",
#' "id2", "2022-10-03",
#' "id3", "2024-01-16",
#' "id4", "2023-03-17",
#' "id5", "2022-02-26",
#' "id6", "2020-08-23")
#' 
#' example_df <- tribble(~patientid, ~recurrence_flag, ~mets_flag,
#' "id2", NA_character_, "Yes", 
#' "id3", "No", NA_character_,
#' "id6", NA_character_, NA_character_)
#' 
#' example_df %>% 
#' join_replace_na(ads = example_ads, column = recurrence_flag) %>% 
#' join_replace_na(ads = example_ads, column = mets_flag)
join_replace_na <- function(.data, ads, column, na_fill = "Unknown") {
  .data %>% 
    full_join(ads %>% 
                select(patientid) %>% 
                distinct(), 
              by = "patientid") %>% 
    mutate({{column}} := tidyr::replace_na({{column}}, na_fill))
}
  
#' Check patient argument
#'
#' This function checks the input data frame for the typical patient columns, `patientid` and
#' `participantid`. It is intended as data validation tool for the common function design pattern of using
#' an argument to specify a data frame of patients to use in the function (e.g., to filter data to).
#'
#' If either `patientid` and/or `participantid` are present, it returns a vector of the present columns,
#' which can be used in joins or other function logic.
#'
#' If neither `patientid` nor `participantid` are present, it returns an error.
#'
#' If `participantid` is not present and `oc_warning = TRUE`, it warns that `participantid` is not present.
#' It is recommended that `oc_warning` is set to `TRUE` when using this function in the context of
#' OpenClinica data, since OpenClinica data is best filtered by a prioritized `participantid` due to
#' duplicated patients with the same `patientid`.
#'
#' @param .data Data frame to check for the existence of `patientid` and/or `participantid` columns.
#' @param oc_warning Warn when `participantid` is missing? Defaults to FALSE. It is recommended that
#'   `oc_warning` is set to `TRUE` when using this function in the context of OpenClinica data, since
#'   OpenClinica data is best filtered by a prioritized `participantid` due to duplicated patients with the
#'   same `patientid`.
#'
#' @return A vector of the present patient columns ("patientid" and/or "participantid").
#' 
#' @examples
#' x <- tribble(~patientid, ~participantid,
#'              1, 1,
#'              2, 2,
#'              3, 3)
#' x %>% check_patient_arg()
#' x %>% select(-patientid) %>% check_patient_arg()
#' x %>% select(-participantid) %>% check_patient_arg(oc_warning = TRUE)
check_patient_arg <- function(.data, oc_warning = FALSE) {
  if (!is.null(.data)) {
    
    .data_name <- deparse(substitute(.data))
    
    if (!is.data.frame(.data)) abort(glue("{.data_name} is not a data frame"))
    
    patient_cols <- intersect(names(.data), c("patientid", "participantid"))
    
    if (length(patient_cols) == 0 ) abort(glue("patientid and participantid are missing from {.data_name} - at least one of them must be present"))
    
    if (!"participantid" %in% patient_cols & oc_warning) warn(glue("participantid is not in {.data_name} - it is best practice to use a prioritized participantid when pulling from OC due to patient duplication!\nUsing patientid only"))
    
    patient_cols
    
  } else {
    
    NULL
    
  }
}
