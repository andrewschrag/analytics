## Functions for building datasets and supporting quick reports.


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


# Make table1
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
  call_col <- enquo(call_col)
  call.map <-
    tbl(spmd_con('clone'), in_schema('ca', 'map_biomarker_call')) %>% collect
  
  df %>%
    mutate(call = ifelse(
      tolower(biomarkertype) == "wild type",
      "NEGATIVE",
      toupper(call)
    )) %>% 
    left_join(call.map, by = c("call")) %>%
    mutate(call = factor(
      call_simple,
      levels = c(
        "Positive",
        "Equivocal",
        "Negative",
        "Low",
        "Indeterminate",
        "QNS",
        "Unknown"
      )
    )) %>%
    select(-call_simple)
}





## Get Last Contact date
get_last_contact <- function(cohort_query, index = 'dx') {
  .con = cohort_query$src$con
  .query = cohort_query %>% select(patientid)
  .cohort = .query %>% collect
  
  tictoc::tic("====>> run time")
  print(
    glue::glue(
      "{timestamp()} - finding last contact for {nrow(.cohort %>% distinct(patientid))} patients..."
    )
  )
  
  suppressWarnings({
    last_contact.encounter <-
      tbl(.con, in_schema('ca', 'all_encounters')) %>%
      inner_join(.query) %>%
      filter(encounter_date < today()) %>%
      distinct(patientid, encounter_date) %>%
      group_by(patientid) %>%
      summarise(
        last_encounter_date =  max(encounter_date, na.rm = T),
        first_encounter_date = min(encounter_date, na.rm = T)
      ) %>%
      ungroup %>%
      collect %>%
      convert_dates
    
    # last_contact.procedure <-
    #   tbl(.con, in_schema('mdr', 'procedure')) %>%
    #   inner_join(.query) %>%
    #   filter(startdts < today()) %>%
    #   distinct(patientid,
    #            startdts,
    #            startdts_partial,
    #            enddts,
    #            enddts_partial) %>%
    #   collect %>%
    #   impute_dates %>%
    #   group_by(patientid) %>%
    #   summarise(
    #     last_proc_date  = pmax(max(enddts, na.rm = T), max(startdts, na.rm = T), na.rm = T),
    #     first_proc_date = min(startdts, na.rm = T)
    #   ) %>%
    #   ungroup %>%
    #   collect %>%
    #   convert_dates
    #
    #
    # last_contact.diagnosis <-
    #   tbl(.con, in_schema('mdr', 'diagnosis')) %>%
    #   inner_join(.query) %>%
    #   filter(coalesce(diagnosisdate, startdts) < today()) %>%
    #   distinct(patientid, diagnosisdate, startdts) %>%
    #   group_by(patientid) %>%
    #   summarise(
    #     last_dx_date  = pmax(
    #       max(diagnosisdate, na.rm = T),
    #       max(startdts, na.rm = T),
    #       na.rm = T
    #     ),
    #     first_dx_date = pmin(
    #       min(diagnosisdate, na.rm = T),
    #       min(startdts, na.rm = T),
    #       na.rm = T
    #     )
    #   ) %>%
    #   ungroup %>%
    #   collect %>%
    #   convert_dates
    #
    # last_contact.tumor <- tbl(.con, in_schema('mdr', 'tumor')) %>%
    #   inner_join(.query) %>%
    #   distinct(patientid, diagnosisdate, diagnosisdate_partial) %>%
    #   collect %>%
    #   impute_dates %>%
    #   group_by(patientid) %>%
    #   summarise(
    #     last_tumor_date  = max(diagnosisdate, na.rm = T),
    #     first_tumor_date = min(diagnosisdate, na.rm = T)
    #   ) %>%
    #   ungroup %>%
    #   convert_dates
    
    # if (!exists('meds.raw')) {
    #   meds.raw <<- tbl(.con, in_schema('ca', "medications")) %>%
    #     inner_join(.query) %>%
    #     filter(startdts <= today()) %>%
    #     mutate(generic_name = coalesce(drugproduct, ordername)) %>%
    #     distinct(patientid,
    #              generic_name,
    #              startdate = startdts,
    #              enddate = enddts) %>%
    #     collect
    # }
    #
    # last_contact.meds <- meds.raw %>%
    #   filter(startdate <= today()) %>%
    #   group_by(patientid) %>%
    #   summarise(
    #     last_tx_date = pmax(max(enddate, na.rm = T), max(startdate, na.rm = T), na.rm = T),
    #     first_tx_date = min(startdate, na.rm = T)
    #   ) %>%
    #   ungroup  %>%
    #   arrange(patientid, desc(last_tx_date)) %>%
    #   group_by(patientid) %>%
    #   filter(row_number() == 1) %>%
    #   ungroup %>%
    #   convert_dates
    last_contact.oc <-
      tbl(.con, in_schema('openclinica', 'follow_up')) %>%
      distinct(patientid = syapse_patient_id, last_contact_date_oc = ma_dateoflastcontact) %>%
      inner_join(.query %>% mutate(patientid = as.character(patientid))) %>%
      collect
    
    death_dates <-
      tbl(.con, in_schema('mdr', 'patient')) %>%
      rename(patientid = id) %>%
      inner_join(.query) %>%
      distinct(patientid, deceaseddate) %>%
      collect
    
    output <- .cohort %>%
      left_join(last_contact.encounter) %>%
      #left_join(last_contact.tumor) %>%
      # left_join(last_contact.meds) %>%
      # left_join(last_contact.procedure) %>%
      left_join(last_contact.oc) %>%
      left_join(death_dates) %>%
      distinct %>%
      convert_dates %>%
      mutate_if(is.Date, ~ if_else(. > today(), NA_Date_, .)) %>%
      group_by(patientid) %>%
      mutate(
        first_contact_date = first_encounter_date,
        #if_else(index == 'dx', first_tumor_date, first_encounter_date),
        last_contact_date = pmax(
          last_encounter_date,
          # last_tx_date,
          last_contact_date_oc,
          na.rm = TRUE
        ),
        deceaseddate = if_else(deceaseddate < last_contact_date, NA_Date_, deceaseddate)
      ) %>%
      summarise(
        last_contact_date = max(last_contact_date, na.rm = TRUE),
        first_contact_date = min(first_contact_date, na.rm = TRUE),
        deceaseddate = max(deceaseddate, na.rm = TRUE)
      ) %>%
      ungroup
  })
  
  print(glue::glue("{timestamp()} - get_last_contact() complete"))
  tictoc::toc()
  return(output)
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





## Functions ====
# NAACCR Map
get_naaccr_map <- function(column, spmd = spmd_con()) {
  tbl(spmd, in_schema("ca", "map_naaccr_dict")) %>%
    filter(tolower(naaccr_item_id) == column)
}

search_naaccr_cols <- function(pattern) {
  naaccr_colnames <-
    tbl(spmd, in_schema('ca', 'registry_naaccr')) %>% colnames
  
  naaccr_colnames[grepl(pattern, naaccr_colnames)]
}

# Search NAACCR
search_naaccr <-
  function(column,
           search_term = NULL,
           spmd = spmd_con(),
           cols = c()) {
    desc_col = sym(paste0(column, '_description'))
    tbl(spmd, in_schema("ca", "registry_naaccr")) %>%
      left_join(get_naaccr_map(column, spmd) %>%
                  select(code, description),
                by = setNames("code", column)) %>%
      rename({{ desc_col }} := description) %>% 
      {
        if (!is.null(search_term))
          filter(., description == search_term)
        else
          .
      } %>%
      #filter(!is.na(description)) %>%
      select(patientid = clientid,
             mrn = patientidnumber,
             sourcename,
             column,
             {{ desc_col }},
             all_of(cols))
  }


## NAACCR Functions ====

search_naaccr_cols <- function(pattern) {
  naaccr_colnames <-
    tbl(spmd, in_schema('registry', 'aurora_naaccr')) %>% head(0) %>% collect %>% names
  
  naaccr_colnames[grepl(pattern, naaccr_colnames)]
}


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








get_confirmed_surgeries <- function (df,
                                     pat_col = patientid,
                                     schema = "mdr",
                                     cancer_type = "breast") {
  tictoc::tic("====>> run time")
  print(
    glue::glue(
      "{timestamp()} - pulling confirmed (related encounter) surgeries for patients with '{cancer_type}'..."
    )
  )
  surgery_regex_values = dplyr::tribble(
    ~ cancer_type,
    ~ regex,
    "breast",
    "lumpectomy|quadrantectomy|mastectomy|excisional biopsy|axillary lymphadenectomy|wide excision",
  )
  surgery_regex = surgery_regex_values %>% filter(cancer_type ==
                                                    cancer_type) %>% pull(regex)
  pat_col <- enquo(pat_col)
  procedures_confirmed <-
    tbl(spmd, in_schema("mdr", "procedure")) %>%
    filter(proceduretype == "Surgery") %>%
    select(
      patientid,
      proceduretype,
      proceduretype_rawvalue,
      procedure_rawvalue,
      startdts,
      sourceschema
    ) %>%
    collect %>%
    inner_join(df %>%
                 distinct(!!pat_col), by = quo_name(pat_col)) %>% mutate(
                   procedure = tolower(procedure_rawvalue),
                   relevant_surgery = ifelse(grepl(surgery_regex, procedure,
                                                   ignore.case = T),
                                             "Yes",
                                             "No")
                 ) %>% filter(relevant_surgery ==
                                "Yes") %>% left_join(tbl(spmd, in_schema("ca", "map_definitive_surgery")) %>%
                                                       as_tibble %>% filter(cancer_type == cancer_type),
                                                     by = "procedure") %>%
    mutate(
      startdts = as.Date(startdts),
      procedure = case_when(
        procedure ==
          "mastectomy, simple" ~ "mastectomy",
        TRUE ~ procedure
      )
    ) %>%
    filter(!is.na(rollup)) %>% left_join(
      tbl(spmd, in_schema(schema,
                          "encounter")) %>% filter(encountertype == "Inpatient") %>%
        select(patientid, startdts) %>% collect %>% mutate(
          startdts = as.Date(startdts),
          related_encounter = "Yes"
        ),
      by = c("patientid", "startdts")
    ) %>%
    filter(related_encounter == "Yes" | (grepl(
      "ascension",
      sourceschema, ignore.case = T
    ))) %>% select(
      -c(
        procedure,
        cancer_type,
        relevant_surgery,
        sourceschema,
        related_encounter,
        proceduretype
      )
    ) %>% distinct
  print(glue::glue("{timestamp()} - get_confirmed_surgeries() complete"))
  tictoc::toc()
  return(procedures_confirmed)
}


build_med_regex <- function (med_regex, spmd = spmd_con())
{
  if (!is.character(med_regex)) {
    stop("med_regex is not a character vector")
  }
  if (length(med_regex) > 1) {
    stop("med_regex is not a character vector of length one (not a string)")
  }
  if (!spmd_list_tables("ca") %>% filter(table_name == "map_atc") %>%
      count() %>% pull(n) > 0) {
    stop(
      "ca.map_atc does not exist! Please rebuild and try again:\n",
      "library(syhelpr)\n",
      "map_atc <- build_atc_table()\n",
      "spmd_write_table(map_atc,'map_atc',overwrite = TRUE)"
    )
  }
  tbl(spmd, in_schema("ca", "map_atc")) %>%
    filter(grepl(med_regex, all_names_regex, ignore.case = TRUE)) %>%
    collect %>%
    mutate(all_names_regex =  gsub('/\\W+/g', '.*', all_names_regex, perl = T)) %>%
    pull(all_names_regex) %>%
    str_split("\\|") %>% c(str_split(med_regex, "\\|"), .) %>%
    unlist() %>% unique() %>% paste0(collapse = "|", sep = "")
}



get_medrio_biomarkers_breast <- function(gene) {
  medrio_biomarker.query <-
    tbl(spmd, in_schema('medrio', 'breast_biomarkers')) %>%
    inner_join(tbl(.$src$con, in_schema('ca', 'map_spmd_subjectid'))) %>%
    filter(!is.na(vargroup1row))
  medrio_biomarker.query %>% head(0) %>% collect %>% names
  brca.ma <- medrio_biomarker.query %>%
    collect %>%
    mutate(
      biomarker_name = coalesce(breast_germline_biomarker, breast_somatic_biomarker),
      genomicsource = case_when(
        !is.na(breast_germline_biomarker) ~ 'Germline',
        !is.na(breast_somatic_biomarker) ~ 'Somatic'
      )
    ) %>%
    filter(grepl(gene, biomarker_name),
           mutation_state == 'Mutated') %>%
    mutate(
      genomicsource = biomarker_class,
      specific_variant = ifelse(
        grepl('not stated|unknown|unk', specific_variant, ignore.case = TRUE),
        'Variant Not Stated',
        specific_variant #gsub(' .*$|\\:|\\(|\\)|\\,|\\;', '', specific_variant)
      ),
      codingchange  = trimws(gsub(
        ' .*$|\\:|\\(|\\)|\\,|\\;',
        ' ',
        gsub(
          'c(\\.){0,1}',
          'c.',
          str_extract(
            gsub('\\(|\\)', ' ', biomarker_raw),
            'c(\\.){0,1}[0-9]{2,4}.*'
          )
        )
      )),
      aminoacidchange = trimws(gsub(
        ' .*$|\\:|\\(|\\)|\\,|\\;',
        ' ',
        str_extract(
          gsub('\\(|\\)', ' ', biomarker_raw),
          'p(\\.){0,1}[A-Z]{1}[a-z]{2,3}.*'
        )
      )),
      variant = trimws(gsub(
        '\\(|\\)', '', paste(
          replace_na(codingchange, ''),
          ifelse(
            grepl('p\\.', specific_variant),
            specific_variant,
            replace_na(aminoacidchange, '')
          )
        )
      ), which = 'both'),
      variant = ifelse(variant == '', specific_variant, variant),
      aminoacidchange = coalesce(aminoacidchange, str_extract(variant, 'p\\..*')),
      gene = gsub(' Mutation', '', biomarker_name)
    ) %>%
    impute_dates %>%
    rename(reportdate = report_date,
           call = mutation_state,
           biomarkertype = biomarker_class) %>%
    process_call %>%
    select(
      patientid,
      reportdate,
      call_clean,
      logical,
      genomicsource,
      gene,
      variant,
      codingchange,
      aminoacidchange
    ) %>%
    distinct %>%
    group_by(patientid) %>%
    mutate_at(.vars = c('variant'), ~ paste(unique(unlist(str_split(
      ., " "
    ))), collapse = " ")) %>%
    mutate(aminoacidchange = trimws(str_remove(
      aminoacidchange, replace_na(codingchange, 'NONE')
    )),
    codingchange = trimws(str_remove(
      codingchange, replace_na(aminoacidchange, 'NONE')
    )))
  return(brca.ma)
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
           realtime = F,
           schema = 'cohorts',
           con = spmd_con()) {
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
  spmd_write = syhelpr::spmd_con("clone", write = TRUE)
  tictoc::tic('====>> run time')
  rlang::inform(as.character(syhelpr::timestamp()))
  
  if (!recreate_mat_view) {
    syhelpr::run_query("REFRESH MATERIALIZED VIEW ca.cohort_counts;", spmd_write)
  } else {
    message('Dropping Existing Materialized View - ca.all_encounters')
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


# source(
#   file.path(
#     '/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils/cohort_func.R'
#   )
# )

source(
  file.path(
    '/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils/swimmer_general.R'
  )
)
source(
  file.path(
    '/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils/swimmer_urology.R'
  )
)




impute_date_oc <-   function (.data,
                              date_col_nm,
                              day = "15",
                              month_day = "06-30",
                              impute_year = T,
                              keep_partial_date_field = F)
{
  day = as.character(day)
  gran_col_nm = str_c(date_col_nm, "_granularity")
  y_col_nm <- str_c(date_col_nm, "_y")
  ym_col_nm <- str_c(date_col_nm, "_ym")
  unk_col_nm <- str_c(date_col_nm, "_unk")
  contains_unk = unk_col_nm %in% (.data %>% colnames)
  .data = .data %>% check_for_column(!!rlang::sym(date_col_nm),
                                     .class = c("character", "Date")) %>% check_for_column(!!rlang::sym(y_col_nm),
                                                                                           .class = c("character", "Date")) %>% check_for_column(!!rlang::sym(ym_col_nm),
                                                                                                                                                 .class = c("character", "Date"))
  .data_transformed = .data %>% mutate(
    `:=`(!!date_col_nm,
         as.character(!!rlang::sym(date_col_nm))),
    `:=`(
      !!gran_col_nm,
      case_when(
        !is.na(!!rlang::sym(date_col_nm)) ~ "DAY",
        !is.na(!!rlang::sym(ym_col_nm)) ~ "MONTH",!is.na(!!rlang::sym(y_col_nm)) ~
          "YEAR",
        TRUE ~ "NONE"
      )
    ),
    `:=`(!!date_col_nm,
         lubridate::as_date(
           case_when(
             !!rlang::sym(gran_col_nm) ==
               "DAY" ~ !!rlang::sym(date_col_nm),!!rlang::sym(gran_col_nm) ==
               "MONTH" ~ str_c(str_sub(!!rlang::sym(ym_col_nm),
                                       1, 7), "-", day),!!rlang::sym(gran_col_nm) == "YEAR" &
               impute_year ~ str_c(str_sub(!!rlang::sym(y_col_nm),
                                           1, 4), "-", month_day),
             TRUE ~ NA_character_
           )
         )),
    .after = !!date_col_nm
  )
  if (!keep_partial_date_field) {
    .data_transformed = .data_transformed %>% select(-all_of(c(y_col_nm,
                                                               ym_col_nm)))
    if (contains_unk) {
      .data_transformed = .data_transformed %>% select(-all_of(unk_col_nm))
    }
  }
  return(.data_transformed)
}


impute_date_oc <-   function (.data,
                              date_col_nm,
                              day = "15",
                              month_day = "06-30",
                              impute_year = T,
                              keep_partial_date_field = F)
{
  day = as.character(day)
  gran_col_nm = str_c(date_col_nm, "_granularity")
  y_col_nm <- str_c(date_col_nm, "_y")
  ym_col_nm <- str_c(date_col_nm, "_ym")
  unk_col_nm <- str_c(date_col_nm, "_unk")
  contains_unk = unk_col_nm %in% (.data %>% colnames)
  .data = .data %>% check_for_column(!!rlang::sym(date_col_nm),
                                     .class = c("character", "Date")) %>% check_for_column(!!rlang::sym(y_col_nm),
                                                                                           .class = c("character", "Date")) %>% check_for_column(!!rlang::sym(ym_col_nm),
                                                                                                                                                 .class = c("character", "Date"))
  .data_transformed = .data %>% mutate(
    `:=`(!!date_col_nm,
         as.character(!!rlang::sym(date_col_nm))),
    `:=`(
      !!gran_col_nm,
      case_when(
        !is.na(!!rlang::sym(date_col_nm)) ~ "DAY",
        !is.na(!!rlang::sym(ym_col_nm)) ~ "MONTH",!is.na(!!rlang::sym(y_col_nm)) ~
          "YEAR",
        TRUE ~ "NONE"
      )
    ),
    `:=`(!!date_col_nm,
         lubridate::as_date(
           case_when(
             !!rlang::sym(gran_col_nm) ==
               "DAY" ~ !!rlang::sym(date_col_nm),!!rlang::sym(gran_col_nm) ==
               "MONTH" ~ str_c(str_sub(!!rlang::sym(ym_col_nm),
                                       1, 7), "-", day),!!rlang::sym(gran_col_nm) == "YEAR" &
               impute_year ~ str_c(str_sub(!!rlang::sym(y_col_nm),
                                           1, 4), "-", month_day),
             TRUE ~ NA_character_
           )
         )),
    .after = !!date_col_nm
  )
  if (!keep_partial_date_field) {
    .data_transformed = .data_transformed %>% select(-all_of(c(y_col_nm,
                                                               ym_col_nm)))
    if (contains_unk) {
      .data_transformed = .data_transformed %>% select(-all_of(unk_col_nm))
    }
  }
  return(.data_transformed)
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
