library(reticulate)
library(botor)


set_aws_token <- function() {
  message(glue::glue("Setting up AWS Session Token"))
  
  Sys.setenv("AWS_DEFAULT_REGION" = "us-west-2")
  session <- boto3$Session(region_name = Sys.getenv("AWS_DEFAULT_REGION"))
  credentials <- session$get_credentials()
  current_credentials <- credentials$get_frozen_credentials()
  
  Sys.setenv("AWS_ACCESS_KEY_ID" = current_credentials['access_key'],
             "AWS_SECRET_ACCESS_KEY" = current_credentials['secret_key'],
             "AWS_SESSION_TOKEN" = current_credentials['token'])
}


setup_iam_auth <- function() {
  
  # path here depends on where miniconda is installed in dockerfile.
  # reticulate::py_install by default looks for miniconda installation and prompts Y/n for installation.
  # In order to avoid the interactive mode we set the python path before installation
  Sys.setenv("RETICULATE_PYTHON"="/root/miniconda3/bin/python")
  reticulate::py_install('boto3')
  set_aws_token()
}



mdx_redux_con <- function() {
  port = 5432
  user = "view"
  host = 'mdx-redux-prod.czsq80p56jgd.us-west-2.rds.amazonaws.com'
  svc <- paws::rds()
  token <-
    svc$build_auth_token(
      endpoint = sprintf("%s:%s", host, port),
      region = "us-west-2",
      user = user
    )
  conn <- DBI::dbConnect(
    odbc::odbc(),
    Driver = "postgresql",
    database = "mdx_service",
    servername = host,
    Port = port,
    UID = user,
    PWD = token,
    sslmode = "require",
    MaxVarChar = 128000
  )
  return(conn)
}

get_mdx_results <- function(mdx = mdx_redux_con()) {
  ## Call Map ====
  call.map <-
    tbl(spmd_con(), in_schema('ca', 'map_biomarker_call')) %>%
    collect
  
  health_systems = tbl(spmd_con(), in_schema('mdr', 'patient')) %>% 
    distinct(sourcename) %>% 
    collect %>% 
    pull(sourcename)
  
  tumor_map = tbl(spmd_con(), dbplyr::in_schema("maps", "tumorhierarchy")) %>%
    select(1:7, tumortype = ort_drt_display_label) %>%
    mutate(
      code_stripped = stringr::str_replace(code, "\\.", ""),
      benign = grepl("^D", code)
    ) %>%
    filter(!benign) %>% 
    filter(
      !tumortype %in% c(
        'Affected Secondary organ',
        'Non-Neoplastic Hematologic',
        'Illdefined Primary Cancer',
        'Unknown if Malignant or Benign',
        'Benign or Carcinoid Tumors'
      )) %>% 
    collect
  
  ## Test Results ====
  test_results <- tbl(mdx_redux_con(), 'test_results') %>%
    select(-c(test_name, created_at, updated_at)) %>%
    filter(
      biomarker_type %in% c(
        'GENOMIC_ALTERATION',
        'PROTEIN_EXPRESSION_ALTERATION',
        'MICROSATELLITE_INSTABILITY',
        'TUMOR_MUTATION_BURDEN',
        'LOSS_OF_HETEROZYGOSITY_ALTERATION',
        'TRANSLOCATION',
        'HRD_ALTERATION'
      )
    ) %>% 
    inner_join(
      tbl(.$src$con, 'patient_information') %>%
        select(
          molecular_data_id,
          patientid = syapse_patient_id,
          mrn,
          last_name,
          first_name,
          dob,
          sex,
          diagnosis,
          health_system = source_name,
          sub_org
        ),
      by = 'molecular_data_id'
    ) %>%
    inner_join(
      tbl(.$src$con, 'test_details') %>%
        distinct(
          molecular_data_id,
          lab_name,
          test_name,
          ordered_date,
          received_date,
          report_date
        )
    )  %>%
    left_join(tbl(.$src$con,'clinical_trials') %>% 
                distinct(molecular_data_id, nct_id) %>% 
                group_by(molecular_data_id) %>% 
                summarise(nct_id = str_flatten(nct_id, collapse = '; '))) %>%  
    left_join(tbl(.$src$con, 'healthcare_organization') %>% 
                select(molecular_data_id, health_system_test = name)) %>% 
    mutate_at(vars(dob, report_date, ordered_date, received_date),
              ~ as.Date(.)) %>%
    mutate(
      genes = str_remove_all(str_replace_all(as.character(genes), ',', '::'), '\\{|\\}'),
      biomarker_name = str_replace_all(biomarker_name, 'PD-L1', 'PDL1'),
      result_raw = coalesce(
        as.character(biomarker_state),
        as.character(result),
        as.character(hrd_status),
        as.character(mmr_call),
        as.character(msi_call),
        as.character(comment),
        as.character(mutation_burden_call),
        ifelse(as.logical(is_expressed), 'Expressed', 'No Expression')
      ),
      result_raw = tolower(ifelse(
        grepl('Not Detected', biomarker_name),
        'Negative',
        result_raw
      )),
      report_diagnosis = str_remove_all(as.character(diagnosis), '\\{|\\}|\\\\|"'),
      report_year = lubridate::year(report_date)
    ) %>%
    # filter(genes != '', !is.na(result_raw)) %>%
    select(
      id,
      molecular_data_id,
      patientid,
      mrn,
      last_name,
      first_name,
      dob,
      sex,
      health_system,
      health_system_test,
      sub_org,
      lab_name,
      test_name,
      ordered_date,
      received_date,
      report_year,
      report_date,
      report_diagnosis,
      platform_technology,
      genomic_source,
      genes,
      biomarker_type,
      result_raw,
      biomarker_name,
      alteration_details,
      matches('hgvs'),
      hrd_score,
      nct_id
      #matches('call|^is|score|hgvs|status')
    ) %>%
    collect %>% 
    mutate(
      hs_ext = str_extract(tolower(health_system_test), paste0(health_systems, collapse='|')),
      health_system = coalesce(health_system, hs_ext),
      sub_org = coalesce(sub_org, health_system_test)) %>% 
    select(-c(hs_ext, health_system_test)) %>% 
    inner_join(call.map, by = c("result_raw" = "call")) %>%
    rename(result = call_simple) %>%
    distinct %>% 
    group_by(molecular_data_id) %>% 
    mutate(icd10 = sort(unlist(strsplit(as.character(report_diagnosis), ",")))[1],
           icd10_stripped = stringr::str_replace(icd10, "\\.", "")) %>% 
    ungroup %>% 
    left_join(
      tumor_map %>%
        filter(code_system == "ICD-10-CM") %>% 
        select(code_stripped, tumortype),
      by = c(icd10_stripped = "code_stripped"))
  
  
  return(test_results)
}



build_ct_cohort <- function(followup = TRUE,
                            age_breaks = c(0, 50, 64, 74, Inf),
                            age_labels = c('<50', '50-64', '65-74', '75+'),
                            con = syhelpr::spmd_con('prod'),
                            write_table = F,
                            ...) {
  tictoc::tic('====>> build_cohort() run time')
  message(glue::glue('{syhelpr::timestamp()} - building cohort'))
  
  ## Variables and Objects ====
  output = list()
  schema = 'mdr'
  
  index_cols =
    c('patientid',
      'diagnosisdate',
      'metastasisdate',
      'advanceddate')
  
  tumor_map = tbl(con, dbplyr::in_schema("maps", "tumorhierarchy")) %>%
    select(1:7, tumortype = ort_drt_display_label) %>%
    mutate(
      code_stripped = stringr::str_replace(code, "\\.", ""),
      benign = grepl("^D", code)
    ) %>%
    filter(
      !tumortype %in% c(
        'Affected Secondary organ',
        'Non-Neoplastic Hematologic',
        'Illdefined Primary Cancer',
        'Unknown if Malignant or Benign',
        'Benign or Carcinoid Tumors'
      )
    )
  
  ## Get MDX Results
  test_results = get_mdx_results() 
  
  #### Cohort ====
  cohort_query = tbl(con, in_schema('mdr', 'patient')) %>%
    filter(id %in% !!unique(test_results$patientid)) %>%
    select(patientid = id,
           mrn,
           lastname,
           firstname,
           sex,
           birthdate,
           isdeceased,
           sourcename,
           suborg)
  
  #### Cancer Data ====
  ## MDR Diagnosis - ICD10
  cancer.mdr = cohort_query %>%
    inner_join(
      tbl(.$src$con, in_schema('mdr', 'tumor')) %>%
        select(
          patientid,
          tumorid = id,
          diagnosisdate,
          tumortype,
          primarysite,
          histology
        )
    ) %>%
    filter(
      !tumortype %in% c(
        'Affected Secondary organ',
        'Non-Neoplastic Hematologic',
        'Illdefined Primary Cancer',
        'Unknown if Malignant or Benign',
        'Benign or Carcinoid Tumors'
      )
    ) %>% 
    collect %>%
    group_by(patientid) %>%
    arrange(desc(diagnosisdate)) %>%
    slice(1) %>%
    ungroup %>%
    mutate(cancer_source = 'mdr.tumor')
  
  cancer.icd10 =
    tbl(con, dbplyr::in_schema('mdr', "diagnosis")) %>%
    inner_join(cohort_query %>% select(patientid)) %>%
    filter(
      grepl('^C[0-9]{2}', diagnosis_rawvalue),
      !patientid %in% !!unique(cancer.mdr$patientid)
    ) %>%
    select(-startdts) %>%
    mutate(icd10_stripped = stringr::str_replace(diagnosis_rawvalue, "\\.", "")) %>%
    inner_join(
      tumor_map %>%
        filter(code_system == "ICD-10-CM"),
      by = c(icd10_stripped = "code_stripped")
    ) %>%
    left_join(tbl(.$src$con, in_schema('mdr', 'encounter')) %>%
                select(encounterid = id, startdts, admitdts)) %>%
    group_by(patientid, tumortype) %>%
    mutate(diagnosisdate = as.Date(pmin(
      noteddate, diagnosisdate, startdts, admitdts, na.rm = T
    ))) %>%
    ungroup %>%
    select(patientid,
           sourcename,
           suborg,
           tumortype,
           icd10_stripped,
           diagnosisdate) %>%
    group_by(patientid, tumortype) %>%
    dbplyr::window_order(patientid, tumortype, diagnosisdate, icd10_stripped) %>%
    filter(row_number() == 1) %>%
    ungroup %>%
    left_join(
      tbl(con, dbplyr::in_schema("ca", "map_icd10_icdo3")) %>%
        mutate(icd10_stripped = stringr::str_replace(icd10, "\\.", "")),
      by = c("icd10_stripped")
    ) %>%
    collect %>%
    group_by(patientid) %>%
    arrange(desc(diagnosisdate)) %>%
    slice(1) %>%
    ungroup %>%
    select(
      patientid,
      diagnosisdate,
      tumortype,
      primarysite = icdo3t_term,
      histology = icdo3m_term
    ) %>%
    mutate(cancer_source = 'mdr.diagnosis')
  
  output$data$cancer = bind_rows(cancer.mdr, cancer.icd10)
  
  ## Staging
  output$data$stage =
    get_stage(cohort_query, cols = c('t', 'n', 'm'))
  
  ## Metastasis
  output$data$metastasis =
    get_metastasis_date(cohort_query, 'first')
  
  output$data$treated_med = dplyr::tbl(con, dbplyr::in_schema("ca", "medications")) %>% 
    filter(drugcategory == 'Antineoplastic') %>% 
    inner_join(cohort_query %>% select(patientid)) %>% 
    collect %>% 
    inner_join(output$data$cancer %>% select(patientid, diagnosisdate)) %>% 
    filter(startdts>=diagnosisdate) %>% 
    mutate(med = coalesce(drugproduct, ordername)) %>% 
    group_by(patientid) %>% 
    summarise(antineoplastics = paste(sort(unique(med)), collapse = '; ')) %>% 
    ungroup
  
  output$data$surgery = dplyr::tbl(con, in_schema('mdr', 'procedure')) %>% 
    filter(proceduretype == 'Surgery') %>% 
    inner_join(cohort_query %>% select(patientid)) %>% 
    collect %>% 
    inner_join(output$data$cancer %>% select(patientid, diagnosisdate)) %>% 
    filter(startdts>=diagnosisdate) %>% 
    group_by(patientid) %>% 
    summarise(surgeries = paste(sort(unique(procedure)), collapse = '; ')) %>% 
    ungroup
  
  ## Followup
  if (followup) {
    output$data$last_contact = get_last_contact(cohort_query)
  }
  
  #### Final Cohort ====
  output$cohort = output$data %>%
    reduce(left_join) %>%
    mutate(
      diagnosisyear = lubridate::year(diagnosisdate),
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
          patientid %in% c(output$data$metastasis$patientid)
      ),
      metastasisdate = if_else(grepl('IV', stage_group),
                               diagnosisdate,
                               metastasisdate),
      index_mets_date = metastasisdate,
      advanced = metastatic |
        stage_category %in% c('Metastatic', 'Locally Advanced'),
      advanceddate = case_when(
        grepl('III[B-C]|IV', stage_group) ~ diagnosisdate,
        metastatic ~ metastasisdate,
        TRUE ~ as.Date(NA)
      ),
      treated_med = ifelse(patientid %in% output$data$treated_med$patientid, 'Yes', 'No'),
      treated_surgery = ifelse(patientid %in% output$data$surgery$patientid, 'Yes', 'No')
    ) %>%
    rename(first_mets_site = mets_site,
           health_system = sourcename) %>%
    rename_at(.vars = c(t, n, m), ~paste0('stage_', .)) %>% 
    mutate(isdeceased = if_else(!is.na(deceaseddate), TRUE, isdeceased)) %>% 
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
  
  test_results.unknown <<- test_results %>%
    filter(!patientid %in% output$cohort$patientid) %>%
    rename(birth_date = dob) %>%
    mutate(cancer_source = 'unmatched')
  
  
  output$data$test_results = test_results
  output$cohort = bind_rows(
    output$cohort %>%
      inner_join(
        output$data$test_results %>%
          select(-any_of(
            c('mrn',
              'last_name',
              'first_name',
              'dob',
              'sex',
              'health_system',
              'tumortype')
          )) %>%
          distinct,
        by = c('patientid')
      ) %>%
      rename(
        last_name  = lastname,
        first_name = firstname,
        birth_date = birthdate
      ),
    test_results.unknown
  ) %>%
    rename(
      codingchange_raw = hgvs_coding_change,
      aminoacidchange_raw = hgvs_protein_change
    )
  
  alteration_details = output$cohort %>%
    filter(call_flag) %>% 
    select(id, molecular_data_id, alteration_details) %>%
    filter(!is.na(alteration_details)) %>%
    as.tbl_json('alteration_details') %>%
    enter_object(transcriptAlterationDetails) %>%
    spread_all %>%
    as_tibble %>%
    distinct(id, molecular_data_id, transcriptid_raw = transcriptID)

  # biomarker_norm = output$cohort %>%
  #   select(patientid, id, molecular_data_id, codingchange_raw, aminoacidchange_raw, genes) %>% 
  #   left_join(alteration_details) %>% 
  #   kms_normalize_variant('codingchange_raw', 'aminoacidchange_raw', 'transcriptid_raw', 'genes') %>%
  #   filter(!is.na(coalesce(aminoacidchange, codingchange))) %>%
  #   distinct(patientid,
  #            genes,
  #            aminoacidchange,
  #            codingchange,
  #            oncokb_key)
  
  output$cohort = output$cohort  %>%
    left_join(alteration_details) #%>%
    # mutate_at(c('treated_med', 'treated_surgery'), ~replace_na(., 'Unknown'))
    # left_join(biomarker_norm, by = c('patientid', 'genes')) %>%
  #select(-c(aminoacidchange_raw, codingchange_raw)) %>%
  #mutate(biomarker_name = coalesce(oncokb_key, biomarker_name))
  
  message(glue::glue('{syhelpr::timestamp()} - build_ct_cohort() complete'))
  tictoc::toc()
  return(output)
}


kms_api_call <- function(api_root, body) {
  output = NULL
  tryCatch({
    response = httr::POST(url = api_root, body = body)
    httr::stop_for_status(response)
    output = RJSONIO::fromJSON(httr::content(response, as = "text", encoding = 'UTF-8'),
                               warn = F)$data %>% map_df(.f = ~ .x)
    ## ...
  }, http_error = function(e) {
    output = NULL
  })
  
  return(output)
}


kms_normalize_variant <- function(df,
                                  codingchange_col = NULL,
                                  aminoacidchange_col = NULL,
                                  transcriptid_col = NULL,
                                  gene_col = 'genes',
                                  cols = c()) {
  api_root = 'https://ci-knowledge.syapse.com/api/V1/normalize_variants_batch'
  .gene_col = sym(gene_col)
  .col_aminoacid = sym(aminoacidchange_col)
  .col_codingchange = sym(codingchange_col)
  .col_transcriptid = sym(transcriptid_col)
  
  # variant_list = .data
  variant_data = output$cohort %>%
    left_join(alteration_details) %>%
    filter(biomarker_type == 'GENOMIC_ALTERATION', call_flag) %>%
    distinct(gene_raw = genes,
             aminoacidchange_raw,
             codingchange_raw,
             transcriptid_raw) %>%
    filter(!is.na(coalesce(aminoacidchange_raw, codingchange_raw))) %>%
    toJSON()
  
  result = kms_api_call(api_root, variant_data)
  
  output =  output$cohort %>%
    left_join(alteration_details) %>% distinct(patientid,
                                               gene_raw = genes,
                                               aminoacidchange_raw,
                                               codingchange_raw,
                                               transcriptid_raw) %>%
    # distinct(
    #   patientid,
    #   {{ .gene_col }} := toupper(gsub(' ', '', {{ .gene_col }})),
    #   {{ .col_codingchange }} := gsub(' ', '', {{ .col_codingchange }}),
    #   {{ .col_aminoacid }} := gsub(' ', '', {{ .col_aminoacid }}),
    #   {{ .col_transcriptid }},
    #   across(any_of(cols))) %>%
    left_join(result,
              by = c(
                setNames('codingchange_raw', codingchange_col),
                setNames('aminoacidchange_raw', aminoacidchange_col),
                setNames('transcriptid_raw', transcriptid_col),
                setNames('gene_raw', gene_col)
              ))
  
  return(output)
}

# spmd_con <-
#   function (write = FALSE,
#             user = Sys.getenv("SPMD_USER"),
#             max_char = 128000)
#   {
#     source = 'prod'
#  
#     if (write) {
#       if (user == "") {
#         user = "ca"
#       }
#     }
#     else {
#       if (user == "") {
#         user = "view"
#       }
#     }
#     port = "5432"
#     region = "us-west-2"
#     if (write) {
#       host = "spmd-prod.cluster-czsq80p56jgd.us-west-2.rds.amazonaws.com"
#     }
#     else {
#       host = "spmd-prod-replica.czsq80p56jgd.us-west-2.rds.amazonaws.com"
#     }
#     
#     svc <- paws::rds()
#     token <- svc$build_auth_token(
#       endpoint = sprintf("%s:%s",
#                          host, port),
#       region = region,
#       user = user
#     )
#     spmd <- DBI::dbConnect(
#       odbc::odbc(),
#       Driver = "postgresql",
#       database = "spmd",
#       servername = host,
#       Port = port,
#       UID = user,
#       PWD = token,
#       sslmode = "require",
#       MaxVarChar = max_char
#     )
#     return(spmd)
#   }



# S3 buckets ----
s3_bucket_path <- list()
s3_bucket_path$id <- glue::glue("s3://syapse-ephemeral/internal_ads/qc/")
write_analytical_dataset <- function(data, dataset_name, type) {
  if (type == "deid") {
    sink <- glue::glue("{s3_bucket_path$deid}{dataset_name}.parquet")
  } else if (type == "id") {
    sink <- glue::glue("{s3_bucket_path$id}{dataset_name}.parquet")
  } else {
    rlang::abort(glue::glue("type = '{type}' not valid option. Use 'deid' or 'id' instead."))
  }
  message(glue::glue("Writing {type} ADS {dataset_name} to S3"))
  # added this step to refresh aws token
  set_aws_token()
  data %>%
    arrow::write_parquet(sink = sink)
  
  message(glue::glue("Finished writing {type} {dataset_name} to S3"))
}

not_all_na <- function(x)
  any(!is.na(x))


make_table <-
  function (df,
            ...,
            sort = c(all_categorical() ~ "frequency"),
            missing = 'ifany',
            add_overall = F) {
    gttable = df %>%
      tbl_summary(
        missing = missing,
        sort = sort,
        ...,
        type = all_continuous() ~ "continuous2",
        statistic = all_continuous() ~ c("{N_nonmiss}",
                                         "{median} ({p25}, {p75})",
                                         "{min}, {max}")
      ) %>%
      bold_labels()
    if (add_overall)
      gttable = gttable %>% add_overall()
    gttable %>%
      as_gt()
  }
