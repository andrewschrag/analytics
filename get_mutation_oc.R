#' Get OpenClinica Biomarker
#' Allows for querying OpenClinica formdata directly to return
get_mutation_oc <- function(gene,
                            normalize = T,
                            limit = 'ALL',
                            patterns = list(p_pattern =
                                              '^((p.){0,1}\\s{0,1}([a-z]{1,3}[0-9]{3}_[a-z]{1,3}[0-9]{3}|[a-z][0-9]{3}[a-z]|[a-z]{3}[0-9]{3}[a-z]{3})(dup|delins|del|ins){0,1}(ala|arg|asn|asp|cys|glu|gln|gly|his|ile|leu|lys|met|phe|pro|ser|thr|trp|tyr|val){0,4}[ACTGactg]{0,12})',
                                            c_pattern =
                                              '^((c.){0,1}\\s{0,1}([0-9]{4}[a-z]{1}.[a-z]{1}|[0-9]{4}(_|\\-)[0-9]{4})(dup|delins|del|ins){0,1}[actg]{0,12})',
                                            e_pattern =
                                              '([Ee]xon)\\s{0,1}(1[0-9]|2[0-8]|[1-9])'))
{
  ## Setup
  .gene = tolower(gene)
  output = list()
  .variant_cols = c('variantnamec', 'variantnamep')#, 'variantnameexon')

  ## Regex patterns
  p_pattern = patterns$p_pattern
  c_pattern = patterns$c_pattern
  e_pattern = patterns$e_pattern

  ## Get Data
  con = spmd_con('prod', max_char = 96000)
  json = run_query(
    glue::glue(
      "SELECT
     patientid,
     formdata::json -> 'biomarker_report' as biomarker_report
   FROM openclinica.formdata,
        jsonb_array_elements(formdata->'biomarker_report') biomarker_report,
        jsonb_array_elements(biomarker_report->'biomarker') biomarker
   WHERE biomarker->>'ma_genenamesomatic' = '{.gene}'
      limit {limit};"
    ),
    connection = con
  )

  ## Wrangle Data
  biomarker_report_raw = json %>%
    as.tbl_json('biomarker_report') %>%
    select(patientid) %>%
    gather_array()

  biomarker_report = biomarker_report_raw %>%
    enter_object(report) %>%
    spread_all() %>%
    as_tibble %>%
    distinct %>%
    select(-any_of(c('reportdate', 'specimendate'))) %>%
    rename_all(.funs = ~ gsub('report\\.|ma_', '', .))

  tested = biomarker_report_raw %>%
    enter_object(biomarker) %>%
    gather_array('array.index.2') %>%
    spread_all %>%
    as_tibble %>%
    distinct %>%
    rename_all(.funs = ~ gsub('report\\.|ma_', '', .)) %>%
    filter(genenamesomatic == .gene, alterationtype == 'mutation') %>%
    syhelpr::map_oc(mutationstate, 'mutation')

  positive <<- tested %>%
    filter(mutationstate == 'Mutated') %>%
    left_join(
      biomarker_report,
      by = c('patientid', 'array.index'),
      relationship = "many-to-many"
    ) %>%
    mutate_at(.variant_cols, ~ trimws(na_if(., ''))) %>%
    distinct(
      patientid,
      gene = coalesce(genenamesomatic, genenamegermline),
      alterationtype,
      variantnamec,
      variantnamep,
      variantnameexon,
      reportcopy
    ) %>%
    mutate(
      reportcopy = gsub(
        '(no t790m mutation|genes tested.*|clinical significance.*|list of egfr exon.*|interpretation:.*|technical note:.*|\\(|\\))',
        ' ',
        tolower(reportcopy)
      ),
      reportcopy = gsub('(p\\.|c\\.)', ' ', reportcopy)
    )

  ## Output ==
  output$results$collected = positive %>%
    # select(-reportcopy) %>%
    filter(if_any(.variant_cols, ~ !is.na(.))) %>%
    distinct() %>%
    mutate(
      source = 'collected',
      variantnameexon = coalesce(str_extract(reportcopy, e_pattern), variantnameexon),
      varianttype = coalesce(
        str_extract(paste(variantnameexon, reportcopy),
                    '(deletion|insertion)'),
        str_extract(variantnamep,
                    'deletion|insertion|delins|del|ins')
      ),
      variantnameexon = str_extract(variantnameexon, '(1[0-9]|2[0-8]|[1-9])')
    )

  extracted_words = positive %>%
    filter(!is.na(reportcopy),
           !patientid %in% output$results$collected$patientid) %>%
    group_by(patientid,
             gene,
             variantnameexon,
             reportcopy) %>%
    unnest_tokens(
      word,
      reportcopy,
      drop = FALSE,
      token = "regex",
      pattern = "[\'\",:;!?=()\\[\\]\\s]"
    ) %>%
    anti_join(get_stopwords())

  output$results$extracted = extracted_words %>%
    mutate(
      variantnamec = gsub("^([^c.].*)", "c.\\1", str_extract(word, c_pattern)),
      variantnamep = gsub("^([^p.].*)", "p.\\1", str_extract(word, p_pattern)),
      variantnameexon = coalesce(variantnameexon, str_extract(reportcopy, e_pattern)),
      codingchangeid = str_extract(word, 'nm_([0-9]){6}(\\.([0-9]){0,4}){0,1}'),
      varianttype = coalesce(
        str_extract(paste(variantnameexon, reportcopy), '(deletion|insertion)'),
        str_extract(variantnamep, 'deletion|insertion|delins|del|ins')
      ),
      variantnameexon = str_extract(variantnameexon, '(1[0-9]|2[0-8]|[1-9])'),
    ) %>%
    select(-word) %>%
    ungroup %>%
    group_by(patientid,
             gene,
             reportcopy) %>%
    summarise(across(
      c(
        variantnamec,
        variantnamep,
        codingchangeid,
        variantnameexon,
        varianttype
      ),
      ~ first(na.omit(.))
    )) %>%
    ungroup %>%
    filter(if_any(c(.variant_cols, 'variantnameexon'), ~ !is.na(.))) %>%
    distinct %>%
    mutate(source = 'extracted')

  output$results$unknown = positive %>%
    filter(if_all(c(.variant_cols, 'variantnameexon'), ~ is.na(.))) %>%
    mutate(source = 'unknown') %>%
    distinct()

  output$patient_lists$tested = unique(tested$patientid)
  output$patient_lists$positive = unique(positive$patientid)

  ## Counts ==
  counts = c(
    'tested' = tested %>% count_patients %>% pull(patients),
    'positive' = positive %>% count_patients %>% pull(patients),
    'collected' = output$results$collected %>% count_patients %>% pull(patients),
    'no_detail' = positive %>% count_patients %>% pull(patients) - output$results$collected %>% count_patients %>% pull(patients),
    'with_report' = positive %>%
      filter(
        !is.na(reportcopy),!patientid %in% output$results$collected$patientid
      ) %>%
      count_patients %>%
      pull(patients),
    'extracted' = output$results$extracted %>% count_patients %>% pull(patients)
  )

  ## Messages ==
  .messages = c(
    'Results for {toupper(gene)}:',
    '{counts[["tested"]]} tested for {toupper(gene)} mutation',
    '{counts[["positive"]]} ({syhelpr::as_percent(counts[["positive"]] / counts[["tested"]], round = 1)}) with positive result',
    '{counts[["collected"]]} of {counts[["positive"]]} ({syhelpr::as_percent(counts[["collected"]] / counts[["positive"]])}) with variant detail collected',
    '{counts[["no_detail"]]} of {counts[["positive"]]} ({syhelpr::as_percent({counts[["no_detail"]]} / {counts[["positive"]]})}) with no detail collected',
    '{counts[["with_report"]]} of {counts[["no_detail"]]} ({syhelpr::as_percent({counts[["with_report"]]} / {counts[["no_detail"]]})}) with no detail collected, but RAW report available',
    '{counts[["extracted"]]} of {counts[["with_report"]]} ({syhelpr::as_percent({counts[["extracted"]]} / {counts[["with_report"]]}, round = 1)}) with variant detail extracted from RAW report'
  )
  lapply(.messages, function(x)
    message(glue::glue(x))) %>% invisible

  return(output)
}
