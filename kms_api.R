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

kms_normalize_variant <-
  function(df,
           codingchange_col = NULL,
           aminoacidchange_col = NULL,
           gene_col = 'gene',
           cols = c()) {

    api_root = 'https://ci-knowledge-sqa.syapse.com/api/V1/variant_normalization?'
    .gene_col = sym(gene_col)
    .col_aminoacid = sym(aminoacidchange_col)
    .col_codingchange = sym(codingchange_col)
    results_list = list()

    .data = df %>%
      distinct(gene = gsub(' ', '', {{ .gene_col }}),
               codingchange = gsub(' ', '', {{ .col_codingchange }}),
               aminoacidchange = gsub(' ', '', {{ .col_aminoacid }})) %>%
      as.list

    variant_list = .data

    for(i in 1:length(.data[[1]])) {
      svMisc::progress(i, length(.data[[1]]))
      params = lapply(list(
        gene = toupper(.data[['gene']][[i]]),
        codingchange = .data[['codingchange']][[i]],
        aminoacidchange = .data[['aminoacidchange']][[i]]
      ), function(x) x[!is.na(x)])

      results_list[[i]] = kms_api_call(api_root, params)
      
      if (i == length(.data[[1]])) message("Done!")
    }

    output = df %>%
      distinct(
        patientid,
        {{ .gene_col }} := toupper(gsub(' ', '', {{ .gene_col }})),
        {{ .col_codingchange }} := gsub(' ', '', {{ .col_codingchange }}),
        {{ .col_aminoacid }} := gsub(' ', '', {{ .col_aminoacid }}),
        across(any_of(cols))) %>%
      left_join(results_list %>% bind_rows,
                by = c(setNames('codingchange_raw', codingchange_col),
                       setNames('aminoacidchange_raw', aminoacidchange_col),
                       setNames('gene_raw', gene_col)))

    return(output)
}
