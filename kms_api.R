kms_api_call <- function(api_root, params = list()){
  res = httr::GET(api_root, query = params)
  output = try(RJSONIO::fromJSON(
    httr::content(res, as="text", encoding = 'UTF-8'), warn = F)$data %>%
      map_df(.f = ~.x))
  return(output)
}

kms_normalize_variant <-
  function(df,
           codingchange_col = NULL,
           aminoacidchange_col = NULL,
           gene_col = 'gene',
           cols = c()) {

    .gene_col = sym(gene_col)
    .col_aminoacid = sym(aminoacidchange_col)
    .col_codingchange = sym(codingchange_col)
    results_list = list()

    .data = df %>%
      distinct(gene = gsub(' ', '', {{ .gene_col }}),
               codingchange = gsub(' ', '', {{ .col_codingchange }}),
               aminoacidchange = gsub(' ', '', {{ .col_aminoacid }})) %>%
      as.list
    api_root = 'https://ci-knowledge-sqa.syapse.com/api/V1/variant_normalization?'

    for(i in 1:length(.data[[1]])) {
      params = lapply(list(
        gene = toupper(.data[['gene']][[i]]),
        codingchange = .data[['codingchange']][[i]],
        aminoacidchange = .data[['aminoacidchange']][[i]]
      ), function(x) x[!is.na(x)])

      results_list[[i]] = kms_api_call(api_root, params)
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
                       setNames('gene', gene_col)))

    return(output)
}
