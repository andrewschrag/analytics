#' SPMD Write/Append/Drop Table
#'
#' Write/Append/Drop a table in the ca schema of the SPMD
#' 
#' NOTE: It is important to understand though that data only flows from production to the clone during the 
#' cloning process - the clone is completely overwritten during the cloning process and data added to/updated 
#' in production after the clone was last refreshed will not appear in the clone until it is next refreshed. 
#' Data added to/updated in the clone will never appear in production unless explicitly copied to it. 
#'
#' @param .data a data frame to write/append to the table
#' @param table_name a character vector representing the name of the table in spmd
#' @param spmd_write_con a DBIConnection write object
#' @param overwrite logical. If TRUE, overwrite the existing table.
#' @param convert_ids logical. If TRUE, convert the spmd 'id$' columns (if any) to a uuid.
#' @param silent logical. If TRUE, do not output any messages
#' @param fail_if_missing logical. If TRUE, rlang::abort if the table_name does not exist.
#' @param slack_notification logical. If TRUE, sends a slack notification to the #ca_schema channel via the `send_slack_message()` function. 
#'
#' @return NULL
#' @export
spmd_write_table <- function(.data, table_name, spmd_write_con = spmd_con('prod', write=TRUE),
                             overwrite=FALSE, convert_ids=TRUE, silent=FALSE, slack_notification=FALSE) {
  
  if(DBI::dbExistsTable(spmd_write_con, table_name) & !overwrite) {
    rlang::abort(glue::glue(
      "'{table_name}' already exists.",
      "`table_name` must not already exist or overwrite must be set to TRUE."
    ))
  }
  message("spmd_write_table - table exists check")
  
  if(as.numeric(lobstr::obj_size(.data)) > 1000000) {
    rlang::warn(glue::glue("{as_label(.data)} is size {lobstr::obj_size(.data)}"))
  }
  
  if (grepl("clone", spmd_write_con@info$servername)) {
    rlang::warn("Data added to clone will not be reflected in production and will be overwritten at the next clone refresh.")
  }
  message("spmd_write_table - copy_to executing")
  copy_to(
    spmd_write_con, .data, dbplyr::in_schema("ca", table_name), temporary = FALSE, overwrite = overwrite
  )
  message("spmd_write_table - copy_to executed")
  
  if(convert_ids) {
    message("spmd_write_table - converting ids")
    # valid_id_cols = uuid id columns + all table_name + 'id' combinations.
    valid_id_cols = tbl(spmd_write_con, dbplyr::in_schema("information_schema","columns")) %>% filter(
      grepl('id$',column_name),
      data_type == 'uuid',
      !grepl('_conceptid$',column_name),
      column_name != 'id'
    ) %>% select(column_name) %>%
      union_all(
        # table_name + column_name where column_name == 'id'
        tbl(spmd_write_con, dbplyr::in_schema("information_schema","columns")) %>%
          filter(column_name == 'id',data_type == 'uuid') %>%
          mutate(column_name = stringr::str_c(table_name,column_name)) %>%
          select(column_name)
      ) %>%
      distinct(column_name) %>% pull
    
    for(col in intersect(colnames(.data),valid_id_cols)) {
      message(glue::glue("spmd_write_table - converting ids, {col}, making it uuid"))
      DBI::dbExecute(
        spmd_write_con,
        glue::glue("ALTER TABLE ca.{table_name} ALTER COLUMN {col} TYPE uuid using {col}::uuid;")
      )
      message(glue::glue("spmd_write_table - converting ids, {col}, making it index"))
      DBI::dbExecute(
        spmd_write_con,
        glue::glue("CREATE INDEX IF NOT EXISTS idx_{table_name}_{col} ON ca.{table_name}({col});")
      )
    }
  }
  message(glue::glue("spmd_write_table - granting select to users"))
  read_users = c('PUBLIC')
  for(read_user in read_users) {
    DBI::dbExecute(spmd_write_con,glue::glue('GRANT SELECT ON ca.{table_name} TO {read_user}'))
  }
  message(glue::glue("spmd_write_table - sending messages"))
  if(!silent) {
    rlang::inform(glue::glue("CREATED TABLE: ca.{table_name}"))
  }
  if(slack_notification) {
    if(!overwrite) { 
      send_slack_message(glue::glue("ca.{table_name} has been created."))
    }
    else {
      send_slack_message(glue::glue("ca.{table_name} has been updated."))
    }
  }
}