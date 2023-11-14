
# Function for writing analytical datasets to S3. Must specify the dataset name and type to ensure it is 
# written to the correct location. type must be one of "deid" (if de-identified dataset is being saved) or 
# "id" (if identified/internal dataset is being saved).
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
    confirm_unique(patientid,
                   error_msg = "patientid is not unique in final data frame! ADS not saved!") %>%
    write_parquet(sink = sink)

  message(glue::glue("Finished writing {type} {dataset_name} to S3"))
}
