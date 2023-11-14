#
# NOTE(ian) 04/21/2023:
# This utility file is intended to hold functions that can
# be used for data validation and other similar quality checks.
#


#' Check that there are no duplicate patient IDs in a given dataset.
#' If the check fails, the function aborts with an error.
#'
#' @param data The dataframe you want to validate. We expect the dataframe to have a "patientid" column.
#' @param group_number A string representing the group number you wish to report on, in the event the check fails.
#' @return NULL
#' @export
#'
#' @examples
#' confirm_unique_patients(lung$ads, "2")
#'
confirm_unique_patients <- function(data, group_number) {
    data %>%
      confirm_unique(patientid,
                     error_msg = "{as_label(.col)} is not unique in dataset after group {group_number} items incorporated.",
                     return_df = FALSE) %>%
      invisible()
}
