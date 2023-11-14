# General functions 

# Adds row to attrition table
add_attrition <- function(data, criterion, attrition = NULL) {
  if (is.null(attrition)) {
    data %>% 
      count_patients(patients = "patientid") %>% 
      mutate(step = criterion)
  } else {
    attrition %>% 
      bind_rows(data %>% 
                  count_patients(patients = "patientid") %>% 
                  mutate(step = criterion))
  }
}

# Retrieves patients and their tumor information associated with a given cohort/cohort_type
# The tumor record associated with the patients' cohort membership is prioritized because 
#   1) the tumor algorithm prioritizes correctly (source ranking of MA > registry > clinical and then the 
#      earliest diagnosis date per patient per tumor type)
#   2) this is the record displayed across products
# We will use the patientid and tumortype from mdr.tumor to trace back to msmdro.tumor in another function.
get_cohort <- function(cohort, cohort_type = "curated", schema = "mdr", spmd = spmd_con()) {
  cohort_type <- tolower(cohort_type)
  cohort <- tolower(cohort)
  # Cohort terminologies are updating: Spectrum -> Structured and Spotlight -> Curated
  if (cohort_type %in% c("spotlight", "spectrum")) {
    lifecycle::deprecate_warn(when = "2023-02-17", 
                              what = "SPMD::get_cohort(cohort_type = 'values should change per the below')", 
                              details = c("'spectrum' should change to 'structured'",
                                          "'spotlight' should change to 'curated'"))
  }
  # Give a warning if the given cohort type is not in the possible set and map between old and new 
  # terminologies depending on input and available options
  # There will be a transition period between the old terms and the new terms across the environments 
  # The new terms will be present in the lower environments first and then the upper environments later 
  # This will handle that transition period without forcing us to constantly change our arguments. We will 
  # just get a warning.
  possible_cohort_types <- tbl(spmd, in_schema("cohorts", "cohort")) %>% 
    pull(cohorttype) %>% 
    tolower() %>% 
    unique()
  if (!cohort_type %in% possible_cohort_types) {
    base_warning <- glue::glue("In an SPMD environment ('{str_extract(spmd@info$servername, 'prod|qc|sqa|clone')}') ")
    if (cohort_type %in% c("structured", "curated")) {
      # We are in an environment that has the old terms still
      old_term <- if (cohort_type == "structured") "spectrum" else "spotlight"
      rlang::warn(paste0(base_warning, glue::glue("that uses '{old_term}' instead of '{cohort_type}', converting '{cohort_type}' to '{old_term}' for the get_cohort query.")))
      cohort_type <- old_term
    } else if (cohort_type %in% c("spectrum", "spotlight")) {
      # We are in an environment that has the new terms
      new_term <- if (cohort_type == "spectrum") "structured" else "curated"
      rlang::warn(paste0(base_warning, glue::glue("that uses '{new_term}' instead of '{cohort_type}', converting '{cohort_type}' to '{new_term}' for the get_cohort query.")))
      cohort_type <- new_term
    } else {
      rlang::warn(paste0(base_warning, glue::glue("where cohort_type = '{cohort_type}' is not present! Possible values for cohort_type are:\n- {paste0(possible_cohort_types, collapse = '\n- ')}")))
    }
  }
  # Warning if the cohort name does not exist in the environment
  possible_names <- tbl(spmd, in_schema("cohorts", "cohort")) %>% 
    pull(name) %>% 
    tolower() %>% 
    unique()
  if (!cohort %in% possible_names) {
    rlang::warn(glue::glue(base_warning, glue::glue("where cohort = '{cohort}' is not present! Possible values for cohort are:\n- {paste0(possible_names, collapse = '\n- ')}")))
  }
  # Query for the patients in the specified cohort
  tbl(spmd, in_schema("cohorts", "cohort")) %>%
    filter(tolower(name) == cohort,
           tolower(cohorttype) == cohort_type) %>%
    select(cohortid = id, cohortname = name, cohorttype) %>%
    inner_join(tbl(spmd, in_schema("cohorts", "currentcohortmember")), 
               by = "cohortid") %>%
    select(patientid, tumorid, cohorttype, cohortname) %>%
    inner_join(tbl(spmd, in_schema(schema, "tumor")) %>%  select(-any_of(c("provenance"))),
               by = c("tumorid" = "id", "patientid"))
}

# Maps OC yes, no, and unknown codes 
convert_oc_yes_no <- function(.data, column) {
  column <- enquo(column)
  .data %>% 
    mutate(!!column := case_when(!!column == "1" ~ "Yes",
                                 !!column == "0" ~ "No", 
                                 !!column == "unk" ~ "Unknown",
                                 !!column == "" ~ NA_character_,
                                 TRUE ~ !!column))
}

# Covnert empty string to NA_character_ (default)
convert_empty_string <- function(.data, column, to = NA_character_) {
  column <- enquo(column)
  .data %>% 
    mutate(!!column := if_else(!!column == "", to, !!column, !!column))
}

# Applies specified OC mapping to the given column. When other = TRUE, replaces "Other" with values in the 
# corresponding "other" column, and errors if "other" column not found or multiple "other" columns found for 
# the given column. Also, converts empty strings ("") to NAs.
map_oc <- function(.data, column, valueset_name, other = FALSE, to = NA_character_, silent = FALSE) {
  # determine column variables
  column_base <- deparse(substitute(column))
  cols <- names(.data)
  column_other <- cols[grepl(column_base, cols) & grepl("other", cols)]
  if (other) {
    if (is_empty(column_other)) {
      stop(str_glue("other = TRUE and no matching 'other' column found for {column_base}"))
    } else if (length(column_other) > 1) {
      stop(paste0(c(str_glue("other = TRUE and multiple 'other' columns found for {column_base}:"), column_other), collapse = "\n\t"))
    }
    column_other <- rlang::sym(column_other)
  }
  column <- enquo(column)
  # clean-up the data and apply the mapping
  .data %>%
    convert_empty_string(!!column) %>%
    {if (other) convert_empty_string(., !!column_other) else .} %>%
    mutate(!!column := coalesce(get_oc_values(!!column, valueset_name, silent = silent), !!column)) %>%
    {if (other) mutate(., !!column := ifelse(!!column == "Other", !!column_other, !!column)) else .} %>% 
    mutate(!!column := str_trim(!!column))
}

# Removes the "ma_" prefix on OC columns when extracted from openclinica.formdata
rename_oc <- function(.data) {
  .data %>% 
    rename_with(~ str_remove(., "^ma_"))
}

# Adds sorting to grade while handling unexpected values
sort_grade <- function(.data, grade) {
  grade <- enquo(grade)
  .data %>% 
    # this implementation of add_sorting handles unexpected values not present in the sort_vector
    # by argument isn't needed
    add_sorting_custom_handle_missing(!!grade, 
                                      c("Grade IV",
                                        "G3: Poorly differentiated, undifferentiated",
                                        "Grade III",
                                        "G2: Moderately differentiated",
                                        "Grade II",
                                        "High grade", # either G2 or G3 and > G1
                                        "G1: Well differentiated",
                                        "Grade I",
                                        "Low grade",
                                        "GB: Borderline Tumor",
                                        "Grade cannot be assessed (GX); Unknown",
                                        "Grade/differentiation unknown, not stated, or not applicable",
                                        # the below do not apply to ovarian and so are placed at the bottom
                                        "B-cell",
                                        "T-cell",
                                        "Null cell",
                                        "NK (natural killer) cell",
                                        "8", # same NK cell -> missing from OC map
                                        # These are thought to be bad data
                                        "M", 
                                        "G",
                                        "D",
                                        "C",
                                        "A"))
} 

# Registry data occasionally has multiple records per patient within the cancer type of interest, we need
# to prioritize amongst by selecting ovarian cancer cases first, then earliest diagnosis. Selection of 
# ovarian cancer cases first is due to the Spotlight definition, which limits tumorttype = "GU: Ovary" as
# of 2022-10-21. It will be fixed sometime in Q4 2022. 
# We "lose" a handful of patients due to NA diagnosisdates, but they are not included in the Spotlight as
# is.
prioritize_registry_ovarian <- function(.data) {
  lifecycle::deprecate_warn("2022-11-15",
                            "prioritize_registry_ovarian()",
                            "get_registry()")
  .data %>% 
    mutate(schema_rank = ifelse(stagingalgorithmschema == "ovary", 1, 2),
           dateofdiagnosis_year = as.numeric(str_extract(dateofdiagnosis, "^\\d{4}")),
           dateofdiagnosis = ymd(dateofdiagnosis, quiet = TRUE)) %>% 
    distinct() %>% 
    group_by(patientid) %>% 
    filter(schema_rank == min(schema_rank),
           dateofdiagnosis == suppressWarnings(min(dateofdiagnosis, na.rm = TRUE)),
           dateofdiagnosis_year == suppressWarnings(min(dateofdiagnosis_year, na.rm = TRUE))) %>% 
    ungroup() %>% 
    select(-stagingalgorithmschema, -schema_rank, -dateofdiagnosis, -dateofdiagnosis_year)
}

# We can supplement the data in mdr with registry data. Patients can have multiple registry records - some 
# which pertain to the cancer they are being included with and others that do not (i.e., related to a 
# different cancer type or from a prior or later version of the cancer). We need a method to extract the 
# data from the registry record(s) that correspond most closely to the view of the patent's cancer that we 
# have in the mdr.tumor table (i.e., that align by primarysite and diagnosis date). We do this by 
# 1) filtering to specific registry cancers qualifying for this cancer type, then 
# 2) picking registry records that match the patient's selected primary site over those that do not, then 
# 3) picking registry records that are closest to the diagnosis date, then 
# 4) picking registry records that are on/before the diagnosis date
# NOTE - After steps 1-3, there can still be ties and all corresponding registry records are returned
get_registry <- function(cohort_query, cancer_type, cols = c()) {
  
  # Pull in the data
  registry_data <- tbl(cohort_query$src$con, in_schema("ca", "registry_naaccr")) %>% 
    # stagingalgorithmschema = Describes the SEER Staging Algorithm (CS TNM or EOD) Schema ID. This specifies 
    # which staging system is used for the staging algorithms. We can use it to filter to specific cancer types
    filter(grepl(cancer_type, stagingalgorithmschema, ignore.case = TRUE)) %>% 
    inner_join(cohort_query %>%
                 select(patientid, primarysite_code, diagnosisdate),
               by = c("patientid")) %>% 
    select(patientid, stagingalgorithmschema, dateofdiagnosis, primarysite, primarysite_code, diagnosisdate, all_of(cols)) %>%
    distinct() %>%
    collect() %>%
    # remove "." from ICD codes to allow us to compare the registry primarysite code to the mdr.tumor.primarysite_code
    mutate(primarysite = str_replace(primarysite_code, "\\.", ""), 
           primarysite_code = str_replace(primarysite_code, "\\.", ""), 
           # extract year of registry dateofdiagnosis in case the date is incomplete
           dateofdiagnosis_year = as.numeric(str_extract(dateofdiagnosis, "^\\d{4}")),
           # convert registry dateofdiagnosis to date object for difference calculation
           dateofdiagnosis = ymd(dateofdiagnosis, quiet = TRUE),
           # calculate the difference between the registry diagnosis date and the mdr.tumor.diagnosisdate
           diff = day_difference(dateofdiagnosis, diagnosisdate))
  # Apply the prioritization
  registry_data_prioritized <- registry_data %>% 
    group_by(patientid) %>% 
    arrange(patientid,
            primarysite != primarysite_code, # primary sites matches over mismatches 
            abs(diff), # closest diagnosis dates
            diff < 0, # registry diagnosis dates on/before over ones after (diff is negative when registry date is after diagnosis date) and F comes before T in arrange
            dateofdiagnosis_year) %>% # in case any just have the year
    slice_head() %>% 
    select(patientid, stagingalgorithmschema, dateofdiagnosis, dateofdiagnosis_year, primarysite) %>% 
    ungroup()
  # Filter to corresponding matches in case there are still ties (user musst decide how to de-dupe specific
  # data elements)
  registry_data %>% 
    semi_join(registry_data_prioritized, 
              by = c("patientid", "stagingalgorithmschema", "dateofdiagnosis", "primarysite", "dateofdiagnosis_year")) %>% 
    select(-primarysite_code, -diagnosisdate, -diff)
}

# Difference in months
month_difference <- function(time1, time2, truncate = FALSE) {
  difference <- interval(time1, time2) / months(1)
  if (truncate) {
    floor(difference)
  } else {
    difference
  }
}

# Converts days to months and leaves number un rounded.
convert_day_to_months <- function(x) {
  x/30.41
}

# Debugging
dupes <- function(.data) {
  .data %>% 
    group_by(patientid) %>% 
    filter(n() > 1) %>% 
    arrange(patientid)
}

# Coverts a vecotr of ICD10 codes to a regex-compatible string
create_icd10_regex <- function(codes) {
  paste0(c(str_remove_all(codes, "\\."), str_replace_all(codes, "\\.", "\\\\.")), collapse = "|")
}

# Extracts year for a given date using all OpenClinica columns associated with that date (e.g., 
# therapystartdate, therapystartdateym, and therapystartdatey). All columns must be present in .data. 
extract_year_oc <- function(.data, date_col) {
  # OC date column names
  date_col_nm <- deparse(substitute(date_col))
  date_col <- enquo(date_col)
  gran_col_nm = rlang::sym(str_c(date_col_nm, "_granularity"))
  y_col_nm <- rlang::sym(str_c(date_col_nm, "y"))
  ym_col_nm <- rlang::sym(str_c(date_col_nm, "ym"))
  extracted_year <- rlang::sym(str_c(date_col_nm, "_year"))
  # Check columns exist
  .data = .data %>% 
    check_for_column(!!date_col, .class = c("character", "Date")) %>%
    check_for_column(!!y_col_nm, .class = c("character", "Date")) %>%
    check_for_column(!!ym_col_nm, .class = c("character", "Date"))
  # Extract year
  .data %>% 
    mutate(!!extracted_year := case_when(!!gran_col_nm == "DAY" ~ year(!!date_col),
                                         !!gran_col_nm == "MONTH" ~ as.numeric(str_extract(!!ym_col_nm, "^\\d{4}")),
                                         !!gran_col_nm == "YEAR" ~ as.numeric(str_extract(!!y_col_nm, "^\\d{4}")),
                                         TRUE ~ NA_real_))
}

# Returns the last of the month for given date
impute_date_last_of_month <- function(date) {
  lubridate::ceiling_date(lubridate::ymd(date), unit = "months") - lubridate::days(1)
}

# Returns the first of the month for given date
impute_date_first_of_month <- function(date) {
  trunc(date, units = "months")
}

# For reporting patient counts from data frames
# Add invisible(.data) call or return_df argument to make pipe compatible? 
report_patient_count <- function(.data, message = "{patients} patients") {
  if (!grepl("\\{patients\\}", message)) rlang::warn("{patients} not included in message, no count displayed!")
  patients <- .data %>% 
    count_patients(na_rm = TRUE) %>% 
    pull()
  print(glue::glue(message))
}

# Replaces -Inf or Inf dates (which are confusingly displayed as NA when data frames are printed to the console)
# with NA_Date_. -Inf and Inf are returned by min and max when na.rm = TRUE and all dates are NA. This is 
# intended to correct that side-effect.
replace_inf_dates <- function(.data) {
  .data %>% 
    mutate(across(where(is.Date), ~ if_else(is.infinite(.x), NA_Date_, .x, .x)))
}

# Full join one data frame to another and replace missing values in column with na_fill 
# Intended for filling NAs "at the site of definition" of ADS variables as opposed to at the end when all 
# data frames are joined together. .data and ads must have patientid.
join_replace_na <- function(.data, ads, column, na_fill = "Unknown") {
  column <- enquo(column)
  .data %>% 
    full_join(ads %>% 
                select(patientid) %>% 
                distinct(), 
              by = "patientid") %>% 
    mutate(!!column := replace_na(!!column, na_fill))
}

# Using sourceschema, prioritizes MA > registry > all others
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

# Converts three letter amino acid abbreviation to single letter abbreviation ("Ile" -> "I")
# Some of the raw OC data comes with three letter abbreviations
convert_three_letter_aa <- function(.data, three_letter_column) {
  three_letter_column <- enquo(three_letter_column)
  .data %>% 
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Ala", "A")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Arg", "R")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Asn", "N")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Asp", "D")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Cys", "C")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Gln", "Q")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Glu", "E")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Gly", "G")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "His", "H")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Ile|lle", "I")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Leu", "L")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Lys", "K")) %>% 
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Met", "M")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Phe|phe", "F")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Pro", "P")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Pyl", "O")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Ser", "S")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Sec", "U")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Thr", "T")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Trp", "W")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Tyr", "Y")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Val", "V")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Asx", "B")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Glx", "Z")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Xaa", "X")) %>%
    mutate(!!three_letter_column := str_replace_all(!!three_letter_column, "Xle", "J"))
}

# Function to summarize biomarker status according to normalized call values and a status map matching those
# values. Selects the earliest test date per patient per normalized call, and then prioritizes the biomarker 
# call per patient according to normalized_call_levels. Finally, it maps the call_normalized values in .data
# to the values in status_map (by names in status_map).
# .data must contain the following columns: patientid, call_normalized, biomarkertype, gene, sourceschema,
# call, and call_rawvalue (the last 5 help determine which biomarker summary produced a warning)
# .data can contain call_normalized values not present in normalized_call_levels or status_map. These will 
# produce warnings to prompt users to fix them, but unexpected and un-normalized values in .data are allowed
# to propagate through to the summary. By default, these are deprioritized in the sorting (i.e., selected 
# last if present).
# Per discussion with Rayna 2022-12-06, we prioritize Positive > Unknown > Negative by default. This is the 
# most conservative method and reduces the false negatives we may have.
summarize_biomarker_status <- function(.data, 
                                       normalized_call_levels = c("Positive", "Unknown", "Negative"), 
                                       status_map = list(Positive = "Positive", Negative = "Negative", Unknown = "Unknown (all tests inconclusive)")) {
  # Check for normalized call column
  .data %>% 
    check_for_column(call_normalized, return_df = FALSE) %>% 
    invisible()
  # Check to see if any normalized values are missing from the expected levels
  normalized_call_levels_additions <- .data %>% 
    filter(!call_normalized %in% normalized_call_levels) %>% 
    pull(call_normalized) %>% 
    sort()
  if (!identical(character(0), normalized_call_levels_additions)) {
    rlang::warn("Unexpected normalized call values in summarize_biomarker_status! Please update the `normalized_call_levels` argument or the normalization to handle:")
    .data %>% 
      filter(!call_normalized %in% normalized_call_levels) %>% 
      count(biomarkertype, gene, sourceschema, call_normalized, call, call_rawvalue) %>% 
      print(n = 50)
  }
  # Check to see if the status map names and expected normalized call levels align
  status_map_names <- names(status_map)
  if (!identical(sort(status_map_names), sort(normalized_call_levels))) {
    types <- .data %>% pull(biomarkertype) %>% unique() %>% sort()
    genes <- .data %>% pull(gene) %>% unique() %>% sort()
    rlang::warn(glue::glue("Names of `status_map` and values in `normalized_call_levels` do not align! Occurs for biomarker type(s) {paste0(types, collapse = ', ')} in gene(s) {paste0(genes, collapse = ', ')}. Please update either or both:
                           ... normalized_call_levels = {paste0(normalized_call_levels, collapse = ', ')}
                           ... status_map names = {paste0(status_map_names, collapse = ', ')}"))
  }
  # Convert status_map to data frame for joining
  status_map <- status_map %>% 
    enframe() %>% 
    unnest(value) %>% 
    rename(call_normalized = name,
           status_value = value)
  # Summarize the data finding the earliest report date for each distinct normalized call per patient and 
  # then selecting the normalized call with the highest level 
  # We include unexpected/un-normalized calls, but put them at the end of the ordering
  normalized_call_levels <- c(normalized_call_levels, normalized_call_levels_additions)
  .data %>% 
    # take the earliest report date per patient per normalized call
    group_by(patientid, call_normalized) %>% 
    summarize(test_date = suppressWarnings(min(reportdate, na.rm = TRUE)),
              .groups = "drop") %>% 
    convert_inf_date(test_date) %>% 
    # prioritize normalized call according to the order in normalized_call_levels
    add_sorting(call_normalized, "custom", normalized_call_levels) %>% 
    group_by(patientid) %>% 
    arrange(call_normalized) %>% 
    slice_head() %>% 
    ungroup() %>% 
    # map normalized call to the values specified in status_map according to names in the status map
    left_join(status_map,
              by = "call_normalized") %>% 
    # if mapping isn't available, we keep the "normalized" call value (this is actually and unexpected and
    # likely un-normalized call)
    mutate(status = coalesce(status_value, as.character(call_normalized))) %>%
    select(patientid, status, test_date)
}

# Functions for handling JSON data
# # Can we create a function to take arbitrary json input and split it into a dataframe?
# # This might take too long maybe takes too long?
# example2 <- run_query("SELECT 
#           patientid,
#           formdata::json -> 'family_personal_history' AS family_personal_history
#           FROM openclinica.formdata 
#           WHERE formdata::json -> 'family_personal_history' IS NOT NULL") %>% 
#   arrange(patientid)
# example2
# example2 %>% 
#   mutate(family_personal_history = purrr::map(family_personal_history, jsonlite::fromJSON)) %>%
#   unnest_wider(family_personal_history) %>% 
#   unnest_wider(group) %>% 
#   unnest_wider(family)
# example2 <- run_query("SELECT 
#           patientid,
#           formdata::json -> 'family_personal_history' AS fh
#           FROM openclinica.formdata 
#           WHERE formdata::json -> 'family_personal_history' -> 'group' IS NOT NULL") %>% 
#   arrange(patientid)
# example2
# example2_p1 <- example2 %>% 
#   mutate(fh = purrr::map(fh, jsonlite::fromJSON))
# example2_p2 <- example2_p1 %>% 
#   unnest_wider(fh)
# names(example2_p2)
# example2_p2 %>% 
#   unnest_wider(group)


get_antineoplastics <- function(patients = NULL, spmd = spmd_con(), silent = TRUE) {
  message("get_antineplastics")
  message("using spmd = ", spmd@info$username)
  antineoplastics_regex <- dplyr::tbl(spmd, dbplyr::in_schema("ca","map_atc")) %>%
    dplyr::filter(class.level2 %in% c("L01: antineoplastic agents","L02: endocrine therapy")) %>%
    dplyr::pull(all_names_regex) %>%
    unique() %>%
    paste0(collapse = "|")
  message(glue::glue("*** Build med regex ***"))
  antineoplastics_regex <- build_med_regex(antineoplastics_regex,spmd=spmd)
  if (!silent) {
    rlang::inform(glue::glue("Using med_regex = '{antineoplastics_regex}'"))
  }

  antineoplastic_ts = str_split(antineoplastics_regex,"\\|")
  antineoplastic_ts = str_c(paste0("(",unlist(antineoplastic_ts),")"),collapse="|")
  antineoplastic_ts = str_replace_all(antineoplastic_ts," ","&") # Space symbols not indexed, thus required
  antineoplastics = dplyr::tbl(spmd, dbplyr::in_schema("ca","medications")) %>%
    dplyr::filter(
      dbplyr::sql(glue::glue("drugsearch_index_col @@ to_tsquery('simple','{antineoplastic_ts}') OR sourceschema = 'ma'"))
    )

  if(!is.null(patients)) {
    antineoplastics = antineoplastics %>%
      dplyr::filter(patientid %in% patients)
  }

  return(antineoplastics)
}

build_med_regex <- function(med_regex, spmd = spmd_con()) {
  message("using spmd = ", spmd@info$username)
  if (!is.character(med_regex)) {
    stop("med_regex is not a character vector")
  }
  if (length(med_regex) > 1) {
    stop("med_regex is not a character vector of length one (not a string)")
  }
  if (!spmd_list_tables("ca", spmd = spmd) %>% filter(table_name == "map_atc") %>% count() %>% pull(n) > 0) {
    stop("ca.map_atc does not exist! Please rebuild and try again:\n",
         "library(syhelpr)\n",
         "map_atc <- build_atc_table()\n",
         "spmd_write_table(map_atc,'map_atc',overwrite = TRUE)")
  }
  tbl(spmd, in_schema("ca", "map_atc")) %>%
    filter(grepl(med_regex, all_names_regex, ignore.case = TRUE)) %>%
    pull(all_names_regex) %>%
    str_split("\\|") %>%
    c(str_split(med_regex, "\\|"), .) %>%
    unlist() %>%
    unique() %>%
    paste0(collapse = "|", sep = "")
}

# Converts Inf and -Inf to NA_Date_
# summarize(date = min(date, na.rm = TRUE)) causes Inf to be returned when there are only NA dates (max 
# causes -Inf to be returned). We need to clean these up since a) writing Inf/-Inf dates to SPMD using 
# spmd_write_table (through dplyr::copy_to) doesn't work and b) filter(is.na(date)) doesn't match date == Inf
convert_inf_date <- function(.data, date_column) {
  date_column <- enquo(date_column)
  .data %>% 
    check_for_column(!!date_column, .class = "Date") %>% 
    mutate(!!date_column := if_else(!!date_column %in% c(Inf, -Inf), 
                                    NA_Date_,
                                    !!date_column,
                                    !!date_column))
}


get_demographics <- function(df, patients, spmd = spmd_con(), schema_name = 'mdr', cols = c(), ma_check = TRUE, inplace = FALSE, rename = TRUE, quiet = TRUE) {
  demo_query <- dplyr::tbl(spmd, in_schema(schema_name, 'patient')) %>%
    dplyr::select(
      patientid = id,
      birthdate,
      deceaseddate,
      isdeceased,  # TRUE, FALSE, or NA
      # vitalsign, # formatted isdeceased
      sex,
      race,
      ethnicity,
      sourceschema,
      tidyselect::all_of(cols)
    )
  
  # dataframe df is supplied
  if (!missing(df)) {
    if (!missing(patients)) {
      # patients is the name of the patientid column
      patient_col <- quo_name(enquo(patients))
    }
    else {
      # the patientid column is taken from several deault options in the data frame
      patient_col <- intersect(
        c("patient_id", "patientid", "patients", "patient", "subjectid"),
        names(df))[1]
    }
    if (is.na(patient_col))
      rlang::abort("missing patient column")
    patients <- pull(df, patient_col)
    dq <- demo_query %>% filter(patientid %in% patients) %>% collect()
  }
  # dataframe df is not provided
  else {
    if (!missing(patients)) {
      # 'patients' is a vector of patient IDs
      dq <- demo_query %>% filter(patientid %in% patients) %>% collect()
    } else {
      # gather demographics for all patients, if
      # there is no dataframe and no patients defined
      dq <- demo_query %>% collect()
    }
    df <- NULL  # this allows us to merge back to an "original" df when `inplace = TRUE`
    patient_col <- sym('patientid')
  }
  
  if (!quiet) print(str_glue('Joining, by = c("patientid" == "{patient_col}")'))
  if(!ma_check) message("ma_check != TRUE: Deceased dates not corrected by MA check.")
  
  if(ma_check & (all(dq$sourceschema != "ma"))) {
    rlang::abort(glue::glue(
      "ma_check requires at least one patient from the patient list be abstracted (MA)"))}
  
  ## MA date check (can be removed after platform fix)
  if (ma_check) {
    followup_ma <- DBI::dbGetQuery(paste0("SELECT
                                          patientid,
                                          formdata::json -> 'follow_up' ->> 'group' AS follow_up
                                          FROM openclinica.formdata
                                          WHERE formdata::json -> 'follow_up' IS NOT NULL
                                          AND patientid IN ('",
                                          paste(dq$patientid, collapse = "','"),
                                          "')"
    ),
    conn = spmd) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(follow_up = purrr::map(follow_up, jsonlite::fromJSON)) %>%
      tidyr::unnest_wider(follow_up) %>%
      dplyr::rename_with(~ str_remove(., "^ma_")) %>%
      dplyr::mutate(last_abstraction_date = as_date(dateoflastabstraction),
                    dateoflastcontact = as_date(dateoflastcontact)) %>%
      select(patientid, last_abstraction_date, vitalstatus, dateoflastcontact) %>%
      # prioritize most recent date
      group_by(patientid) %>%
      arrange(is.na(last_abstraction_date),
              desc(last_abstraction_date)) %>%
      slice_head() %>%
      ungroup() %>%
      distinct(patientid, vitalstatus, dateoflastcontact)
    
    dq <- dq %>%
      dplyr::left_join(followup_ma, by = "patientid") %>%
      dplyr::mutate(
        isdeceased = case_when(vitalstatus == 2 ~ TRUE,
                               vitalstatus == 1 & deceaseddate < dateoflastcontact ~ FALSE,
                               vitalstatus == 1 & deceaseddate > dateoflastcontact  ~ TRUE,
                               vitalstatus == 1 & is.na(deceaseddate) ~ FALSE,
                               TRUE ~ isdeceased),
        deceaseddate = case_when(vitalstatus == 2 ~ dateoflastcontact,
                                 vitalstatus == 1 & deceaseddate < dateoflastcontact ~ NA_Date_,
                                 vitalstatus == 1 & deceaseddate > dateoflastcontact ~ deceaseddate,
                                 TRUE ~ deceaseddate)) %>%
      dplyr::distinct()
    if ("vitalsign" %in% cols) {
      dq <- dq %>%
        dplyr::mutate(vitalsign = case_when(isdeceased == TRUE ~ "Deceased",
                                            isdeceased == FALSE ~ "Alive",
                                            is.na(isdeceased) ~ "Unknown"))
    }
    dq <- dq %>%
      dplyr::select(-vitalstatus, -dateoflastcontact)
  }
  
  if (!"sourceschema" %in% cols) {
    dq <- dq %>%
      select(-sourceschema)}
  
  dq %>%
    dplyr::mutate(
      sex = replace_na(sex, 'Unknown'),
      race = replace_na(race, 'Not Provided'),
      ethnicity = replace_na(ethnicity, 'Unknown'),
      # add race_ethnicity column
      race_ethnicity = case_when(
        ethnicity == "Hispanic/Latino" ~ "Hispanic/Latino",
        race == "Black or African American" ~ "Black",
        race %in% c("Asian","Native Hawaiian or Other Pacific Islander") ~ "Asian/Hawaiian/Pacific Islander",
        race == "White" ~ "Non-Hispanic White",
        race == "American Indian or Alaska Native" ~ "American Indian or Alaska Native",
        TRUE ~ "Other or Unknown"
      )) %>%
    dplyr::distinct()  %>%
    {`if`(inplace == TRUE & !is.null(df), dplyr::right_join(x=., y=df, by = c('patientid' = patient_col)), .)} %>%
    {`if`(rename == FALSE, rename(., !!patient_col := patientid), .)}
  
  
}


# Function to ID individuals by tumor number
id_ind_by_tumor_n <- function(.data, .tumor_n) {
  .data %>%
    group_by(patientid) %>% filter(n()==.tumor_n) %>%
    arrange(patientid)
}

# Function to check that there are no patientids with multiple entries when they should only have 1 entry
check_multiple_entries <- function(.data) {
  .data %>%
    group_by(patientid) %>%
    add_count(patientid, name = "n_patientids") %>%
    filter(n_patientids>1) %>%
    {if (dim(.)[1] > 0) rlang::abort("There is at least 1 patientid currently marked as resolved that currently have >1 tumorid.") else invisible()}
}

# An implementation of add_sorting(.data, .field, by = "custom", sort_vector = c("some", "values")) that 
# handles values in .field missing from the sort_vector. Missing values are added to the end of the 
# sort_vector. See summarize_biomarker_status for first implementation in the context of biomarker calls. 
# Implemented as part of https://syapse.atlassian.net/browse/CA-3727
# THIS SHOULD BE USED FOR ALL FIELDS THAT REQUIRE SOME SORT OF DE-DUPLICATION BY PATIENT (PRIORITIZATION) 
# BUT FOR WHICH UNDERLYING SOURCE DATA IS NOT CONSTRAINED TO A SPECIFIC VALUE SET (e.g., registry).
add_sorting_custom_handle_missing <- function(.data, 
                                              .field,
                                              sort_vector) {
  .field <- enquo(.field)
  # Check for normalized call column
  .data %>% 
    check_for_column(!!.field, return_df = FALSE) %>% 
    invisible()
  # Check for missing sort vector
  if (is.null(sort_vector)) {
    rlang::abort("sort_vector required")
  }
  # Check for missing values from the sort_vector, warn if any missing, and add them to the sort vector
  sort_additions <- .data %>% 
    filter(!(!!.field %in% sort_vector)) %>% 
    pull(!!.field) %>% 
    unique() %>% 
    sort() %>% 
    remove_na()
  if (length(sort_additions)) {
    missing_vals = paste(paste0("'", sort_additions, "'"), collapse = ",")
    rlang::warn(glue::glue("add_sorting_custom_handle_new: sort_vector does not contain all values in data. Missing {missing_vals}\n
                           These values will be appended to the end of sort_vector!"))
    
  }
  sort_vector <- c(sort_vector, sort_additions)

  .data %>% 
    add_sorting(!!.field, 
                by = "custom", 
                sort_vector = sort_vector)
}

# Wrapper of jsonlite::fromJSON that handles NA/NULL inputs
fromJSON_na <- function(x) {
  if (is.na(x) | is.null(x)) {
    data.frame()
  } else {
    jsonlite::fromJSON(x)
  }
}

# Converts all empty strings ("") in character columns to the to argument (default is NA)
convert_empty_strings <- function(.data, to = NA_character_) {
  .data %>% 
    mutate_if(is.character, list(~ ifelse(. == "", to, .))) 
}

# An updated version of tidyr is in prod but not QC. Therefore, separate_longer_delim currently works in prod but not QC. 
# Therefore, below is parts of the original function code from https://github.com/tidyverse/tidyr/blob/HEAD/R/separate-wider.R
# TO DO: Review with Connor when he returns from vacation.
#' Split a string into rows
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Each of these functions takes a string and splits it into multiple rows:
#'
#' * `separate_longer_delim()` splits by a delimiter.
#' * `separate_longer_position()` splits by a fixed width.
#'
#' @export
#' @param delim For `separate_longer_delim()`, a string giving the delimiter
#'   between values. By default, it is interpreted as a fixed string; use
#'   [stringr::regex()] and friends to split in other ways.
#' @inheritParams separate_wider_delim
#' @return A data frame based on `data`. It has the same columns, but different
#'   rows.
#' @examples
#' df <- tibble(id = 1:4, x = c("x", "x y", "x y z", NA))
#' df %>% separate_longer_delim(x, delim = " ")
#'
#' # You can separate multiple columns at once if they have the same structure
#' df <- tibble(id = 1:3, x = c("x", "x y", "x y z"), y = c("a", "a b", "a b c"))
#' df %>% separate_longer_delim(c(x, y), delim = " ")
#'
#' # Or instead split by a fixed length
#' df <- tibble(id = 1:3, x = c("ab", "def", ""))
#' df %>% separate_longer_position(x, 1)
#' df %>% separate_longer_position(x, 2)
#' df %>% separate_longer_position(x, 2, keep_empty = TRUE)
separate_longer_delim_r <- function(data, cols, delim, ...) {
#  check_installed("stringr")
#  check_data_frame(data)
#  check_required(cols)
#  check_string(delim)
#  check_dots_empty()
  
#  if (is_bare_string(delim)) {
#    delim <- stringr::fixed(delim)
#  }
  
  map_unchop(data, {{ cols }}, stringr::str_split, pattern = delim)
}


# helpers -----------------------------------------------------------------

map_unchop <- function(data, cols, fun, ...) {
  cols <- tidyselect::eval_select(
    enquo(cols),
    data = data,
    allow_rename = FALSE,
    allow_empty = FALSE,
#    error_call = .error_call
  )
  col_names <- names(cols)
  
  for (col in col_names) {
    data[[col]] <- fun(data[[col]], ...)
  }
  
  unchop(
    data = data,
    cols = all_of(col_names),
#    keep_empty = .keep_empty,
    ptype = character(),
#    error_call = .error_call
  )
}

# Registry date columns are stored as string and support partial dates. This function puts partial dates in 
# registry date columns into their own column (date_partial). The output of this function can then be piped 
# into impute_dates for processing. 
create_partial_date_registry <- function(.data, date_col) {
  partial_col <- rlang::sym(str_c(deparse(substitute(date_col)), "_partial"))
  date_col <- rlang::enquo(date_col)
  .data %>% 
    check_for_column(!!date_col, .class = "character", return_df = FALSE)
  .data %>% 
    mutate(!!partial_col := case_when(nchar(!!date_col) == 6 ~ paste0(str_extract(!!date_col, "^\\d{4}"), "-", str_extract(!!date_col, "\\d{2}$"), "-UN"),
                                      nchar(!!date_col) == 4 ~ paste0(!!date_col, "-UN-UN"),
                                      TRUE ~ NA_character_))
  
}


## Bucket stop reasons
## categorizes stop reasons into roll ups
## adds sorting to allow for easy prioritzation
bucket_therapy_stop_reasons <- function(.data, stop_reason, therapy_ongoing) {
  stop_reason <- rlang::enquo(stop_reason)
  therapy_ongoing <- rlang::enquo(therapy_ongoing)
  .data %>% 
    check_for_column(!!stop_reason, .class = "character") %>% 
    check_for_column(!!therapy_ongoing, .class = "character", return_df = FALSE)
  .data %>% 
    mutate(stop_reason_category = case_when(!!stop_reason %in% c("Undergoing definitive surgery/radiation, without evidence of progression/worsening cancer",
                                                               "End of planned therapy, without evidence of progression/worsening cancer",
                                                               "Treatment gap of > 365 days where there is no other reason for discontinuation available"
    ) ~ "End of planned therapy",
    ## other 
    !!stop_reason %in% c("Biomarker testing found new actionable mutation,without evidence of progression/worsening cancer",
                       "Treatment for other diseases including hospitalization, without evidence of progression/worsening cancer",
                       "Insurance/cost any changes in, without evidence of progression/worsening cancer",
                       "Changes in insurance, without evidence of progression/worsening cancer",
                       "Patient’s choice not captured in other options, without evidence of progression/worsening cancer",
                       "Physician’s choice not captured in other options, without evidence of progression/worsening cancer"
    ) ~ "Other",
    ## lost to follow up
    !!stop_reason %in% c("Transferred care outside of health system, without evidence of progression/worsening cancer",
                       "Hospice referral, without evidence of progression/worsening cancer",
                       "Patient lost to follow up, without evidence of progression/worsening cancer"
    ) ~ "Lost to follow up",
    ## ongoing
    !!therapy_ongoing == "Yes" ~ "Treatment ongoing",
    is.na(!!stop_reason) ~ "Unknown",
    TRUE ~ !!stop_reason)
    ) %>% 
    add_sorting(stop_reason_category, by = "custom", sort_vector = c("Death",
                                                                     "Progression/worsening cancer",
                                                                     "Intolerance/Toxicity, without evidence of progression/worsening cancer",
                                                                     "Other",
                                                                     "Treatment ongoing",
                                                                     "Lost to follow up",
                                                                     "End of planned therapy",
                                                                     "Unknown"))
}

# Function checks input data frame (typically the argument patients in a function) for patientid and participantid
# If neither are present, errors. If oc_warning is set to TRUE, will warn if participantid is not present.
check_patient_arg <- function(.data, oc_warning = FALSE) {
  if (!is.null(.data)) {
    .data_name <- deparse(substitute(.data))
    if (!is.data.frame(.data)) rlang::abort(glue::glue("{.data_name} is not a data frame"))
    patient_cols <- intersect(names(.data), c("patientid", "participantid"))
    if (length(patient_cols) == 0 ) rlang::abort(glue::glue("patientid and participantid are missing from {.data_name} - at least one of them must be present"))
    if (!"participantid" %in% patient_cols & oc_warning) rlang::warn(glue::glue("participantid is not in {.data_name} - it is best practice to use a prioritized participantid when pulling from OC due to patient duplication!\nUsing patientid only"))
    patient_cols
  } else {
    NULL
  }
}
