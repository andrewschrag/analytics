# Staging functions

# Generics (pan-ADS) ----
# Maps registry TNM codes to standard values
# Could eventually write a wrapper for this like syhelpr::impute_dates
map_tnm <- function(.data, column, mapping, spmd = spmd_con()) {
  # # We could do the below to allow for tbl_lazy input, however in the pipeline there are a lot of joins 
  # # already and I think it is faster to do this in R (i.e., not a tbl_lazy)
  # if ("tbl_lazy" %in% class(.data)) {
  #   mappingdf <- tbl(.data$src$con, in_schema("maps", "concept")) %>% 
  #     filter(vocabulary == "NAACCR" & domain == "SEER_DI" & class == mapping) %>% 
  #     select(code, name)
  # } else if (is.data.frame(.data)) {
  #   mappingdf <- tbl(spmd, in_schema("maps", "concept")) %>% 
  #     filter(vocabulary == "NAACCR" & domain == "SEER_DI" & class == mapping) %>% 
  #     select(code, name) %>% 
  #     collect()
  # } else {
  #   rlang::abort(".data is not tbl_lazy, data.frame, or tibble")
  # }
  mappingdf <- tbl(spmd, in_schema("maps", "concept")) %>%
    filter(vocabulary == "NAACCR" & domain == "SEER_DI" & class == mapping) %>%
    select(code, name) %>%
    collect()
  column <- enquo(column)
  .data %>% 
    mutate(code = !!column) %>% 
    left_join(mappingdf, by = "code") %>% 
    mutate(!!column := coalesce(name, !!column)) %>% 
    # "88" = "Not applicable, no code assigned for this case in the current AJCC Staging Manual."
    # if_else necessary to preserve data type as character if the column is all NA
    mutate(!!column := if_else(!!column == "88", NA_character_, !!column)) %>% 
    select(-name, -code)
}

# Maps FIGO stage returned from registry data. This field holds a mix of the old codes, new codes, and an 
# edge case we need to handle. 
# stage should be a string, still need to figure out how to implement using non-standard evaluation
# May be ovarian specific - would any others use this?
map_figo_stage <- function(.data, stage, spmd = spmd_con()) {
  map_figo_old <- tbl(spmd, in_schema("maps", "concept")) %>%
    filter(vocabulary == "NAACCR" & domain == "SEER_DI" & class == "FIGOStage") %>% 
    select(code, name) %>% 
    collect()
  .data %>% 
    # Apply the "new" mapping
    coalesce_map_naaccr(stage, spmd) %>%
    # Apply the "old" mapping
    left_join(map_figo_old, by = setNames("code", stage)) %>% 
    mutate(!!rlang::sym(stage) := coalesce(name, !!rlang::sym(stage))) %>% 
    select(-name) %>% 
    # Fix edge case: somehow the NAACCR API returned the wrong code for this mapping (returned code = "IC2" 
    # instead of "1C2")
    # if_else necessary to preserve data type as character if the column is all NA
    mutate(!!rlang::sym(stage) := if_else(!!rlang::sym(stage) == "1C2", 
                                          "FIGO Stage IC2", 
                                          !!rlang::sym(stage)))
}

# Prioritizes stage by taking non-missing staging, then CTR corrected > pathological/FIGO > clinical, then 
# highest available. Requires `source` column in `.data`.
prioritize_stage <- function(.data, stage_column = stagegroup) {
  stage_column <- enquo(stage_column)
  .data %>% 
    # gives us an more easily sortable stage
    mutate(!!stage_column := str_trim(str_remove_all(!!stage_column, "FIGO|Stage"))) %>% 
    mutate(stage = bucket_stage(!!stage_column, subdivide = FALSE),
           # CTR correct pre-populated data, so we prioritize this field
           source_rank = case_when(source == "CTR" ~ 1, 
                                   source %in% c("pathological", "FIGO") ~ 2, # we treat these as equivalent
                                   source == "clinical" ~ 3)) %>% 
    group_by(patientid) %>% 
    # prioritize non-missing staging, then CTR corrected > pathological/FIGO > clinical, then highest
    # bucket_stage (above) returns NA for "Not found", "Unable to stage", "Unknown", "Not applicable", etc.
    arrange(is.na(stage), 
            source_rank,
            desc(stage), 
            !!stage_column) %>% 
    slice_head() %>% 
    ungroup()
}

#' Bundle TNM
#'
#' This function rolls up T, N, and M fields into succinct buckets.
#'
#' @param x A character vector
#' 
#' @export
bucket_tnm <- function(.data, .vars = c(...)) {
  .data %>%
    mutate(across(all_of(.vars),
                  list(detailed = ~ ., # copy current, detailed into column with _detailed name
                       # some special cases, unsure if breast specific
                       clean = ~case_when(. == "Tis" ~ ., 
                                          . == "N1mi" ~ .,
                                          grepl("X", .) ~ .,
                                          TRUE ~ stringr::str_extract(pattern = '^.\\d',# removes everything after the stage digit 
                                                                      string = .)
                       )
                  ),
                  .names = "{.col}_{.fn}")
    ) %>% 
    ## now remove orig. cols so cleaned version can replace
    select(-all_of(.vars)) %>% 
    rename_with(.cols = ends_with("_clean"), 
                ~stringr::str_remove(.,"_clean"))
}

# Bladder ----

# Function for prioritizing stage for bladder. If can take a grouped df as input. Please see bladder.R for 
# an example of using. 
# Prioritization applied (by grouping):
# 1) Non-missing stage, T, N, or M over records missing these
# 2) clinical > pathological
# 3) Highest stage
# Non-missing stage is defined as at least one of stage group, T, N, or M is present.
# Highest stage is defined as
# 1) highest non-missing stage group first (IV > III > II > I > 0is > 0a)
# 2) then highest non-missing M (M1 > M0)
# 3) then highest non-missing T (T4b > T4a > ... > T1 > Tis > Ta)
# Notes:
# - TX, NX, and MX are counted as missing and are therefore the same as NA, however we want to keep these 
#   values if that is all the patient has so please create NA versions of TX, NX, and MX to make sorting 
#   easier
prioritize_stage_bladder <- function(.data, 
                                     stagegroup = stagegroup, 
                                     stagegroup_raw = stagegroup_raw, 
                                     t = t_clean, 
                                     n = n_clean, 
                                     m = m_clean, 
                                     source = source) {
  if (!is_grouped_df(.data)) {
    rlang::warn("Stage data frame supplied ungrouped to prioritize_stage_bladder! Stage will be prioritized by patient! 
                If performing multiple tumor prioritization, please group by patientid and tumorid.")
    .data <- .data %>% 
      group_by(patientid)
  }
  stagegroup <- enquo(stagegroup)
  stagegroup_raw <- enquo(stagegroup_raw)
  t <- enquo(t)
  n <- enquo(n)
  m <- enquo(m)
  source <- enquo(source)
  # adding source and T prioritization
  .data <- .data %>% 
    add_sorting(!!source, 
                sort_vector = c("clinical", "pathological")) %>% 
    add_sorting_custom_handle_missing(!!t, 
                                      c("T4b",
                                        "T4a",
                                        "T4",
                                        "T3b",
                                        "T3a",
                                        "T3",
                                        "T2b",
                                        "T2a",
                                        "T2",
                                        "T1",
                                        "Tis",
                                        "Ta",
                                        "T0"))
  .data %>% 
    # at least one or more non-missing: if one or more is present, this becomes FALSE and FALSE comes before TRUE
    arrange(is.na(!!stagegroup) & is.na(!!t) & is.na(!!n) & is.na(!!m), 
            # clinical > pathological
            !!source,
            # highest stage (IVC to 0a)
            desc(!!stagegroup),
            # highest M (assuming stage is missing)
            !grepl("1", !!m), 
            # highest T (assuming stage is missing)
            !!t, 
            # this last one breaks any ties when stage is one of the flavors of unknown by using the raw 
            # version
            desc(!!stagegroup_raw)
    ) %>% 
    slice_head() %>% 
    ungroup() %>% 
    mutate(!!source := as.character(!!source))
}

# Function to take cleaned (normalized, imputed, prioritized) staging data for bladder patients and configure
# it into the structure and value set expected in the ADS. Can produce the overall, clinical, and pathological
# versions of the variables (set the `type` parameter).
wrangle_staging_bladder <- function(.data, type = NULL) {
  
  if (!is.null(type)) {
    if (!type %in% c("pathological", "clinical")) {
      rlang::abort(glue::glue("type = '{type}' is not a valid option. Use one of 'pathological' or 'clinical'"))
    }
  }
  
  .data <- .data %>% 
    # replace the artifical NAs in stagegroup
    mutate(stage_group = coalesce(stagegroup, stagegroup_raw)) %>% 
    # Added additional mapping
    mutate(stage_group = case_when(stage_group == "IVC" ~ "IV",
                                   stage_group == "IIIC" ~ "III",
                                   stage_group %in% c("IIB", "IIC") ~ "II", 
                                   stage_group == "IC" ~ "I", 
                                   stage_group == "Not applicable - Registry" ~ "Unknown",
                                   stage_group == "Unknown, not staged" ~ "Not staged",
                                   stage_group == "Unknown, unstaged" ~ "Not staged",
                                   TRUE ~ stage_group)) %>% 
    replace_na(list(stage_group = "Unknown")) %>% 
    # not using t_clean etc. because we want TX etc.
    rename_with(~ str_c("stage_", .), c(t, n, m, mismatched_staging_flag, mismatched_staging_highest)) %>% 
    rename(stage_group_type = source) %>% 
    select(patientid, 
           stage_t, 
           stage_n, 
           stage_m, 
           stage_group, 
           stage_group_type, 
           stage_imputation_flag, 
           stage_mismatched_staging_flag, 
           stage_mismatched_staging_highest)
  
  filter_rename_type <- function(.data, type) {
    .data %>% 
      # if patient is missing a given staging type, we want to leave the unknowns, but prefer known over unknown (NA)
      filter(stage_group_type == type | is.na(stage_group_type)) %>% 
      group_by(patientid) %>% 
      arrange(stage_group_type) %>% 
      slice_head() %>% 
      ungroup() %>% 
      # currently not reporting
      select(-c(stage_group_type, stage_imputation_flag, stage_mismatched_staging_flag, stage_mismatched_staging_highest)) %>% 
      rename_with(~ str_replace(., "^stage_", glue::glue("{type}_stage_"))) %>% 
      rename_with(~ str_c(., "_dx"), .cols = !patientid) 
      
  }
  
  if (is.null(type)) {
    .data %>% 
      rename_with(~ str_replace(., "^stage_", "prioritized_stage_")) %>% 
      rename_with(~ str_replace(., "stage_group", "stage_group_dx")) %>%
      rename_with(~ str_replace(., "stage", "stage_dx"), .cols = matches("imputation|mismatched")) %>%
      rename_with(~ str_c(., "_dx"), .cols = matches("_(t|n|m)$")) %>% 
      rename(prioritized_stage_group_detailed_dx = prioritized_stage_group_dx)
  } else if (type == "pathological") {
    .data %>% 
      filter_rename_type("pathological")
  } else if (type == "clinical") {
    .data %>% 
      filter_rename_type("clinical")
  }
}

# Breast ----

# Function to prioritize breast cancer stages. Can take a grouped df as input. 
# Prioritization applied (by grouping):
# 1) Non-missing stage, T, N, or M over records missing these
# 2) CTR corrected > pathological > clinical
# 3) Highest stage
# Non-missing stage is defined as at least one of stage group, T, N, or M is present.
# Highest stage is defined as
# 1) highest non-missing stage group first (IV > III > II > I > 0is > 0a)
# 2) then highest non-missing M (M1 > M0(i+) > M0)
# 3) then highest non-missing T (T4b > T4a > ... > T1 > Tis > Ta > T0 > TX)
# 4) then highest non-missing N (N4b > N4a > ... > N1 > N0 > NX)
# Notes:
# - [not implemented] If TX or T0 or Tis we might wish to flag the patient to check manually. 
prioritize_stage_breast <- function(.data, 
                                    stagegroup = stagegroup, 
                                    t = t, 
                                    n = n, 
                                    m = m, 
                                    type = type) {
  if (!is_grouped_df(.data)) {
    rlang::warn("Stage data frame supplied ungrouped to prioritize_stage_breast! Stage will be prioritized by patient!")
    .data <- .data %>% 
      group_by(patientid)
  }
  
  stagegroup <- enquo(stagegroup)
  t <- enquo(t)
  n <- enquo(n)
  m <- enquo(m)
  type <- enquo(type)
  # adding prioritizations
  .data <- .data %>% 
    add_sorting(!!type, 
                sort_vector = c("MA",  
                                "AJCC pathological",
                                "AJCC clinical")) %>% 
    add_sorting_custom_handle_missing(!!t, 
                                      c("T4d", "T4c", "T4b", "T4a", "T4", 
                                        "T3c", "T3b", "T3a", "T3mi", "T3",
                                        "T2c", "T2b", "T2a", "T2mi", "T2",
                                        "T1c", "T1b", "T1a", "T1mi", "T1",
                                        "Tis", "Ta", "T0", "TX", "Unknown, not staged")) %>% 
    add_sorting_custom_handle_missing(!!n, 
                                      c("N3c", "N3b", "N3a", "N3", 
                                        "N2b", "N2a", "N2", 
                                        "N1c", "N1b", "N1a", "N1mi", "N1", 
                                        "N0(mol+)", "N0(i+)", "N0(mol-)", "N0(i-)", "N0", 
                                        "NX", "Unknown, not staged")) %>% 
    add_sorting_custom_handle_missing(!!m, 
                                      c("M1", "M0(i+)", "M0", 
                                        "MX", "Unknown, not staged"))
  
  .data %>% 
    # gives us an more easily sortable stage
    mutate(stage = bucket_stage(!!stagegroup, subdivide = FALSE)) %>% 
    # bucket_stage (above) returns NA for "Not found", "Unable to stage", "Unknown", "Not applicable", etc.
    # at least one or more non-missing: if one or more is present, this becomes FALSE and FALSE comes before TRUE
    arrange(is.na(stage), # & is.na(!!t) & is.na(!!n) & is.na(!!m),  # removing in the absence of mismatch control
            # MA > pathological > clinical
            desc(!!type),
            # highest stage (IVC to 0a)
            desc(stage), # arranges the unknown flavors to last
            desc(!!stagegroup), # ensures IIIc > IIIb, etc
            # highest M 
            desc(!!m), 
            # highest T 
            desc(!!t), 
            # highest N 
            desc(!!n)
    ) %>% 
    slice_head() %>% 
    ungroup() %>% 
    mutate(stage = coalesce(stage, !!stagegroup), ## pull more detailed unknown back to stage
           !!type := as.character(!!type)) %>% 
    rename(stage_detailed = !!stagegroup)
}

