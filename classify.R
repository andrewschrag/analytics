# ---------- helpers ----------

.site_major <- function(x) {
  x %>%
    str_replace("\\s*\\(.*\\)$", "") %>%      # drop parenthetical
    str_replace("\\s*\\bNOS\\b.*$", "") %>%   # drop NOS and anything after
    str_replace(",\\s*.*$", "") %>%           # drop trailing ", ..."
    str_trim() %>%
    str_to_lower()
}

# Convert various time types to numeric days
.to_days_numeric <- function(x) {
  if (is.null(x)) return(Inf)
  if (is.numeric(x)) return(as.numeric(x))
  if (inherits(x, "Duration")) return(as.numeric(x, units = "days"))
  if (inherits(x, "difftime")) return(as.numeric(x, units = "days"))
  if (inherits(x, "Period")) return(lubridate::time_length(x, "days"))
  Inf
}

# ---------- main function ----------

#' Classify primary cancers vs recurrences and link concurrent primaries
#'
#' @param df tibble with one row per diagnosis episode
#' @param patient_id_col patient id column (sym), default npm_patient_id
#' @param dx_date_col diagnosis date (Date or parseable), default diagnosis_date
#' @param site_desc_col site description (e.g., icdo3_site_desc)
#' @param laterality_col laterality ("Left","Right", etc.), can be NA
#' @param histology_col optional; set split_by_histology=TRUE to split primaries by histology
#' @param paired_sites character vector of paired organs (lowercase)
#' @param concurrent_window_days diagnoses within this gap (days) form same concurrent group
#' @param force_new_primary_after_days numeric days or lubridate period/duration/difftime
#' @param split_by_histology logical; default FALSE
#' @param quiet logical; print summary if FALSE
#'
#' @return tibble with added columns:
#'   site_major, dx_date, primary_flag, recurrence_flag,
#'   primary_id, primary_sequence, recurrence_number, concurrent_primary_id
#'
classify_primary <- function(
    df,
    patient_id_col = npm_patient_id,
    dx_date_col = diagnosis_date,
    site_desc_col = site_group,
    laterality_col = ca_laterality,
    histology_col = NULL,
    paired_sites = c("lung","breast","kidney","ovary","testis","eye","parotid","submandibular"),
    concurrent_window_days = 60L,
    force_new_primary_after_days = lubridate::years(5),
    split_by_histology = FALSE,
    quiet = TRUE
) {
  # NSE capture
  patient_id_col <- enquo(patient_id_col)
  dx_date_col    <- enquo(dx_date_col)
  site_desc_col  <- enquo(site_desc_col)
  laterality_col <- enquo(laterality_col)
  histology_col  <- enquo(histology_col)

  # numeric threshold days (robust to Period/Duration/difftime)
  threshold_days <- .to_days_numeric(force_new_primary_after_days)
  if (!is.finite(threshold_days) || is.na(threshold_days)) threshold_days <- Inf

  has_site_group <- "site_group" %in% names(df)

  # Base normalization
  x <- df %>%
    mutate(
      .patient_id   = as.character(!!patient_id_col),
      .dx_date_raw  = !!dx_date_col,
      .dx_date      = lubridate::as_date(.dx_date_raw),
      .site_desc    = as.character(!!site_desc_col),
      .site_major   = .site_major(.site_desc),
      .site_desc_lc = str_to_lower(.site_desc),
      .laterality_raw = as.character(!!laterality_col)
    )

  if (!quo_is_null(histology_col)) {
    x <- x %>% mutate(.histology = as.character(!!histology_col))
  } else {
    x <- x %>% mutate(.histology = NA_character_)
  }

  # Laterality normalization (handles " R", "Rt", "L", etc.)
  x <- x %>%
    mutate(
      .lat_lc = .laterality_raw %>% str_to_lower() %>% str_trim(),
      .laterality_norm = case_when(
        is.na(.lat_lc) | .lat_lc == "" ~ "unk",
        .lat_lc %in% c("left","l","lt","lft") | str_detect(.lat_lc, "^l\\b") ~ "left",
        .lat_lc %in% c("right","r","rt","rgt") | str_detect(.lat_lc, "^r\\b") ~ "right",
        str_detect(.lat_lc, "bilat|both|two") ~ "bilateral",
        TRUE ~ "unk"
      )
    )

  # Derive organ for primary/paired logic
  if (has_site_group) {
    x <- x %>%
      mutate(.organ = .data$site_group %>% tolower() %>% str_trim())
  } else {
    x <- x %>%
      mutate(
        .organ = case_when(
          str_detect(.site_desc_lc, "\\blung\\b") ~ "lung",
          str_detect(.site_desc_lc, "\\bbreast\\b") ~ "breast",
          str_detect(.site_desc_lc, "\\bkidney\\b") ~ "kidney",
          str_detect(.site_desc_lc, "\\bovary|ovarian\\b") ~ "ovary",
          str_detect(.site_desc_lc, "\\btestis|testicle\\b") ~ "testis",
          str_detect(.site_desc_lc, "\\beye\\b") ~ "eye",
          str_detect(.site_desc_lc, "\\bparotid\\b") ~ "parotid",
          str_detect(.site_desc_lc, "\\bsubmandibular\\b") ~ "submandibular",
          # if pattern like "Upper lobe, lung", capture organ after comma
          str_detect(.site_desc_lc, ",\\s*(lung|breast|kidney|ovary|testis|eye|parotid|submandibular)\\b") ~
            str_match(.site_desc_lc, ",\\s*(lung|breast|kidney|ovary|testis|eye|parotid|submandibular)\\b")[,2],
          TRUE ~ .site_major
        )
      )
  }

  paired_sites <- str_to_lower(paired_sites)

  x <- x %>%
    mutate(
      .is_paired  = .organ %in% paired_sites,
      .paired_key = if_else(.is_paired, .laterality_norm, "nonpaired")
    )

  # Primary key (optionally include histology)
  x <- x %>%
    mutate(
      .primary_key = if (isTRUE(split_by_histology)) {
        paste(.organ, .paired_key, coalesce(.histology, "na"), sep = "::")
      } else {
        paste(.organ, .paired_key, sep = "::")
      }
    ) %>%
    arrange(.patient_id, .dx_date, .organ)

  # Split into blocks if gap > threshold_days within same primary_key
  x <- x %>%
    group_by(.patient_id, .primary_key) %>%
    mutate(
      .gap_days = as.integer(.dx_date - lag(.dx_date)),
      .force_split = !is.na(.gap_days) & (.gap_days > threshold_days),
      .primary_block = cumsum(coalesce(.force_split, FALSE)) + 1L
    ) %>%
    ungroup()

  # Mark index primaries table (the FIRST row by date in each block)
  primaries <- x %>%
    group_by(.patient_id, .primary_key, .primary_block) %>%
    slice_min(.dx_date, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      .patient_id, .primary_key, .primary_block,
      .primary_index_date = .dx_date,
      .is_primary = TRUE
    )

  # Join primary flags back
  x <- x %>%
    left_join(primaries, by = c(".patient_id", ".primary_key", ".primary_block")) %>%
    mutate(
      primary_flag = coalesce(.is_primary, FALSE),
      recurrence_flag = !primary_flag
    )

  # Assign primary_id & sequence per patient using primaries ONLY
  primaries_with_ids <- primaries %>%
    distinct(.patient_id, .primary_key, .primary_block, .primary_index_date) %>%
    arrange(.patient_id, .primary_index_date) %>%
    group_by(.patient_id) %>%
    mutate(
      primary_sequence = dplyr::row_number(),
      primary_id = sprintf("%s-P%03d", .patient_id, primary_sequence)
    ) %>%
    ungroup()

  # Bring primary_id to all rows
  x <- x %>%
    left_join(primaries_with_ids,
              by = c(".patient_id", ".primary_key", ".primary_block"))

  # Number recurrences within each primary_id
  x <- x %>%
    group_by(.patient_id, primary_id) %>%
    arrange(.dx_date, .by_group = TRUE) %>%
    mutate(
      recurrence_number = if_else(recurrence_flag,
                                  cumsum(recurrence_flag),
                                  NA_integer_)
    ) %>%
    ungroup()

  # Compute concurrent groups *from primaries_with_ids* (has .primary_index_date!)
  prim_concurrent <- primaries_with_ids %>%
    arrange(.patient_id, .primary_index_date) %>%
    group_by(.patient_id) %>%
    mutate(
      .delta = as.integer(.primary_index_date - lag(.primary_index_date)),
      .new_grp = if_else(is.na(.delta) | (.delta > concurrent_window_days), 1L, 0L),
      .grp_index = cumsum(.new_grp),
      concurrent_primary_id = sprintf("%s-C%02d", .patient_id, .grp_index)
    ) %>%
    ungroup() %>%
    select(.patient_id, primary_id, concurrent_primary_id)

  # Join concurrent group id to all rows by primary_id
  x <- x %>%
    left_join(prim_concurrent, by = c(".patient_id", "primary_id"))

  # Final shape: preserve user columns first, then outputs
  out <- x %>%
    arrange(.patient_id, .dx_date, dplyr::desc(primary_flag)) %>%
    mutate(
      site_major = .organ,
      dx_date = .dx_date
    ) %>%
    select(
      any_of(names(df)),
      site_major, dx_date,
      primary_flag, recurrence_flag,
      primary_id, primary_sequence, recurrence_number,
      concurrent_primary_id
    )

  # ---- summary without referencing internal .patient_id column ----
  if (!quiet) {
    pid_name <- rlang::as_name(enquo(patient_id_col))  # original column name
    patients <- dplyr::n_distinct(out[[pid_name]])
    primaries_n <- sum(out$primary_flag, na.rm = TRUE)
    recurrences_n <- sum(out$recurrence_flag, na.rm = TRUE)
    message(sprintf("Summary: patients=%d, primaries=%d, recurrences=%d",
                    patients, primaries_n, recurrences_n))
  }

  out
}




#' Classify Lung Cancer Histology as Squamous / Non-Squamous / Other
#'
#' @description
#' Adds a `squamous_group` column that classifies each row into
#' **"Squamous"**, **"Non-Squamous"**, or **"Other"** using available
#' histology signals from free text and/or ICD-O-3 morphology codes.
#' The logic targets the usual NSCLC split (squamous vs non-squamous),
#' treating small cell as **"Other"** and including several rare/sarcomatoid
#' variants in **"Non-Squamous"**.
#'
#' @details
#' **Signals used (when present):**
#' - **ICD-O-3 morphology**:  
#'   - Squamous: `8050–8084`, `8123`  
#'   - Adenocarcinoma (adeno-like): `8140–8389`  
#'   - Large cell: `8012–8014`  
#'   - NSCLC-NOS: `8046`  
#'   - Adenosquamous: `8560`  
#'   - Small cell: `8041–8045` (classified as **"Other"**)  
#'   - LCNEC: `8013` (classified as **"Non-Squamous"** by default)  
#'   - Extra non-squamous variants: `8550`, `8551`, `8490`, `8576`,
#'     `8030–8034`, `8022`, `8004`, `8980`, `8972`, `8023`
#'
#' - **Histology text regex** (case-insensitive):
#'   - "squamous" → **"Squamous"**
#'   - "non\\s*small\\s*cell" or "\\bnsclc\\b" → **"Non-Squamous"**
#'   - "adeno|adenoca|adenocarcinoma" → **"Non-Squamous"**
#'   - "large\\s*cell" or "large\\s*cell\\s*neuroendocrine|\\blcnec\\b" → **"Non-Squamous"**
#'   - "adenosquamous" → **"Non-Squamous"**
#'   - "small\\s*cell" → **"Other"**
#'   - Rare/variants mapped to **"Non-Squamous"**: "signet\\s*ring",
#'     "sarcomatoid|carcinosarcoma|spindle|pleomorphic|giant\\s*cell|pseudosarcom|polygonal",
#'     "blastoma", "\\bnut[-\\s]*associated|nuclear\\s*protein\\s*in\\s*testis"
#'
#' **Lung context heuristic:**  
#' Rows are considered in lung context if any lung-ish signal is present
#' (e.g., NSCLC, squamous, adeno, large cell, LCNEC, small cell, the extra
#' variant set above) or the cancer type/histology text explicitly contains
#' "lung". Non-lung rows default to **"Other"** when clearly non-lung;
#' otherwise remain `NA`.
#'
#' @section Output:
#' Returns the input `df` with an added character column `squamous_group`
#' taking values in `c("Squamous", "Non-Squamous", "Other", NA)`.
#'
#' @param df A data frame containing histology signals.
#' @param cancer_type_col Optional. Name of a column with cancer/site type
#'   information (e.g., `"ca_type"`, `"cancer_type"`). Used only for context
#'   (detecting "lung" or "small cell"/"NSCLC" phrases).
#' @param histology_text_col Optional. Name of a column containing free-text
#'   histology labels (e.g., `"histology"`, `"histology_label"`).
#' @param icdo3_morph_col Optional. Name of a column containing ICD-O-3
#'   morphology codes (e.g., `"icdo3_morph"`, `"icdo3_histology"`). Codes may
#'   include behavior suffixes like `"/3"`; these are normalized internally.
#'
#' @return A tibble/data frame with a new column `squamous_group`.
#'
#' @examples
#' \dontrun{
#' # Minimal example with ICD-O-3 and text
#' library(dplyr)
#' df <- tibble::tibble(
#'   icdo3_morph = c("8070/3","8140/3","8041/3","8013/3","8032/3","8550/3","8046/3"),
#'   histology   = c("Squamous cell carcinoma, NOS",
#'                   "Adenocarcinoma, NOS",
#'                   "Small cell carcinoma",
#'                   "Large cell neuroendocrine carcinoma",
#'                   "Spindle cell carcinoma, NOS",
#'                   "Acinar cell carcinoma",
#'                   "Non-small cell carcinoma")
#' )
#' classify_squamous_group(df) %>% dplyr::count(squamous_group)
#'
#' # With your dataset:
#' # df2 <- read.delim("temp (17).tsv", sep = "\t")
#' # out <- classify_squamous_group(df2)  # auto-detects icdo3_morph + histology
#' # dplyr::count(out, squamous_group)
#' }
#'
#' @seealso
#' - ICD-O-3 coding guides for morphology categories
#'
#' @export

classify_squamous_group <- function(
    df,
    cancer_type_col = NULL,      # optional: e.g., "ca_type"
    histology_text_col = NULL,   # optional: e.g., "histology", "histology_label"
    icdo3_morph_col = NULL       # optional: e.g., "icdo3_morph", "icdo3_histology"
) {
  # --- Auto-detect columns if not provided ---
  if (is.null(histology_text_col)) {
    histology_text_col <- intersect(c("histology","histology_label","histology_text"), names(df)) %>% head(1)
    if (length(histology_text_col) == 0) histology_text_col <- NULL
  }
  if (is.null(icdo3_morph_col)) {
    icdo3_morph_col <- intersect(c("icdo3_morph","icdo3_histology","morphology","icdo_morph"), names(df)) %>% head(1)
    if (length(icdo3_morph_col) == 0) icdo3_morph_col <- NULL
  }
  if (is.null(cancer_type_col)) {
    cancer_type_col <- intersect(c("ca_type","cancer_type","primary_site","site_type"), names(df)) %>% head(1)
    if (length(cancer_type_col) == 0) cancer_type_col <- NULL
  }

  pull_or_na <- function(dat, col) if (!is.null(col) && col %in% names(dat)) dat[[col]] else NA

  ct <- pull_or_na(df, cancer_type_col)
  hx <- pull_or_na(df, histology_text_col)
  ic <- pull_or_na(df, icdo3_morph_col)

  ct_str <- tolower(as.character(ct))
  hx_str <- tolower(as.character(hx))

  # Normalize ICD-O-3 morphology to 4-digit numeric string (strip suffix like "/3")
  ic_str <- if (all(is.na(ic))) {
    rep(NA_character_, nrow(df))
  } else {
    x <- as.character(ic)
    x <- gsub("[^0-9]", "", x)
    x[nchar(x) >= 4] <- substr(x[nchar(x) >= 4], 1, 4)
    x[nchar(x) > 0 & nchar(x) < 4] <- sprintf("%04d", suppressWarnings(as.integer(x[nchar(x) > 0 & nchar(x) < 4])))
    x[nchar(x) == 0] <- NA_character_
    x
  }

  in_code_range <- function(code_chr, lo, hi) {
    suppressWarnings({
      num <- as.integer(code_chr)
      !is.na(num) & num >= lo & num <= hi
    })
  }

  # --- Core ICD-O groupings ---
  ic_squamous       <- in_code_range(ic_str, 8050, 8084) | ic_str == "8123"
  ic_adeno          <- in_code_range(ic_str, 8140, 8389)
  ic_large_cell     <- in_code_range(ic_str, 8012, 8014)
  ic_adenosquamous  <- ic_str == "8560"
  ic_nsclc_nos      <- ic_str == "8046"
  ic_small_cell     <- in_code_range(ic_str, 8041, 8045)
  ic_lcnec          <- ic_str == "8013"   # Large cell neuroendocrine carcinoma (treat as Non-Squamous)

  # --- Explicit Non-Squamous extras (rare/sarcomatoid/variants outside simple ranges) ---
  ic_nonsq_extras <- ic_str %in% c(
    "8550","8551",        # Acinar (adeno variant)
    "8490",               # Signet ring cell carcinoma (adeno variant)
    "8576",               # Hepatoid adenocarcinoma
    "8030","8031","8032","8033","8034",  # Giant/spindle/pseudosarcomatoid/polygonal
    "8022","8004",        # Pleomorphic, malignant spindle cell
    "8980","8972",        # Carcinosarcoma, Pulmonary blastoma
    "8023"                # NUT carcinoma
  )

  # --- Text patterns ---
  pat_squamous       <- "squamous"
  pat_small_cell     <- "small\\s*cell"
  pat_nsclc_text     <- "non\\s*small\\s*cell|\\bnsclc\\b"
  pat_adeno_like     <- "adeno|adenoca|adenocarcinoma"
  pat_large_cell     <- "large\\s*cell"
  pat_adenosquamous  <- "adenosquamous"
  pat_lcnec          <- "large\\s*cell\\s*neuroendocrine|\\blcnec\\b"

  # Rare/variant textual cues
  pat_signet         <- "signet\\s*ring"
  pat_sarcomatoid    <- "sarcomatoid|carcinosarcoma|spindle|pleomorphic|giant\\s*cell|pseudosarcom|polygonal"
  pat_blastoma       <- "blastoma"
  pat_nut            <- "\\bnut[-\\s]*associated|nuclear\\s*protein\\s*in\\s*testis"

  hx_has_squamous       <- !is.na(hx_str) & str_detect(hx_str, pat_squamous)
  hx_has_small_cell     <- !is.na(hx_str) & str_detect(hx_str, pat_small_cell)
  hx_has_nsclc          <- !is.na(hx_str) & str_detect(hx_str, pat_nsclc_text)
  hx_has_adeno_like     <- !is.na(hx_str) & str_detect(hx_str, pat_adeno_like)
  hx_has_large_cell     <- !is.na(hx_str) & str_detect(hx_str, pat_large_cell)
  hx_has_adenosquamous  <- !is.na(hx_str) & str_detect(hx_str, pat_adenosquamous)
  hx_has_lcnec          <- !is.na(hx_str) & str_detect(hx_str, pat_lcnec)

  hx_has_signet         <- !is.na(hx_str) & str_detect(hx_str, pat_signet)
  hx_has_sarcomatoid    <- !is.na(hx_str) & str_detect(hx_str, pat_sarcomatoid)
  hx_has_blastoma       <- !is.na(hx_str) & str_detect(hx_str, pat_blastoma)
  hx_has_nut            <- !is.na(hx_str) & str_detect(hx_str, pat_nut)

  ct_has_small_cell     <- !is.na(ct_str) & str_detect(ct_str, pat_small_cell)
  ct_has_nsclc          <- !is.na(ct_str) & str_detect(ct_str, pat_nsclc_text)

  # --- Lung context heuristic ---
  is_lungish_signal <- (
    ic_nsclc_nos | ic_squamous | ic_adeno | ic_large_cell | ic_adenosquamous | ic_small_cell | ic_lcnec | ic_nonsq_extras |
      hx_has_nsclc | hx_has_squamous | hx_has_adeno_like | hx_has_large_cell | hx_has_adenosquamous | hx_has_lcnec |
      hx_has_signet | hx_has_sarcomatoid | hx_has_blastoma | hx_has_nut |
      ct_has_nsclc | ct_has_small_cell
  )
  is_lung_word <- (!is.na(ct_str) & str_detect(ct_str, "lung")) | (!is.na(hx_str) & str_detect(hx_str, "lung"))
  lung_context <- is_lungish_signal | is_lung_word

  out <- rep(NA_character_, nrow(df))

  # Non-lung entries (if clearly non-lung) -> "Other"
  non_lung_clear <- !lung_context & !is.na(ct_str) & !str_detect(ct_str, "lung")
  out[non_lung_clear] <- "Other"

  # Priority A: Squamous
  squam_signal <- (ic_squamous | hx_has_squamous) & lung_context
  out[squam_signal] <- "Squamous"

  # Priority B: Small cell -> Other (override)
  small_cell_signal <- (ic_small_cell | hx_has_small_cell | ct_has_small_cell) & lung_context
  out[small_cell_signal] <- "Other"

  # Priority C: Non-squamous NSCLC (adeno, large cell, adenosquamous, NSCLC-NOS, LCNEC, extras)
  nonsq_signal <- (
    (ic_adeno | ic_large_cell | ic_adenosquamous | ic_nsclc_nos | ic_lcnec | ic_nonsq_extras |
       hx_has_adeno_like | hx_has_large_cell | hx_has_adenosquamous | hx_has_nsclc |
       hx_has_signet | hx_has_sarcomatoid | hx_has_blastoma | hx_has_nut | ct_has_nsclc) &
      !squam_signal & !small_cell_signal & lung_context
  )
  out[nonsq_signal] <- "Non-Squamous"

  df %>% mutate(squamous_group = out)
}
