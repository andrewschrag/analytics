# Functions ----
#' Calculates the difference in days between two dates
#'
#' Wrapper for using difftime to get difference in days between two dates.
#'
#' @param time1 A date or timestamp object that time2 is being compared to to find the difference in days.
#' @param time2 A date or timestamp object that's difference, in days, from time1 will be calculated.
#'
#' @return A numeric representing the difference in days between time1 and time2. Can be negative.
day_difference <- function(time1,time2) {
  return(as.double(difftime(time2,time1,units="days")))
}

#' Wrangle medication data
#'
#' Wrangles medication data into a shape useful for generating lines of therapy.
#'
#' @param meds A data frame consisting of the base view of a cohort's `medication` table (i.e., that returned from calling `<cohort_connection>$medication`).
#'
#' @return A data frame of a simplified view of the medication data.
wrangle_medications <- function(meds, cols = c()) {
  meds %>%
    rename_with(tolower) %>%
    mutate(drugproduct = tolower(drugproduct),
           startdate = as.Date(pmin(earlieststartdts, startdts, orderdate, na.rm = T)),
           enddate = if_else(source == "ma",
                             as.Date(pmin(latestenddts, enddts, na.rm = T)),
                             as.Date(pmin(latestenddts, enddts, lateststartdts, orderdate, na.rm = T)))) %>%
    mutate(startdate_year = year(startdate),
           enddate_year = year(enddate)) %>%
    distinct(patientId=patientid, drugproduct, startdate, startdate_year, enddate, enddate_year, across(any_of(cols)))
}

#' Categorize medication timing
#'
#' Categorizes medication data representing medication time spans (i.e., start and end of therapies) relative
#' to a given index date. Divides input data into a named list of four data frames representing four
#' possible categories:
#' \itemize{
#'   \item `after`: Medications that unambiguously start after index date, or that unambiguously end after
#'   the index date plus some input cutoff (`min_days_after_index`).
#'   \item `before`: Medications that unambiguously end before index date, or that unambiguously start before
#'   index and unambiguously end before index date plus some input cutoff (`min_days_after_index`).
#'   \item `unknown`: Medications that have entirely unknown dates, partially unknown dates in the same year
#'   as index (and didn't unambiguously end before index), have a start date before index and no end date
#'   (i.e., potentially straddle the index), have an end date before the index plus some input cutoff
#'   (`min_days_after_index`) and no start date (i.e., potentially before or after), or otherwise ambiguous
#'   timing relative to index (i.e., did not qualify for `after` or `before`).
#'   \item `start_after_end`: Medications that have a start date occurring after an end date. If a
#'   medication starts after it ends (start_date > end_date), then we suspect bad data and possibly cannot
#'   accurately categorize the patient's medications relative to the index. There are scenarios in which
#'   flipped dates would push (inaccurately) a medication from one category to another so this is its own
#'   category.
#' }
#'
#' `after` medications were definitely given in the setting of interest (after the index date), `before`
#' medications were definitely not given in the setting of interest, and for the `unknown` medications it is
#' impossible to know whether or not it was given in the setting of interest. `start_after_end` medications
#' are bad data and we cannot confidently assign them one way or another - these are categorized first.
#'
#' In order to accurately relate medications to the index, `index_date` must be granular to the day (i.e.,
#' not NA). How would we determine if a drug starting "2020-11-09" was before or after index_date =
#' "2020-UN-UN"?
#'
#' It is suggested that users derive the start_date_year and end_date_year columns themselves, especially when
#' working with medication data that could have partial dates in it (e.g., data sourced from MA). This is
#' because, in certain circumstances, we can determine whether or not a medication occurs `before` or `after`
#' index using a partial date. For example, when start_date = "2015-07-09", end_date = "NA", end_date_partial =
#' "2015-UN-UN" (unimputable), and index_date = "2016-03-31" the medication clearly ends before the index_date
#' (2015 < 2016). If the user derived end_date_year from end_date and end_date_partial themselves (getting 2015),
#' then it would be categorized as such (in `before`). BUT, if the user did not do this,
#' \code{categorize_med_timing()} derives end_date_year from end_date (getting NA), and the medication ends up
#' in `unknown`. Patients with `unknown` medications cannot have accurate lines of therapy calculated and
#' therefore this patient would be removed by the line of therapy algorithm, even if all their other
#' medications were eligible.
#'
#' Considering the above two points, it is strongly recommended \code{categorize_med_timing()} is used in
#' conjunction with \code{\link[=wrangle_antineoplastics]{wrangle_antineoplastics()}}.
#'
#' Medications that end at least `min_days_after_index` days after the index date are counted as being given
#' in the setting of interest (and therefore are in `after`) as long as they also clearly started before
#' index. This value defaults to 30, meaning that any medication starting before index and ending at least
#' 30 days after index is included in `after`. Medications with a year-only start date in the same year as
#' index will end up in `unknown` as it is not clear whether or not they occurred in the setting of interest
#' or not. In practice:
#' - picking a large, positive number (10000 = 27 years) makes it so that no drugs starting before index
#'   count as occurring in the setting of interest
#' - picking a large, negative number (-10000 = -27 years), makes it so that almost drugs starting before index
#'   count as occurring in the setting of interest (exceptions being undated drugs, year-only dates in the
#'   same year as index, drug start after drug end)
#'
#' @param data A data frame containing the index_date, start_date, and end_date columns at a minimum. It is
#' suggested that users derive the start_date_year and end_date_year columns themselves (see Details).
#' @param index_date Date column in `data` holding the index date to reference.
#' @param start_date Date column in `data` holding the start date of the medication. Defaults to start_date.
#' @param end_date Date column in `data` holding the end date of the medication. Defaults to end_date.
#' @param min_days_after_index Number of days after the index that medications starting before the index must
#' end on or after to qualify the medication as occurring after index. Defaults to 30.
#' @param impute_index_start Boolean for whether or not to impute the start date to be the index date for
#' therapies starting before index and continuing for min_days_after_index. Defaults to TRUE. If FALSE, time
#' from index to 1L will be negative.
#'
#' @return A named list consisting of three data frames dividing `data` by when each medication occurred
#' relative to the index date (`before`, `after`, `unknown`, or `start_after_end`).
#' @export
#'
categorize_med_timing <- function(data,
                                  index_date,
                                  start_date = start_date,
                                  end_date = end_date,
                                  min_days_after_index = 30,
                                  impute_index_start = TRUE) {

  if (missing(index_date)) {
    rlang::abort('argument "index_date" is missing')
  }

  index_date_year <- rlang::sym(paste0(substitute(index_date), "_year"))
  index_date <- rlang::enquo(index_date)
  start_date_year <- rlang::sym(paste0(substitute(start_date), "_year"))
  start_date <- rlang::enquo(start_date)
  end_date_year <- rlang::sym(paste0(substitute(end_date), "_year"))
  end_date <- rlang::enquo(end_date)

  data <- data %>%
    check_for_column(!!index_date, .class = "Date") %>%
    check_for_column(!!start_date, .class = "Date") %>%
    check_for_column(!!end_date, .class = "Date")
  # Patients without missing index dates can't have their medications categorized
  index_date_missing <- data %>%
    dplyr::filter(is.na(!!index_date)) %>%
    dplyr::summarize(n_distinct(patientid)) %>%
    dplyr::pull()
  if (index_date_missing > 0) {
    rlang::abort(glue::glue("Patients missing index date in {rlang::as_label(index_date)}, N = {index_date_missing}."))
  }
  # We expect users to supply the start_date_year and end_date_year columns, since _year can be derived from
  # partial and full dates, but it is not required
  data <- tryCatch(
    data %>%
      check_for_column(!!start_date_year, .class = "numeric"),
    error = function(cond) {
      message(glue::glue("Can't find {rlang::as_label(start_date_year)} in data frame, calculating from {rlang::as_label(start_date)}."))
      data %>% dplyr::mutate(!!start_date_year := lubridate::year(!!start_date))
    })
  data <- tryCatch(
    data %>%
      check_for_column(!!end_date_year, .class = "numeric"),
    error = function(cond) {
      message(glue::glue("Can't find {rlang::as_label(end_date_year)} in data frame, calculating from {rlang::as_label(end_date)}."))
      data %>% dplyr::mutate(!!end_date_year := lubridate::year(!!end_date))
    })
  # Since index date needs to be populated, we don't message so as not too produce extraneous warnings in
  data <- tryCatch(
    data %>%
      check_for_column(!!index_date_year, .class = "numeric"),
    error = function(cond) {
      data %>% dplyr::mutate(!!index_date_year := lubridate::year(!!index_date))
    })

  # If start_date is after end_date, then we flag this separately as bad data since we cannot accurately
  # categorize the patient's medications. There are scenarios in which flipped dates would inaccurately push
  # a medication from before/after to the opposite time frame.
  start_after_end <- data %>%
    dplyr::filter(!!start_date > !!end_date |
                    !!start_date_year > !!end_date_year)
  data <- data %>%
    dplyr::anti_join(start_after_end, by = colnames(data))

  # Finally, categorize relative to index date
  # Medications that start after the index or that end after the cutoff, count
  # as occurring after the index
  # - starts after index
  after_start <- data %>%
    dplyr::filter(!!start_date >= !!index_date |
                    !!start_date_year > !!index_date_year)
  # - starts before index or no start date, but ends after index + cutoff
  after_end <- data %>%
    filter(!!start_date < !!index_date |
             !!start_date_year < !!index_date_year |
             is.na(!!start_date_year)) %>%
    # lubridate::days(min_days_after_index) / lubridate::years() converts days to years
    filter(!!end_date >= !!index_date + lubridate::days(min_days_after_index) |
             # the extra is.na() check is here to prevent a fully dated end date at the beginning of the year from
             # sneaking in if it occurs before the min_days_after_index
             is.na(!!end_date) & !!end_date_year > !!index_date_year + lubridate::days(min_days_after_index) / lubridate::years()) %>%
    ungroup()
  # Imputation of index date will affect this anti join, exclude after therapies before checking
  data <- data %>%
    dplyr::anti_join(bind_rows(after_start, after_end) %>%
                       distinct(),
                     by = colnames(data))
  # Impute start date for treatments starting before index date if flagged for it
  if (impute_index_start) {
    after_end <- after_end %>%
      # if therapy starts before index, but ends after index + cutoff, then we can impute the start date to
      # be the index date (note, we do not impute when start date is unknown)
      mutate(!!start_date := if_else(!!start_date < !!index_date |
                                       !!start_date_year < !!index_date_year,
                                     !!index_date,
                                     !!start_date,
                                     !!start_date),
             !!start_date_year := if_else(!!start_date_year < !!index_date_year, # since the above modifies start date, do not check the start date
                                          !!index_date_year,
                                          !!start_date_year,
                                          !!start_date_year))
  }
  after <- bind_rows(after_start,
                     after_end) %>%
    distinct()
  # Medications that ending before index or that start before the index and end before the index + cutoff,
  # count as occurring before the index
  # - ends before index
  before_end <- data %>%
    dplyr::filter(!!end_date < !!index_date |
                    !!end_date_year < !!index_date_year)
  # - starts before index and end before index + cutoff, note this of course catches some of the cases as
  # caught directly above
  before_start <- data %>%
    filter(!!start_date < !!index_date |
             !!start_date_year < !!index_date_year) %>%
    dplyr::filter(!!end_date < !!index_date + lubridate::days(min_days_after_index) |
                    # floor is necessary to ensure end dates with year only in same year as index aren't counted
                    # e.g., 2020 < 2020 + days(2) / years(), and we don't want that to count
                    !!end_date_year < !!index_date_year + floor(lubridate::days(min_days_after_index) / lubridate::years()))
  before <- bind_rows(before_end,
                      before_start) %>%
    distinct()
  # medications that have entirely unknown dates, partially unknown dates in the same year as index, have
  # a start date before index and no end date, have no start date and an end date after index (before or
  # after the cutoff), or otherwise ambiguous timing relative to index (that do not qualify for after/before)
  unknown <- data %>%
    dplyr::anti_join(before, by = colnames(data))

  list("after" = after, "before" = before, "unknown" = unknown, "start_after_end" = start_after_end)
}

# Exclude patients ineligible for LoT due to their therapy data
# Filters out patients from .data$after that have therapy data that makes it possible we would compute
# inaccurate LoTs:
# - therapies with start dates after end dates (assuming bad data)
# - therapies that have year only dates (cannot determine when exactly it occurred relative to other dates)
# - therapies missing dates (cannot determine when exactly it occurred relative to other dates)
# - (optionally) therapies occurring before the index date - this only intended for use when the index date
#   is the diagnosis date as it can be difficult to determine whether or not these therapies should be
#   included in the LoT algorithm (are these 1L therapies? are they treatment for something else? should all
#   of them be included? are they bad data?) and there shouldn't be too many cases of this. If there are a
#   lot of these cases, then you should investigate further.
exclude_ineligible_therapy_patients <- function(.data, start_date = start_date, end_date = end_date, exclude_before = FALSE, silent = FALSE) {
  # Set columns and check
  start_date <- rlang::enquo(start_date)
  end_date <- rlang::enquo(end_date)
  .data$after %>%
    check_for_column(!!start_date, .class = "Date") %>%
    check_for_column(!!end_date, .class = "Date", return_df = FALSE)
  # Report mutually-exclusive exclusions if desired
  # Order is intentional to remove cases we cannot programmatically resolve (start after end requires CTR
  # review), then cases we aren't likely to resolve programmatically (incomplete dates), and then cases we
  # may be able to programmatically resolve but which a) is optional and b) is probably least impactful.
  if (!silent) {
    print(glue::glue("Excluding ineligible therapy patients for LoT: Initial N = {bind_rows(.data) %>% count_patients() %>% pull()}"))
    .data$start_after_end %>%
      report_patient_count("... {patients} patients with therapy start after therapy end")
    exclude_unknown <- .data$unknown %>%
      bind_rows(.data$after %>%
                  filter(is.na(!!start_date) | is.na(!!end_date)))
    exclude_unknown %>%
      anti_join(.data$start_after_end, by = "patientid") %>%
      report_patient_count("... {patients} additional patients with incomplete therapy dates")
    if (exclude_before) {
      # if we are excluding the before patients, we simply remove them
      .data$before %>%
        anti_join(.data$start_after_end, by = "patientid") %>%
        anti_join(exclude_unknown, by = "patientid") %>%
        report_patient_count("... {patients} additional patients with therapy before index")
    } else {
      # if we are not, we need to account for the fact that some patients only have before drugs and report
      # this number, otherwise the sums will not make sense and the "exclusion" is hidden
      .data$before %>%
        anti_join(.data$start_after_end, by = "patientid") %>%
        anti_join(exclude_unknown, by = "patientid") %>%
        anti_join(.data$after, by = "patientid") %>%
        report_patient_count("... {patients} additional patients that only have therapies before index")
    }
  }
  # Perform exclusions on data
  before_only <- .data$before %>%
    anti_join(.data$start_after_end, by = "patientid") %>%
    anti_join(.data$unknown, by = "patientid") %>%
    anti_join(.data$after, by = "patientid")
  excluded <- .data$start_after_end %>%
    bind_rows(.data$unknown) %>%
    bind_rows(.data$after %>%
                filter(is.na(!!start_date) | is.na(!!end_date))) %>%
    {if (exclude_before) bind_rows(., .data$before) else bind_rows(., before_only)}
  final <- .data$after %>%
    anti_join(excluded, by = "patientid")
  # Report final remaining and combined exclusions if desired
  if (!silent) {
    final %>%
      report_patient_count(glue::glue("{{patients}} patients remaining for LoT (N excluded = {excluded %>% count_patients() %>% pull()})"))
  }
  final
}

# Function to create a treatment line of therapy using the specific regimen window, treatment gap, overlap,
# switch, and maintenance rules.
# Improvements: Find a way to parameterize therapyname (participates in joins and I haven't had success with
# this before...).
create_lot_tx <- function(.data, start_date = start_date, end_date = end_date, regimen_window = 30, treatment_gap = 120, overlaps_allowed = FALSE, switches_allowed = FALSE, maintenance_rules = "none") {
  start_date <- enquo(start_date)
  end_date <- enquo(end_date)
  # Identify therapy records in the regimen window - this is the initial regimen for this lot
  lot <- .data %>%
    group_by(patientid) %>%
    filter(!!start_date <= suppressWarnings(min(!!start_date)) + days(regimen_window)) %>%
    ungroup()
  # Identify therapy records starting after the regimen window - these are therapies that could be in future
  # lots or the current lot
  lot_after <- .data %>%
    group_by(patientid) %>%
    filter(!!start_date > suppressWarnings(min(!!start_date)) + days(regimen_window)) %>%
    ungroup()
  # Addition of any new antineoplastic after the first regimen_days days of a line will trigger the start of
  # a new treatment line, unless that therapy matches one of the cases below:
  # 1) One or more drugs from the initial treatment line continues (even if one or more of the drugs from
  #    the initial treatment line is dropped). Subsequent rules (switch and treatment gap), also apply to
  #    these overlapping drugs.
  #    Optional, defaults to FALSE. Introduced during ovarian expansion.
  #    NOTE: This and the switch therapy rule can and will capture the same therapies.
  if (overlaps_allowed) {
    lot_overlaps <- lot_after %>%
      inner_join(lot %>%
                   group_by(patientid) %>%
                   summarize(line_end = suppressWarnings(max(!!end_date))),
                 by = "patientid",
                 multiple = 'all') %>%
      filter(!!start_date < line_end) %>%
      select(-line_end)
    lot_after <- lot_after %>%
      anti_join(lot_overlaps, by = colnames(lot_overlaps))
    # this makes it so that subsequent rules apply to the overlapping therapies
    lot <- bind_rows(lot, lot_overlaps)
  }
  # 2) The new antineoplastic is an allowed switch therapy for one of the therapies in the line.
  #    Treatment gap rule also applies to switch therapies.
  if (switches_allowed) {
    lot <- check_switch_therapy(lot_after = lot_after,
                                lot = lot,
                                start_date = rlang::quo_name(start_date),
                                end_date = rlang::quo_name(end_date),
                                treatment_gap = treatment_gap,
                                maintenance_rules = maintenance_rules) %>%
      distinct()
    lot_after <- lot_after %>%
      anti_join(lot, by = colnames(lot))
  }
  # 3) Treatment gaps of more than treatment_gap days trigger a new treatment line.
  lot <- check_treatment_gap(lot_after = lot_after,
                             lot = lot,
                             start_date = rlang::quo_name(start_date),
                             end_date = rlang::quo_name(end_date),
                             treatment_gap = treatment_gap,
                             maintenance_rules = maintenance_rules) %>%
    distinct()
  lot_after <- lot_after %>%
    anti_join(lot, by = colnames(lot))
  # 4) Second overlap check now that treatment gaps and switch therapies have been added to account for
  #    therapies that did not overlap prior to a gap or switch in the initial treatment. We need to run this
  #    twice (once at the beginning and once at the end) for this to work - we can't run it once at the end
  #    since the therapies that would overlap the initial treatment won't be counted as such and will instead
  #    take precedence in the switch and treatment gap checks, causing them to return nothing.
  if (overlaps_allowed) {
    lot_overlaps <- lot_after %>%
      left_join(lot %>%
                  group_by(patientid) %>%
                  summarize(line_end = suppressWarnings(max(!!end_date))),
                by = "patientid") %>%
      filter(!!start_date < line_end) %>%
      select(-line_end)
    lot_after <- lot_after %>%
      anti_join(lot_overlaps, by = colnames(lot_overlaps))
    # this makes it so that subsequent rules apply to the overlapping therapies
    lot <- bind_rows(lot, lot_overlaps)
  }
  lot
}

# Function that checks for allowed switch therapies between current lot and theoretical next lot
# Allowed switch therapies are added to the current lot to prevent them from triggering a new lot
# Treatment gaps between switch therapies are allowed
# Specific rules implemented for ovarian maintenance therapies, see check_switch_therapy_ovarian
check_switch_therapy <- function(lot_after, lot, start_date, end_date, treatment_gap = 120, maintenance_rules = "none") {
  if (!maintenance_rules %in% c("ovarian", "none")) {
    rlang::abort(glue::glue("maintenance_rules = '{maintenance_rules}' not available in check_switch_therapy!"))
  }
  start_date <- rlang::sym(start_date)
  end_date <- rlang::sym(end_date)
  # identify therapies in the lot (for preventing therapies from switching with themselves)
  lot_drugs <- lot %>%
    select(patientid, therapyname) %>%
    distinct()
  # identify potential switches
  lot_switch <- lot_after %>%
    # Decided to remove this filter since we will see "interrupted" treatment records: patient receives drug
    # A for a month, then receives drug B outside of the regimen window for a short while, then receives
    # drug A' (allowed switch for A) a month after drug B. Drug B would therefore always be the only drug
    # captured by this filter and the only drug to compare against drug A for allowed switches/treatment
    # gaps, failing (rightfully) to pass through and be added to the lot. Drug A' would never be checked and
    # included. Drug B would also never count as overlapping with drug A because drug A ended before drug B
    # started. Therefore, the lot would only capture drug A, drug B would become the next lot, and drug A'
    # the lot after that.
    # Note: if overlap was set to false in the above scenario, then drug B wouldn't be captured in the lot
    # and therefore eventually end the lot when we check in step 4 of create_lots (even if drug A' was
    # temporarily added to the lot).
  # group_by(patientid) %>%
  # filter(!!start_date == suppressWarnings(min(!!start_date))) %>%
  # ungroup() %>%
  # ensure that drugs do not exchange with themselves (treatment gap rule checks this)
  anti_join(lot_drugs, by = colnames(lot_drugs)) %>%
    # check for switch therapy matches between current lot and next lot starting therapies
    inner_join(lot %>%
                 group_by(patientid) %>%
                 mutate(line_end = suppressWarnings(max(!!end_date))) %>%
                 select(patientid, switch_therapy, line_end) %>%
                 ungroup(),
               # join on switch_therapy for allowed switches
               by = c("patientid", "switch_therapy"),
               multiple = 'all') %>%
    # therapies without allowed switches have switch_therapy left as NA, remove these
    filter(!is.na(switch_therapy))
  # check to see if there is overlap or an allowed treatment gap in switch therapies (treatment gap rule
  # also applies to switch therapies)
  if (maintenance_rules == "none") {
    # we check between line_end and start_date because as long as one of the therapies in the line continues
    # (even if it's not the therapy that switches), we consider the line to be ongoing and a switch can occur
    # and continue the line
    # i.e., there can be a gap of treatment_gap days between the end of a line and the start of an allowed
    # switch therapy
    lot_switch <- lot_switch %>%
      filter(!!start_date <= line_end + days(treatment_gap))
  } else if (maintenance_rules == "ovarian") {
    # ovarian cancer has it's own rules due to the maintenance therapy rules
    lot_switch <- lot_switch %>%
      check_switch_therapy_ovarian(start_date = rlang::quo_name(start_date),
                                   end_date = rlang::quo_name(end_date),
                                   treatment_gap = treatment_gap)
  }
  # Remove/rename columns so we can call this function again (and join the data again)
  lot_switch <- lot_switch %>%
    select(-line_end) %>%
    distinct()
  # if there are switch eligible therapies, remove from after therapies, add to current lot, and check again
  if (dim(lot_switch)[1] > 0) {
    lot_after <- lot_after %>%
      anti_join(lot_switch, by = colnames(lot_switch))
    lot <- bind_rows(lot, lot_switch)
    lot <- lot %>%
      bind_rows(check_switch_therapy(lot_after,
                                     lot,
                                     rlang::quo_name(start_date),
                                     rlang::quo_name(end_date),
                                     treatment_gap,
                                     maintenance_rules))
  }
  # else, end of recursion - return lot
  lot %>%
    distinct()
}

# Function to filter to allowed switch therapies for ovarian cancer
check_switch_therapy_ovarian <- function(.data, start_date, end_date, treatment_gap = 120) {
  start_date <- rlang::sym(start_date)
  end_date <- rlang::sym(end_date)
  # we check between line_end and start_date because as long as one of the therapies in the line continues
  # (even if it's not the therapy that switches), we consider the line to be ongoing and a switch can occur
  # and continue the line
  # i.e., there can be a gap of treatment_gap days between the end of a line and the start of an allowed
  # switch therapy
  .data %>%
    # treatment gaps > 90 days trigger a new line in the 3L+ setting
    filter((tx_lines > 2 & !!start_date <= line_end + days(90)) |
             # treatments gaps > treatment_gap days trigger a new line for all other cases
             (!!start_date <= line_end + days(treatment_gap)))
}

# Function that checks the treatment gap between current lot and theoretical next lot
# Treatment gaps of more than treatment_gap days trigger a new treatment line, therefore if therapy occurs
# prior to the treatment gap limit, then we add it to the current lot to prevent triggering of new lot.
# Specific rules implemented for ovarian maintenance therapies - more can be added to the case_when.
check_treatment_gap <- function(lot_after, lot, start_date, end_date, treatment_gap = 120, maintenance_rules = "none") {
  if (!maintenance_rules %in% c("ovarian", "none")) {
    rlang::warn(glue::glue("maintenance_rules = '{maintenance_rules}' not available in check_treatment_gap!"))
  }
  start_date <- rlang::sym(start_date)
  end_date <- rlang::sym(end_date)
  # identify allowed treatment gaps
  lot_gap <- lot_after %>%
    # Decided to remove this filter since we will see "interrupted" treatment records: patient receives drug
    # A for a month, then receives drug B outside of the regimen window for a short while, then receives
    # drug A' (allowed switch for A) a month after drug B. Drug B would therefore always be the only drug
    # captured by this filter and the only drug to compare against drug A for allowed switches/treatment
    # gaps, failing (rightfully) to pass through and be added to the lot. Drug A' would never be checked and
    # included. Drug B would also never count as overlapping with drug A because drug A ended before drug B
    # started. Therefore, the lot would only capture drug A, drug B would become the next lot, and drug A'
    # the lot after that.
    # Note: if overlap was set to false in the above scenario, then drug B wouldn't be captured in the lot
    # and therefore eventually end the lot when we check in step 4 of create_lots (even if drug A' was
    # temporarily added to the lot).
  # # determine therapies that start the next lot (in theory)
  # group_by(patientid) %>%
  # filter(!!start_date == suppressWarnings(min(!!start_date))) %>%
  # ungroup() %>%
  inner_join(lot %>%
               group_by(patientid) %>%
               mutate(line_end = suppressWarnings(max(!!end_date))) %>%
               select(patientid, therapyname, line_end),
             # join on therapyname for gaps within a therapy
             by = c("patientid", "therapyname"),
             multiple = 'all') %>%
    # these drugs are too ambiguous to allow in the treatment gap rule
    filter(!therapyname %in% c("other", "unknown"))
  # check to see if there is allowed treatment gap in a therapy
  if (maintenance_rules == "none") {
    # we check between line_end and start_date because as long as one of the therapies in the line continues,
    # we consider the line to be ongoing and a therapy can be restarted and continue the line
    # i.e., there can be a gap of treatment_gap days between the end of a line and the restart of a therapy
    lot_gap <- lot_gap %>%
      filter(!!start_date <= line_end + days(treatment_gap))
  } else if (maintenance_rules == "ovarian") {
    # ovarian cancer has it's own rules due to the maintenance therapy rules
    lot_gap <- lot_gap %>%
      check_treatment_gap_ovarian(start_date = rlang::quo_name(start_date),
                                  end_date = rlang::quo_name(end_date),
                                  treatment_gap = treatment_gap)

  }
  # Remove/rename columns so we can call this function again (and join the data again)
  lot_gap <- lot_gap %>%
    select(-line_end) %>%
    distinct()
  # if there are treatment gap eligible therapies, remove from after therapies, add to current lot, and check again
  if (dim(lot_gap)[1] > 0) {
    lot_after <- lot_after %>%
      anti_join(lot_gap, by = colnames(lot_gap))
    lot <- bind_rows(lot, lot_gap)
    lot <- lot %>%
      bind_rows(check_treatment_gap(lot_after,
                                    lot,
                                    rlang::quo_name(start_date),
                                    rlang::quo_name(end_date),
                                    treatment_gap,
                                    maintenance_rules))
  }
  # else, end of recursion - return lot
  lot %>%
    distinct()
}

# Function to filter to allowed switch therapies for ovarian cancer - same as check_switch_therapy_ovarian
# but leaving as separate for now in case any downstream changes
check_treatment_gap_ovarian <- function(.data, start_date, end_date, treatment_gap = 120) {
  start_date <- rlang::sym(start_date)
  end_date <- rlang::sym(end_date)
  # we check between line_end and start_date because as long as one of the therapies in the line continues,
  # we consider the line to be ongoing and a therapy can be restarted and continue the line
  # i.e., there can be a gap of treatment_gap days between the end of a line and the restart of a therapy
  .data %>%
    # treatment gaps > 90 days trigger a new line in the 3L+ setting
    filter((tx_lines > 2 & !!start_date <= line_end + days(90)) |
             # treatments gaps > treatment_gap days trigger a new line for all other cases
             (!!start_date <= line_end + days(treatment_gap)))
}

# Function that wraps specific maintenance rule implementations (functions)
create_lot_mtx <- function(.data, maintenance_rules, start_date = start_date, end_date = end_date, treatment_gap = 120, switches_allowed = FALSE, patient_lines = NULL) {
  if (!maintenance_rules %in% c("ovarian")) {
    rlang::abort(glue::glue("maintenance_rules = '{maintenance_rules}' not available in create_lot_mtx!"))
  }
  if (dim(.data)[1] == 0) {
    # we don't have any maintenance lines since this is the first iteration
    return(.data)
  }
  start_date <- enquo(start_date)
  end_date <- enquo(end_date)
  if (maintenance_rules == "ovarian") {
    .data %>%
      create_lot_mtx_ovarian(patient_lines = patient_lines,
                             start_date = !!start_date,
                             end_date = !!end_date,
                             treatment_gap = treatment_gap,
                             switches_allowed = switches_allowed)
  }
}

# Function to implement ovarian cancer specific maintenance line of therapy rules
create_lot_mtx_ovarian <- function(.data, patient_lines, start_date = start_date, end_date = end_date, treatment_gap = 120, switches_allowed = TRUE) {
  if (is.null(patient_lines)) {
    rlang::abort("patient_lines = NULL in create_lot_mtx_ovarian!")
  }
  start_date <- enquo(start_date)
  end_date <- enquo(end_date)
  # Identify initial maintenance therapies
  lot <- .data %>%
    group_by(patientid) %>%
    filter(!!start_date == suppressWarnings(min(!!start_date))) %>%
    ungroup()
  # Identify therapies after the initial maintenance therapies that could be in the maintenance lot or in
  # future lots
  lot_after <- .data %>%
    group_by(patientid) %>%
    filter(!!start_date > suppressWarnings(min(!!start_date))) %>%
    ungroup()
  # 1) Add PARPi that qualify for combination therapy to the current lot if they started after bevacizumab
  #    - bevacizumab + PARPi combination = bevacizumab and PARPi start on the same day, or PARPi is added to
  #      bevacizumab monotherapy within 120 (1L, 2L) or 90 (3L+) days of treatment line ending and they
  #      overlap for > 1 day (we do not need to check the overlap at this step, we will check it later when
  #      assigning regimens)
  #    - Note: PARPi to bevacizumab switches trigger new treatment line
  #    - Note: We don't have to check the contents of the current lot since we know it's maintenance and
  #            either starts with bevacizumab or PARPi or both and in all cases, the addition of a(nother)
  #            PARPi within the appropriate time window is allowed
  combination_parpi <- lot_after %>%
    inner_join(patient_lines %>%
                 # by definition for maintenance therapy, the prior line must be a treatment line so this
                 # isn't necessary at the moment, but if the rules change we want this in place to force us
                 # to account for it
                 filter(line_type == "Treatment") %>%
                 group_by(patientid) %>%
                 summarize(line_end = suppressWarnings(max(line_end)),
                           .groups = "drop") %>%
                 ungroup(),
               by = "patientid",
               multiple = 'all') %>%
    # we care about the prior treatment line number for determining the PARPi combination window
    filter((switch_therapy == "PARPi" & tx_lines <= 2 & !!start_date <= line_end + days(120)) |
             switch_therapy == "PARPi" & tx_lines > 2 & !!start_date <= line_end + days(90)) %>%
    select(-line_end)
  lot_after <- lot_after %>%
    anti_join(combination_parpi, by = colnames(lot_after))
  lot <- bind_rows(lot, combination_parpi)
  # 2) Check for switch therapies for one of the therapies in the line.
  #    Treatment gap rule also applies to switch therapies.
  if (switches_allowed) {
    lot <- check_switch_therapy(lot_after = lot_after,
                                lot = lot,
                                start_date = rlang::quo_name(start_date),
                                end_date = rlang::quo_name(end_date),
                                treatment_gap = treatment_gap,
                                maintenance_rules = "ovarian") %>%
      distinct()
    lot_after <- lot_after %>%
      anti_join(lot, by = colnames(lot))
  }
  # 3) Check for allowed treatment gaps of maintenance therapies.
  lot <- check_treatment_gap(lot_after = lot_after,
                             lot = lot,
                             start_date = rlang::quo_name(start_date),
                             end_date = rlang::quo_name(end_date),
                             treatment_gap = treatment_gap,
                             maintenance_rules = "ovarian") %>%
    distinct()
  lot_after <- lot_after %>%
    anti_join(lot, by = colnames(lot))
  # 4) Add PARPi that qualify for bevacizumab to PARPi switch to the current lot
  #    - bevacizumab to PARPi switch = PARPi is added to a bevacizumab monotherapy within 120 (1L, 2L) or
  #      90 (3L+) days of bevacizumab ending
  #    - Note: PARPi to bevacizumab switches trigger new treatment line
  #    - Note: This is basically the same thing as #2, but we check it after accounting for switching and
  #            treatment gaps in case the bevacizumab was continued for a long time prior to switching to
  #            PARPi. This also allows us to keep the stricter rule above and modify this rule as necessary.
  bevacizumab_switch_parpi <- lot_after %>%
    inner_join(lot %>%
                 filter(switch_therapy == "bevacizumab") %>%
                 group_by(patientid) %>%
                 summarize(bevacizumab_end = suppressWarnings(max(!!end_date)),
                           .groups = "drop"),
               by = "patientid",
               multiple = 'all') %>%
    filter((switch_therapy == "PARPi" & tx_lines <= 2 & !!start_date <= bevacizumab_end + days(120)) |
             switch_therapy == "PARPi" & tx_lines > 2 & !!start_date <= bevacizumab_end + days(90)) %>%
    select(-bevacizumab_end)
  lot <- bind_rows(lot, bevacizumab_switch_parpi)
  lot
}

# Assigns the line number and type for the input therapies using the patients' history of lots and maintenance rules
# As more maintenance rules are added, we will probably want to split them into individual functions
# Rules are current only available for ovarian cancer
# lot_after = data frame of therapies remaining to be assigned to a line of therapy (i.e., therapies that
#             did not qualify for the current lot)
# patient_lines = "running" data frame of patient line information as lots are derived (patientid, line,
#                 line_type, line_end)
# maintenance_rules = maintenance line of therapy rules to implement (cannot be left unspecified)
# start_date = column containing therapy start dates
# end_date = column containing therapy end dates
assign_line_and_type <- function(lot_after, patient_lines, maintenance_rules, start_date = start_date, end_date = end_date) {
  if (!maintenance_rules %in% c("ovarian")) {
    rlang::abort(glue::glue("maintenance_rules = '{maintenance_rules}' not available in assign_line_type!"))
  }
  start_date <- enquo(start_date)
  end_date <- enquo(end_date)
  # Rules for ovarian cancer are defined at present
  if (maintenance_rules == "ovarian") {
    next_lot_line_type <- lot_after %>%
      # Remove these columns if they exist since we will be deriving these fresh and we don't want the joins
      # to duplicate them
      select(-matches("^line$|line_end|line_type|tx_lines")) %>%
      group_by(patientid) %>%
      # Filter to therapies that start the next lot - these help determine what the next line will be
      filter(!!start_date == suppressWarnings(min(!!start_date))) %>%
      # Next line's type is determined by the number of treatment lines the patient has received and the
      # most recent line's type
      # Get the maximum treatment line number
      inner_join(patient_lines %>%
                   filter(line_type == "Treatment") %>%
                   group_by(patientid) %>%
                   summarize(tx_lines = suppressWarnings(max(line)),
                             .groups = "drop"),
                 by = "patientid",
                 multiple = 'all') %>%
      # Get the most recent line type to determine the next line's type
      inner_join(patient_lines %>%
                   group_by(patientid) %>%
                   filter(line_end == suppressWarnings(max(line_end))) %>%
                   select(patientid, line_type, line_end) %>%
                   ungroup(),
                 by = "patientid") %>%
      # If any are not maintenance therapies, then it is a treatment line
      mutate(line_type = case_when(any(!maintenance_therapy) ~ "Treatment",
                                   # Bevacizumab monotherapy is maintenance when it overlaps with a prior
                                   # treatment line or starts within 90 days of a treatment line ending
                                   all(line_type == "Treatment" & switch_therapy == "bevacizumab" & !!start_date <= line_end + days(90)) ~ "Maintenance",
                                   # PARPi monotherapy is maintenance when it overlaps with a prior
                                   # treatment line or starts within 120 days of a treatment line ending in
                                   # the 1L or 2L setting
                                   all(line_type == "Treatment" & switch_therapy == "PARPi" & tx_lines <= 2 & !!start_date <= line_end + days(120)) ~ "Maintenance",
                                   # PARPi monotherapy is maintenance when it overlaps with a prior
                                   # treatment line or starts within 90 days of a treatment line ending in
                                   # the 3L+ setting
                                   all(line_type == "Treatment" & switch_therapy == "PARPi" & tx_lines > 2 & !!start_date <= line_end + days(90)) ~ "Maintenance",
                                   # Combination therapy is maintenance when bevacizumab overlaps with a
                                   # prior treatment line and PARPi starts within 120 days or bevacizumab
                                   # starts within 90 days of the prior treatment line and PARPi starts
                                   # within 120 days in the 1L or 2L setting (or 90 days in the 3L+ setting)
                                   # Note: We only have to check bevacizumab in this scenario because it has
                                   # the stricter criterion (<= 90 days) in all settings (1L, 2L, 3L+)
                                   # Note 2: This could be used on it's own for bevacizumab monotherapy
                                   # second entry in the case_when) and combination, but both are included
                                   # for clarity
                                   any(line_type == "Treatment" & switch_therapy == "bevacizumab" & !!start_date <= line_end + days(90)) ~ "Maintenance",
                                   # All other cases are handled as treatment lines, including:
                                   # - bevacizumab to PARPi switch that is >120 (1L, 2L) or >90 days (3L+)
                                   #   after the bevacizumab end date
                                   # - PARPi starting >120 (1L, 2L) >90 days (3L+) days after the end of
                                   #   another PARPi or a treatment line
                                   # - bevacizumab treatment gap of > 120 days
                                   # - PARPi to bevacizumab switch or the addition of bevacizumab to a PARPi
                                   #   line (regardless of the time)
                                   TRUE ~ "Treatment")) %>%
      ungroup() %>%
      # Derive the line number using the most advanced treatment line. For a new treatment line, this is
      # calculated as prior line number + 1 (e.g., 1L -> 2L). For a new maintenance line, the maintenance
      # line number is the same as the most advanced treatment line number, i.e., maintenance line name
      # should follow the corresponding treatment line (e.g., MTx following 2L should be 2L MTx even if it
      # the first maintenance therapy a patient receives).
      mutate(line = case_when(line_type == "Treatment" ~ coalesce(tx_lines + 1, 1),
                              line_type == "Maintenance" ~ coalesce(tx_lines, 1),
                              # Default behavior is to increment by 1 from prior treatment line number
                              TRUE ~ coalesce(tx_lines + 1, 1))) %>%
      distinct(patientid, line, line_type, tx_lines)
  }
  # Attach the assigned line and line type for next lot to all remaining therapies (allows us to avoid
  # having assign it again to therapies that don't trigger the next lot but are included)
  lot_after %>%
    select(-matches("^line$|line_end|line_type|tx_lines")) %>%
    left_join(next_lot_line_type, by = "patientid")
}

# Function to create lots
create_lots <- function(.data, start_date = start_date, end_date = end_date, regimen_window = 30, treatment_gap = 120, overlaps_allowed = FALSE, switches_allowed = FALSE, maintenance_rules = "none", patient_lines = NULL) {
  # 1) End of recursion check
  if (dim(.data)[1] == 0) {
    # If there are no more therapies to assign to lines, then we are done
    return(.data)
  }
  # 2) data checks (column checks) -> move into wrapper function? otherwise the checks execute on every iteration
  start_date <- enquo(start_date)
  end_date <- enquo(end_date)
  .data %>%
    check_for_column(!!start_date, .class = "Date") %>%
    check_for_column(!!end_date, .class = "Date") %>%
    check_for_column(patientid, .class = c("character", "numeric")) %>%
    check_for_column(therapyname, .class = c("character"), return_df = FALSE)
  if (switches_allowed) {
    .data %>%
      check_for_column(switch_therapy, return_df = FALSE)
  }
  if (maintenance_rules != "none") {
    .data %>%
      check_for_column(maintenance_therapy, return_df = FALSE)
  }
  # 3) Determine the current lot
  # Assign line = 1 if this is the first iteration through
  if (!"line" %in% names(.data)) {
    .data <- .data %>%
      mutate(line = 1)
  }
  if (maintenance_rules == "none") {
    # Derive the current lot according to the treatment line rules
    lot <- .data %>%
      create_lot_tx(start_date = !!start_date,
                    end_date = !!end_date,
                    regimen_window = regimen_window,
                    treatment_gap = treatment_gap,
                    overlaps_allowed = overlaps_allowed,
                    switches_allowed = switches_allowed,
                    maintenance_rules = maintenance_rules) %>%
      distinct()
  } else {
    # We derive the current lot depending on whether it is a treatment line or a maintenance line and
    # according to the number of prior treatment lines received. We use line_type and tx_lines for these
    # respectively. These columns will be used in downstream functions.
    # If this is the first iteration, by definition it is the first line and therefore line_type won't exist
    # as a column yet so we add it and set it to Treatment.
    if (!"line_type" %in% names(.data)) {
      .data <- .data %>%
        mutate(line_type = "Treatment")
    }
    # If this is the first iteration, by definition there have been no prior treatment lines and therefore
    # tx_lines won't exist as a column yet so we add it and set it to 0.
    if (!"tx_lines" %in% names(.data)) {
      .data <- .data %>%
        mutate(tx_lines = 0)
    }
    # Derive the treatment lots
    lot_tx <- .data %>%
      filter(line_type == "Treatment") %>%
      create_lot_tx(start_date = !!start_date,
                    end_date = !!end_date,
                    regimen_window = regimen_window,
                    treatment_gap = treatment_gap,
                    overlaps_allowed = overlaps_allowed,
                    switches_allowed = switches_allowed,
                    maintenance_rules = maintenance_rules) %>%
      distinct()
    # Derive the MTx lots
    lot_mtx <- .data %>%
      filter(line_type == "Maintenance") %>%
      create_lot_mtx(maintenance_rules = maintenance_rules,
                     start_date = !!start_date,
                     end_date = !!end_date,
                     treatment_gap = treatment_gap,
                     switches_allowed = switches_allowed,
                     patient_lines = patient_lines) %>%
      distinct()
    # Combine
    lot <- bind_rows(lot_tx, lot_mtx)
  }
  # 4) Clean-up
  lot_after <- .data %>%
    anti_join(lot, by = colnames(.data))
  # Remove therapies in the current lot that start on/after therapies that would trigger a new lot and add to
  # therapies that will create future lots -> move to the individual create_lot_tx and create_lot_mtx
  # functions?
  lot_after_additions <- lot %>%
    inner_join(lot_after %>%
                 group_by(patientid) %>%
                 summarize(next_line_start = suppressWarnings(min(!!start_date))),
               by = "patientid",
               multiple = 'all') %>%
    filter(!!start_date >= next_line_start) %>%
    select(-next_line_start)
  lot <- lot %>%
    anti_join(lot_after_additions, by = colnames(lot))
  lot_after <- lot_after %>%
    bind_rows(lot_after_additions)
  # For maintenance therapies that were included by the overlap, treatment gap, or treatment switch rules,
  # but that actually start on or after the end of all the non-maintenance therapies in the lot, remove
  # these from the lot and add to the next lot. Single-day maintenance therapies that start and end on the
  # end date of non-maintenance therapies do not count.
  # When there are overlaps, the start date of the maintenance line is the current line end date + 1 day.
  if (maintenance_rules != "none") {
    lot_mtx_start_after <- lot %>%
      filter(maintenance_therapy) %>%
      inner_join(lot %>%
                   filter(!maintenance_therapy) %>%
                   group_by(patientid) %>%
                   summarize(line_end = max(!!end_date)),
                 by = "patientid",
                 multiple = 'all') %>%
      # starts on or after the end of the non-maintenance therapies
      filter(!!start_date >= line_end &
               # and is not a single-day therapy with the same end date as the non-maintenance therapies
               line_end != !!end_date)
    lot <- lot %>%
      anti_join(lot_mtx_start_after, by = colnames(lot))
    lot_mtx_start_after <- lot_mtx_start_after %>%
      # When there are overlaps, the start date of the maintenance line is the current line end date + 1 day
      # In this case, overlap would be when the start date of the maintenance therapy is the end date of the
      # non-maintenance therapies in the current line. The maintenance therapy should be updated to reflect
      # it starting one day after the end of the current line.
      mutate(!!start_date := if_else(!!start_date == line_end,
                                     line_end + days(1),
                                     !!start_date,
                                     !!start_date)) %>%
      select(-line_end)
    lot_after <- lot_after %>%
      bind_rows(lot_mtx_start_after)
  }
  # For maintenance therapies that started in the current lot and continue beyond the rest of the
  # non-maintenance therapies, create two versions: one that ends at the end of the non-maintenance therapies
  # in the current lot and one that starts the day after the end of the non-maintenance therapies in the
  # current lot.
  # This allows us to keep the maintenance therapies in the current lot and also assign them as a
  # maintenance line.
  # When there are overlaps, the start date of the maintenance line is the current line end date + 1 day.
  if (maintenance_rules != "none") {
    lot_mtx_continues <- lot %>%
      filter(maintenance_therapy) %>%
      inner_join(lot %>%
                   filter(!maintenance_therapy) %>%
                   group_by(patientid) %>%
                   summarize(line_end = max(!!end_date)),
                 by = "patientid",
                 multiple = 'all') %>%
      filter(!!end_date > line_end)
    lot <- lot %>%
      anti_join(lot_mtx_continues, by = colnames(lot))
    # Version of the maintenance therapies to keep in the current lot
    lot_mtx <- lot_mtx_continues %>%
      mutate(!!end_date := line_end) %>%
      # if the therapy continues from this LoT to the next, we want to reset the therapy ongoing and stop
      # reason flags so that they only affect the next LoT
      update_continuing_therapy_flags_end()
    lot <- bind_rows(lot, lot_mtx) %>%
      select(-line_end)
    # Version of the maintenance therapies to add to future lot creation
    lot_after_mtx <- lot_mtx_continues %>%
      mutate(!!start_date := line_end + days(1)) %>%
      # if the therapy was continued into the next lot, we want to reset the therapy start reasons so they
      # only affect the current LoT
      update_continuing_therapy_flags_start()
    lot_after <- bind_rows(lot_after_mtx, lot_after) %>%
      select(-line_end)
  }
  # For therapies that overlap with the next lot, create two versions of the therapy: one in the current lot
  # and one in the next lot.
  # When there are overlaps, we keep the line start date but end it at the start date of the next lot - 1 day.
  # We do this after the maintenance continuation check since a maintenance line could continue between two
  # treatment lines.
  lot_continues <- lot %>%
    inner_join(lot_after %>%
                 group_by(patientid) %>%
                 summarize(next_line_start = suppressWarnings(min(!!start_date))),
               by = "patientid",
               multiple = 'all') %>%
    # therapies that end on the same day as the next line start don't count as continuing
    filter(!!end_date > next_line_start)
  lot <- lot %>%
    anti_join(lot_continues, by = colnames(lot))
  # Version of the continuing therapies to keep in the current lot
  lot_continues_keep <- lot_continues %>%
    mutate(!!end_date := next_line_start - days(1)) %>%
    # if the therapy continues from this LoT to the next, we want to reset the therapy ongoing and stop
    # reason flags so that they only affect the next LoT
    update_continuing_therapy_flags_end()
  lot <- bind_rows(lot, lot_continues_keep) %>%
    select(-next_line_start)
  # Version of the continuing therapies to add to future lot creation
  lot_continues_after <- lot_continues %>%
    mutate(!!start_date := next_line_start) %>%
    # if the therapy was continued into the next lot, we want to reset the therapy start reasons so they
    # only affect the current LoT
    update_continuing_therapy_flags_start()
  lot_after <- bind_rows(lot_continues_after, lot_after)  %>%
    select(-next_line_start)
  # Update line number and assign line type for the next lot
  if (maintenance_rules == "none") {
    # When there aren't maintenance rules, the next lot line number = current line number + 1
    lot_after <- lot_after %>%
      mutate(line = line + 1)
  } else {
    # When there are maintenance rules, we need to account for prior line_type and line numbers when
    # assigning next lot line_type and line number
    # Extract current line number and type to add it to a "running" data frame of patient line information
    patient_lines <- lot %>%
      group_by(patientid) %>%
      mutate(line_end = suppressWarnings(max(!!end_date))) %>%
      distinct(patientid, line, line_type, line_end) %>%
      ungroup() %>%
      bind_rows(patient_lines)
    # Assign the line number and type for the next lot to all remaining therapies (assigning to all
    # remaining therapies allows us to avoid having to do it again after we derive the next lot)
    lot_after <- lot_after %>%
      assign_line_and_type(patient_lines = patient_lines,
                           maintenance_rules = maintenance_rules,
                           start_date = !!start_date,
                           end_date = !!end_date)
  }
  # 5) Recursion
  lot %>%
    bind_rows(create_lots(lot_after,
                          start_date = !!start_date,
                          end_date = !!end_date,
                          regimen_window = regimen_window,
                          treatment_gap = treatment_gap,
                          overlaps_allowed = overlaps_allowed,
                          switches_allowed = switches_allowed,
                          maintenance_rules = maintenance_rules,
                          patient_lines = patient_lines))
}

# Helper for investigating data
get_random_patient <- function(.data) {
  patient <- .data %>%
    distinct(patientid) %>%
    sample_n(1)
  .data %>%
    semi_join(patient, by = "patientid")
}
# Helper for investigating data
plot_tx_journey <- function(.data, include_lots = TRUE) {
  .data <- .data %>%
    get_random_patient() %>%
    group_by(patientid) %>%
    mutate(tx_start = min(therapystartdate)) %>%
    ungroup() %>%
    mutate(therapystart = day_difference(tx_start, therapystartdate),
           therapyend = day_difference(tx_start, therapyenddate))
  if (include_lots) {
    .data_lines <- .data %>%
      select(patientid, tx_start, line_name, line_start_date, line_end_date) %>%
      distinct() %>%
      mutate(therapystart = day_difference(tx_start, line_start_date),
             therapyend = day_difference(tx_start, line_end_date),
             therapyname = line_name)
    .data <- .data %>%
      bind_rows(.data_lines)
  }
  patient_selected <- .data %>% pull(patientid) %>% unique()
  .data %>%
    ggplot() +
    geom_segment(aes(x = therapystart, xend = therapyend, y = therapyname, yend = therapyname, color = therapyname)) +
    # plots single day drugs
    geom_point(data = function(x) subset(x, therapystart == therapyend),
               aes(x = therapystart, y = therapyname, color = therapyname)) +
    ggtitle(glue::glue("patientid = {patient_selected}")) +
    labs(x = "Treatment journey (days since start of therapy)",
         y = "Therapy",
         color = "Therapy")
}

# Helper function to update therapyongoing and stopreason when a therapy continues from one LoT to another
# We want to reset these flags to "No" and <NA> respectively in the LoT they originally started in because
# we truncate the end date of the therapy to the start of the next LoT. This would inadvertently make it
# appear as if the therapyongoing and stop reason applied to that end date, when in reality they apply to
# the true end date (which is now in the next LoT).
update_continuing_therapy_flags_end <- function(.data) {
  # If the therapy has an ongoing status (Yes, No, Unknown, NA) and continues into the next LoT, set the
  # ongoing statuses for the current lot to No (regardless of the original status, it is now not ongoing
  # in the context of the current LoT). The next LoT will have it's original status.
  if (any(c("therapyongoing", "therapy_ongoing", "ongoing") %in% names(.data))) {
    therapy_ongoing <- intersect(c("therapyongoing", "therapy_ongoing", "ongoing"), names(.data))[1]
    therapy_ongoing <- rlang::sym(therapy_ongoing)
    .data <- .data %>%
      # if_else is required to preserve data types
      mutate(!!therapy_ongoing := if_else(!is.na(!!therapy_ongoing),
                                          "No",
                                          !!therapy_ongoing))
  }
  # If the therapy has a non-null stop reason (there are many) and continues into the next LoT, set the
  # stop reason statuses for the current lot to NA (regardless of the original status, it is now not
  # stopping for this reason in the context of the current LoT). The next LoT will have it's original status.
  if (any(c("stopreason", "stop_reason") %in% names(.data))) {
    stop_reason <- intersect(c("stopreason", "stop_reason"), names(.data))[1]
    stop_reason <- rlang::sym(stop_reason)
    .data <- .data %>%
      # if_else is required to preserve data types
      mutate(!!stop_reason := if_else(!is.na(!!stop_reason),
                                      NA_character_,
                                      !!stop_reason))
  }
  .data
}

# Helper function to update startreason when a therapy continues from one LoT to another
# We want to reset this flag to "No" in the subsequent LoT because we "add" a therapy continued into the
# next LoT by creating a copy of the therapy that has its start date shifted to be the start of the next
# LoT. This would inadvertently make it appear as if the startreason to that start date, when in reality
# it applies to the true start date (which is now in the prior LoT).
update_continuing_therapy_flags_start <- function(.data) {
  # If the therapy has a start reason (Yes, No, Unknown, NA) and was continued into the next LoT, set the
  # start reason for the next lot to "No" (regardless of the original status, it is now not starting for
  # that reason in the context of the next LoT. The prior LoT will have it's original status.
  if (any(c("startreason", "start_reason") %in% names(.data))) {
    start_reason <- intersect(c("startreason", "start_reason"), names(.data))[1]
    start_reason <- rlang::sym(start_reason)
    .data <- .data %>%
      # if_else is required to preserve data types
      mutate(!!start_reason := if_else(!is.na(!!start_reason),
                                       "No",
                                       !!start_reason))
  }
  .data
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
