patient_journey <- swimmer.general <- function(patient, x_min = 6) {
  #patient = 'D501C5AF-872E-4DBF-8340-628128ECFD93'
  #patient='47D6A0E0-B924-4BC7-952F-31CE32DF4968'
  patient_tbl <- tribble( ~ patientid, patient)
  spmd_con()
  
  events.cancer_diagnosis <-
    events.diagnosis <-
    events.molecular_reports <-
    events.encounters <-
    events.oncology <-
    events.therapy <-
    events.procedures <-
    events.labs <-
    events.imaging <-
    events.end_date <-
    tribble(
      ~ patientid,
      ~ index_date,
      ~ event_date,
      ~ event_end_date,
      ~ n_days,
      ~ n_months,
      ~ months_since_index,
      ~ event,
      ~ event_text
    )
  
  ## Tooltips ====
  tooltip <- list(
    bordercolor = "transparent",
    font = list(
      color = "#ffffff",
      size = 14,
      family = 'Spectral'
    ),
    align = 'left'
  )
  ## Mapping of SNOMED Result Codes ====
  snomed_map <-
    list(
      snomed_code = c(260373001,
                      260415000,
                      419984006,
                      125154007,
                      10828004),
      result_desc = c(
        'Detected',
        'Not Detected',
        'Inconclusive',
        'Specimen Unsatisfactory',
        'Detected'
      )
    ) %>%
    as_tibble
  
  ## Systemic Therapies ====
  meds.raw <- get_ca_meds(patients = patient_tbl$patientid) %>%
    filter(!is.na(drugcategory)) %>%
    mutate(
      generic_name = coalesce(drugproduct, ordername),
      tx_start = startdts,
      tx_end = enddts
    ) %>%
    collect
  
  ## Index Dates and Cancer Diagnoses
  # cohort.dataset <- cohort_profile(cohort = patient_tbl)
  patient_query <-  tbl(spmd, in_schema('ca', 'cancer')) %>% 
    filter(patientid %in% !!patient_tbl$patientid)
  cancer <-
    patient_query %>% 
    collect %>% 
    mutate(cancer_type = tumortype)
  index_date.data <- cancer %>%
    filter(!is.na(diagnosisdate),
           !grepl('Unknown|Illdefined', cancer_type)) %>% #, cancer_type == 'GU: Prostate') %>%
    distinct(patientid, index_date = min(diagnosisdate, na.rm = T)) %>%
    filter(patientid == patient)
  
  
  
  ## Diagnosis Events
  events.cancer_diagnosis <- cancer %>%
    filter(!is.na(diagnosisdate),
           !grepl('Unknown|Illdefined', cancer_type)) %>%
    inner_join(index_date.data) %>%
    # left_join(
    #   tbl(spmd, in_schema('mdr', 'stage')) %>%
    #     filter(patientid == patient) %>%
    #     select(patientid, stagegroup) %>%
    #     collect %>%
    #     rename(stage = stagegroup) %>%
    #     clean_stage(stage) %>%
    #     group_by(patientid, tumorid) %>%
    #     arrange(patientid, stage) %>%
    #     filter(row_number() == 1) %>%
    #     ungroup,
  #   by = c('patientid', 'tumorid')
  # ) %>%
  mutate(cancer_type = gsub('.*\\: ', '', cancer_type)) %>%
    group_by(patientid, index_date, diagnosisdate, stagegroup) %>%
    summarise(cancer_type = paste0(sort(unique(
      stringr::str_to_title(cancer_type)
    )), collapse = ', ')) %>%
    ungroup %>%
    distinct(
      patientid,
      index_date,
      event_date = diagnosisdate,
      event_end_date = NA_Date_,
      event = 'Cancer',
      event_text = paste0(cancer_type, ifelse(
        is.na(stagegroup), '', paste0(' :: Stage ', toupper(as.character(stagegroup)))
      )),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    )
  
  
  ## Last Contact
  last_contact <- suppressWarnings(get_last_contact(patient_query))
  
  
  meds <- meds.raw %>%
    inner_join(index_date.data %>% distinct(patientid, index_date)) %>%
    filter(!is.na(generic_name)) %>%
    left_join(last_contact) %>%
    convert_dates('start|end|date') %>%
    mutate(last_contact_date = coalesce(deceaseddate, last_contact_date)) %>%
    filter(tx_start < last_contact_date) %>%
    filter(tx_start >= index_date - days(60), !is.na(generic_name)) %>%
    arrange(tx_start)
  
  ## Generate Tx History for LoT
  tx_hist <-
    construct_tx_history(meds, 'tx_start', 'tx_end', 30) %>%
    mutate(drug = generic_name) %>%
    arrange(tx_start, tx_end)
  tx_hist$generic_name <-
    factor(tx_hist$generic_name, levels = rev(unique(tx_hist$generic_name)))
  
  ## Calculate LoTs ====
  lots <- calculate_lot(tx_hist, 30, 1)
  
  
  ## Combine Meds in LoT to Regimens
  regimens <<- lots %>%
    group_by(patientid, line, window_start, line_end) %>%
    summarise(regimen = paste0(sort(unique(
      stringr::str_to_title(generic_name)
    )), collapse = ', '),
    line = as.character(line)) %>%
    ungroup %>%
    distinct %>%
    left_join(last_contact) %>%
    convert_dates('start|end|date') %>%
    mutate(
      tx_end = line_end,
      tx_start = window_start,
      last_contact_date = coalesce(deceaseddate, last_contact_date),
      tx_end = if_else(tx_end > last_contact_date, last_contact_date, tx_end),
      n_days = syhelpr::day_difference(tx_start, tx_end),
      n_months = round(days(n_days) / months(1), 1),
      line_end = if_else(
        line == max(line) & is.na(line_end),
        last_contact_date,
        line_end
      )
    ) %>%
    suppressWarnings()
  
  
  ## Molecular Reports
  events.molecular_reports <- index_date.data %>%
    inner_join(
      tbl(spmd, in_schema('mdr', 'molecularreport')) %>%
        filter(patientid == patient) %>%
        rename(molecularreportid = id) %>%
        collect %>%
        left_join(
          tbl(spmd, in_schema('mdr', 'molecularbiomarker')) %>%
            filter(patientid == patient) %>%
            distinct(
              molecularreportid,
              patientid,
              platformtechnology,
              biomarkername,
              call
            ) %>%
            collect %>%
            process_call %>%
            filter(logical) %>%
            group_by(molecularreportid, patientid, platformtechnology) %>%
            summarise(biomarkername = paste0(biomarkername, collapse = '<br>')),
          by = c('patientid', 'molecularreportid')
        ) %>%
        distinct(
          patientid,
          labname,
          event_date = reportdate,
          event_end_date = NA_Date_,
          testname = str_to_title(testname),
          biomarkername = replace_na(biomarkername, 'No Mutations Found')
        )
    ) %>%
    mutate(
      months_since_index = round(interval(index_date, event_date) / months(1), 1),
      event = 'Molecular Report',
      # event_text = glue::glue('{labname}: {testname} <br><b>{biomarkername}</b>'),
      event_text = glue::glue('{testname} <br><span style="font-size: .8em;">{biomarkername}</span>'),
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date = NA_Date_,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    )
  
  
  ## Therapy Events ====
  events.therapy <- index_date.data %>%
    inner_join(
      regimens %>% distinct(
        patientid,
        regimen,
        event = 'Systemic Therapy',
        event_text = regimen,
        event_text_lot = paste0(
          'Line ',
          line,
          ': ',
          regimen,
          '<br>',
          "duration of ",
          round(n_months, 1),
          ' months'
        ),
        event_date = tx_start,
        event_end_date = tx_end,
        n_days,
        n_months
      )
    ) %>%
    mutate(months_since_index = round(interval(index_date, event_date) / months(1), 1)) %>% #,
    # event_text = glue::glue(
    #   '<b>{regimen}<br>{months_since_index}</b> months after diagnosis') %>%
    select(-regimen) %>%
    arrange(patientid, event) %>%
    distinct %>%
    select(-matches('tx'))
  
  
  ## Procedures ====
  events.procedures <- index_date.data %>%
    inner_join(tbl(spmd, in_schema('mdr', 'procedure')) %>%
                 filter(patientid == patient) %>%
                 collect) %>%
    filter(suppressWarnings(is.na(as.numeric(procedure)))) %>% 
    convert_dates %>%
    mutate(months_since_index = round(interval(index_date, enddts) / months(1), 1)) %>%
    distinct(
      patientid,
      event_date = startdts,
      event_end_date = enddts,
      event = 'Procedures',
      event_text = str_to_title(procedure),
      # event_text = glue::glue(
      #   '<b>{str_to_title(procedure)}<br>{months_since_index}</b> months after diagnosis'
      # ),
      n_days = day_difference(event_date, event_end_date),
      n_months = NA,
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    arrange(months_since_index)
  
  
  
  ## Encounters
  #nppes <- tbl(spmd, in_schema('ca', 'providerspecialty')) %>% collect
  encounters <- 
    tbl(spmd, in_schema('mdr', 'encounter')) %>%
    filter(
      patientid == patient,
      !is.na(encountertype)
      #encountertype %in% c('Ambulatory', 'Outpatient', 'Inpatient', 'Emergency')
    ) %>%
    distinct(patientid, encounterid = id, encountertype, startdts, enddts) %>%
    left_join(
      tbl(.$src$con, in_schema('mdr', 'providerrelationship')) %>%
        distinct(encounterid, providerid),
      by = c('encounterid')
    ) %>%
    left_join(tbl(.$src$con, in_schema('ca', 'providerspecialty')), by = c('providerid')) %>%
    collect %>%
    convert_dates %>%
    janitor::remove_empty() %>% 
    mutate(
      encountertype = case_when(
        encountertype %in% c('Ambulatory', 'Outpatient', 'Other', NA_character_) ~ 'Outpatient',
        encountertype %in% c('Inpatient') ~ 'Inpatient',
        encountertype %in% c('Emergency', 'Emergency Medicine') ~ 'Emergency',
        TRUE ~ 'Outpatient'
      ))
  
  
  events.encounters <- index_date.data %>%
    inner_join(encounters %>% filter(!grepl(
      'urology|oncology', specialty, ignore.case = T
    )), by = c('patientid')) %>%
    convert_dates() %>%
    mutate(
      event_date = startdts,
      event_end_date = enddts,
      event = encountertype,
      event_text = glue::glue(
        'Encounter {ifelse(!is.na(specialty), glue::glue(": {specialty}"), "")}'
      ),
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    ) %>%
    arrange(months_since_index)
  
  
  events.oncology <- index_date.data %>%
    inner_join(encounters %>% filter(grepl('oncology', specialty, ignore.case = TRUE)), by = c('patientid')) %>%
    convert_dates() %>%
    mutate(
      event_date = startdts,
      event_end_date = enddts,
      event = 'Oncology',
      event_text = encountertype,
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    ) %>%
    arrange(months_since_index)
  
  
  ## Diagnosis ====
  diagnosis.raw <-
    tbl(spmd, in_schema('mdr', 'diagnosis')) %>%
    filter(patientid == patient) %>%
    left_join(
      tbl(spmd, in_schema('mdr', 'procedure')) %>%
        distinct(procedureid = id, encounterid_p = encounterid),
      by = c('procedureid')
    ) %>%
    mutate(encounterid = coalesce(encounterid, encounterid_p)) %>%
    left_join(tbl(spmd, in_schema('mdr', 'providerrelationship')) %>%
                distinct(encounterid, providerid),
              by = c('encounterid')) %>%
    left_join(
      tbl(spmd, in_schema('mdr', 'provider')) %>%
        distinct(providerid = id, npi, specialty_rawvalue),
      by = c('providerid')
    ) %>%
    collect %>%
    rename(icd10 = diagnosis_rawvalue) %>%
    janitor::remove_empty()
  
  
  comorbidities <- diagnosis.raw %>%
    mutate(
      icd10 = gsub('\\.', '', icd10),
      diagnosisdate = coalesce(diagnosisdate, noteddate)
    ) %>%
    filter(!grepl('^C', icd10)) %>%
    inner_join(icd::icd10_map_quan_elix %>%
                 unlist() %>%
                 as_tibble() %>%
                 rename('icd10' = 'value')) %>%
    left_join(icd10()) %>%
    distinct(patientid,
             encounterid,
             icd10,
             icd_description,
             diagnosisdate) %>%
    group_by(patientid, icd_description, diagnosisdate) %>%
    arrange(patientid, diagnosisdate, icd_description) %>%
    filter(row_number() == 1) %>%
    ungroup
  
  
  diagnosis <- diagnosis.raw %>%
    left_join(encounters %>% distinct(encounterid, specialty)) %>%
    mutate(
      icd10 = gsub('\\.', '', icd10),
      specialty = coalesce(specialty, specialty_rawvalue),
      diagnosisdate = coalesce(diagnosisdate, noteddate)
    ) %>%
    filter(!grepl('^C', icd10)) %>%
    left_join(icd10()) %>%
    filter(
      !grepl(
        'encounter|encntr|Gingivitis|dental|malignant|tumors',
        icd_description,
        ignore.case = T
      ),!icd10 %in% comorbidities$icd10
    ) %>%
    distinct(patientid,
             encounterid,
             specialty,
             icd_description,
             diagnosisdate) %>%
    group_by(patientid, icd_description, diagnosisdate) %>%
    arrange(patientid, diagnosisdate, icd_description, specialty) %>%
    filter(row_number() == 1) %>%
    ungroup
  
  
  
  
  events.diagnosis <- bind_rows(
    index_date.data %>%
      inner_join(diagnosis,
                 by = 'patientid') %>%
      mutate(
        event_date = diagnosisdate,
        event_end_date = NA_Date_,
        event = 'Non-Cancer',
        n_days = day_difference(event_date, event_end_date),
        n_months = days(n_days) / months(1),
        months_since_index = round(interval(index_date, event_date) / months(1), 1)
      ) %>%
      distinct(
        patientid,
        index_date,
        event_date,
        event_end_date,
        n_days,
        n_months,
        months_since_index,
        event,
        event_text = str_to_title(icd_description)
      ),
    events.cancer_diagnosis
  ) %>%
    filter(!is.na(event_text)) %>%
    group_by(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event
    ) %>%
    summarise(event_text = paste0(sort(unique(event_text)), collapse = '<br>')) %>%
    ungroup
  
  
  events.comorbid <- index_date.data %>%
    inner_join(comorbidities,
               by = 'patientid') %>%
    mutate(
      event_date = diagnosisdate,
      event_end_date = NA_Date_,
      event = 'Non-Cancer',
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text = str_to_title(icd_description)
    )
  
  ## Mets
  events.mets <- index_date.data %>%
    inner_join(tbl(spmd, in_schema('mdr', 'metastasis')) %>%
                 filter(patientid == patient) %>%
                 collect) %>%
    convert_dates() %>%
    group_by(patientid, index_date, bodysite) %>%
    arrange(patientid, metastasisdate) %>%
    filter(metastasisdate >= index_date, row_number() == 1) %>%
    ungroup %>%
    mutate(
      event_date = metastasisdate,
      event_end_date = NA_Date_,
      event = 'Cancer',
      event_text = coalesce(str_to_title(bodysite), ''),
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    )
  
  
  ## Clinical Labs
  events.labs <- index_date.data %>%
    inner_join(
      tbl(spmd, in_schema('mdr', 'clinicallab')) %>%
        filter(patientid == patient) %>%
        collect %>%
        janitor::remove_empty()
    ) %>%
    convert_dates() %>%
    distinct(
      patientid,
      index_date,
      orderpanel,
      testcomponent,
      resultdts,
      resultvalue,
      resultunit
    ) %>%
    left_join(snomed_map, by = c('resultvalue' = 'snomed_code')) %>%
    mutate(resultvalue = coalesce(result_desc, as.character(resultvalue))) %>%
    group_by(patientid, index_date, resultdts) %>%
    summarise(event_text = paste0(
      str_to_title(testcomponent),
      '  <b>',
      ifelse(!is.na(resultvalue), resultvalue, ''),
      ' ',
      ifelse(!is.na(resultunit), resultunit, ''),
      '</b>',
      collapse = '<br>'
    )) %>%
    ungroup %>%
    mutate(
      event_date = resultdts,
      event_end_date = NA_Date_,
      event = 'Clinical Lab',
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    ) %>%
    arrange(months_since_index)
  
  
  ## Imaging ====
  events.imaging <- index_date.data %>%
    inner_join(
      tbl(spmd, in_schema('mdr', 'imaging')) %>%
        filter(patientid == patient, !is.na(startdts)) %>%
        collect %>%
        janitor::remove_empty()
    ) %>%
    convert_dates() %>%
    mutate(
      event_date = startdts,
      event_end_date = enddts,
      event = 'Imaging',
      event_text = str_to_title(imagingprocedure),
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days,
      n_months,
      months_since_index,
      event,
      event_text
    )
  
  
  
  ## End Dates
  events.end_date <-
    index_date.data %>%
    left_join(last_contact) %>%
    convert_dates() %>%
    mutate(
      event_date = last_contact_date,
      event_end_date = NA_Date_,
      event = 'Last Contact',
      event_text = glue::glue('{ifelse(!is.na(deceaseddate), "Deceased", "")}'),
      n_days = day_difference(event_date, event_end_date),
      n_months = days(n_days) / months(1),
      months_since_index = round(interval(index_date, event_date) / months(1), 1)
    ) %>%
    distinct(
      patientid,
      index_date,
      event_date,
      event_end_date,
      n_days = NA,
      n_months = NA,
      months_since_index,
      event,
      event_text
    )
  
  
  ## All Events ====
  events.data <- bind_rows(
    events.cancer_diagnosis,
    events.diagnosis,
    # %>% filter(event != 'Cancer'),
    events.molecular_reports,
    events.encounters,
    events.oncology,
    #events.urology,
    events.therapy,
    events.procedures,
    events.labs,
    events.imaging,
    events.end_date
  ) %>%
    ungroup %>%
    filter(!is.na(event),!is.na(months_since_index)) %>%
    mutate(
      event_text = glue::glue(
        "{event_text}<br><b><span style='font-size:14; float:right !important;'>{months_since_index}</b></span> <span style='font-size:12; float:right !important;'>months after index</span>"
      )
    ) %>%
    convert_dates %>%
    distinct %>%
    arrange(months_since_index)
  
  
  ####################################################################################################
  #### Define how/which events show up in the plot (only singular events; durations handled below)
  ####################################################################################################
  #events = c(unique(events.data$event))
  events <-
    c(
      "<b>Diagnosis            </b>",
      "Cancer",
      "Molecular Report",
      "Clinical Lab",
      "Imaging",
      "Non-Cancer",
      "<b>Treatment            </b>",
      "Systemic Therapy",
      "Procedures",
      "<b>Encounters          </b>",
      "Oncology",
      #"Urology",
      "Emergency",
      "Inpatient",
      "Outpatient",
      "<b>Followup        </b>",
      "Last Contact"
    )
  
  #scales::show_col(sypalette$hex)
  color_palette <-
    c(setNames(
      c(
        # Diagnosis
        '#FFFFFF',
        '#D00000',
        sypalette$hex[11],
        sypalette$hex[1],
        sypalette$hex[2],
        '#8C99A6',
        #sypalette$hex[1],
        # Treatment
        '#FFFFFF',
        sypalette$hex[8],
        colorspace::lighten(sypalette$hex[8], .2),
        # sypalette$hex[4],
        #'#561D9A',
        # Encounters
        '#FFFFFF',
        '#B5179E',
       # colorspace::darken(sypalette$hex[8], .5),
        '#FF500A',
        #Int'l Orange
        '#FF6D33',
        '#8C99A6',
        #sypalette$hex[4],
        # Last Contact
        '#FFFFFF',
        '#000000',
        ifelse(
          events.end_date$event_text == 'Deceased',
          '#000000',
          sypalette$hex[3]
        ),
        '#840000'
      ),
      c(events, 'Metastasis')
    ))
  
  #scales::show_col(sypalette$hex)
  shape_palette <- c(setNames(
    c(
      # Diagnosis
      'line-ns-open',
      'circle-open',
      'star-triangle-up-open-dot',
      'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      # Treatment
      'triangle-down',
      #'line-ns-open',
      'circle', #'triangle-up-open',
      'circle', #'triangle-up',
      # Encounters
      #'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      'line-ns-open',
      # Last Contact
      #'line-ns-open',
      ifelse(
        events.end_date$event_text == 'Deceased',
        'octagon',
        'triangle-right'
      ),
      'circle-open-dot'
    ),
    c(events, 'Metastasis')
  ))
  events.data$event <-
    factor(events.data$event, levels = events)
  
  ## Output ====
  
  # The initial plot just creates the singular events of interest
  plot.swimmer <- plot_ly(data = events.data,
                          y = ~ event,
                          x = ~ months_since_index)  %>%
    add_markers(
      data = events.data %>% filter(
        !event %in% c('Systemic Therapy', 'Procedures', 'Last Contact')
      ),
      y = ~ event,
      x = ~ months_since_index,
      name = ~ event,
      hovertemplate = ~ glue::glue("{replace_na(event_text, '')}"),
      hoverlabel = tooltip,
      type = 'scatter',
      mode = 'markers',
      symbol = ~ event,
      symbols = shape_palette,
      color = ~ event,
      colors = color_palette,
      marker = list(
        size = 22,
        opacity = 0.65,
        line = list(width = 2.5)
      ),
      showlegend = TRUE
    )
  
  if (nrow(events.comorbid) > 0) {
    plot.swimmer <- plot.swimmer %>%
      add_markers(
        data = events.comorbid,
        y = ~ event,
        x = ~ months_since_index,
        name = 'Comorbidity',
        hovertemplate = ~ glue::glue("{replace_na(event_text, '')}"),
        hoverlabel = tooltip,
        type = 'scatter',
        mode = 'markers',
        symbol = ~ event,
        symbols = shape_palette,
        color = ~ event,
        colors = color_palette,
        marker = list(
          size = 22,
          opacity = 1,
          color = sypalette$hex[9],
          line = list(width = 3)
        ),
        showlegend = TRUE
      )
  }
  
  if (nrow(events.mets) > 0) {
    plot.swimmer <- plot.swimmer %>%
      add_markers(
        data = events.mets,
        y = ~ event,
        x = ~ months_since_index,
        name = 'Metastasis',
        hovertemplate = ~ glue::glue("{replace_na(event_text, '')}"),
        hoverlabel = tooltip,
        type = 'scatter',
        mode = 'markers',
        symbol = 'Metastasis',
        symbols = shape_palette,
        marker = list(
          size = 22,
          opacity = 0.85,
          color = colorspace::darken('#D00000', 0.5),
          line = list(width = 3)
        ),
        showlegend = TRUE
      )
  }
  
  
  
  
  # if (!is.na(sum(events.therapy$n_days, na.rm = TRUE))) {
  #   # events.therapy$event <-
  #   #   factor(events.therapy$event, levels = events)
  #   lot.palette <-
  #     colorRampPalette(c(
  #       colorspace::lighten('#8C99A6', .25),
  #       '#8C99A6',
  #       colorspace::darken('#8C99A6', .25),
  #       colorspace::darken('#8C99A6', .5)
  #     ))(length(events.therapy$event))
  #   for (i in 1:(nrow(events.therapy))) {
  #     plot.swimmer <- plot.swimmer %>%
  #       add_markers(
  #         mode = "lines",
  #         data = events.therapy,
  #         x = c(
  #           events.therapy$months_since_index[i],
  #           events.therapy$months_since_index[i] + events.therapy$n_months[i] -
  #             .25
  #         ),
  #         y = 'Systemic Therapy',
  #         name = 'Line of Therapy',
  #         #hovertemplate = ~ glue::glue('{events.therapy$event_text_lot[i]}'),
  #         hoverlabel = tooltip,
  #         hoverinfo = "none",
  #         opacity = 0.25,
  #         line = list(color = lot.palette[i], width = 22),
  #         marker = list(color = 'rgba(0,0,0,0)'),
  #         showlegend = F
  #       )
  #   }
  # }
  
  plot.swimmer <- plot.swimmer %>%
    add_markers(
      data = events.data %>% filter(event %in% c(
        'Systemic Therapy', 'Procedures', 'Last Contact'
      )),
      y = ~ event,
      x = ~ months_since_index,
      name = ~ event,
      hovertemplate = ~ glue::glue("{replace_na(event_text, '')}"),
      hoverlabel = tooltip,
      type = 'scatter',
      mode = 'markers',
      symbol = ~ event,
      symbols = shape_palette,
      color = ~ event,
      colors = color_palette,
      marker = list(
        size = 18,
        opacity = 0.90,
        line = list(width = 2.5)
      ),
      showlegend = TRUE
    )
  
  min_event <-
    min(events.data$months_since_index, na.rm = TRUE) - 2
  max_event <-
    max(events.data$months_since_index, na.rm = TRUE) + 2
  
  plot.final <- plot.swimmer %>%
    layout(
      #title = c(list(text=patient), font_medium),
      hovermode = 'closest',
      uniformtext = list(minsize = 18, mode = 'hide'),
      margin = plot_margin,
      legend = list(
        font = font_large,
        orientation = "h",
        xanchor = "center",
        x = 0.5,
        y = -0.2
      ),
      yaxis = list(
        title = "",
        showgrid = FALSE,
        tickfont =  list(font = list(size = 16, family = 'Spectral')),
        categoryorder = "array",
        categoryarray = rev(events)
      ),
      xaxis = list(
        title = list(
          text = "Months from Intial Cancer Diagnosis",
          font = list(size = 16, family = 'Spectral')
        ),
        dtick = 6,
        tickfont =  list(font = list(size = 16, family = 'Spectral')),
        showgrid = FALSE,
        range = c(pmax(x_min * -1, min_event, na.rm = TRUE), max_event)
      )
    )
  return(plot.final)
  
}
