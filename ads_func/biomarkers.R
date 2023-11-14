# Function that retrieves biomarkers from SPMD with optional arguments for `patients` to filter to,
# whether or not to exclude MA biomarkers (default TRUE since need some updates made to OC biomarker ingestion), 
# and the spmd connection and schema to use. 
# Use this with a set of patients, otherwise it will attempt to return the whole mdr.molecularbiomarker table
# Joins against the specimen table to bring in specimen information
get_biomarkers_spmd <- function(patients = NULL, exclude_ma = TRUE, spmd = spmd_con(), schema = "mdr") {
  # Check patients argument
  if (is.null(patients)) {
    rlang::warn("patients not specified! Getting all biomarker data...")
  } else if (is.data.frame(patients)) {
    patients <- patients %>% 
      check_for_column(patientid) %>% 
      pull(patientid)
  }
  # Biomarker data
  biomarkers <- tbl(spmd, in_schema(schema, "molecularbiomarker")) %>%
    {if (exclude_ma) filter(., sourceschema != "ma") else .} %>%
    {if (!is.null(patients)) filter(., patientid %in% patients) else .} %>% 
    select(-contains("provenance")) %>% 
    collect()
  # Specimen data
  specimens <- tbl(spmd, in_schema(schema, "specimen")) %>%
    {if (exclude_ma) filter(., sourceschema != "ma") else .} %>%
    {if (!is.null(patients)) filter(., patientid %in% patients) else .} %>%   
    # specimens can be associated with other events, want ones directly tied to biomarkers
    filter(!is.na(molecularreportid)) %>%
    mutate(specimendate = coalesce(receiveddate, collectiondate)) %>%
    select(molecularreportid, specimentype, bodysite, specimendate) %>%
    collect()
  # Report data
  reports <- tbl(spmd, in_schema(schema, "molecularreport")) %>% 
    {if (exclude_ma) filter(., sourceschema != "ma") else .} %>%
    {if (!is.null(patients)) filter(., patientid %in% patients) else .} %>% 
    # certain biomarker sources (registry) do not create report ids, do not want to join on NA
    filter(!is.na(id)) %>%
    # at the moment, the only time labname is missing and labname_rawvalue is not, the data is unusable
    # https://syapse.atlassian.net/browse/DQL-245
    mutate(testname = coalesce(testname, testname_rawvalue)) %>% 
    mutate(testname = case_when(testname == "unspecified" ~ "Unknown",
                                testname == "ALK Lung FDA/NY" ~ "ALK Lung",
                                testname == "New York Lung Targeted Profile" ~ "Lung Targeted Profile",
                                TRUE ~ testname)) %>%
    select(molecularreportid = id, labname, testname) %>% 
    collect()
  # Join 
  biomarkers <- biomarkers %>%
    left_join(specimens,
              by = "molecularreportid") %>% 
    left_join(reports,
              by = "molecularreportid") %>%
    left_join(labname_masking_map, by = "labname") %>%
    mutate(labname_normalized = if_else(is.na(labname), "Unknown", labname_normalized, NA_character_))
  
  labname_unmapped <- biomarkers %>% filter(is.na(labname_normalized))
  
  if(nrow(labname_unmapped) > 0){
    rlang::warn(glue::glue("{nrow(labname_unmapped)} biomarker results contain an unmapped lab name. Unmapped values are below:"))
    labname_unmapped %>% distinct(labname) %>% print(n=Inf)
  }  
  
  # Use normalized labname
  biomarkers <-
    biomarkers %>%
    mutate(labname_normalized = coalesce(labname_normalized, labname))
  
  # Impute dates - imputation of report date with specimen date is left up to the user
  biomarkers <- biomarkers %>%
    # partial specimen dates don't exist (no columns for them), convenience for impute_dates
    mutate(specimendate_partial = NA_character_) %>% 
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>%
    mutate(specimendate_year = year(specimendate),
           reportdate_year = case_when(reportdate_granularity %in% c("DAY", "MONTH") ~ year(reportdate),
                                       reportdate_granularity == "YEAR" ~ as.numeric(str_extract(reportdate_partial, "\\d{4}"))))
  # Cleaning up genes to make a singlular gene field for all results other than fusions/rearrangements
  # Going to want to include this step in the documentation
  biomarkers %>% 
    mutate(genes = str_remove_all(genes, "\\{|\\}")) %>% 
    mutate(gene = ifelse(biomarkertype %in% c("Fusion", "Rearrangement"),
                         genes, 
                         str_split(genes, ","))) %>% 
    unnest(gene) %>%
    convert_empty_string(gene) %>%
    distinct()
}

# Function to pull biomarker results from OC for the specified `cancer` site. Optionally filters to a set of
# patients as well. 
get_biomarkers_oc <- function(cancer, patients = NULL, spmd = spmd_con()) {
  # Check patients argument
  if (!is.null(patients) & is.data.frame(patients)) {
    patients <- patients %>% 
      check_for_column(patientid) %>% 
      pull(patientid)
  }
  # The JSON extraction maintains the references (i.e., the report and biomarker correspond to each other)
  # json_array_elements is okay in this scenario since we do not need any of the flagging data and are purely 
  # interested in biomarkers.
  biomarkers <- run_query(glue::glue("SELECT
                                        patientid,
                                        json_array_elements(formdata::json -> 'biomarker_report') -> 'report' AS report,
                                        json_array_elements(json_array_elements(formdata::json -> 'biomarker_report') -> 'biomarker') AS biomarker
                                      FROM openclinica.formdata
                                      WHERE formdata::json -> 'primary_cancer_diagnosis_information' -> 'group' ->> 'ma_site' = '{cancer}'
                                      AND formdata::json -> 'biomarker_report' IS NOT NULL"),
                          connection = spmd) %>%
    {if (!is.null(patients)) filter(., patientid %in% patients) else .} %>% 
    mutate(report = purrr::map(report, fromJSON_na),
           biomarker = purrr::map(biomarker, fromJSON_na)) %>%
    unnest_wider(report) %>%
    unnest_wider(biomarker) %>%
    # reportdate and specimendate are extraneous OC columns telling you the format of the respective dates
    # they will collide with rename_oc when it attempts to convert ma_reportdate to reportdate, producing an
    # error
    select(-reportdate, -specimendate) %>%
    rename_oc()
  
  # Convert empty strings
  biomarkers <- biomarkers %>%
    mutate(across(.cols = where(is.character),
                  .fns = ~ na_if(., "")))
  
  # Map common OC biomarker data items that will be used immediately
  biomarkers <- biomarkers %>% 
    map_oc(specimentype, "specimenType") %>%
    map_oc(specimensite, "specimenSite") %>%
    map_oc(genomicsource, "genomicSource") %>%
    map_oc(genenamesomatic, "geneNameSomatic (r)") %>%
    map_oc(genenamegermline, "geneNameGermline (r)") %>%
    map_oc(alterationtype, "alterationType (r)") %>%
    map_oc(mutationstate, "mutationState (r)") %>%
    map_oc(unknownsignificance, "unknownSignificance (r)") %>%
    map_oc(biomarkerinterp, "interpretation") %>%
    map_oc(expressionstate, "expressionState (r)") %>%
    map_oc(cnvstate, "cnvState (r)") %>%
    map_oc(rearrangementstate, "rearrangementState (r)") %>%
    map_oc(fusionstate, "fusionState (r)") %>%
    map_oc(msistate, "msiState (r)") %>%
    map_oc(tmbstate, "tmbState (r)") %>%
    map_oc(mmrstate, "mmrState (r)") %>%
    map_oc(lohstage, "lohStage (r)") %>%
    map_oc(methylation, "methylation (r)") %>%
    map_oc(hrdresults, "hrdResults (r)") %>% 
    map_oc(labname, valueset_name = "labName") %>% 
    map_oc(testname, valueset_name = "testName") %>% 
    # not in the map
    mutate(mmrstate = case_when(mmrstate == "notchanged" ~ "Not changed", 
                                mmrstate == "abnormal" ~ "Abnormal",
                                TRUE ~ mmrstate)) %>% 
    map_oc(alterationnotstated, "alterationNotStated (r)") %>% 
    # combine genes
    mutate(gene = coalesce(genenamesomatic, genenamegermline)) %>% 
    # combine lab and test name
    mutate(labname = coalesce(labnameshow, labname),
           testname = coalesce(testnameshow, testname)) %>% 
    # standardize unknowns
    mutate(labname = ifelse(labname == "unspecified", "Unknown", labname),
           testname = case_when(testname == "unspecified" ~ "Unknown",
                                testname == "New York Lung Targeted Profile" ~ "Lung Targeted Profile",
                                TRUE ~ testname)) %>%
    left_join(labname_masking_map, by = "labname") %>%
    mutate(labname_normalized = if_else(is.na(labname), "Unknown", labname_normalized, NA_character_))
  
  labname_unmapped <- biomarkers %>% filter(is.na(labname_normalized))
  
  if(nrow(labname_unmapped) > 0){
    rlang::warn(glue::glue("{nrow(labname_unmapped)} biomarker results contain an unmapped lab name. Unmapped values are below:"))
    labname_unmapped %>% distinct(labname) %>% print(n=Inf)
  }  
  # Columns not mapped at the moment since not heavily used and will not appear in the core set of biomarker
  # dfs
  # not mapping these since not going into their own dfs at the moment
  # map_oc(genomicinstability, "genomicInstability (r)") %>%
  # map_oc(tairesults, "taiResults (r)") %>%
  # map_oc(lstresults, "lstResults (r)") %>%
  # map_oc(ca125results, "ca125Results (r)")
  
  biomarkers <- biomarkers %>% 
    mutate(labname_normalized = coalesce(labname_normalized, labname),
           platformtechnology = case_when(platformtechnology == "x" ~ tolower(platformtechnologyother),
                                          TRUE ~ toupper(platformtechnology)),
           platformtechnology = case_when(
             alterationtype == "Expression" & grepl("ihc", platformtechnology, ignore.case=TRUE) ~ "IHC",
             alterationtype == "Expression" & is.na(platformtechnology) ~ "IHC",
             alterationtype == "CNV" & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype == "CNV" & grepl("cish", platformtechnology, ignore.case=TRUE) ~ "CISH",
             alterationtype == "CNV" & grepl("ish", platformtechnology, ignore.case=TRUE) ~ "ISH",
             alterationtype == "Mutation" & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype %in% c("Fusion", "Rearrangement") & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype %in% c("Fusion", "Rearrangement") & grepl("ish", platformtechnology, ignore.case=TRUE) ~ "ISH",
             alterationtype == "MSI" & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype == "TMB" & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype == "LOH" & grepl("seq|pcr|ngs", platformtechnology, ignore.case=TRUE) ~ "NGS",
             alterationtype == "MMR" & grepl("ihc", platformtechnology, ignore.case=TRUE) ~ "IHC" ,
             TRUE ~ platformtechnology)
    )
  
  # Impute dates - note, we are not back-filling missing report dates with specimen date
  biomarkers %>%
    impute_dates(impute_year = FALSE, keep_partial_date_field = TRUE) %>%
    extract_year_oc(reportdate) %>% 
    extract_year_oc(specimendate)
}

# Helper function to select the OC biomarker columns we want to include in all biomarker dfs and rename them
# ... allows you to pass an arbitrary number of columns to additional selection from .data
# Also sets source = "ma"
select_biomarker_columns_oc <- function(.data, ...) {
  .data %>% 
    select(patientid, 
           lab_name = labname_normalized,
           test_name = testname,
           platform_technology = platformtechnology, 
           report_date = reportdate, 
           report_date_granularity = reportdate_granularity, 
           report_date_year = reportdate_year, 
           specimen_date = specimendate, 
           specimen_date_granularity = specimendate_granularity, 
           specimen_date_year = specimendate_year, 
           specimen_type = specimentype,
           specimen_site = specimensite, 
           genomic_source = genomicsource,
           gene,
           ...) %>% 
    mutate(source = "ma")
}

# Helper function to select the SPMD biomarker columns we want to include in all biomarker dfs and rename them
# ... allows you to pass an arbitrary number of columns to additional selection from .data
# Also sets source according to the sourceschema
select_biomarker_columns_spmd <- function(.data, ...) {
  .data %>% 
    mutate(source = case_when(sourceschema == "ma" ~ "ma",
                              sourceschema %in% c("mdx", "nlp") ~ "Structured sources: Third party lab integrations",
                              grepl("registry", sourceschema) ~ "Structured sources: Registry",
                              TRUE ~ "Other")) %>% 
    select(patientid, 
           # molecularreportid, # skipping for now
           lab_name = labname_normalized,
           test_name = testname,
           platform_technology = platformtechnology, 
           report_date = reportdate, 
           report_date_granularity = reportdate_granularity, 
           report_date_year = reportdate_year, 
           specimen_date = specimendate, 
           specimen_date_granularity = specimendate_granularity, 
           specimen_date_year = specimendate_year, 
           specimen_type = specimentype,
           specimen_site = bodysite, 
           genomic_source = genomicsource,
           gene,
           call,
           ...,
           source)
}

# Function that takes a biomarker data frame from SPMD and OC, filters to sequence variants, normalizes amino
# acid and coding changes if possible/desired, converts three letter amino acid representations to single 
# letter versions if possible/desired, and combines into a single data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_sequence_variant_df <- function(bmk_spmd, bmk_oc, normalize_variants = TRUE, convert_aa = FALSE, spmd = spmd_con(), schema = "mdr") {
  
  variants_spmd <- bmk_spmd %>% 
    filter(biomarkertype %in% c("Sequence Variant", "Genetic Alteration", "Wild Type")) %>% 
    # There is an issue wherein some labs send us Wild Type Absent and others send Wild Type Present
    # In both cases, they intend to convey that no variants were found for this gene. So we set to negative.
    # https://syapse.atlassian.net/browse/DQL-102
    mutate(call = ifelse(biomarkertype == "Wild Type", 
                         "Negative", 
                         call)) %>% 
    select_biomarker_columns_spmd(id, aminoacidchange, codingchange, exons, clinicalsignificance) %>% 
    rename(amino_acid_change = aminoacidchange,
           coding_change = codingchange,
           exon = exons,
           clinical_significance = clinicalsignificance,
           mutation_state = call)
  
  variants_oc <- bmk_oc %>% 
    filter(alterationtype == "Mutation") %>% 
    filter(!is.na(gene)) %>% 
    select_biomarker_columns_oc(mutationstate, variantnamep, variantnamec, variantnameexon, biomarkerinterp, unknownsignificance) %>% 
    rename(amino_acid_change = variantnamep,
           coding_change = variantnamec,
           exon = variantnameexon,
           clinical_significance = biomarkerinterp,
           unknown_significance = unknownsignificance,
           mutation_state = mutationstate)
  
  # This introduces "p." and "c." notation where the joins work
  if (normalize_variants) {
    # Pull normalized sequence variants (will be used as a mapping table)
    normalized_variants <- tbl(spmd, in_schema(schema, "sequencevariant")) %>%
      filter(is_normalized) %>%
      filter(!is.na(codingchange_rawvalue) | !is.na(aminoacidchange_rawvalue)) %>%
      select(molecularbiomarkerid, genes, genes_rawvalue, codingchange, codingchange_rawvalue, aminoacidchange, aminoacidchange_rawvalue) %>%
      collect() %>%
      mutate(genes = str_remove_all(genes, "\\{|\\}"),
             genes_rawvalue = str_remove_all(genes, "\\{|\\}")) %>%
      distinct()
    # Join normalized variants onto raw MA variants to see if we can normalize any of them
    variants_oc <- variants_oc %>%
      left_join(normalized_variants %>%
                  select(-genes, -molecularbiomarkerid) %>%
                  distinct(),
                # could do gene = genes, but this didn't producemore normalization in bladder
                by = c("gene" = "genes_rawvalue", "coding_change" = "codingchange_rawvalue", "amino_acid_change" = "aminoacidchange_rawvalue")) %>% 
      mutate(amino_acid_change = coalesce(aminoacidchange, amino_acid_change),
             coding_change = coalesce(codingchange, coding_change)) %>% 
      select(-codingchange, -aminoacidchange)
    variants_spmd <- variants_spmd %>%
      left_join(normalized_variants %>%
                  select(-genes, -genes_rawvalue),
                by = c("id" = "molecularbiomarkerid")) %>%
      mutate(amino_acid_change = coalesce(aminoacidchange, amino_acid_change),
             coding_change = coalesce(codingchange, coding_change)) %>%
      select(-id, -aminoacidchange, -aminoacidchange_rawvalue, -codingchange, -codingchange_rawvalue)
  }
  # Combine, convert AA if specified, and return
  bind_rows(variants_spmd,
            variants_oc) %>% 
    select(-matches("^id$")) %>% 
    relocate(source, .after = last_col()) %>% 
    {if (convert_aa) convert_three_letter_aa(., amino_acid_change) else .} %>% 
    distinct()
}

# Function that takes a biomarker data frame from SPMD and OC, filters to copy number variants and combines 
# into a single data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_copy_number_df <- function(bmk_spmd, bmk_oc) {
  
  cnv_spmd <- bmk_spmd %>% 
    filter(biomarkertype == "Copy Number Variant") %>% 
    select_biomarker_columns_spmd(copynumber) %>% 
    rename(cnv_state = call,
           copy_number = copynumber)
  
  cnv_oc <- bmk_oc %>% 
    filter(alterationtype == "CNV") %>% 
    filter(!is.na(gene)) %>% 
    select_biomarker_columns_oc(cnvstate) %>%
    rename(cnv_state = cnvstate)
  
  bind_rows(cnv_spmd,
            cnv_oc) %>% 
    distinct()
}

# Function that takes a biomarker data frame from SPMD and OC, filters to fusions and rearrangements, 
# creates a fusion_genes column, and and combines into a single data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_fusion_rearrangement_df <- function(bmk_spmd, bmk_oc) {
  
  rfs_spmd <- bmk_spmd %>% 
    filter(biomarkertype %in% c("Rearrangement", "Fusion")) %>% 
    # Converting "," to "::" to denote fusion genes, but not checking if rearrangements or not
    # Then splitting the "," genes (which were left comma-separated for SPMD data) into their own rows
    # This creates a row per gene involved in the fusion/rearrangement
    mutate(fusion_genes = ifelse(grepl("\\,", gene),
                                 str_replace_all(gene, ",", "::"),
                                 NA_character_),
           gene = str_split(gene, ",")) %>% 
    unnest(gene) %>% 
    select_biomarker_columns_spmd(biomarkertype, fusion_genes, clinicalsignificance) %>% 
    rename(alteration_type = biomarkertype,
           state = call,
           clinical_significance = clinicalsignificance)
  
  rfs_oc <- bmk_oc %>% 
    filter(alterationtype %in% c("Rearrangement", "Fusion")) %>% 
    filter(!is.na(gene) | !is.na(fusiongenes)) %>% 
    mutate(state = coalesce(fusionstate, rearrangementstate)) %>% 
    select_biomarker_columns_oc(alterationtype, state, fusiongenes, biomarkerinterp, unknownsignificance) %>%
    # fusion_genes field left as-is becuase the data is highly irregular and hard to normalize
    rename(alteration_type = alterationtype,
           clinical_significance = biomarkerinterp,
           unknown_significance = unknownsignificance,
           fusion_genes = fusiongenes)
  
  bind_rows(rfs_spmd,
            rfs_oc) %>% 
    relocate(source, .after = last_col()) %>% 
    distinct()
  
}

# Function that takes a biomarker data frame from SPMD and OC, filters to expressions, and combines into a single data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_expression_df <- function(bmk_spmd, bmk_oc) {
  
  # Not sure if any of these columns can be consolidated
  expression_spmd <- bmk_spmd %>% 
    filter(biomarkertype == "Protein Expression") %>% 
    select_biomarker_columns_spmd(combinedpositivescore, tumorproportionscore, stainpercent) %>% 
    rename(expression_state = call,
           combined_positive_score = combinedpositivescore,
           tumor_proportion_score = tumorproportionscore,
           stain_percent = stainpercent)
  
  expression_oc <- bmk_oc %>% 
    filter(alterationtype == "Expression") %>% 
    filter(!is.na(gene)) %>% 
    mutate(stain_percent = coalesce(pdl1scorepositivity, score)) %>% 
    select_biomarker_columns_oc(expressionstate, pdl1scorecps, pdl1scoretps, stain_percent) %>%
    rename(expression_state = expressionstate,
           combined_positive_score = pdl1scorecps,
           tumor_proportion_score = pdl1scoretps)
  
  bind_rows(expression_spmd,
            expression_oc) %>% 
    distinct()
}

# Function that takes a biomarker data frame from SPMD and OC, filters to LOH, MSI, MMR, TMB, HRD, and 
# methylation results, and combines into a single data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_tumor_marker_df <- function(bmk_spmd, bmk_oc) {
  
  markers_spmd <- bmk_spmd %>% 
    # methylation results not present in MDR, filtering for it just in case
    filter(biomarkertype %in% c("Loss of Heterozygosity",
                                "Microsatellite Instability",
                                "Mismatch Repair",
                                "Tumor Mutation Burden",
                                "HRD Alteration",
                                "Methylation")) %>% 
    select_biomarker_columns_spmd(biomarkertype, tumormutationburdenscore, tumormutationburdenunit, hrdscore, hrdstatus) %>% 
    mutate(msi_state = ifelse(biomarkertype == "Microsatellite Instability", call, NA_character_),
           tmb_state = ifelse(biomarkertype == "Tumor Mutation Burden", call, NA_character_),
           mmr_state = ifelse(biomarkertype == "Mismatch Repair", call, NA_character_),
           loh_state = ifelse(biomarkertype == "Loss of Heterozygosity", call, NA_character_),
           methylation_state = ifelse(biomarkertype == "Methylation", call, NA_character_), 
           hrd_state = ifelse(biomarkertype == "HRD Alteration", call, NA_character_),
           # hrdstatus and call are duplicates 100% of the time currently, just in case
           hrd_state = coalesce(hrd_state, hrdstatus)) %>% 
    select(-gene, -hrdstatus, -call) %>% 
    rename(alteration_type = biomarkertype,
           tmb_value = tumormutationburdenscore,
           tmb_unit = tumormutationburdenunit,
           hrd_value = hrdscore)
  
  markers_oc <- bmk_oc %>%
    filter(is.na(alterationtype)) %>% 
    filter(alterationtype == "Methylation" |
             gene %in% c("MMR - Mismatch Repair", 
                         "MSI - Microsatellite Instability", 
                         "TMB",
                         "LOH - Loss of Heterozygosity",
                         "HRD")) %>% 
    mutate(alterationtype = case_when(alterationtype == "Methylation" ~ alterationtype,
                                      gene == "LOH - Loss of Heterozygosity" ~ "Loss of Heterozygosity",
                                      gene == "MSI - Microsatellite Instability" ~ "Microsatellite Instability",
                                      gene == "MMR - Mismatch Repair" ~ "Mismatch Repair",
                                      gene == "TMB" ~ "Tumor Mutation Burden",
                                      gene == "HRD" ~ "HRD Alteration")) %>% 
    select_biomarker_columns_oc(alterationtype, msistate, tmbstate, tmbvalue, mmrstate, lohstage, lohvalue, methylation, hrdresults, hrdvalue) %>%
    mutate(tmbvalue = as.numeric(tmbvalue)) %>% 
    select(-gene) %>% 
    rename(alteration_type = alterationtype,
           msi_state = msistate,
           tmb_state = tmbstate,
           tmb_value = tmbvalue,
           mmr_state = mmrstate,
           loh_state = lohstage,
           loh_value = lohvalue,
           methylation_state = methylation,
           hrd_state = hrdresults,
           hrd_value = hrdvalue)
  
  bind_rows(markers_spmd,
            markers_oc) %>% 
    select(patientid, 
           # molecularreportid, # skipping for now
           lab_name, test_name, 
           platform_technology, 
           report_date, report_date_granularity, report_date_year, 
           specimen_date, specimen_date_granularity, specimen_date_year, specimen_type, specimen_site, 
           genomic_source, 
           alteration_type, 
           msi_state, 
           tmb_state, tmb_value, tmb_unit, 
           mmr_state, 
           loh_state, loh_value,
           # methylation_state, #removing for now, as 100% NA
           hrd_state, hrd_value,
           source) %>%
    distinct()
}

# Function that takes a biomarker data frame from SPMD and OC, filters to NOS, and combines into a single 
# data frame
# UPDATES: 
# - Make both data frames optional so that user could only specify one if desired? 
# - Add column/data checks
create_nos_biomarker_df <- function(bmk_spmd, bmk_oc) {
  
  nos_spmd <- bmk_spmd %>% 
    # currently, nos, Not stated, and entirely missing biomarker type, only come from MA so 
    # until MDR holds all the relevant information from MA, this won't return anything
    filter(biomarkertype_rawvalue %in% c("nos", "Not stated") |
             (is.na(biomarkertype_rawvalue) & !is.na(genes))) %>% 
    select_biomarker_columns_spmd() %>% 
    rename(state = call)
  
  nos_oc <- bmk_oc %>% 
    filter(alterationtype %in% c("Not stated", "nos") |
             (is.na(alterationtype) & 
                !gene %in% c("MMR - Mismatch Repair", 
                             "MSI - Microsatellite Instability", 
                             "TMB",
                             "LOH - Loss of Heterozygosity",
                             "HRD",
                             "GENOMICINSTABILITY") &
                !is.na(gene))) %>% 
    select_biomarker_columns_oc(alterationtype, alterationnotstated, cnvstate, expressionstate, mutationstate, alterationnotstated, rearrangementstate, fusionstate) %>% 
    mutate(state = ifelse(alterationtype %in% c("Not stated", "nos"), 
                          alterationnotstated, 
                          coalesce(alterationnotstated, cnvstate, expressionstate, mutationstate, alterationnotstated, rearrangementstate, fusionstate))) %>% 
    select(-c(alterationtype, cnvstate, expressionstate, mutationstate, alterationnotstated, rearrangementstate, fusionstate))
  
  
  bind_rows(nos_spmd,
            nos_oc) %>% 
    distinct()
}
