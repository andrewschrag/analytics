# SPMD connections ----
# Based on the environment variable ENV which is set at the infrastructure level
env <- Sys.getenv("ENV")
if (env == "") {
  # on RSSP, we do not have ENV set, so we set it here depending on the environment we want to run in
  env <- "prod"
}
spmd <- spmd_con(env, max_char = 32784)
spmd_write <- spmd_con(env, TRUE, user = Sys.getenv("SPMD_USER"))
spmd_oc <- spmd_con(env, max_char = 32784) # special connection to allow us to read in large OC JSON

# S3 buckets ----
s3_bucket_path <- list()
# De-identified
if (env == "prod") {
  s3_bucket_path$deid <- "s3://syapse-deidentify-emr-data/deidentify/ads/"
} else {
  s3_bucket_path$deid <- glue::glue("s3://syapse-deidentify-{env}-emr-data/deidentify/ads/")
}
# Identified (internal)
s3_bucket_path$id <- glue::glue("s3://syapse-ephemeral/internal_ads/{env}/")

# OC Value map ----
# TEMPORARY workaround until syhelpr v2.18.0 released (or as needed)
map_oc_value <- tbl(spmd, in_schema("ca", "map_oc_value")) %>%
  collect()

# Reference date ----
reference_date <- as.Date("1970-01-01")

# Lab Name Health System Masking Map ----
# googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1WCAFPvjp3Kv-3Dt_vSjQ4vNE8fZEzjVjEk91rOqLx5I/edit#gid=0",
#                            sheet = "OC Lab Names") %>%
#   bind_rows(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1WCAFPvjp3Kv-3Dt_vSjQ4vNE8fZEzjVjEk91rOqLx5I/edit#gid=0",
#                                       sheet = "SPMD Lab Names")) %>%
#   distinct() %>%
#   write_rds(file = "/var/lib/rstudio-server/rstudio-users/syapse-shared/ADS/maps/labname_hs_masking.rds")
#
#  labname_masking_map <- read_rds(file = "/var/lib/rstudio-server/rstudio-users/syapse-shared/ADS/maps/labname_hs_masking.rds") %>%
#    mutate(labname = case_when(labname == "NA" ~ NA,
#                               is.na(labname) ~ "",
#                              TRUE ~ labname))
#
#  labname_masking_map %>%
#    spmd_write_table("map_labname_hs_masking", spmd_write_con = spmd_con("prod", write = TRUE), overwrite = TRUE)

labname_masking_map <- tbl(spmd, in_schema("ca", "map_labname_hs_masking")) %>%
  collect()

# Primary site codes ----
# Put in each individual script or define collectively here?
# Define in a tibble?
# Source definition from somewhere?
primary_site_codes <- list()
primary_site_codes$ovarian <- c("C56","C56.9", "C57.0", "C48.1", "C48.2")
primary_site_codes$breast <- c("C50")
primary_site_regexes <- list()
primary_site_regexes$ovarian <- create_icd10_regex(primary_site_codes$ovarian)
primary_site_regexes$breast <- create_icd10_regex(primary_site_codes$breast)

# Surgery regex ----
surgery_regex<- list()
surgery_regex$ovarian <- "colectomy|lobectomy|lumpectomy|oophorectomy|salpingo|mastectomy|hysterectomy|resect|debulking|cytoreductive|salpingectomy"

# Bladder Primary Codes ----
primary_site_codes$bladder_cohort <- c("C65.9","C66.9", "C67.0", "C67.1", "C67.2", "C67.3",
                                       "C67.4", "C67.5", "C67.6", "C67.7", "C67.8", "C67.9",
                                       "C68.0", "C68.1", "C68.2", "C68.3",
                                       "C68.4", "C68.5", "C68.6", "C68.7", "C68.8", "C68.9")
primary_site_codes$bladder_only <- c("C67.3", "C67.5", "C67.9", "C67.1", "C67.2", "C67.8", "C67.4", "C67.0", "C67.7")
## MIBC related definitions
primary_site_codes$mibc_unknown_sites <- c("C68.8", "C68.9")
primary_site_codes$mibc_irrelvant_sites <- c("C65.9", "C66.9", "C68.1", "C67.6", "C68.0")

## Primary regexes ----
primary_site_regexes$bladder_cohort <- create_icd10_regex(primary_site_codes$bladder_cohort)
primary_site_regexes$bladder_only <- create_icd10_regex(primary_site_codes$bladder_only)


# Lung Primary Codes
primary_site_codes$lung_cohort <- c("C34.0","C34.1","C34.2","C34.3","C34.4","C34.5","C34.6","C34.7","C34.8","C34.9")
primary_site_regexes$lung_cohort <- create_icd10_regex(primary_site_codes$lung_cohort)


# Bladder-specific definitions
bladder_cohort_name <- "curated"
diagnosis_days_apart <- 45

## primary surgery sites for radical cystectomy
radical_cystectomy_site = c("50bl", "60bl", "61bl", "62bl", "63bl", "64bl", "70bl", "71bl", "72bl" , "73bl", "80bl")
radical_cystectomy_site_display <- map_oc_value %>%
  filter(code %in% radical_cystectomy_site &
           valueset_name == "rxSummSurgPrimSite (r)") %>%
  pull(display) %>%
  unique()
partial_cystectomy_site = c("30bl")
partial_cystectomy_site_display <- map_oc_value %>%
  filter(code %in% partial_cystectomy_site &
           valueset_name == "rxSummSurgPrimSite (r)") %>%
  pull(display) %>%
  unique()

mibc_t_stage <- c( "T2", "T2a", "T2b",
                   "T3", "T3a", "T3b",
                   "T4", "T4a", "T4b")
mibc_stage <- c("II", "IIA", "IIB",
                "III", "IIIA", "IIIB",
                "IV", "IVA", "IVB")
nmibc_t_stage <- c("Tis", "Ta", "T1")
nmibc_stage <- c("0is", "0a", "I")


# Bladder ICD-O-3 tumor histology
# Spreadsheet from https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit?usp=sharing
bladder_cat_histology_map <- tibble(histology_code = c("8000/3",
                                                       "8002/3",
                                                       "8010/2",
                                                       "8010/3",
                                                       "8012/3",
                                                       "8013/3",
                                                       "8020/3",
                                                       "8021/3",
                                                       "8022/3",
                                                       "8031/3",
                                                       "8032/3",
                                                       "8033/3",
                                                       "8035/3",
                                                       "8041/3",
                                                       "8045/3",
                                                       "8046/3",
                                                       "8050/2",
                                                       "8050/3",
                                                       "8052/2",
                                                       "8052/3",
                                                       "8070/2",
                                                       "8070/3",
                                                       "8071/3",
                                                       "8072/3",
                                                       "8074/3",
                                                       "8082/3",
                                                       "8120/2",
                                                       "8120/3",
                                                       "8121/3",
                                                       "8122/3",
                                                       "8130/1",
                                                       "8130/2",
                                                       "8130/3",
                                                       "8131/2",
                                                       "8131/3",
                                                       "8140/2",
                                                       "8140/3",
                                                       "8144/3",
                                                       "8230/3",
                                                       "8246/3",
                                                       "8249/3",
                                                       "8255/3",
                                                       "8260/2",
                                                       "8260/3",
                                                       "8265/3",
                                                       "8310/3",
                                                       "8312/3",
                                                       "8317/3",
                                                       "8323/3",
                                                       "8380/3",
                                                       "8480/3",
                                                       "8490/3",
                                                       "8560/3",
                                                       "8574/3",
                                                       "8680/3",
                                                       "8700/3",
                                                       "8714/3",
                                                       "8720/3",
                                                       "8800/3",
                                                       "8805/3",
                                                       "8858/3",
                                                       "8890/3",
                                                       "8894/3",
                                                       "8900/3",
                                                       "8902/3",
                                                       "8910/3",
                                                       "8963/3",
                                                       "8980/3",
                                                       "9180/3",
                                                       NA_character_),
                                    histology_group = c("Unknown",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Non-urothelial",
                                                        "Unknown"))

# Registry staging algorithm schema regexes ----
# regexes for searching registry records to identify those that match specific cohorts - obtained by
# filtering registry records to cohort-specific primary sites from the curated cohort definitions and
# generating a frequency count of stagingalgorithmschema and primarysite matching those primary sites codes
# and then deciding upon a regex to match the appropriate stagingalgorithmschema values
# NOTE: primarysite in the registry is a period-stripped ICD-O-3 topography code (e.g., "C560")
staging_algorithm_regex <- list()
staging_algorithm_regex$bladder <- "^bladder$|kidney_renal_pelvis|urethra|urinary_other" # (^ and $ wrap bladder so we don't also get "gallbladder")
staging_algorithm_regex$breast <- "breast"
staging_algorithm_regex$lung <- "lung"
staging_algorithm_regex$ovarian <- "fallopian|ovary|peritoneum"

# Breast-specific definitions ----
## Histology map and group updates ----
### This map has been curated by Giovanna Cruz in response to errors found in
### the ca schema breast dictionary. Further work is to be done to update this
### mapping to be year of dx specific as is appropriate, and will replace current
### ca schema mapping once fully curated. For now, Breast ADS uses the mapping below
# breast_histology_map <- tbl(spmd, in_schema("ca", "map_breast_histology_dict")) %>%
#   select(histology_code, histology) %>% # reinclude this once updated # , histology_group = subcategory
#   filter(histology != "Papillary carcinoma, encapsulated of thyroid",
#          !is.na(histology_code)) %>%
#   mutate(histology_map_code = substr(histology_code, 1, 4) %>%  as.numeric()) %>%
# mutate(histology_group = case_when(histology_map_code %in% c(8500, 8501, 8035, 8201, 8022, 8230, 8401, 8502, 8508, 8010) ~ "Ductal carcinoma",
#                                    histology_map_code %in% c(8520, 8519) ~ "Lobular carcinoma",
#                                    histology_map_code %in% c(8530) ~ "Inflammatory carcinoma",
#                                    histology_map_code %in% c(8510, 8512, 8513) ~ "Medullary carcinoma",
#                                    histology_map_code %in% c(8480, 8481) ~ "Mucinous carcinoma",
#                                    histology_map_code %in% c(8503, 8504, 8507, 8509, 8050) ~ "Papillary carcinoma",
#                                    histology_map_code %in% c(8211) ~ "Tubular",
#                                    histology_map_code %in% c(8523, 8524, 8255) ~ "Mixed histologies",
#                                    histology_map_code %in% c(8522) ~ "Duct and lobular carcinoma",
#                                    TRUE ~ "Other"
#                                    )) %>%
#   distinct(histology_code, histology, histology_group) %>%
#   collect()

# Bladder drug map from Bladder Working Data Dictionary: Bladder Systemic Therapies
# https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit#gid=1736968236&range=A73:B104
bladder_drug_map <- tibble(drug_name = c("atezolizumab",
                                         "avelumab",
                                         "bacille calmette-guerin (bcg)",
                                         "bevacizumab",
                                         "cabozantinib",
                                         "carboplatin",
                                         "chemotherapy, nos",
                                         "cisplatin",
                                         "docetaxel",
                                         "doxorubicin",
                                         "durvalumab",
                                         "enfortumab vedotin",
                                         "erdafitinib",
                                         "fluorouracil",
                                         "gemcitabine",
                                         "immunotherapy, nos",
                                         "ipilimumab",
                                         "methotrexate",
                                         "mitomycin",
                                         "nab-paclitaxel",
                                         "nivolumab",
                                         "other",
                                         "paclitaxel",
                                         "palbociclib",
                                         "pembrolizumab",
                                         "pemetrexed",
                                         "sacituzumab govitecan-hziy",
                                         "trial drug",
                                         "unknown",
                                         "valrubicin",
                                         "vinblastine"),
                           drug_cat = c("Biologic therapy (BRM, immunotherapy)",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Targeted therapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Chemotherapy",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Other",
                                        "Chemotherapy",
                                        "Targeted therapy",
                                        "Biologic therapy (BRM, immunotherapy)",
                                        "Chemotherapy",
                                        "Targeted therapy",
                                        "Trial drug",
                                        "Unknown",
                                        "Chemotherapy",
                                        "Chemotherapy"))

# Bladder neoadjuvant chemotherapy regimen map from Bladder Working Data Dictionary: MIBC --> RC neoadjuvant chemo regimens
# https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit#gid=56396634
bladder_nac_map <- tibble(neoadjuvant_analysis_regimen = c("cisplatin + gemcitabine",
                                                           "cisplatin + doxorubicin + methotrexate + vinblastine",
                                                           "carboplatin + gemcitabine",
                                                           "cisplatin",
                                                           "cisplatin + gemcitabine + gemocitabine",
                                                           "gemcitabine",
                                                           "carboplatin + cisplatin + gemcitabine",
                                                           "cisplatin + doxorubicin + vinblastine",
                                                           "cisplatin + doxorubicin + gemcitabine + methotrexate + vinblastine",
                                                           "carboplatin",
                                                           "carboplatin + gemcitabine + gemcitabine",
                                                           "chemotherapy, nos",
                                                           "cisplatin + cisplatin + doxorubicin + doxorubicin + gemcitabine + gemcitabine + methotrexate + vinblastine",
                                                           "cisplatin + cisplatin + doxorubicin + gemcitabine + methotrexate + vinblastine",
                                                           "cisplatin + doxorubicin + gemcitabine",
                                                           "cisplatin + doxorubicin + methotrexate",
                                                           "docetaxel + gemcitabine",
                                                           "doxorubicin + gemcitabine",
                                                           "fluorouracil",
                                                           "cisplatin + gemcitabine + gemcitabine",
                                                           "cisplatin + doxorubicin + fluorouracil + methotrexate + vinblastine"),
                          is_neoadjuvant_chemo = c("Yes",
                                                   "Yes",
                                                   "Yes",
                                                   "No",
                                                   "Yes",
                                                   "No",
                                                   "Yes",
                                                   "Yes",
                                                   "Yes",
                                                   "No",
                                                   "Yes",
                                                   "Unknown",
                                                   "Yes",
                                                   "Yes",
                                                   "Unknown",
                                                   "Yes",
                                                   "No",
                                                   "No",
                                                   "No",
                                                   "Yes",
                                                   "No"))

# Bladder neoadjuvant therapy regimen map from "CA-3991:Final neoadjuvant flag/regimen distribution" tab in
# bladder data dictionary: https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit#gid=1805844252
bladder_non_neoadjuvant_therapy_map <- tribble(~neoadjuvant_regimen_therapies, ~neoadjuvant_therapy_flag,
                                               "cisplatin", "No",
                                               "cisplatin + other", "Unknown",
                                               "carboplatin + other", "Unknown",
                                               "chemotherapy, nos", "Unknown",
                                               "fluorouracil + other", "No",
                                               "gemcitabine", "No")

# Bladder trimodal therapy regimens from CA-3717 ticket
# https://syapse.atlassian.net/browse/CA-3717
# # also includes taxanes only
bladder_trimodal_systemic_therapy_regimens <- c("cisplatin",
                                                "gemcitabine",
                                                "fluorouracil + mitomycin",
                                                "cisplatin + fluorouracil",
                                                "cisplatin + paclitaxel",
                                                "fluorouracil",
                                                "capecitabine")
# + taxane only (use therapy_sublcass)

# bladder TURBT from "primary site surgery map" tab of
# https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit#gid=1372439317
bladder_turbt <- as_tibble(c("02bl",
                             "20bl",
                             "20bl",
                             "21bl",
                             "22bl",
                             "23bl",
                             "24bl",
                             "25bl",
                             "26bl",
                             "27bl")) %>%
  rename(surgery_code = value) %>%
  mutate(surgery = surgery_code) %>%
  map_oc(surgery, "rxSummSurgPrimSite (r)")

# Bladder non-adjuvant therapies
bladder_adjuvant_therapy_other <- c("mitomycin", "fluorouracil")

ignore_grades_bladder <- c("Grade cannot be assessed (GX); Unknown",
                           "Grade/differentiation unknown, not stated, or not applicable",
                           "B-cell",
                           "T-cell",
                           "Null cell",
                           "NK (natural killer) cell",
                           "8",
                           "M",
                           "G",
                           "D",
                           "C",
                           "A")

# Tumor Size sorting schema ----
sort_tumor_size <- c("989 millimeters or larger",
                     "2 mm to 988 mm",
                     "1 mm or described as less than 1 mm",
                     "Microscopic focus or foci only and no size of focus is given",
                     "Diffuse",
                     "Unknown; size not stated; Not documented in patient record; Not applicable",
                     "No mass/tumor found")



# Breast-specific definitions ----
## Histology map and group updates ----
### This map has been curated by Giovanna Cruz in response to errors found in
### the ca schema breast dictionary. Further work is to be done to update this
### mapping to be year of dx specific as is appropriate, and will replace current
### ca schema mapping once fully curated. For now, Breast ADS uses the mapping below
# breast_histology_map <- tbl(spmd, in_schema("ca", "map_breast_histology_dict")) %>%
#   select(histology_code, histology) %>% # reinclude this once updated # , histology_group = subcategory
#   filter(histology != "Papillary carcinoma, encapsulated of thyroid",
#          !is.na(histology_code)) %>%
#   mutate(histology_map_code = substr(histology_code, 1, 4) %>%  as.numeric()) %>%
# mutate(histology_group = case_when(histology_map_code %in% c(8500, 8501, 8035, 8201, 8022, 8230, 8401, 8502, 8508, 8010) ~ "Ductal carcinoma",
#                                    histology_map_code %in% c(8520, 8519) ~ "Lobular carcinoma",
#                                    histology_map_code %in% c(8530) ~ "Inflammatory carcinoma",
#                                    histology_map_code %in% c(8510, 8512, 8513) ~ "Medullary carcinoma",
#                                    histology_map_code %in% c(8480, 8481) ~ "Mucinous carcinoma",
#                                    histology_map_code %in% c(8503, 8504, 8507, 8509, 8050) ~ "Papillary carcinoma",
#                                    histology_map_code %in% c(8211) ~ "Tubular",
#                                    histology_map_code %in% c(8523, 8524, 8255) ~ "Mixed histologies",
#                                    histology_map_code %in% c(8522) ~ "Duct and lobular carcinoma",
#                                    TRUE ~ "Other"
#                                    )) %>%
#   distinct(histology_code, histology, histology_group) %>%
#   collect()

# # Last run 10/21/2022 - Don't need to run again
# # Ovarian Managing Physician Specialty groups ----
# googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1qqlhh3rYMANkwS_1d_B63GK8kPK6S6h0S513mGiIf3w/edit#gid=1505625152",
#                           sheet = "Managing Physician Specialties") %>%
#   select(`Final groupings`, `Specialties to include`) %>%
#   rename(specialty_group = `Final groupings`,
#          specialty = `Specialties to include`) %>%
#   spmd_write_table("map_managing_physician_specialty_ovarian")

# # Last run 11/21/22 - Don't need to run again
# # OncoKB BRCA Map ----
# oncokb_brca1 <- read_csv(paste0(syapse_shared_path, "/data/az-turin/oncokb_brca1.csv"),
#                          col_types = cols()) %>%
#   mutate(gene = "BRCA1")
# oncokb_brca2 <- read_csv(paste0(syapse_shared_path, "/data/az-turin/oncokb_brca2.csv"),
#                          col_types = cols()) %>%
#   mutate(gene = "BRCA2")
# oncokb <- bind_rows(oncokb_brca1,
#                     oncokb_brca2) %>%
#   mutate(alteration = str_replace(alteration, "\\*", "\\\\*"))
# oncokb %>%
#   select(gene, alteration, oncogenic, mutation_effect) %>%
#   distinct() %>%
#   spmd_write_table("map_brca_oncokb", overwrite = TRUE, spmd_write_con = spmd_write)

# # Run on 4/15/2023
# # Bladder Physician Specialty groups ----
# googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1cxVgNUt4Z23ZQoPcACKDmCDZBRcVx8yPJMuBPOEsxxo/edit#gid=1054024096",
#                           sheet = "Surgery/managing physician NPI specialty") %>%
#   select(specialty_group_bladder, specialty) %>%
#   rename(specialty_group = specialty_group_bladder) %>%
#   filter(!is.na(specialty)) %>%
#   spmd_write_table("map_physician_specialty_bladder", spmd_write_con = spmd_con("prod", TRUE), overwrite = TRUE)

# Therapy class and subclass mapping ----
# therapy_class_subclass_map <-
#   googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1kcKPkiTOjnpLFTdRad6Rc2lZiATBcw9OIqGv7OSoUqg/edit#gid=100652015",
#                           sheet = "Therapy Mapping",
#                           skip = 1) %>%
#   select(-c(`ADS class/subclass standardization`, `...6`))
#
# therapy_class_subclass_map %>%
#   spmd_write_table("map_therapy_class_subclass", spmd_write_con = spmd_con("prod", write = TRUE), overwrite = TRUE)

therapy_class_subclass_map <- tbl(spmd, in_schema("ca", "map_therapy_class_subclass")) %>% collect()

# # Lung histology map ----
# # There was a mapping used called ca.map_ads_lung_histology, which was a near-exact recreation of the original
# # ca.map_lung_histology_dict. We could never track down the source of ca.map_ads_lung_histology, as it was
# # indicated that https://docs.google.com/spreadsheets/d/1CoqLNbeYainV1XScUtpO_EILAiwmEl5jb6EZ30Nwuyg/edit#gid=940587027
# # was the source, but that is the exact source for ca.map_lung_histology_dict...
# # So, maps were consolidated to the Lung tab of https://docs.google.com/spreadsheets/d/1e-mvlJd3i0rbO7XfByiz3-KHFo3gtK1MwNO8dxidx2U/edit#gid=769482253
# # and will be maintained from there going forwards
# googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1e-mvlJd3i0rbO7XfByiz3-KHFo3gtK1MwNO8dxidx2U/edit#gid=769482253",
#                           sheet = "Lung") %>%
#   spmd_write_table("map_lung_histology_dict",
#                    spmd_write_con = spmd_con("prod", TRUE),
#                    overwrite = TRUE)
