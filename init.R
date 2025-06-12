## Initialize Environment
reticulate::use_condaenv('/var/lib/rstudio-server/rstudio-users/aschrag/envs/r-miniconda/')
reticulate::py_require("boto3")

## Load Packages ====
pacman::p_load(
  data.table,
  dtplyr,
  tidyverse,
  dbplyr,
  syhelpr,
  DT,
  DBI,
  RPostgres,
  glue,
  janitor,
  vroom,
  tidyjson,
  tidytext,
  gt,
  gtsummary,
  plotly,
  botor,
  labelled,
  httr,
  jsonlite,
  syrnapi,
  tictoc
)

options(width = 140)


## Set standard plot margins ====
plot_margin <-
  list(
    l = 5,
    r = 5,
    b = 40,
    t = 40,
    pad = 4
  )


#This is the connection config
rds_client = botor_client("rds", region_name = "us-west-2")
host = host = 'spmd-prod.cluster-czsq80p56jgd.us-west-2.rds.amazonaws.com' #'spmd-prod-clone.cluster-czsq80p56jgd.us-west-2.rds.amazonaws.com'
port = '5432'
user = 'view'
spmd <- DBI::dbConnect(
  RPostgres::Postgres(),
  #Driver = "postgresql",
  #database = "spmd",
  dbname = 'spmd', 
  host = host,
  port = port,
  user = user,
  password = rds_client$generate_db_auth_token(host, port, user),
  sslmode = "require"
)

shared_root = 'https://raw.githubusercontent.com/andrewschrag/analytics/main'
source(file.path(shared_root, 'func.R'))
source(file.path(shared_root, 'cohort_func.R'))
source(file.path(shared_root, 'lot_func.R'))
source(file.path(shared_root, 'ads_func.R'))
#source(file.path(shared_root, 'plots.R'))
#source(file.path(shared_root, 'plotly_style.R'))
#source(file.path(shared_root, 'her2_func.R'))
#source(file.path(shared_root, 'get_mutation_oc.R'))
source(file.path(shared_root, 'kms_api.R'))
#source(file.path(shared_root, 'swimmer_general.R'))
#source(file.path(shared_root, 'swimmer_urology.R'))
