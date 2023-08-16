## Initialize Environment
## Load Packages ====
pacman::p_load(
  plotly,
  dtplyr,
  tidyverse,
  dbplyr,
  syhelpr,
  table1,
  glue,
  janitor,
  vroom,
  tidyjson,
  tidytext,
  gt,
  gtsummary,
  botor,
  labelled,
  httr,
  jsonlite,
  syrnapi
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
  odbc::odbc(),
  Driver = "postgresql",
  database = "spmd",
  servername = host,
  Port = port,
  UID = user,
  PWD = rds_client$generate_db_auth_token(host, port, user),
  sslmode = "require",
  MaxVarChar = 65568
)

shared_root = 'https://raw.githubusercontent.com/andrewschrag/analytics/main/'
source(file.path(shared_root, 'func.R'))
source(file.path(shared_root, 'cohort_func.R'))
source(file.path(shared_root, 'lot_func.R'))
#source(file.path(shared_root, 'plots.R'))
#source(file.path(shared_root, 'plotly_style.R'))
source(file.path(shared_root, 'her2_func.R'))
source(file.path(shared_root, 'get_mutation_oc.R'))
source(file.path(shared_root, 'kms_api.R'))
source(file.path(shared_root, 'swimmer_general.R'))
source(file.path(shared_root, 'swimmer_urology.R'))
