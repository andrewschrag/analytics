## Initialize Environment
## Load Packages ====
pacman::p_load(
  syhelpr,
  table1,
  glue,
  plotly,
  survival,
  survminer,
  janitor,
  vroom,
  tidyverse,
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

options(width = 160)


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
host = 'spmd-prod-clone.cluster-czsq80p56jgd.us-west-2.rds.amazonaws.com'
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


source(file.path('/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils', 'func.R'))
source(file.path('/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils', 'cohort_func.R'))
source(file.path('/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils', 'plots.R'))
source(file.path('/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils', 'plotly_style.R'))
source(file.path('/var/lib/rstudio-server/rstudio-users/syapse-shared/aschrag/utils', 'her2_func.R'))
