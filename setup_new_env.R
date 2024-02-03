## Script for installing packages into new renv
options(install.packages.compile.from.source = "always")
install.packages(c('pacman', 'devtools'))
devtools::install_github('https://github.com/cran/quantreg')
devtools::install_github('https://github.com/jackwasey/icd')
pacman::p_load(
  'devtools',
  'syhelpr',
  'tidyverse',
  'dbplyr',
  'dtplyr',
  'plotly',
  'glue',
  'janitor',
  'vroom',
  'tidyjson',
  'tidytext',
  'gt',
  'gtsummary',
  'botor',
  'labelled',
  'httr',
  'jsonlite',
  'syrnapi',
  'formatR',
  'rvest',
  'RJSONIO',
  'pdftools',
  'arrow',
  'odbc',
  'paws',
  'reticulate'
)

if (!reticulate::condaenv_exists(here::here('renv/r-miniconda/envs/r-reticulate/'))) {
  reticulate::install_miniconda(path = here::here('renv/r-miniconda/'))
}
reticulate::use_miniconda(here::here('renv/r-miniconda/envs/r-reticulate/bin/python'))
reticulate::py_install('boto3')
