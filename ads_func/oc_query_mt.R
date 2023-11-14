oc_query_mt <- function(queries, cohort, con = spmd_con){
  library(doParallel)
  library(parallel)
  cores = parallel::detectCores()
  doParallel::registerDoParallel(cl <- parallel::makeCluster(cores))
  
  results_list <- foreach(i = 1:length(queries), .packages = c("syhelpr", "glue")) %dopar% {
    spmd = con()
    cohort = cohort
    run_query(glue::glue(queries[i]), spmd)
  }
  parallel::stopCluster(cl)
  names(results_list) <- names(oc_queries)
  return(results_list)
}


# oc_query_mt <- function(queries, cohort, con = spmd_con){
# library(doParallel)
# library(parallel)
# cores = parallel::detectCores()
# doParallel::registerDoParallel(cl <- parallel::makeCluster(cores))
# 
# results_list <- foreach(i = 1:length(queries), .packages = c("syhelpr", "glue")) %dopar% {
#   spmd = con()
#   cohort = cohort
#   run_query(glue::glue(queries[i]), spmd)
# }
# parallel::stopCluster(cl)
# names(results_list) <- names(oc_queries)
# return(results_list)
# 
# # foreach(m=isplitCols(X2, chunks=ncores), .combine='cbind',
# #         .packages='missForest') %dopar% {
# }
# 
# 

source(file.path('~', 'Projects', 'utils', 'init.R'))
spmd <- spmd_con('clone', max_char = 128000)

refresh_cohort_counts(T)   
list_cohorts('breast')
