library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim2_fn.R")

set.seed(524879)

nsim <- 500
ncores <- 16

#### Simulation 2.1: Number of sites ####
cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  nPA = c(25, 50, 100, 150, 200), 
  S = c(50, 200, 800)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim2_fn.R")
})


results_list <- parLapply(cl, 1:nrow(specs_df), 
                          run_many_sim2_parallel_wrapper, 
                          nsim = nsim, specs_df = specs_df)

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
for (i in 1:nrow(specs_df)) {
  this_result <- results_list[[i]]
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
}


saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim2_1.RDS")

stopCluster(cl)
rm(cl)


#### Simulation 2.2: Information content ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  alpha0 = c(-1, -0.5, 0, 0.5, 1), # D1 information content: alpha0
  theta0 =c(-6, -4, -2) # D2 information content: theta0
)) %>%
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim2_fn.R")
})


results_list <- parLapply(cl, 1:nrow(specs_df), 
                          run_many_sim2_parallel_wrapper, 
                          nsim = nsim, specs_df = specs_df)

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
for (i in 1:nrow(specs_df)) {
  this_result <- results_list[[i]]
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
}


saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim2_2.RDS")

stopCluster(cl)
rm(cl)


#### Simulation 2.3: Coverage of x ####


cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = 1:5/5
)) %>%
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim2_fn.R")
})


results_list <- parLapply(cl, 1:nrow(specs_df), 
                          run_many_sim2_parallel_wrapper, 
                          nsim = nsim, specs_df = specs_df)

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
for (i in 1:nrow(specs_df)) {
  this_result <- results_list[[i]]
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
}


saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim2_3.RDS")

stopCluster(cl)
rm(cl)



#### Simulation 2.4: Effect of noise in secondary dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  sigma = c(0.1, 0.5, 1, 2, 4), # noise in D2
  nPA = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim2_fn.R")
})


results_list <- parLapply(cl, 1:nrow(specs_df), 
                          run_many_sim2_parallel_wrapper, 
                          nsim = nsim, specs_df = specs_df)

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
for (i in 1:nrow(specs_df)) {
  this_result <- results_list[[i]]
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
}


saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim2_4.RDS")

stopCluster(cl)
rm(cl)



#### Simulation 2.5: Effect of bias on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  zeta = -1 * 0:8 / 4, # bias in D2
  nPA = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim2_fn.R")
})


results_list <- parLapply(cl, 1:nrow(specs_df), 
                          run_many_sim2_parallel_wrapper, 
                          nsim = nsim, specs_df = specs_df)

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
for (i in 1:nrow(specs_df)) {
  this_result <- results_list[[i]]
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df =specs_df),
        "sim_output/sim2_5.RDS")

stopCluster(cl)

rm(cl)


