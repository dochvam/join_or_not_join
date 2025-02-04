library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim2_fn.R")

set.seed(524879)

nsim <- 200
ncores <- 16

#### Simulation 2.1: decreasing marginal benefits of 2nd dataset ####
cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  nPA = c(30, 60, 120, 250), 
  S = c(60, 120, 600)
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


#### Simulation 2.2: Effect of data dimensions on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.4, 0.7, 1), # coverage
  alpha0 = c(-1, -1/3, 1/3, 1), # D1 information content: alpha0
  theta0 =c(-6, -4, -2), # D2 information content: theta0
  J = 3, eta = 100
)) %>% 
  bind_rows(as.data.frame(expand.grid(
    xi = c(0.1, 0.4, 0.7, 1), # coverage
    J = c(3, 4, 5, 6), # D1 information content: J
    eta = c(10, 100, 1000), # D2 information content: eta
    alpha0 = 0, theta0 = -3
  ))) %>% 
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

#### Simulation 2.3: Effect of noise, bias in secondary dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.5, 1), # coverage
  sigma = c(0.1, 0.5, 1, 2, 4), # noise in D2
  nPA = c(60, 90, 120)
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



#### Simulation 2.4: Effect of bias on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  zeta = c(0.1, 0.5, 1), # bias in D2
  xi = c(0.1, 0.5, 1),   # coverage
  alpha0 = c(-1, -1/3, 1/3, 1),
  nPA = 60
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
        "sim_output/sim2_4.RDS")

stopCluster(cl)

rm(cl)


