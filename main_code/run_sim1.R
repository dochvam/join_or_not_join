library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim1_fn.R")

nsim <- 500
ncores <- 32


#### Simulation 1.1: decreasing marginal benefits of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
    n1 = c(30, 60, 120, 250), n2 = c(250, 500, 1000)
  )) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim1_1.RDS")

stopCluster(cl)

#### Simulation 1.2: Effect of data dimensions on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.4, 0.7, 1),       # coverage
  alpha0_1 = c(-1, -1/3, 1/3, 1), # D1 information content: alpha0_1
  alpha0_2 = c(-1, 0, 1),         # D2 information content: alpha0_2
  sigma = c(0.1, 1)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim1_2.RDS")

stopCluster(cl)

#### Simulation 1.3: Effect of noise on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.5, 1), # coverage
  sigma = c(0.1, 0.5, 1, 2, 4), # noise in D2
  n1 = c(30, 60, 120), n2 = 120
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df =specs_df),
        "sim_output/sim1_3.RDS")

stopCluster(cl)

rm(cl)


#### Simulation 1.4: Effect of bias on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  zeta = c(0.1, 0.5, 1),        # bias in D2
  xi = c(0.1, 0.5, 1), # coverage
  alpha0_1 = c(-1, -1/3, 1/3, 1),
  n1 = 60, 
  n2 = 120
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             specs_df =specs_df),
        "sim_output/sim1_4.RDS")

stopCluster(cl)

rm(cl)




