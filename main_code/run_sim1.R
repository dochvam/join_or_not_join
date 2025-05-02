library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim1_fn.R")

nsim <- 500
ncores <- 32


#### Simulation 1.1: data dimensions, computational cost ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
    n1 = c(25, 50, 100, 150, 200), n2 = c(50, 200, 800)
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

#### Simulation 1.2: Effect of detection rate on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  alpha0_1 = c(-1, -0.5, 0, 0.5, 1), # D1 information content: alpha0_1
  alpha0_2 = c(-1, 0, 1)             # D2 information content: alpha0_2
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

#### Simulation 1.3: Effect of x coverage ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = 1:5/5
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
        "sim_output/sim1_3.RDS")

stopCluster(cl)

#### Simulation 1.4: Effect of noise on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  sigma = c(0.1, 0.5, 1, 2, 4), # noise in D2
  n1 = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
ppc_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl, ppc = TRUE)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  ppc_results[[i]]        <- this_result$ppc_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             ppc_results = ppc_results,
             specs_df =specs_df),
        "sim_output/sim1_4.RDS")

stopCluster(cl)

rm(cl)


#### Simulation 1.5: Effect of bias on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  zeta = -1 * 0:8 / 4,        # bias in D2
  n1 = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

estimation_results <- list()
cv_results <- list()
runtime_results <- list()
ppc_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim1(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl, ppc = TRUE)
  estimation_results[[i]] <- this_result$estimation_result
  cv_results[[i]]         <- this_result$cv_result
  runtime_results[[i]]    <- this_result$runtime_result
  ppc_results[[i]]        <- this_result$ppc_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             cv_results = cv_results, 
             runtime_results = runtime_results,
             ppc_results = ppc_results,
             specs_df = specs_df),
        "sim_output/sim1_5.RDS")

stopCluster(cl)

rm(cl)


#### Simulation 1.6: Effect of baseline occupancy on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  beta0 = c(-1, -0.5, 0, 0.5, 1)
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
        "sim_output/sim1_6.RDS")

stopCluster(cl)