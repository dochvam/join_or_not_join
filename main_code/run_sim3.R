library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
library(PointedSDMs)
library(sp)
library(spatstat)
library(INLA)
library(inlabru)
library(raster)
library(terra)
library(gstat)

source("support_fn/sim3_fn.R")

nsim <- 500
ncores <- 10


#### Simulation 3.1: decreasing marginal benefits of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
    n1 = c(25, 50, 100, 150, 200),
    n2 = c(50, 100, 200)
  )) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim3_fn.R")
})

estimation_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim3(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             runtime_results = runtime_results,
             specs_df = specs_df),
        "sim_output/sim3_1.RDS")

stopCluster(cl)
rm(cl)

#### Simulation 3.3: Effect of noise on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  sigma = c(0.1, 0.5, 1, 2, 4), # noise in D2
  n1 = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim3_fn.R")
})

estimation_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim3(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             runtime_results = runtime_results,
             specs_df =specs_df),
        "sim_output/sim3_3.RDS")

stopCluster(cl)

rm(cl)


#### Simulation 3.4: Effect of bias on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  zeta = -1 * 0:8 / 4, # bias in D2
  n1 = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 100000))

capture <- clusterEvalQ(cl, {
  source("support_fn/sim3_fn.R")
})

estimation_results <- list()
runtime_results <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  this_result <- run_many_sim3(specs_df_onerow = specs_df[i, ], 
                               nsim = nsim, cl)
  estimation_results[[i]] <- this_result$estimation_result
  runtime_results[[i]]    <- this_result$runtime_result
  pb$tick()
}

saveRDS(list(estimation_results = estimation_results,
             runtime_results = runtime_results,
             specs_df =specs_df),
        "sim_output/sim3_4.RDS")

stopCluster(cl)

rm(cl)




