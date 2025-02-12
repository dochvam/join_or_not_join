library(tidyverse)
library(parallel)

capture <- lapply(list.files("support_fn", full.names = T), source)

nsim <- 100
ncores <- 4

#### Compare a priori for sim 1 ####


cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  n1 = c(30, 60, 120, 250), 
  zeta = -6:0 / 4,
  sigma = 0:6/3
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 1000000))

capture <- clusterEvalQ(cl, {
  source("support_fn/apriori_fn.R")
  source("support_fn/sim2_fn.R")
})

apriori_results_list <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  apriori_results_list[[i]] <- compare_apriori_many_sim1(
    specs_df_onerow = specs_df[i, ], 
    nsim = nsim, 
    cl = NULL
  )
  pb$tick()
}
apriori_results <- bind_rows(apriori_results_list) %>% 
  left_join(specs_df, by = "scenario")

stopCluster(cl)
rm(cl)
write_csv(apriori_results, "sim_output/apriori_1_1.csv")

#### Compare a priori for sim 2 ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  nPA = c(30, 60, 120, 250), 
  zeta = -6:0 / 4,
  sigma = 0:6/3
)) %>% 
  mutate(scenario = row_number(),
         seed = 1 + floor(runif(n()) * 1000000))

capture <- clusterEvalQ(cl, {
  source("support_fn/apriori_fn.R")
  source("support_fn/sim2_fn.R")
})

apriori_results_list <- list()
pb <- progress::progress_bar$new(total = nrow(specs_df))
for (i in 1:nrow(specs_df)) {
  apriori_results_list[[i]] <- compare_apriori_many_sim2(
                                specs_df_onerow = specs_df[i, ], 
                                nsim = nsim, 
                                cl
                          )
  pb$tick()
}
apriori_results <- bind_rows(apriori_results_list) %>% 
  left_join(specs_df, by = "scenario")

stopCluster(cl)
rm(cl)
write_csv(apriori_results, "sim_output/apriori_2_1.csv")


