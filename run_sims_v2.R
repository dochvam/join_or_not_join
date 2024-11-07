set.seed(23663962)
library(tidyverse)

source("sim_explore_fn.R")

nrep <- 50

make_condition <- function(ncam = 30, cam_pct_xcover = 1, PO_avg_effort = 50,
                           target = "") {
  data.frame(ncam = ncam, 
             cam_pct_xcover = cam_pct_xcover,
             PO_avg_effort = PO_avg_effort,
             target = target)
}

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

condition_df <-
  expand.grid.df(
    bind_rows(
      make_condition(ncam = c(10, 30, 50), target = "PA effort"),
      make_condition(cam_pct_xcover = c(1, 0.5, 0.1), target = "PA coverage"),
      make_condition(PO_avg_effort = c(5, 25, 50, 250), target = "PO effort")
    ),
    expand.grid(
      PO_noise = c(0, 0.5, 1),
      PO_bias = c(0, -0.5, -1)
    )
  ) %>% 
  mutate(est_group = as.numeric(as.factor(ncam)))

dataset_list <- list()
for (i in 1:length(unique(condition_df$est_group))) {
  dataset_list[[i]] <- list()
}

ct <- 0
for (cond in 1:nrow(condition_df)) {
  group <- condition_df$est_group[cond]
  for (i in 1:nrep) {
    dataset_list[[group]][[length(dataset_list[[group]]) + 1]] <- 
      simulate_joint_data_v2(sim_scenario = cond, 
                             ncamera = condition_df$ncam[cond],
                             cam_pct_xcover = condition_df$cam_pct_xcover[cond],
                             PO_avg_effort = condition_df$PO_avg_effort[cond],
                             PO_error_sd = condition_df$PO_noise[cond],
                             PO_bias = condition_df$PO_bias[cond], 
                             PO_zi = 0.3,
                             dataset_ID = i)
  }
}

saveRDS(dataset_list, paste0("output_v2/sim_datasets.RDS"))
  
  results <- fit_nimblemodel_MLE(datlist = dataset_list, nsim = nrep, sim_scenario = sim,
                                 modtype = "joint", holdout_frac = 0.2)
  saveRDS(results, paste0("output/scenario", sim, "_joint_results.RDS"))
  
  results_cameraOnly <- fit_nimblemodel_MLE(dataset_list, nsim = nrep, sim_scenario = sim,
                                            modtype = "camera_only", holdout_frac = 0.2)
  saveRDS(results_cameraOnly, paste0("output/scenario", sim, "_cameraOnly_results.RDS"))
  
  results_POOnly <- fit_nimblemodel_MLE(datlist = dataset_list, nsim = nrep, sim_scenario = sim,
                                        modtype = "PO_only", holdout_frac = 0.2)
  saveRDS(results_POOnly, paste0("output/scenario", sim, "_POOnly_results.RDS"))
}

