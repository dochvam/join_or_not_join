set.seed(23663962)
library(tidyverse)
library(parallel)

source("sim_explore_fn.R")

nrep <- 40

# make_condition <- function(ncam = 30, cam_pct_xcover = 1, PO_avg_effort = 50,
#                            target = "") {
#   data.frame(ncam = ncam, 
#              # cam_pct_xcover = cam_pct_xcover,
#              PO_avg_effort = PO_avg_effort,
#              target = target)
# }

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

condition_df <-
    expand.grid(
      PO_avg_effort = 50,
      ncam = c(30, 60, 120),
      cam_pct_xcover = c(1, 0.5, 0.1),
      PO_noise = c(0, 0.25, 0.5, 1, 1.5, 2),
      PO_bias = c(0, -0.25, -0.5, -1, -1.5, -2),
      intensity_beta = c(1, 2, 3)
    ) %>% 
    mutate(scenario = row_number(),
           estimation_group = as.numeric(as.factor(paste0(ncam, "_", cam_pct_xcover))))

write_csv(condition_df, "output_v2/conditions.csv")

dataset_list <- list()
ct <- 0

one_dataset <- function(i, cond, condition_df) {
  simulate_joint_data_v2(sim_scenario = cond, 
                         ncamera = condition_df$ncam[cond],
                         cam_pct_xcover = condition_df$cam_pct_xcover[cond],
                         PO_avg_effort = condition_df$PO_avg_effort[cond],
                         PO_error_sd = condition_df$PO_noise[cond],
                         PO_bias = condition_df$PO_bias[cond], 
                         intensity_b1 = condition_df$intensity_beta[cond], 
                         PO_zi = 0.3,
                         dataset_ID = i)
}

cl <- makeCluster(16)
capture <- clusterEvalQ(cl, {
  source("sim_explore_fn.R")
  library(tidyverse)
})

pb <- progress::progress_bar$new(total = nrow(condition_df))
for (cond in 1:nrow(condition_df)) {
  pb$tick()
  dataset_list[[cond]] <- parLapply(cl, 1:nrep, one_dataset, 
                                    cond = cond, condition_df = condition_df)
}

# saveRDS(dataset_list, paste0("output_v2/sim_datasets.RDS"))


dataset_list_tofit <- list()

for (i in 1:max(condition_df$estimation_group)) {
  indices <- which(condition_df$estimation_group == i)
  
  dataset_list_tofit[[i]] <- do.call("c", dataset_list[indices])
}

all_results_joint <- parLapply(cl,
                         X = dataset_list_tofit,
                         fun = fit_nimblemodel_MLE,
                         modtype = "joint",
                         progress = FALSE)
saveRDS(all_results_joint, "output_v2/all_results_joint.RDS")

all_results_camera <- parLapply(cl,
                         X = dataset_list_tofit,
                         fun = fit_nimblemodel_MLE,
                         modtype = "camera_only",
                         progress = FALSE)
saveRDS(all_results_camera, "output_v2/all_results_camera.RDS")

stopCluster(cl)





