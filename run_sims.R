set.seed(23663962)

source("sim_explore_fn.R")

nrep <- 100

for (sim in 3:5) {
  dataset_list <- list()
  for (i in 1:nrep) {
    dataset_list[[i]] <- simulate_joint_data(sim_scenario = sim, 
                                             ncamera = 50, PO_zi = 0.3,
                                             dataset_ID = i, holdout_frac = 0.2)
  }
  
  saveRDS(dataset_list, paste0("output/scenario", sim, "_datasets.RDS"))
  
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

