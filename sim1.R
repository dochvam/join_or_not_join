
source("sim_explore_fn.R")





dataset_list <- list() 

for (i in 1:10) {
  dataset_list[[i]] <- simulate_joint_data(sim_scenario = 1, dataset_ID = i)
}

results <- fit_nimblemodel(dataset_list, 10, modtype = "joint", 
                           ni = 10000, nb = 5000, nc = 3, nt = 1)

true_values <- dataset_list[[1]]$true_params %>% 
  mutate(param = c("intensity_int", "intensity_beta[1]",
                   "intensity_beta[2]", "intensity_beta[3]",
                   "cam_p_int", "det_beta[2]", "theta0", "theta1",
                   "overdisp"),
         mean = value)

ggplot(results, aes(param, mean)) + 
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_jitter(width = 0.2)) +
  geom_point(data = true_values, col = "red") + 
  coord_flip() +
  theme_minimal()



