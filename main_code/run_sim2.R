library(spOccupancy)
library(tidyverse)
library(parallel)
source("support_fn/sim2_fn.R")

nsim <- 80
ncores <- 20

#### Simulation 2.1: decreasing marginal benefits of 2nd dataset ####
cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  nPA = c(30, 60, 120, 250), S = c(60, 120, 600)
)) %>% 
  mutate(scenario = row_number())

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
             runtime_results = runtime_results),
        "sim_output/sim2_1.RDS")

stopCluster(cl)

#### Visualizing 2.1 ####

# Plot: improvement to precision of Beta0
b0_prec_plot_2.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  mutate(CI_width = `97.5%` - `2.5%`) %>% 
  select(type, scenario, iter, CI_width) %>% 
  pivot_wider(names_from = type, values_from = CI_width) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(nPA), mean, ymin = lb, ymax = ub, group = S, col = as.factor(S))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0")

# Plot: improvement to MSE of Beta0
b0_mse_plot_2.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(nPA), mean, ymin = lb, ymax = ub, group = S, col = as.factor(S))) +
  geom_pointrange(position = position_dodge(width = 0.3), show.legend = F) +
  geom_line(position = position_dodge(width = 0.3), show.legend = F) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Improvement in MSE beta0 estimate\ndue to data integration") + 
  ggtitle("(A) Improvement in MSE")



# Plot: Improvement in power to detect beta1. I don't like this plot as much
#         bc it's totally contingent on the value of beta1 that we pick
beta1_power_plot_2.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta1") %>% 
  mutate(nonzero = sign(`2.5%`) > 0) %>% 
  group_by(scenario, type) %>% 
  summarize(power = mean(nonzero)) %>% 
  select(type, scenario, power) %>% 
  pivot_wider(names_from = type, values_from = power) %>% 
  mutate(improvement = joint - one) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(nPA), improvement, group = S, col = as.factor(S))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Improvement in power to detect\nnonzero effect of beta1") 



# Plot: Increasing computational cost as data size increases
runtime_plot_2.1 <- bind_rows(runtime_results) %>% 
  group_by(scenario, type) %>% 
  summarize(value = mean(vaule)) %>% 
  select(type, scenario, value) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(improvement = joint - one) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(nPA), improvement, group = S, col = as.factor(S))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Increase in runtime (s)") + 
  ggtitle("(B) Computational cost")


# Save fig 2.1
p2.1 <- arrangeGrob(b0_mse_plot_2.1, runtime_plot_2.1, widths = c(1, 0.8), nrow = 1)
ggsave("plots/Fig2-1.jpg", p2.1, width = 12, height = 4, dpi = 600)


#### Simulation 2.2: Effect of data dimensions on usefulness of 2nd dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.4, 0.7, 1), # coverage
  alpha0 = c(-1, -1/3, 1/3, 1), # D1 information content: alpha0
  theta0 =c(-4, -3, -2), # D2 information content: theta0
  J = 3, eta = 100
)) %>% 
  bind_rows(as.data.frame(expand.grid(
    xi = c(0.1, 0.4, 0.7, 1), # coverage
    J = c(3, 5, 10), # D1 information content: J
    eta = c(20, 100, 500, 2500), # D2 information content: eta
    alpha0 = 0, theta0 = -3
  ))) %>% 
  mutate(scenario = row_number())

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
             runtime_results = runtime_results),
        "sim_output/sim2_2.RDS")

stopCluster(cl)

#### Visualizing 2.2 ####

# Plot: improvement to precision of Beta0
bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  filter(scenario %in% 1:48) %>% 
  mutate(CI_width = `97.5%` - `2.5%`) %>% 
  select(type, scenario, iter, CI_width) %>% 
  pivot_wider(names_from = type, values_from = CI_width) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(alpha0, mean, ymin = lb, ymax = ub, group = theta0, col = as.factor(theta0))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  facet_wrap(~paste0("PA-1 coverage: ", xi), nrow = 1) +
  geom_hline(yintercept = 0) +
  xlab("Avg. logit-scale detection prob. in primary dataset") +
  scale_color_discrete("Avg. detection rate in secondary dataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0")



# Plot: Improvement in power to detect beta1. I don't like this plot as much
#         bc it's totally contingent on the value of beta1 that we pick
bind_rows(estimation_results) %>% 
  filter(param == "beta1") %>% 
  mutate(nonzero = sign(`2.5%`) > 0) %>% 
  group_by(scenario, type) %>% 
  summarize(power = mean(nonzero)) %>% 
  select(type, scenario, power) %>% 
  pivot_wider(names_from = type, values_from = power) %>% 
  mutate(improvement = joint - one) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(nPA), improvement, group = S, col = as.factor(S))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Improvement in power to detect\nnonzero effect of beta1") 






#### Simulation 2.3: Effect of noise, bias in secondary dataset ####

cl <- makeCluster(ncores)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.5, 1), # coverage
  sigma = c(0.1, 0.5, 1, 2), # noise in D2
  zeta = c(-2, -1.5, -1, -0.5, -0.1), # error in D2
  nPA = c(25, 50, 100)
)) %>% 
  mutate(scenario = row_number())

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
             runtime_results = runtime_results),
        "sim_output/sim2_3.RDS")

stopCluster(cl)