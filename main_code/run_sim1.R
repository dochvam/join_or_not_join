library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim1_fn.R")
nsim <- 20

#### Simulation 1.1: decreasing marginal benefits of 2nd dataset ####

cl <- makeCluster(10)

specs_df <- as.data.frame(expand.grid(
    n1 = c(30, 60, 120, 250), n2 = c(250, 500, 1000)
  )) %>% 
  mutate(scenario = row_number())

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
             runtime_results = runtime_results),
        "sim_output/sim1_1.RDS")

stopCluster(cl)

#### Visualizing 1.1 ####

# Plot: improvement to precision of Beta0
b0_prec_plot_1.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  mutate(CI_width = Q975 - Q025) %>% 
  select(type, scenario, iter, CI_width) %>% 
  pivot_wider(names_from = type, values_from = CI_width) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(n1), mean, ymin = lb, ymax = ub, group = n2, col = as.factor(n2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0")

# Plot: improvement to MSE of Beta0
b0_mse_plot_1.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>%
  mutate(MSE = sd^2 + (mean - truth)^2,
         CI_width = Q975 - Q025) %>%
  select(type, scenario, iter, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(n1), mean, ymin = lb, ymax = ub, group = n2,
             col = as.factor(n2))) +
  geom_pointrange(position = position_dodge(width = 0.3), show.legend = F) +
  geom_line(position = position_dodge(width = 0.3), show.legend = F) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Improvement in MSE beta0 estimate\ndue to data integration") + 
  ggtitle("(A) Improvement in MSE")

# Plot: improvement to MSE of Beta1
b1_mse_plot_1.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta1") %>%
  mutate(MSE = sd^2 + (mean - truth)^2,
         CI_width = Q975 - Q025) %>%
  select(type, scenario, iter, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(n1), mean, ymin = lb, ymax = ub, group = n2, col = as.factor(n2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Absolute improvement in\nMSE of estimate of beta1")


# Plot: Improvement in power to detect beta1. I don't like this plot as much
#         bc it's totally contingent on the value of beta1 that we pick
b1_power_plot_1.1 <- bind_rows(estimation_results) %>% 
  filter(param == "beta1") %>% 
  mutate(nonzero = sign(Q025) > 0) %>% 
  group_by(scenario, type) %>% 
  summarize(power = mean(nonzero)) %>% 
  select(type, scenario, power) %>% 
  pivot_wider(names_from = type, values_from = power) %>% 
  mutate(improvement = joint - one) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(n1), improvement, group = n2, col = as.factor(n2))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Improvement in power to detect\nnonzero effect of beta1") 


# Plot: Increasing computational cost as data size increases
runtime_plot_1.1 <- bind_rows(runtime_results) %>% 
  group_by(scenario, type) %>% 
  summarize(value = mean(vaule)) %>% 
  select(type, scenario, value) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(improvement = joint - one) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(n1), improvement, group = n2, col = as.factor(n2))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Increase in runtime (s)") + 
  ggtitle("(B) Computational cost")

# Save fig 1.1
p1.1 <- arrangeGrob(b0_mse_plot_1.1, runtime_plot_1.1, widths = c(1, 0.8), nrow = 1)
ggsave("plots/Fig1-1.jpg", p1.1, width = 12, height = 4, dpi = 600)


#### Simulation 1.2: Effect of data dimensions on usefulness of 2nd dataset ####

cl <- makeCluster(10)

specs_df <- as.data.frame(expand.grid(
  xi = c(0.1, 0.4, 0.7, 1),       # coverage
  alpha0_1 = c(-1, -1/3, 1/3, 1), # D1 information content: alpha0_1
  alpha0_2 = c(-1, 0, 1), # D2 information content: alpha0_2
  sigma = c(0.1, 1)
)) %>% 
  mutate(scenario = row_number())

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
             runtime_results = runtime_results),
        "sim_output/sim1_2.RDS")

stopCluster(cl)


#### Visualizing 1.2 ####


# Plot: improvement to precision of Beta0
bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  mutate(CI_width = Q975 - Q025) %>% 
  select(type, scenario, iter, CI_width) %>% 
  pivot_wider(names_from = type, values_from = CI_width) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob. of \nsecondary\ndataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0") +
  facet_grid(paste0("PA-2 noise:", sigma)~paste0("PA-1 coverage: ", xi))


# Plot: improvement to MSE of Beta0
bind_rows(estimation_results) %>% 
  filter(param == "beta0") %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob. of \nsecondary\ndataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0") +
  facet_grid(paste0("PA-2 noise:", sigma)~paste0("PA-1 coverage: ", xi))





