library(spOccupancy)
library(tidyverse)
library(parallel)
source("support_fn/sim2_fn.R")

cl <- makeCluster(10)

#### Illustration 1: decreasing marginal benefits of 2nd dataset ####

specs_df <- as.data.frame(expand.grid(
  n1 = c(40, 80, 160, 320, 640), eta = c(250, 500, 1000)
)) %>% 
  mutate(scenario = row_number())
nsim <- 100

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
  ggplot(aes(as.factor(n1), mean, ymin = lb, ymax = ub, group = n2, col = as.factor(n2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Absolute improvement in\nprecision of estimate of beta0")



# Plot: Improvement in power to detect beta1. I don't like this plot as much
#         bc it's totally contingent on the value of beta1 that we pick
bind_rows(estimation_results) %>% 
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
bind_rows(runtime_results) %>% 
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
  ylab("Increase in runtime (s)")

