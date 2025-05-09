library(tidyverse)
library(gridExtra)

#' Four panels in main figure:
#'  A: Data dimensions:
#'  B: Data dimensions:
#'  C: Steady computational cost
#'  D: Bias/noise effects
#'  
#' In general: use improvement value for cases where it's always good (marginal
#' benefit), use improvement rate when sometimes it's not good


#### Summarize ####
this_res <- readRDS("sim_output/sim3_1.RDS")
estimation_results <- this_res$estimation_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


summary3.1 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(median = median(improvement, na.rm = T), ub = quantile(improvement, 0.75, na.rm = T),
            lb = quantile(improvement, 0.25, na.rm = T)) %>%
  left_join(specs_df, by = "scenario")


runtime3.1 <- bind_rows(runtime_results) %>% 
  mutate(value = vaule) %>% 
  select(type, scenario, iter, value) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(improvement = joint - one) %>% 
  group_by(scenario) %>% 
  summarize(median = median(improvement, na.rm = T), ub = quantile(improvement, 0.75, na.rm = T),
            lb = quantile(improvement, 0.25, na.rm = T)) %>% 
  left_join(specs_df, by = "scenario")


datavol3.1 <- bind_rows(runtime_results) %>% 
  select(type, scenario, iter, num_points ) %>% 
  pivot_wider(names_from = type, values_from = num_points ) %>% 
  group_by(scenario) %>% 
  summarize(mean_joint = mean(joint, na.rm = T), mean_one = mean(one, na.rm = T),
            sd_joint = sd(joint, na.rm = T), sd_one = sd(one, na.rm = T)
            
            ) %>% 
  left_join(specs_df, by = "scenario")



this_res <- readRDS("sim_output/sim3_3.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df

summary3.3 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(median = median(improvement, na.rm = T), ub = quantile(improvement, 0.75, na.rm = T),
            lb = quantile(improvement, 0.25, na.rm = T)) %>%
  left_join(specs_df, by = "scenario")




this_res <- readRDS("sim_output/sim3_4.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df

summary3.4 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(rate = mean(improvement > 0, na.rm = T)) %>%
  left_join(specs_df, by = "scenario")


#### Vis

panelA <- ggplot(summary3.1,
                 aes(beta0_1, median, ymax = ub, ymin = lb, col = as.factor(beta0_2),
                     group = as.factor(beta0_2))) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nsecondary\ndataset", begin = 0, end = 0.7) +
  xlab("Num. sites in primary dataset") +
  ylab("Improvement in MSE\nof occupancy intercept") +
  ggtitle("(A) Data dimensions")


panelB <- ggplot(runtime3.1,
                 aes(beta0_1, median, ymax = ub, ymin = lb, col = as.factor(beta0_2),
                     group = as.factor(beta0_2))) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nsecondary\ndataset", begin = 0, end = 0.7) +
  xlab("Num. sites in primary dataset") +
  ylab("Increase in runtime (s)") +
  ggtitle("(B) Computational cost")




panelE <- ggplot(summary3.3,
                 aes(sigma, rate, col = as.factor(nPA),
                     group = as.factor(nPA))) +
  geom_point(position = position_dodge(width = 0)) +
  geom_line(position = position_dodge(width = 0)) +
  geom_hline(yintercept = 0.5) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites\nin primary\ndataset", begin = 0, end = 0.7,
                        option = "plasma") +
  xlab("Logit-scale random noise in secondary dataset") +
  ylab("Rate at which joint model outperforms\nsingle-dataset model, by MSE") +
  ggtitle("(E) Noise in secondary data")
