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
this_res <- readRDS("sim_output/sim2_1.RDS")
estimation_results <- this_res$estimation_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df

summary2.1 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(median = median(improvement), ub = quantile(improvement, 0.75),
            lb = quantile(improvement, 0.25)) %>%
  left_join(specs_df, by = "scenario")

runtime2.1 <- bind_rows(runtime_results) %>% 
  select(type, scenario, iter, value) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(improvement = joint - one) %>% 
  group_by(scenario) %>% 
  summarize(median = median(improvement), ub = quantile(improvement, 0.75),
            lb = quantile(improvement, 0.25)) %>% 
  left_join(specs_df, by = "scenario")


this_res <- readRDS("sim_output/sim2_2.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df

summary2.2 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(median = median(improvement), ub = quantile(improvement, 0.75),
            lb = quantile(improvement, 0.25)) %>%
  left_join(specs_df, by = "scenario")

this_res <- readRDS("sim_output/sim2_3.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df


summary2.3 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(median = median(improvement), ub = quantile(improvement, 0.75),
            lb = quantile(improvement, 0.25)) %>%
  left_join(specs_df, by = "scenario")


this_res <- readRDS("sim_output/sim2_4.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df

summary2.4 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(rate = mean(improvement > 0)) %>%
  left_join(specs_df, by = "scenario")



this_res <- readRDS("sim_output/sim2_5.RDS")
estimation_results <- this_res$estimation_results
specs_df <- this_res$specs_df

summary2.5 <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario) %>% 
  summarize(rate = mean(improvement > 0)) %>%
  left_join(specs_df, by = "scenario")

# 
# this_res <- readRDS("sim_output/sim2_6.RDS")
# estimation_results <- this_res$estimation_results
# specs_df <- this_res$specs_df
# 
# summary2.6 <- estimation_results %>% 
#   bind_rows() %>% 
#   filter(param %in% c("beta0", "beta1")) %>% 
#   mutate(MSE = sd^2 + (mean - truth)^2) %>% 
#   select(type, scenario, iter, param, MSE) %>% 
#   pivot_wider(names_from = type, values_from = MSE) %>% 
#   mutate(improvement = one - joint) %>% 
#   group_by(scenario) %>% 
#   summarize(mean = mean(improvement), ub = quantile(improvement, 0.75),
#             lb = quantile(improvement, 0.25)) %>%
#   left_join(specs_df, by = "scenario")




#### Sim Fig A - num sites ####
#' x-axis: D2 size
#' y-axis: MSE improvement VALUE
#' color: D1 size

panelA <- ggplot(summary2.1,
                 aes(nPA, median, ymax = ub, ymin = lb, col = as.factor(S),
                     group = as.factor(S))) +
  geom_pointrange(position = position_dodge(width = 5)) +
  geom_line(position = position_dodge(width = 5)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nsecondary\ndataset", begin = 0, end = 0.7) +
  xlab("Num. sites in primary dataset") +
  ylab("Improvement in MSE\nof occupancy intercept") +
  ggtitle("(A) Data dimensions")


#### Sim Fig B - computational cost ####

panelB <- ggplot(runtime2.1,
                 aes(nPA, median, ymax = ub, ymin = lb, col = as.factor(S),
                     group = as.factor(S))) +
  geom_pointrange(position = position_dodge(width = 5)) +
  geom_line(position = position_dodge(width = 5)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nsecondary\ndataset", begin = 0, end = 0.7) +
  xlab("Num. sites in primary dataset") +
  ylab("Increase in runtime (s)") +
  ggtitle("(B) Computational cost")

#### Sim Fig C - information content ####

panelC <- ggplot(summary2.2,
                 aes(alpha0, median, ymax = ub, ymin = lb, col = as.factor(theta0),
                     group = as.factor(theta0))) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  geom_line(position = position_dodge(width = 0.1)) +
  theme_minimal() +
  scale_color_viridis_d("Logit-scale\ndetection\nprobability in \nsecondary\ndataset", begin = 0, end = 0.7) +
  xlab("Logit-scale detection probability in primary dataset") +
  ylab("Improvement in MSE\nof occupancy intercept") +
  ggtitle("(C) Detection rates")


#### Sim Fig D - covariate coverage ####

panelD <- summary2.3 %>% 
  mutate(param = recode(param, "beta0" = "Intercept", "beta1" = "Effect of x")) %>% 
  ggplot(aes(xi, median, ymax = ub, ymin = lb, col = as.factor(param),
             group = as.factor(param))) +
  geom_pointrange(position = position_dodge(width = 0)) +
  geom_line(position = position_dodge(width = 0)) +
  theme_minimal() +
  scale_color_viridis_d("Parameter", 
                        begin = 0.33, end = 0.66, option = "mako") +
  xlab("Pct. of domain of covariate x covered by primary dataset") +
  ylab("Improvement in MSE\nof occupancy intercept") +
  ggtitle("(D) Covariate coverage")


#### Sim Fig E - noise ####
panelE <- ggplot(summary2.4,
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

#### Sim Fig F - bias ####
panelF <- ggplot(summary2.5,
                 aes(zeta, rate, col = as.factor(nPA),
                     group = as.factor(nPA))) +
  geom_point(position = position_dodge(width = 0)) +
  geom_hline(yintercept = 0.5) +
  geom_line(position = position_dodge(width = 0)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites\nin primary\ndataset", begin = 0, end = 0.7,
                        option = "plasma") +
  xlab("Logit-scale bias in the effect of x on secondary dataset") +
  ylab("Rate at which joint model outperforms\nsingle-dataset model, by MSE") +
  ggtitle("(F) Bias in secondary data")

#### Put the panels together ####
sim_results_fig <- arrangeGrob(panelA, panelB, panelC, panelD, panelE, panelF,
                               nrow = 3)
ggsave("sim2_results_suppl.jpg", sim_results_fig, width = 9, height = 10, dpi = 300)

