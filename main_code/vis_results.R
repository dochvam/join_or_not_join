library(tidyverse)
library(gridExtra)

#### Visualizing 1.1 ####

this_res <- readRDS("sim_output/sim1_1.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df

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
  theme_bw() +
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
  theme_bw() +
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
  theme_bw() +
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
  theme_bw() +
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
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Increase in runtime (s)") + 
  ggtitle("(B) Computational cost")

# Save fig 1.1
p1.1 <- arrangeGrob(b0_mse_plot_1.1, runtime_plot_1.1, widths = c(1, 0.8), nrow = 1)
ggsave("plots/Fig1-1.jpg", p1.1, width = 12, height = 4, dpi = 600)


#### Visualizing 1.2 ####

this_res <- readRDS("sim_output/sim1_2.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df

# Plot: improvement to MSE
plot_1.2A <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, MSE, param) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(sigma == 0.1) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob.\nof secondary\ndataset") +
  ylab("Absolute improvement in\nMSE of parameter estimate") +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi)) +
  ggtitle("A. Improvement to MSE due to joint modeling")


# Plot: improvement to abs. error
plot_1.2B <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(error = abs(mean - truth)) %>% 
  select(type, scenario, iter, error, param) %>% 
  pivot_wider(names_from = type, values_from = error) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(sigma == 0.1) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob.\nof secondary\ndataset") +
  ylab("Absolute improvement in\nabs. error of parameter estimate") +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi))

# Plot: improvement to OOS CV
plot_1.2C <- bind_rows(cv_results) %>% 
  pivot_wider(names_from = type, values_from = brier_score) %>% 
  mutate(improvement = 100 * (one - joint) / one) %>% 
  group_by(scenario) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(sigma == 0.1) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob. of \nsecondary\ndataset") +
  ylab("Percent improvement in\nOOS cross validation Brier score") +
  facet_grid(~paste0("PA-1 coverage: ", xi))


# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_1.2D <- bind_rows(
  # Metric 1: OOS CV
    bind_rows(cv_results) %>% 
      pivot_wider(names_from = type, values_from = brier_score) %>% 
      mutate(improvement = 100 * (one - joint) / one) %>% 
      group_by(scenario) %>% 
      summarize(mean = mean(improvement > 0)) %>% 
      mutate(metric = "OOS CV"),
  # Metric 2: MSE
    bind_rows(estimation_results) %>% 
      filter(param %in% c("beta0", "beta1")) %>% 
      mutate(MSE = sd^2 + (mean - truth)^2) %>% 
      select(type, scenario, iter, MSE, param) %>% 
      pivot_wider(names_from = type, values_from = MSE) %>% 
      mutate(improvement = one - joint) %>% 
      group_by(scenario) %>% 
      summarize(mean = mean(improvement > 0)) %>% 
      mutate(metric = "MSE"), 
  # Metric 3: error
    bind_rows(estimation_results) %>% 
      filter(param %in% c("beta0", "beta1")) %>% 
      mutate(error = abs(mean - truth)) %>% 
      select(type, scenario, iter, error, param) %>% 
      pivot_wider(names_from = type, values_from = error) %>% 
      mutate(improvement = one - joint) %>% 
      group_by(scenario) %>% 
      summarize(mean = mean(improvement > 0)) %>% 
      mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = Q975 - Q025) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
  ) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(sigma == 0.1) %>% 
  ggplot(aes(alpha0_1, mean, group = alpha0_2, col = as.factor(alpha0_2))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob.\nof secondary\ndataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

# p1.2 <- arrangeGrob(plot_1.2A, plot_1.2D, nrow = 2)
ggsave("plots/Fig1-2.jpg", plot_1.2D, width = 8, height = 8, dpi = 600)

#### Visualizing 1.3 ####

this_res <- readRDS("sim_output/sim1_3.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


# Improvement in MSE
plot_1.3A <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(n1 = factor(n1, levels = unique(specs_df$n1))) %>% 
  ggplot(aes(sigma, mean, ymin = lb, ymax = ub, group = n1, col = n1)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Num. sites in\nprimary dataset") +
  ylab("Improvement in MSE of estimate")

# Improvement in point error
plot_1.3B <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(error = abs(mean - truth)) %>% 
  select(type, scenario, iter, param, error) %>% 
  pivot_wider(names_from = type, values_from = error) %>% 
  mutate(improvement = (one - joint)) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(n1 = factor(n1, levels = unique(specs_df$n1))) %>% 
  ggplot(aes(sigma, mean, ymin = lb, ymax = ub, group = n1, col = n1)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Num. sites in\nprimary dataset") +
  ylab("Improvement in precision of estimate")

# Improvement in precision
plot_1.3C <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(precision = Q975 - Q025) %>% 
  select(type, scenario, iter, param, precision) %>% 
  pivot_wider(names_from = type, values_from = precision) %>% 
  mutate(improvement = (one - joint)/one) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(n1 = factor(n1, levels = unique(specs_df$n1))) %>% 
  ggplot(aes(sigma, mean, ymin = lb, ymax = ub, group = n1, col = n1)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Num. sites in\nprimary dataset") +
  ylab("Improvement in precision of estimate")


# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_1.3D <- bind_rows(
  # Metric 1: OOS CV
  bind_rows(cv_results) %>% 
    pivot_wider(names_from = type, values_from = brier_score) %>% 
    mutate(improvement = 100 * (one - joint) / one) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "OOS CV"),
  # Metric 2: MSE
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = sd^2 + (mean - truth)^2) %>% 
    select(type, scenario, iter, MSE, param) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "MSE"), 
  # Metric 3: error
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(error = abs(mean - truth)) %>% 
    select(type, scenario, iter, error, param) %>% 
    pivot_wider(names_from = type, values_from = error) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = Q975 - Q025) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(n1 = factor(n1, levels = unique(specs_df$n1))) %>% 
  ggplot(aes(sigma, mean, group = n1, col = n1)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Num. sites in\nprimary dataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

ggsave("plots/Fig1-3.jpg", plot_1.3D, width = 8, height = 8, dpi = 600)

#### Visualizing 1.4 ####

this_res <- readRDS("sim_output/sim1_4.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


# Improvement in MSE
plot_1.4A <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(zeta = factor(zeta, levels = unique(specs_df$zeta))) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = zeta, col = zeta)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Magnitude of bias in\nprimary dataset") +
  ylab("Improvement in MSE of estimate")

# Improvement in point error
plot_1.4B <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(error = abs(mean - truth)) %>% 
  select(type, scenario, iter, param, error) %>% 
  pivot_wider(names_from = type, values_from = error) %>% 
  mutate(improvement = (one - joint)) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(zeta = factor(zeta, levels = unique(specs_df$zeta))) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = zeta, col = zeta)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Magnitude of bias in\nprimary dataset") +
  ylab("Improvement in precision of estimate")

# Improvement in precision
plot_1.4C <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  mutate(precision = Q975 - Q025) %>% 
  select(type, scenario, iter, param, precision) %>% 
  pivot_wider(names_from = type, values_from = precision) %>% 
  mutate(improvement = (one - joint)/one) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(zeta = factor(zeta, levels = unique(specs_df$zeta))) %>% 
  ggplot(aes(alpha0_1, mean, ymin = lb, ymax = ub, group = zeta, col = zeta)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Degree of noise in secondary dataset") +
  scale_color_discrete("Magnitude of bias in\nprimary dataset") +
  ylab("Improvement in precision of estimate")


# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_1.4D <- bind_rows(
  # Metric 1: OOS CV
  bind_rows(cv_results) %>% 
    pivot_wider(names_from = type, values_from = brier_score) %>% 
    mutate(improvement = 100 * (one - joint) / one) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "OOS CV"),
  # Metric 2: MSE
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = sd^2 + (mean - truth)^2) %>% 
    select(type, scenario, iter, MSE, param) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "MSE"), 
  # Metric 3: error
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(error = abs(mean - truth)) %>% 
    select(type, scenario, iter, error, param) %>% 
    pivot_wider(names_from = type, values_from = error) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = Q975 - Q025) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(zeta = factor(zeta, levels = unique(specs_df$zeta))) %>% 
  ggplot(aes(alpha0_1, mean, group = zeta, col = zeta)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Magnitude of bias in\nprimary dataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

ggsave("plots/Fig1-4.jpg", plot_1.4D, width = 8, height = 8, dpi = 600)



#### Visualizing 2.1 ####

this_res <- readRDS("sim_output/sim2_1.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


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
  theme_bw() +
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
  theme_bw() +
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
  theme_bw() +
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
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("Num. sites in primary dataset") +
  scale_color_discrete("Num. sites in\nsecondary\ndataset") +
  ylab("Increase in runtime (s)") + 
  ggtitle("(B) Computational cost")


# Save fig 2.1
p2.1 <- arrangeGrob(b0_mse_plot_2.1, runtime_plot_2.1, widths = c(1, 0.8), nrow = 1)
ggsave("plots/Fig2-1.jpg", p2.1, width = 12, height = 4, dpi = 600)


#### Visualizing 2.2 ####


this_res <- readRDS("sim_output/sim2_2.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


### 2.2A -- changing alpha0 and theta0, effects on MSE

# Plot: improvement to MSE of Beta0 and Beta1 with alpha0, theta0
plot_2.2A <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  filter(scenario %in% 1:48) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(alpha0, mean, ymin = lb, ymax = ub, group = theta0, col = as.factor(theta0))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Avg. logit-scale detection prob. in primary dataset") +
  scale_color_discrete("Avg. detection rate\nin secondary dataset") +
  ylab("Improvement in MSE of estimate")


# Plot: improvement to MSE of Beta0 and Beta1 with J, eta
plot_2.2B <- bind_rows(estimation_results) %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  filter(scenario %in% 49:96) %>% 
  mutate(MSE = sd^2 + (mean - truth)^2) %>% 
  select(type, scenario, iter, param, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  group_by(scenario, param) %>% 
  summarize(mean = mean(improvement), ub = quantile(improvement, 0.975),
            lb = quantile(improvement, 0.025)) %>% 
  left_join(specs_df, by = "scenario") %>% 
  ggplot(aes(as.factor(J), mean, ymin = lb, ymax = ub, group = eta, col = as.factor(eta))) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  facet_grid(paste0("Parameter: ", param)~paste0("PA-1 coverage: ", xi),
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab("Num. replicate visits in primary dataset") +
  scale_color_discrete("Avg. effort per grid cell\nin secondary dataset") +
  ylab("Improvement in MSE of estimate")


# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_2.2A <- bind_rows(
  # Metric 1: OOS CV
  bind_rows(cv_results) %>% 
    pivot_wider(names_from = type, values_from = brier_score) %>% 
    mutate(improvement = 100 * (one - joint) / one) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "OOS CV"),
  # Metric 2: MSE
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = sd^2 + (mean - truth)^2) %>% 
    select(type, scenario, iter, MSE, param) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "MSE"), 
  # Metric 3: error
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(error = abs(mean - truth)) %>% 
    select(type, scenario, iter, error, param) %>% 
    pivot_wider(names_from = type, values_from = error) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = `97.5%` - `2.5%`) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(scenario %in% 1:48) %>% 
  ggplot(aes(alpha0, mean, group = theta0, col = as.factor(theta0))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Mean det. prob.\nof secondary\ndataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

# p1.2 <- arrangeGrob(plot_1.2A, plot_1.2D, nrow = 2)
ggsave("plots/Fig2-2A.jpg", plot_2.2A, width = 8, height = 8, dpi = 600)

# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_2.2B <- bind_rows(
  # Metric 1: OOS CV
  bind_rows(cv_results) %>% 
    pivot_wider(names_from = type, values_from = brier_score) %>% 
    mutate(improvement = 100 * (one - joint) / one) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "OOS CV"),
  # Metric 2: MSE
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = sd^2 + (mean - truth)^2) %>% 
    select(type, scenario, iter, MSE, param) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "MSE"), 
  # Metric 3: error
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(error = abs(mean - truth)) %>% 
    select(type, scenario, iter, error, param) %>% 
    pivot_wider(names_from = type, values_from = error) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = `97.5%` - `2.5%`) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  filter(scenario %in% 49:96) %>% 
  ggplot((aes(as.factor(J), mean, group = eta, col = as.factor(eta)))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Num. replicates per site in primary dataset") +
  scale_color_discrete("Mean per-cell effort\nin secondary\ndataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

# p1.2 <- arrangeGrob(plot_1.2A, plot_1.2D, nrow = 2)
ggsave("plots/Fig2-2B.jpg", plot_2.2B, width = 8, height = 8, dpi = 600)


#### Visualizing 2.3 ####


this_res <- readRDS("sim_output/sim2_3.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
runtime_results <- this_res$runtime_results
specs_df <- this_res$specs_df


# Plot: Pct. of cases where joint outperforms single by all 3 metrics
plot_2.3D <- bind_rows(
  # Metric 1: OOS CV
  bind_rows(cv_results) %>% 
    pivot_wider(names_from = type, values_from = brier_score) %>% 
    mutate(improvement = 100 * (one - joint) / one) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "OOS CV"),
  # Metric 2: MSE
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = sd^2 + (mean - truth)^2) %>% 
    select(type, scenario, iter, MSE, param) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "MSE"), 
  # Metric 3: error
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(error = abs(mean - truth)) %>% 
    select(type, scenario, iter, error, param) %>% 
    pivot_wider(names_from = type, values_from = error) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Abs. error"),
  # Metric 4: precision
  bind_rows(estimation_results) %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(precision = `97.5%` - `2.5%`) %>% 
    select(type, scenario, iter, precision, param) %>% 
    pivot_wider(names_from = type, values_from = precision) %>% 
    mutate(improvement = one - joint) %>% 
    group_by(scenario) %>% 
    summarize(mean = mean(improvement > 0)) %>% 
    mutate(metric = "Precision")
) %>% 
  mutate(metric = factor(metric, levels = c("Abs. error", "Precision",
                                            "MSE", "OOS CV"))) %>% 
  left_join(specs_df, by = "scenario") %>% 
  mutate(nPA = factor(nPA, levels = unique(nPA))) %>% 
  ggplot(aes(sigma, mean, group = nPA, col = nPA)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  xlab("Mean det. prob. of primary dataset") +
  scale_color_discrete("Num. sites in\nprimary dataset") +
  ylab("Percent of cases where joint\noutperforms single-dataset model") +
  facet_grid(metric~paste0("PA-1 coverage: ", xi)) +
  ggtitle("Rate at which joint model outperforms PA-1-only model")

ggsave("plots/Fig2-3.jpg", plot_2.3D, width = 8, height = 8, dpi = 600)

