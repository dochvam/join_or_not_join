library(tidyverse)
library(gridExtra)

ests_list <- list()
rmse_list <- list()
cv_list <- list()

for (i in 1:5) {
  results_joint <- readRDS(paste0("output/scenario", i, "_joint_results.RDS"))
  results_cam <- readRDS(paste0("output/scenario", i, "_cameraOnly_results.RDS"))
  results_PO <- readRDS(paste0("output/scenario", i, "_POOnly_results.RDS"))
  
  ests_list[[i]] <- bind_rows(
    results_joint$estimates,
    results_cam$estimates,
    results_PO$estimates
  ) %>% 
    mutate(scenario = i)
  
  rmse_list[[i]] <- bind_rows(
    results_joint$rmse,
    results_cam$rmse,
    results_PO$rmse
  ) %>% 
    mutate(scenario = i)
  cv_list[[i]] <- bind_rows(
    results_joint$CV,
    results_cam$CV,
    results_PO$CV
  ) %>% 
    mutate(scenario = i)
}

ests_df <- bind_rows(ests_list) %>% 
  mutate(modtype = factor(recode(modtype, "PO_only" = "PO-only",
                                 "camera_only" = "PA-only",
                                 "joint" = "Joint"),
                          levels = c("PA-only", "PO-only", "Joint")))
rmse_df <- bind_rows(rmse_list)
cv_df <- bind_rows(cv_list)

#### Compare intensity parameter estimates ####

### Summarize all info
# Coverage plot
coverage_plot <- ests_df %>% 
  filter(grepl("intensity_beta", param), !is.na(se)) %>% 
  filter(!is.na(covered)) %>% 
  ggplot() + 
  geom_bar(aes(modtype, fill = covered), position = "fill",
           show.legend = F) + 
  coord_flip() +
  geom_hline(yintercept = 0.95) +
  theme_bw() +
  facet_wrap(~paste0("Coverage - Scenario ", scenario), ncol = 1) +
  xlab("") + ylab("Pct. of intensity covar effs.\nw/ 95% CIs covering the true val.") +
  scale_fill_manual("Coverage", values = c("darkgrey", "darkred"))

# Error plot
error_plot <- ests_df %>% 
  filter(grepl("intensity_beta", param)) %>% 
  filter(!is.na(est)) %>% 
  mutate(error = est - true_value) %>% 
  ggplot() + 
  geom_boxplot(aes(modtype, error)) + 
  ylim(c(-7.5, 7.5)) +
  facet_wrap(~paste0("Error - Scenario ", scenario), ncol = 1) +
  coord_flip() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  xlab("") + ylab("Abs. error in estimated intensity covar effects")

precision_plot <- ests_df %>% 
  filter(grepl("intensity_beta", param)) %>% 
  filter(!is.na(se)) %>% 
  ggplot() + 
  geom_boxplot(aes(modtype, se)) + 
  coord_flip() + ylim(c(0, 2.5)) +
  facet_wrap(~paste0("Precision - Scenario ", scenario), ncol = 1) +
  theme_bw() +
  xlab("") + ylab("SE of intensity param. estimate")

all_plots <- arrangeGrob(coverage_plot, error_plot, precision_plot, nrow = 1)
ggsave("plots/inference.jpg", all_plots, width = 10, height = 8)

# Estimates of theta1, by scenario
# stop("Modify below: color points to indicate whether they're nonzero")
theta1_plot <- ests_df %>% 
  filter(grepl("theta1", param), !is.na(se)) %>% 
  mutate(nonzero = sign(est + 1.96*se) == sign(est - 1.96*se)) %>% 
  ggplot() + 
  geom_hline(yintercept = 1) +
  geom_boxplot(aes(as.factor(scenario), est), outliers = F) + 
  geom_point(aes(as.factor(scenario), est, col = nonzero),
             position = position_jitter(width = 0.2),
             alpha = 0.6) + 
  coord_flip() +
  scale_color_manual(values = c("darkgray", "darkred")) +
  theme_minimal() + 
  ylab("Estimate of theta1") + xlab("Scenario")

ggsave("plots/theta1.jpg", theta1_plot, width = 5, height = 3)

#### Compare RMSE ####

rmse_plot <- rmse_df %>% 
    select(dataset_ID, modtype, scenario, RMSE) %>% 
    pivot_wider(names_from = modtype, values_from = RMSE) %>% 
    mutate(improvement = camera_only - joint) %>% 
    ggplot() + 
    geom_hline(yintercept = 0) +
    theme_minimal() + coord_flip() +
    ylim(c(-2, 5)) +
    xlab("Simulation scenario") + 
    ylab("Improvement in RMSE (joint vs. PA-only)") +
    geom_boxplot(aes(as.factor(scenario), group = as.factor(scenario), improvement))

ggsave("plots/prediction.jpg", rmse_plot, width = 4, height = 3)
#### Compare cross-validation ####

# Camera-based Brier score
brier_score_plot <- cv_df %>% 
  filter(score == "brier", type == "camera") %>% 
  select(dataset_ID, modtype, scenario, value) %>% 
  pivot_wider(names_from = modtype, values_from = value) %>% 
  mutate(improvement = camera_only - joint) %>% 
  ggplot() + 
  coord_flip() +
  geom_hline(yintercept = 0) +
  xlab("Simulation scenario") + 
  ylab("Change in OOS CV Brier score (PA-only - joint)") +
  geom_boxplot(aes(as.factor(scenario), group = as.factor(scenario), improvement)) +
  theme_minimal()

# Camera-based Brier score
brier_score_plot <- cv_df %>% 
  filter(score == "brier", type == "camera") %>% 
  select(dataset_ID, modtype, scenario, value) %>% 
  pivot_wider(names_from = modtype, values_from = value) %>% 
  mutate(improvement = camera_only - joint) %>% 
  ggplot() + 
  coord_flip() +
  geom_hline(yintercept = 0) +
  xlab("Simulation scenario") + 
  ylab("Change in OOS CV Brier score (PA-only - joint)") +
  geom_boxplot(aes(as.factor(scenario), group = as.factor(scenario), improvement)) +
  theme_minimal()

ggsave("plots/oos_cv.jpg", brier_score_plot, width = 4, height = 3)


