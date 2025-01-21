library(tidyverse)
library(gridExtra)
library(MuMIn)

ests_list <- list()
rmse_list <- list()
cv_list <- list()

condition_df <- read_csv("output_v2/conditions.csv")
all_results_joint <- readRDS("output_v2/all_results_joint.RDS")
all_results_camera <- readRDS("output_v2/all_results_camera.RDS")


for (i in 1:length(all_results_joint)) {
  ests_list[[i]] <- bind_rows(
    all_results_joint[[i]]$estimates,
    all_results_camera[[i]]$estimates,
  )
  rmse_list[[i]] <- bind_rows(
    all_results_joint[[i]]$rmse,
    all_results_camera[[i]]$rmse
  )
  cv_list[[i]] <- bind_rows(
    all_results_joint[[i]]$CV,
    all_results_camera[[i]]$CV
  )
}

ests_df <- bind_rows(ests_list) %>% 
  mutate(modtype = factor(recode(modtype, "PO_only" = "PO-only",
                                 "camera_only" = "PA-only",
                                 "joint" = "Joint"),
                          levels = c("PA-only", "PO-only", "Joint")))
rmse_df <- bind_rows(rmse_list) %>% 
  mutate(modtype = factor(recode(modtype, "PO_only" = "PO-only",
                                 "camera_only" = "PA-only",
                                 "joint" = "Joint"),
                          levels = c("PA-only", "PO-only", "Joint")))
cv_df <- bind_rows(cv_list) %>% 
  mutate(modtype = factor(recode(modtype, "PO_only" = "PO-only",
                                 "camera_only" = "PA-only",
                                 "joint" = "Joint"),
                          levels = c("PA-only", "PO-only", "Joint")))



#### Compare intensity parameter estimates ####

plot_one <- function(target, outcome, result_df, condition_df, 
                     as_diff = TRUE, hline_yint = NULL, incl_CI = FALSE) {
  
  outcome_codes <- c("coverage", "error", "precision", "rmse", "cv", "theta1")
  outcome_display_names <- c("95% CI coverage", "Error", "Precision",
                             "RMSE (prediction)", "PA cross validation score", 
                             "Join parameter theta1")
  
  which_outcome <- which(outcome_codes == outcome)
  
  stopifnot(outcome %in% outcome_codes)
  
  if (target == "PA effort") {
    groupcol <- "ncam"
    groupname <- "Num. PA sites"
  } else if (target == "PA coverage") {
    groupcol <- "cam_pct_xcover"
    groupname <- "PA pct. coverage of covar. space"
  } else if (target == "PO effort") {
    groupcol <- "PO_avg_effort"
    groupname <- "Avg. PO obs. per cell"
  } else {
    stop(paste0('Invalid target: "', target, '"'))
  }
  
  all_df <- left_join(result_df, condition_df, by = c("scenario")) %>% 
    rename(target_col = target) %>% 
    filter(target_col == target)
  
  if (outcome == "coverage") {
    
    plot_df <- all_df %>% 
      filter(param == "intensity_beta") %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = mean(covered, na.rm = TRUE), 
                CI_LB = mean(covered, na.rm = TRUE),
                CI_UB = mean(covered, na.rm = TRUE),
                .groups = "drop"
                )
    
  } else if (outcome == "precision") {
    
    plot_df <- all_df %>% 
      filter(param == "intensity_beta") %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = median(se, na.rm = TRUE), .groups = "drop")
    
  } else if (outcome == "error") {
    
    plot_df <- all_df %>% 
      filter(param == "intensity_beta") %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = median(est - true_value, na.rm = TRUE), .groups = "drop")
    
  } else if (outcome == "cv") {
    
    plot_df <- all_df %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = median(value, na.rm = TRUE), .groups = "drop")
    
  } else if (outcome == "rmse") {
    
    plot_df <- all_df %>% 
      filter(goal == "prediction") %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = median(RMSE, na.rm = TRUE), .groups = "drop")
    
  } else if (outcome == "theta1") {
    
    plot_df <- all_df %>% 
      filter(param == "theta1") %>% 
      group_by_at(c("scenario", "modtype")) %>% 
      summarize(outcome = median(est, na.rm = TRUE), .groups = "drop")
    
  }
  
  if (as_diff) {
    compare_df <- plot_df %>%
      pivot_wider(names_from = modtype, values_from = outcome) %>% 
      mutate(outcome_diff = `Joint` - `PA-only`) %>% 
      left_join(condition_df, by = "scenario") %>% 
      mutate(bias_name = factor(paste0("PO bias = ", PO_bias),
                                levels = paste0("PO bias = ", unique(sort(condition_df$PO_bias)))))
    colnames(compare_df)[colnames(compare_df) == groupcol] <- "plotgroup"
    
    p <- ggplot(compare_df) +
      geom_line(aes(PO_noise, outcome_diff, group = plotgroup, col = as.factor(plotgroup))) +
      facet_grid(~bias_name) +
      xlab("PO noise") +
      theme_minimal() +
      ylab(paste0(outcome_display_names[which_outcome], " diff (joint - PA-only)"))
  } else {
    plot_df <- plot_df %>% left_join(condition_df, by = "scenario") %>% 
      mutate(bias_name = factor(paste0("PO bias = ", PO_bias),
                                levels = paste0("PO bias = ", unique(sort(condition_df$PO_bias)))))
    
    colnames(plot_df)[colnames(plot_df) == groupcol] <- "plotgroup"
    
    p <- ggplot(plot_df) +
      geom_line(aes(PO_noise, outcome, group = plotgroup, col = as.factor(plotgroup))) +
      facet_grid(modtype~bias_name) +
      xlab("PO noise") +
      ylab(outcome_display_names[which_outcome])
  }
  
  # Handle theming
  p <- p + 
    theme_bw() +
    scale_color_viridis_d(groupname, end = 0.75)
  
  # Add hline
  if (!is.null(hline_yint)) p <- p + geom_hline(yintercept = hline_yint)
  
  return(p)
}

# 
# # Plot coverage
# plot_one("PA coverage", "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)
# plot_one("PA effort",   "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)
# plot_one("PO effort",   "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)
# 
# # Plot precision
# plot_one("PA coverage", "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PA effort",   "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PO effort",   "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)
# 
# # Plot error
# plot_one("PA coverage", "error", ests_df, condition_df, as_diff = F, hline_yint = 0)
# plot_one("PA effort",   "error", ests_df, condition_df, as_diff = F, hline_yint = 0)
# plot_one("PO effort",   "error", ests_df, condition_df, as_diff = F, hline_yint = 0)
# 
# # Plot CV
# plot_one("PA coverage", "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PA effort",   "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PO effort",   "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)
# 
# 
# # Plot CV
# plot_one("PA coverage", "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PA effort",   "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)
# plot_one("PO effort",   "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)



#### What explains CV difference? ####

cv_diff <- cv_df %>% 
  pivot_wider(names_from = "modtype", values_from = "value") %>% 
  mutate(best_model = ifelse(Joint < `PA-only`, "Joint", "PA-only"),
         diff = Joint - `PA-only`,
         joint_preferred = as.numeric(best_model == "Joint"))%>% 
  left_join(condition_df) %>% 
  na.omit()

cv_diff_summary <- cv_diff %>% 
  group_by(scenario) %>% 
  summarize(
    joint_rate = mean(best_model == "Joint"),
    med_diff = median(diff),
    .groups = "drop"
  ) %>% 
  left_join(condition_df)


cv_plot <- cv_diff_summary %>% 
  ggplot() +
  geom_hline(yintercept = 0.5) +
  geom_point(aes(as.factor(PO_bias), joint_rate, col = PO_noise),
             position = position_jitter(width = 0.2)) +
  facet_grid(factor(paste0("PA sites: ", ncam), levels = paste0("PA sites: ", sort(unique(ncam)))) ~
               paste0("Pct coverage: ", cam_pct_xcover)) +
  theme_bw() +
  scale_color_viridis_c(end = 0.9) +
  xlab("Bias in PO beta") +
  ylab("Rate at which the joint model is better, by OOS CV") +
  ggtitle("Cross-validation error")

ggsave("plots/CV_plot_v2.jpg", cv_plot, width = 8, height = 6)

#### What explains precision difference? ####

ests_diff <- ests_df %>% 
  filter(param == "intensity_beta") %>% 
  select(dataset_ID, scenario, modtype, se) %>% 
  pivot_wider(names_from = "modtype", values_from = "se") %>% 
  mutate(best_model = ifelse(Joint < `PA-only`, "Joint", "PA-only"),
         diff = Joint - `PA-only`,
         joint_preferred = as.numeric(best_model == "Joint"))%>% 
  left_join(condition_df) %>% 
  na.omit()

ests_diff_summary <- ests_diff %>% 
  group_by(scenario) %>% 
  summarize(
    joint_rate = mean(best_model == "Joint"),
    med_diff = median(diff),
    .groups = "drop"
  ) %>% 
  left_join(condition_df)

# # Create the full formula with main effects and first-order interactions
# formula_full <- as.formula(joint_preferred ~ PO_bias + cam_pct_xcover + PO_noise + 
#                              ncam + PO_avg_effort + PO_bias:cam_pct_xcover +
#                              PO_bias:PO_noise + PO_bias:ncam + PO_bias:PO_avg_effort +
#                              cam_pct_xcover:PO_noise + PO_noise:ncam + PO_noise:PO_avg_effort +
#                              PO_bias:PO_noise:cam_pct_xcover)
# 
# # Fit the full logistic regression model
# full_model <- glm(formula_full, data = ests_diff, family = binomial,
#                   na.action = na.fail)
# 
# # Perform comprehensive model selection using dredge
# model_set <- dredge(full_model, rank = "AIC")
# 
# # Select the best model based on AIC
# best_model <- get.models(model_set, 1)[[1]]
# summary(best_model)

precision_plot <- ests_diff_summary %>% 
  ggplot() +
  geom_hline(yintercept = 0.5) +
  geom_point(aes(as.factor(PO_bias), joint_rate, col = PO_noise),
             position = position_jitter(width = 0.2)) +
  facet_grid(factor(paste0("PA sites: ", ncam), levels = paste0("PA sites: ", sort(unique(ncam)))) ~
               paste0("Pct coverage: ", cam_pct_xcover)) +
  theme_bw() +
  scale_color_viridis_c(end = 0.9) +
  xlab("Bias in PO beta") +
  ylab("Rate at which the joint model is more precise") +
  ggtitle("Precision")

ggsave("plots/prec_plot_v2.jpg", precision_plot, width = 8, height = 6)


# p1 <- ests_diff_summary %>% 
#   ggplot() +
#   geom_boxplot(aes(as.factor(cam_pct_xcover), joint_rate)) +
#   theme_bw() +
#   scale_color_viridis_c(end = 0.9) + 
#   ylab("Rate at which the joint model beta est. is more precise") + 
#   xlab("Pct. coverage of x domain by PA data")
# 
# p2 <- ests_diff_summary %>% 
#   ggplot() +
#   geom_boxplot(aes(as.factor(ncam), joint_rate)) +
#   theme_bw() +
#   scale_color_viridis_c(end = 0.9) + 
#   ylab("") + 
#   xlab("Num. PA sites")
# 
# p3 <- ests_diff_summary %>% 
#   ggplot() +
#   geom_boxplot(aes(as.factor(PO_avg_effort), joint_rate)) +
#   theme_bw() +
#   scale_color_viridis_c(end = 0.9) + 
#   ylab("") + 
#   xlab("Average num. PO obs. per grid")

# prec_plots <- gridExtra::arrangeGrob(p1, p2, p3, nrow = 1)
# plot(prec_plots)

#### What explains error difference? ####

err_diff <- ests_df %>% 
  mutate(abs_error = abs(true_value - est)) %>% 
  filter(param == "intensity_beta") %>% 
  select(dataset_ID, scenario, modtype, abs_error) %>% 
  pivot_wider(names_from = "modtype", values_from = "abs_error") %>% 
  mutate(best_model = ifelse(Joint < `PA-only`, "Joint", "PA-only"),
         diff = Joint - `PA-only`,
         joint_preferred = as.numeric(best_model == "Joint"))%>% 
  left_join(condition_df) %>% 
  na.omit()

err_diff_summary <- err_diff %>% 
  group_by(scenario) %>% 
  summarize(
    joint_rate = mean(best_model == "Joint"),
    med_diff = median(diff),
    .groups = "drop"
  ) %>% 
  left_join(condition_df)

# library(MuMIn)
# 
# # Create the full formula with main effects and first-order interactions
# formula_full <- as.formula(diff ~ PO_bias + cam_pct_xcover + PO_noise + 
#                              ncam + PO_avg_effort + PO_bias:cam_pct_xcover +
#                              PO_bias:PO_noise + PO_bias:ncam + PO_bias:PO_avg_effort +
#                              cam_pct_xcover:PO_noise + PO_noise:ncam + PO_noise:PO_avg_effort +
#                              PO_bias:PO_noise:cam_pct_xcover)
# 
# # Fit the full logistic regression model
# full_model <- glm(formula_full, data = cv_diff, family = gaussian,
#                   na.action = na.fail)
# 
# # Perform comprehensive model selection using dredge
# model_set <- dredge(full_model, rank = "AIC")
# 
# # Display the top models based on AIC
# print(model_set)
# 
# # Select the best model based on AIC
# best_model <- get.models(model_set, 1)[[1]]
# summary(best_model)



error_plot <- err_diff_summary %>% 
  ggplot() +
  geom_hline(yintercept = 0.5) +
  facet_grid(factor(paste0("PA sites: ", ncam), levels = paste0("PA sites: ", sort(unique(ncam)))) ~
               paste0("Pct coverage: ", cam_pct_xcover)) +
  
  geom_point(aes(as.factor(PO_bias), joint_rate, col = PO_noise),
             position = position_jitter(width = 0.2)) +
  theme_bw() +
  scale_color_viridis_c(end = 0.9) +
  xlab("Bias in PO beta") +
  ylab("Rate at which the joint model is more accurate") +
  ggtitle("Absolute error")

ggsave("plots/error_plot_v2.jpg", error_plot, width = 8, height = 6)


theta1_df <- ests_df %>% 
  filter(param == "theta1") %>% 
  mutate(joint_preferred = sign(est - 1.96*se) == sign(est + 1.96*se))

theta1_summary <- theta1_df %>% 
  group_by(scenario) %>% 
  summarize(
    joint_rate = mean(joint_preferred, na.rm = T),
    .groups = "drop"
  ) %>% 
  left_join(condition_df)

theta1_plot <- theta1_summary %>% 
  ggplot() +
  geom_hline(yintercept = 0.5) +
  facet_grid(paste0("Effect size: beta =", intensity_beta)~
    #factor(paste0("PA sites: ", ncam), levels = paste0("PA sites: ", sort(unique(ncam)))) ~
               paste0("Pct coverage: ", cam_pct_xcover)) +
  
  geom_point(aes(as.factor(PO_bias), joint_rate, col = PO_noise),
             position = position_jitter(width = 0.2)) +
  theme_bw() +
  scale_color_viridis_c(end = 0.9) +
  xlab("Bias in PO beta") +
  ylab("Rate at which theta1 is different from 0") +
  ggtitle("Theta1")

ggsave("plots/theta1_plot_v2.jpg", theta1_plot, width = 8, height = 6)

# 
# 
# err_diff_summary %>% 
#   ggplot() +
#   geom_point(aes(as.factor(ncam), joint_rate, col = PO_noise),
#              position = position_jitter(width = 0.2)) +
#   facet_wrap(~PO_bias) +
#   theme_bw() +
#   scale_color_viridis_c(end = 0.9) +
#   xlab("Pct. coverage of x domain by PA data") +
#   ylab("Rate at which the joint model is better, by OOS CV")
