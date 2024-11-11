library(tidyverse)
library(gridExtra)

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


# Plot coverage
plot_one("PA coverage", "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)
plot_one("PA effort",   "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)
plot_one("PO effort",   "coverage", ests_df, condition_df, as_diff = F, hline_yint = 0.95)

# Plot precision
plot_one("PA coverage", "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PA effort",   "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PO effort",   "precision", ests_df, condition_df, as_diff = T, hline_yint = 0)

# Plot error
plot_one("PA coverage", "error", ests_df, condition_df, as_diff = F, hline_yint = 0)
plot_one("PA effort",   "error", ests_df, condition_df, as_diff = F, hline_yint = 0)
plot_one("PO effort",   "error", ests_df, condition_df, as_diff = F, hline_yint = 0)

# Plot CV
plot_one("PA coverage", "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PA effort",   "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PO effort",   "cv", cv_df, condition_df, as_diff = T, hline_yint = 0)


# Plot CV
plot_one("PA coverage", "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PA effort",   "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)
plot_one("PO effort",   "rmse", rmse_df, condition_df, as_diff = T, hline_yint = 0)


