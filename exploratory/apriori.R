library(tidyverse)

test_corr <- function(this_dat) {
  cam_df <- this_dat$camera_df %>% 
    mutate(cam_ID = row_number()) %>% 
    select(cam_ID, PO_agg_x, PO_agg_y, all_of(this_dat$visit_cols)) %>%
    pivot_longer(cols = this_dat$visit_cols) %>% 
    group_by(cam_ID, PO_agg_x, PO_agg_y) %>% 
    summarize(n = n(), obs = sum(value) > 0, .groups = "drop")
  compare_df <- left_join(
    cam_df, this_dat$po_df[, c("PO_agg_x", "PO_agg_y", "effort", "y")], 
    by = join_by(PO_agg_x, PO_agg_y)
  ) %>% 
    mutate(PO_rate = y / effort) %>% 
    filter(!is.na(PO_rate), effort > 0)
  
  result <- NA
  fit <- NA
  catch <- tryCatch({
    fit <- glm(obs ~ PO_rate, family = binomial, data = compare_df)

    s <- summary(fit)
    
    result <- data.frame(
      dataset_ID = this_dat$dataset_ID,
      slope = s$coefficients[2, 1], 
      slope_pval = s$coefficients[2, 4],
      cor = cor(compare_df$obs, compare_df$PO_rate),
      deviance = fit$deviance,
      null_deviance = fit$null.deviance
    )
  }, error = function(e) { })
  
  if (!is.data.frame(result) && is.na(result)) {
    return(data.frame(
      dataset_ID = this_dat$dataset_ID,
      slope = NA, 
      slope_pval = NA,
      cor = NA,
      deviance = NA, null_deviance = NA
    ))
  }
  
  return(result)
}

slopes_df_list <- list()
for (i in 1:5) {
  sim_data <- readRDS(paste0("output/scenario", i, "_datasets.RDS"))
  
  slopes <- lapply(sim_data, test_corr)
  # for (d in 1:length(sim_data)) {
  #   slopes[[d]] <- test_corr(this_dat = sim_data[[d]])
  # }
  
  slopes_df_list[[i]] <- bind_rows(slopes) %>% 
    mutate(scenario = i)
}

slopes_df_list %>% 
  bind_rows() %>% 
  filter(abs(slope) <= 10) %>% 
  ggplot() +
  geom_boxplot(aes(slope, scenario, group = scenario)) +
  geom_vline(xintercept = 0, linewidth = 1) +
  theme_minimal()

slopes_df_list %>% 
  bind_rows() %>% 
  ggplot() +
  geom_boxplot(aes(cor, scenario, group = scenario)) +
  geom_vline(xintercept = 0, linewidth = 1) +
  theme_minimal()

slopes_df_list %>% 
  bind_rows() %>% 
  mutate(deviance_expl = (null_deviance - deviance) / null_deviance) %>% 
  ggplot() +
  geom_boxplot(aes(deviance_expl, scenario, group = scenario)) +
  geom_vline(xintercept = 0, linewidth = 1) +
  theme_minimal()



