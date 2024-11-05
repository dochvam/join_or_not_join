library(tidyverse)
library(gridExtra)

ests_list <- list()
rmse_list <- list()
cv_list <- list()

# For each sim x dataset, I want to choose btw cam-only and joint based on:
# 1. Unobserved, correct measures of performance
#   1a. Prediction RMSE
#   1b. Inference RMSE (betas only)
#   1c. Inference RMSE (all params) <- can't do this as written bc num. of params is different
# 2. Observed measures of performance
#   2a. OOS CV of camera data
#   2b. theta1

decisions_list <- list()

for (i in 1:5) {
  results_joint <- readRDS(paste0("output/scenario", i, "_joint_results.RDS"))
  results_cam <- readRDS(paste0("output/scenario", i, "_cameraOnly_results.RDS"))
  
  decisions_list[[i]] <- data.frame(
    dataset = sort(unique(results_joint$estimates$dataset_ID)),
    scenario = i,
    d1a = NA,
    d1b = NA,
    d2a = NA,
    d2b = NA,
    RMSE_camera_1a = NA,
    RMSE_joint_1a = NA,
    RMSE_camera_1b = NA,
    RMSE_joint_1b = NA
  )
  
  # Compare RMSE
  RMSE_all <- bind_rows(results_joint$rmse, results_cam$rmse) %>% 
    pivot_wider(names_from = modtype, values_from = RMSE) %>% 
    arrange(dataset_ID)
  
  this_RMSE <- RMSE_all %>% filter(goal == "prediction")
  decisions_list[[i]]$d1a <- ifelse(this_RMSE$joint < this_RMSE$camera_only,
                                    "Joint", "PA-only")
  decisions_list[[i]]$RMSE_joint_1a <- this_RMSE$joint
  decisions_list[[i]]$RMSE_camera_1a <- this_RMSE$camera_only
  
  this_RMSE <- RMSE_all %>% filter(goal == "inference_betas")
  decisions_list[[i]]$d1b <- ifelse(this_RMSE$joint < this_RMSE$camera_only,
                                    "Joint", "PA-only")
  decisions_list[[i]]$RMSE_joint_1b <- this_RMSE$joint
  decisions_list[[i]]$RMSE_camera_1b <- this_RMSE$camera_only
  
  # Compare CV
  CV_all <- bind_rows(results_joint$CV, results_cam$CV) %>% 
    filter(type == "camera") %>% 
    pivot_wider(names_from = modtype, values_from = value) %>% 
    arrange(dataset_ID)
  decisions_list[[i]]$d2a <- ifelse(CV_all$joint < CV_all$camera_only,
                                    "Joint", "PA-only")
  
  # Look at theta1 estimate
  theta1_est <- results_joint$estimates %>% 
    filter(param == "theta1") %>% 
    mutate(positive = sign(est - 1.96*se > 0)) %>% 
    arrange(dataset_ID)
  decisions_list[[i]]$d2b <- ifelse(theta1_est$positive,
                                    "Joint", "PA-only")
}

decisions <- bind_rows(decisions_list)

# What was the decision?

decisions %>% 
  select(-RMSE_camera_1a, -RMSE_camera_1b,
         -RMSE_joint_1a, -RMSE_joint_1b,) %>% 
  rename(prediction = d1a, inference = d1b, 
         CV = d2a, corr = d2b) %>%
  pivot_longer(cols = c("CV", "corr")) %>% 
  rename(obs = value, criterion = name) %>% 
  filter(!is.na(obs)) %>% 
  ggplot() + 
  geom_bar(aes(criterion, fill = obs)) +
  coord_flip() +
  facet_grid(scenario~"Decision") +
  scale_fill_manual("What was the\ndecision?",
                    values = c("#202C59", "#D95D39"), na.value = "white") +
  theme_minimal()


# Was the decision correct?
decisions %>% 
  select(-RMSE_camera_1a, -RMSE_camera_1b,
         -RMSE_joint_1a, -RMSE_joint_1b,) %>% 
  rename(prediction = d1a, inference = d1b, 
         CV = d2a, theta1 = d2b) %>%
  pivot_longer(cols = c("CV", "theta1")) %>% 
  rename(obs = value, criterion = name) %>% 
  mutate(pred = obs == prediction,
         inf  = obs == inference) %>%
  select(-prediction, -inference) %>% 
  pivot_longer(cols = c("pred", "inf")) %>% 
  rename(agreed = value, goal = name) %>% 
  filter(!is.na(agreed)) %>% 
  mutate(agreed = ifelse(agreed, "Yes", "No")) %>% 
  ggplot() + 
  geom_bar(aes(criterion, fill = agreed)) +
  coord_flip() +
  facet_grid(scenario~paste0("Goal: ", goal)) +
  scale_fill_manual("Was the decision\ncorrect?",
                    values = c("#cc3333", "#333333"), na.value = "white") +
  theme_bw()


# How much did we lose by making the decision we made?
decisions %>% 
  rename(prediction = d1a, inference = d1b, 
       CV = d2a, corr = d2b) %>%
  pivot_longer(cols = c("CV", "corr")) %>% 
  rename(obs = value, criterion = name) %>% # obs = our decision 
  mutate(
    loss_pred = ifelse(obs == "Joint", RMSE_joint_1a - RMSE_camera_1a,
                       RMSE_camera_1a - RMSE_joint_1a),
    loss_inf = ifelse(obs == "Joint", RMSE_joint_1b - RMSE_camera_1b,
                      RMSE_camera_1b - RMSE_joint_1b),
  ) %>%
  mutate(correct_pred = loss_pred < 0, correct_inf = loss_inf < 0) %>% 
  mutate(
    loss_pred = abs(loss_pred),
    loss_inf = abs(loss_inf)
  ) %>% 
  filter(!is.na(correct_inf)) %>% 
  ggplot() + 
  geom_boxplot(aes(loss_inf, correct_inf, group = correct_inf, col = correct_inf), boundary = 0) +
  facet_grid(scenario~paste0(criterion)) +
  scale_color_manual("Correct choice?",
                    values = c("#cc3333", "#333333"), na.value = "white") +
  theme_bw() +
  xlim(c(0, 2.5)) +
  xlab("Diff. in RMSE btw. best and worst model") +
  ggtitle("Goal: inference") +
  ylab("") + theme(axis.ticks = element_blank(), 
                   axis.text.y = element_blank())

decisions %>% 
  rename(prediction = d1a, inference = d1b, 
         CV = d2a, corr = d2b) %>%
  pivot_longer(cols = c("CV", "corr")) %>% 
  rename(obs = value, criterion = name) %>% # obs = our decision 
  mutate(
    loss_pred = ifelse(obs == "Joint", RMSE_joint_1a - RMSE_camera_1a,
                       RMSE_camera_1a - RMSE_joint_1a),
    loss_inf = ifelse(obs == "Joint", RMSE_joint_1b - RMSE_camera_1b,
                      RMSE_camera_1b - RMSE_joint_1b),
  ) %>%
  mutate(correct_pred = loss_pred < 0, correct_inf = loss_inf < 0) %>% 
  mutate(
    loss_pred = abs(loss_pred),
    loss_inf = abs(loss_inf)
  ) %>% 
  filter(!is.na(correct_pred)) %>% 
  ggplot() + 
  geom_boxplot(aes(loss_pred, correct_pred, group = correct_pred, col = correct_pred), boundary = 0) +
  facet_grid(scenario~paste0(criterion)) +
  scale_color_manual("Correct choice?",
                     values = c("#cc3333", "#333333"), na.value = "white") +
  theme_bw() +
  xlim(c(0, 5)) +
  xlab("Diff. in RMSE btw. best and worst model") +
  ggtitle("Goal: prediction") +
  ylab("") + theme(axis.ticks = element_blank(), 
                   axis.text.y = element_blank())
