library(tidyverse)
library(spOccupancy)
library(gridExtra)

this_res <- readRDS("sim_output/sim1_5.RDS")
estimation_results <- this_res$estimation_results
cv_results <- this_res$cv_results
ppc_results <- this_res$ppc_results
specs_df <- this_res$specs_df

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# How different are point estimates of beta0 and beta1 btw. single-dataset and
# joint models?

comparison_df <- estimation_results %>% 
  bind_rows() %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  select(param, iter, scenario, type, mean) %>% 
  pivot_wider(names_from = "type", values_from = "mean") %>% 
  mutate(diff = joint - one) %>% 
  group_by(scenario) %>% 
  summarize(median = median(diff), ub = quantile(diff, 0.75),
            lb = quantile(diff, 0.25)) %>%
  left_join(specs_df, by = "scenario")

cv_df <- cv_results %>% 
  bind_rows() %>% 
  select(iter, scenario, type, brier_score) %>% 
  pivot_wider(names_from = "type", values_from = "brier_score") %>% 
  mutate(diff = joint - one) %>% 
  group_by(scenario) %>% 
  summarize(median = median(diff), ub = quantile(diff, 0.75),
            lb = quantile(diff, 0.25)) %>%
  left_join(specs_df, by = "scenario")

ppc_dev_df <- ppc_results %>% 
  bind_rows() %>% 
  select(iter, scenario, type, ppc_pval_dev) %>% 
  filter(type == "D2") %>% 
  group_by(scenario) %>% 
  summarize(median = median(ppc_pval_dev), ub = quantile(ppc_pval_dev, 0.75),
            lb = quantile(ppc_pval_dev, 0.25)) %>%
  left_join(specs_df, by = "scenario") %>% 
  mutate(metric = "Deviance")
ppc_cs_df <- ppc_results %>% 
  bind_rows() %>% 
  select(iter, scenario, type, ppc_pval_cs) %>% 
  filter(type == "D2") %>% 
  group_by(scenario) %>% 
  summarize(median = median(ppc_pval_cs), ub = quantile(ppc_pval_cs, 0.75),
            lb = quantile(ppc_pval_cs, 0.25)) %>%
  left_join(specs_df, by = "scenario") %>% 
  mutate(metric = "Chi-squared")


#### Make the panels
slopes_plot <- comparison_df %>% 
  ggplot(aes(zeta, median, ymax = ub, ymin = lb, col = as.factor(n1),
             group = as.factor(n1))) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  geom_line(position = position_dodge(width = 0.1)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nprimary\ndataset", begin = 0, end = 0.7) +
  xlab("Logit-scale bias in the effect of x on PA-2") +
  ylab("Difference in point estimates") +
  geom_hline(yintercept = 0) +
  ggtitle("(3) Comparing joint vs. single-dataset ests.")


# Cross-validation
cv_plot <- cv_df %>% 
  ggplot(aes(zeta, median, ymax = ub, ymin = lb, col = as.factor(n1),
                     group = as.factor(n1))) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  geom_line(position = position_dodge(width = 0.1)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nprimary\ndataset", begin = 0, end = 0.7) +
  xlab("Logit-scale bias in the effect of x on PA-2") +
  ylab("Difference in OOS Brier score") +
  geom_hline(yintercept = 0) +
  ggtitle("(2) Cross validation")
  

# Posterior predictive checks
ppc_cs_plot <- ppc_cs_df %>% 
  ggplot(aes(zeta, median, ymax = ub, ymin = lb, col = as.factor(n1),
             group = as.factor(n1))) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  geom_line(position = position_dodge(width = 0.1)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nprimary\ndataset", begin = 0, end = 0.7) +
  xlab("Logit-scale bias in the effect of x on PA-2") +
  ylab("Bayesian p-value") +
  geom_hline(yintercept = c(0.05, 0.95)) +
  ggtitle("(1a) Cross validation (chi-squared)")

ppc_dev_plot <- ppc_dev_df %>% 
  ggplot(aes(zeta, median, ymax = ub, ymin = lb, col = as.factor(n1),
             group = as.factor(n1))) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  geom_line(position = position_dodge(width = 0.1)) +
  theme_minimal() +
  scale_color_viridis_d("Num. sites in\nprimary\ndataset", begin = 0, end = 0.7) +
  xlab("Logit-scale bias in the effect of x on PA-2") +
  ylab("Bayesian p-value") +
  geom_hline(yintercept = c(0.05, 0.95)) +
  ggtitle("(1b) Cross validation (deviance)")


legend <- get_legend(ppc_dev_plot)


plot_all <- arrangeGrob(
  ppc_cs_plot + theme(legend.position = "None"),
  ppc_dev_plot + theme(legend.position = "None"),
  cv_plot + theme(legend.position = "None"),
  slopes_plot + theme(legend.position = "None"),
  right = legend
)
ggsave("aposteriori_fig.jpg", plot_all, width = 9, height = 6)
