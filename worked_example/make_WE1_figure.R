library(tidyverse)
library(gridExtra)

source('worked_example/preprocess_data.R')

results <- read_csv("worked_example/joint_model_results_theta1fixed.csv") 
results2 <- read_csv("worked_example/joint_model_results_wNDVIdet.csv") 

covar_names <- data.frame(
  param = c("beta0", "beta1", "alpha0", "alpha1", "alpha2"),
  name = c("State intercept", "Effect of NDVI on state", "Det. intercept", "Effect of road dist. on det.", "Effect of NDVI on det.")
) %>% 
  mutate(name = factor(name, levels = rev(name)))

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

data_panel <- grid_cell_df %>% 
  left_join(as.data.frame(main_grid, xy = TRUE), by = c("grid_cell" = "lyr.1")) %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = n_inat)) +
  geom_point(aes(x, y, col = n_ct > 0, alpha = n_ct > 0)) +
  scale_color_manual("Cell contains CT", values = c(NA, "#dd3344")) +
  scale_alpha_manual("Cell contains CT", values = c(0.1, 1)) +
  scale_fill_viridis_c("Num. iNat obs.",
                       trans = "log", na.value = "#dddddd", breaks = c(1, 8, 64, 400)) +
  theme_minimal() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("(A) Sampling effort") +
  guides(colour = guide_legend(override.aes = list(alpha=1)))


results_panel <- results %>% 
  filter(param %in% covar_names$param) %>% 
  left_join(covar_names) %>% 
  mutate(type = recode(type, joint = "Joint", one = "CT only")) %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(name, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.3)) +
  theme_bw() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~species, nrow = 2) + coord_flip() + xlab("") +
  scale_color_manual("", values = c(`Joint` = "black", `CT only` = "#dd3344")) +
  ylab("Param. estimate (95%CI)") +
  ggtitle("(B) Results of Model 1")


results2_panel <- results2 %>% 
  filter(param %in% covar_names$param) %>% 
  left_join(covar_names) %>% 
  mutate(type = recode(type, joint = "Joint", one = "CT only")) %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(name, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.3)) +
  theme_bw() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~species, nrow = 2) + coord_flip() + xlab("") +
  scale_color_manual("", values = c(`Joint` = "black", `CT only` = "#dd3344")) +
  ylab("Param. estimate (95%CI)") +
  ggtitle("(C) Results of Model 2")

legend <- get_legend(results_panel)



fig_combined <- arrangeGrob(data_panel, 
                            results_panel + theme(legend.position = "None"), 
                            results2_panel + theme(legend.position = "None"), 
                            right = legend,
                            nrow = 1, widths = c(1, 0.7, 0.7))
ggsave("WE1_fig.jpg", fig_combined, width = 15, height = 5)
