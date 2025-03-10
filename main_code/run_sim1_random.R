library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
source("support_fn/sim1_fn.R")

nsim <- 32 * 1000
ncores <- 32

cl <- makeCluster(ncores)

capture <- clusterEvalQ(cl, {
  source("support_fn/sim1_fn.R")
})

this_result <- run_many_sim_randomspecs(nsim = nsim, cl)

stopCluster(cl)
rm(cl)

saveRDS(this_result, "sim_output/sim1_random_specs.RDS")



# Look at MSE of beta0
mod_results <- this_result$estimation_result %>% 
  filter(param == "beta1") %>%
  select(param, iter, type, MSE) %>% 
  pivot_wider(names_from = type, values_from = MSE) %>% 
  mutate(improvement = one - joint) %>% 
  left_join(this_result$specs_df, by = "iter")

levers <- c("n1", "n2", "J", "xi", "beta0", "beta1", "alpha0_1",
            "alpha0_2", "alpha1_1", "alpha1_2", "abs(zeta)", "sigma")
interactions <- t(combn(levers, 2)) %>% apply(1, function(x) paste0(x, collapse = ":"))

fmla <- as.formula(paste0("improvement ~ ", paste(c(levers, interactions), collapse = "+")))

fit <- lm(fmla, data = mod_results)
summary(fit)


vis_one <- function(dim1, dim2, fit, levers) {
  
  vec1 <- switch(dim1, 
                 n1 = seq(20, 200, length.out = 10),
                 n2 = seq(200, 2000, length.out = 10),
                 J = 2:6,
                 xi = seq(0.1, 1, length.out = 9),
                 beta0 = seq(-2, 2, length.out = 11),
                 beta1 = seq(-2, 2, length.out = 11),
                 alpha0_1 = seq(-2, 2, length.out = 11),
                 alpha0_2 = seq(-2, 2, length.out = 11),
                 alpha1_1 = seq(-2, 2, length.out = 11),
                 alpha1_2 = seq(-2, 2, length.out = 11),
                 zeta = seq(0, 2, length.out = 11),
                 sigma = seq(0, 5, length.out = 10)
                 )
  vec2 <- switch(dim2, 
                 n1 = seq(20, 200, length.out = 10),
                 n2 = seq(200, 2000, length.out = 10),
                 J = 2:6,
                 xi = seq(0.1, 1, length.out = 9),
                 beta0 = seq(-2, 2, length.out = 11),
                 beta1 = seq(-2, 2, length.out = 11),
                 alpha0_1 = seq(-2, 2, length.out = 11),
                 alpha0_2 = seq(-2, 2, length.out = 11),
                 alpha1_1 = seq(-2, 2, length.out = 11),
                 alpha1_2 = seq(-2, 2, length.out = 11),
                 zeta = seq(0, 2, length.out = 11),
                 sigma = seq(0, 5, length.out = 11)
                 )
  
  newdata <- as.data.frame(expand.grid(vec1, vec2))
  colnames(newdata) <- c(dim1, dim2)
  for (i in 1:length(sim1_params)) {
    if (!sim1_params[i] %in% colnames(newdata)) {
      newdata[[sim1_params[i]]] <- sim1_defaults[sim1_params[i]]
    }
  }
  
  newdata$predicted <- predict(fit, newdata = newdata)
  
  ggplot(newdata) +
    geom_tile(aes(.data[[dim1]], .data[[dim2]], fill = predicted)) +
    scale_fill_viridis_c(option = "plasma", end = 0.9) +
    theme_minimal()
}


