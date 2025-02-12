

compare_apriori_many_sim1 <- function(specs_df_onerow, nsim, cl = NULL) {
  stopifnot(all(colnames(specs_df_onerow) %in% c(names(sim1_defaults), "scenario", "seed")))
  
  sim1_values <- sim1_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim1_params) {
      sim1_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  if (!is.null(cl)) {
    result_list <- parLapply(cl, X = 1:nsim, fun = compare_apriori_one_sim1,
                             scenario = specs_df_onerow$scenario,
                             seed = specs_df_onerow$seed,
                             n1 = sim1_values["n1"], 
                             n2 = sim1_values["n2"], 
                             J = sim1_values["J"], 
                             xi = sim1_values["xi"],
                             beta0 = sim1_values["beta0"],
                             beta1 = sim1_values["beta1"], 
                             alpha0_1 = sim1_values["alpha0_1"], 
                             alpha1_1 = sim1_values["alpha1_1"],
                             alpha0_2 = sim1_values["alpha0_2"], 
                             alpha1_2 = sim1_values["alpha1_2"],
                             zeta = sim1_values["zeta"], 
                             sigma = sim1_values["sigma"])
    
  } else {
    result_list <- list()
    for (i in 1:nsim) {
      result_list[[i]] <- compare_apriori_one_sim1(iter = i,
                                       n1 = sim1_values["n1"], 
                                       n2 = sim1_values["n2"], 
                                       J = sim1_values["J"], 
                                       xi = sim1_values["xi"],
                                       beta0 = sim1_values["beta0"],
                                       beta1 = sim1_values["beta1"], 
                                       alpha0_1 = sim1_values["alpha0_1"], 
                                       alpha1_1 = sim1_values["alpha1_1"],
                                       alpha0_2 = sim1_values["alpha0_2"], 
                                       alpha1_2 = sim1_values["alpha1_2"],
                                       zeta = sim1_values["zeta"], 
                                       sigma = sim1_values["sigma"]) 
    }
  }
  
  return(result_list)
}

compare_apriori_one_sim1 <- function(iter, 
                                     scenario = 0,
                                     n1 = 30, 
                                     n2 = 30, 
                                     J = 5, 
                                     xi = 1,
                                     beta0 = 0.5,
                                     beta1 = 1, 
                                     alpha0_1 = 0, 
                                     alpha1_1 = 1,
                                     alpha0_2 = -0.5, 
                                     alpha1_2 = 0.5,
                                     zeta = 0, 
                                     sigma = 0,
                                     seed) {
  set.seed(seed + (iter * 17))
  
  dat <- simulate_data_sim1(n1, n2, J, xi, beta0, beta1, alpha0_1, alpha1_1, 
                            alpha0_2, alpha1_2, zeta, sigma)
  
  # Comparison one: difference in slopes of I(y1) ~ x, I(y2) ~ x
  c1_y1 <- glm(rowSums(dat$y$y_1) > 0 ~ dat$occ.covs$x[dat$sites[[1]]], family = binomial())
  c1_y2 <- glm(rowSums(dat$y$y_2) > 0 ~ dat$occ.covs$x[dat$sites[[2]]], family = binomial())
  
  # Comparison one: difference in slopes of sum(y1) ~ x, sum(y2) ~ x
  c2_y1 <- glm(cbind(rowSums(dat$y$y_1 > 0), rowSums(dat$y$y_1 == 0)) ~ dat$occ.covs$x[dat$sites[[1]]], family = binomial())
  c2_y2 <- glm(cbind(rowSums(dat$y$y_2 > 0), rowSums(dat$y$y_2 == 0)) ~ dat$occ.covs$x[dat$sites[[2]]], family = binomial())
  
}