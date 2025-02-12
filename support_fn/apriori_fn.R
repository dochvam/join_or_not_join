

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
    }
  }
  
  return(bind_rows(result_list))
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
  
  # browser()
  # Comparison one: difference in slopes of I(y1) ~ x, I(y2) ~ x
  c1_y1 <- glm(rowSums(dat$y$y_1) > 0 ~ dat$occ.covs$x[dat$sites[[1]]], family = binomial())
  c1_y2 <- glm(rowSums(dat$y$y_2) > 0 ~ dat$occ.covs$x[dat$sites[[2]]], family = binomial())
  
  # Comparison one: difference in slopes of sum(y1) ~ x, sum(y2) ~ x
  c2_y1 <- glm(cbind(rowSums(dat$y$y_1 > 0), rowSums(dat$y$y_1 == 0)) ~ dat$occ.covs$x[dat$sites[[1]]], family = binomial())
  c2_y2 <- glm(cbind(rowSums(dat$y$y_2 > 0), rowSums(dat$y$y_2 == 0)) ~ dat$occ.covs$x[dat$sites[[2]]], family = binomial())
  
  metrics_df <- data.frame(
    d1 = abs(summary(c1_y1)$coefficients[2, 1] - summary(c1_y2)$coefficients[2, 1]),
    d2 = abs(summary(c2_y1)$coefficients[2, 1] - summary(c2_y2)$coefficients[2, 1]),
    iter = iter,
    scenario = scenario
  )
  
  metrics_df
  
}

compare_apriori_many_sim2 <- function(specs_df_onerow, nsim, cl = NULL) {
  stopifnot(all(colnames(specs_df_onerow) %in% c(names(sim2_defaults), "scenario", "seed")))
  
  sim2_values <- sim2_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim2_params) {
      sim2_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  if (!is.null(cl)) {
    # browser()
    result_list <- parLapply(cl, X = 1:nsim, fun = compare_apriori_one_sim2,
                             scenario = specs_df_onerow$scenario,
                             seed = specs_df_onerow$seed,
                             S = sim2_values["S"],
                             nPA = sim2_values["nPA"],
                             J = sim2_values["J"],
                             xi = sim2_values["xi"],
                             eta = sim2_values["eta"],
                             q = sim2_values["q"],
                             beta0 = sim2_values["beta0"],
                             beta1 = sim2_values["beta1"],
                             alpha0 = sim2_values["alpha0"],
                             alpha1 = sim2_values["alpha1"],
                             theta0 = sim2_values["theta0"],
                             theta1 = sim2_values["theta1"],
                             phi = sim2_values["phi"],
                             zeta = sim2_values["zeta"],
                             sigma = sim2_values["sigma"])
    
  } else {
    result_list <- list()
    for (i in 1:nsim) {
      result_list[[i]] <- compare_apriori_one_sim2(scenario = specs_df_onerow$scenario,
                                                   seed = specs_df_onerow$seed,
                                                   iter = i,
                                                   S = sim2_values["S"],
                                                   nPA = sim2_values["nPA"],
                                                   J = sim2_values["J"],
                                                   xi = sim2_values["xi"],
                                                   eta = sim2_values["eta"],
                                                   q = sim2_values["q"],
                                                   beta0 = sim2_values["beta0"],
                                                   beta1 = sim2_values["beta1"],
                                                   alpha0 = sim2_values["alpha0"],
                                                   alpha1 = sim2_values["alpha1"],
                                                   theta0 = sim2_values["theta0"],
                                                   theta1 = sim2_values["theta1"],
                                                   phi = sim2_values["phi"],
                                                   zeta = sim2_values["zeta"],
                                                   sigma = sim2_values["sigma"]) 
    }
  }
  
  return(bind_rows(result_list))
}


compare_apriori_one_sim2 <- function(iter, 
                                     scenario = 0,
                                     S = 100,
                                     nPA = 50,
                                     J = 3,
                                     xi = 1,
                                     eta = 100,
                                     q = 0.7,
                                     beta0 = 0.5,
                                     beta1 = 1,
                                     alpha0 = 0,
                                     alpha1 = 1,
                                     theta0 = -3,
                                     theta1 = 1,
                                     phi = 0.25,
                                     zeta = 0,
                                     sigma = 0,
                                     seed) {
  set.seed(seed + (iter * 19))
  
  dat <- simulate_data_sim2(S, nPA, J, xi, eta, q, beta0, beta1, alpha0,
                            alpha1, theta0, theta1, phi, zeta, sigma)
  
  unified_df <- data.frame(
    x = dat$data$x,
    y_PO = dat$data$y_PO,
    E_PO = dat$data$E_PO
  ) %>% 
    mutate(r_PO = y_PO / E_PO) %>% 
    left_join(data.frame(y_PA_any = rowSums(dat$data$y_PA) > 0,
                         x = dat$data$x_PA), by = "x") %>% 
    filter(!is.na(y_PA_any), !is.nan(r_PO)) %>% 
    mutate(y_PA_any = as.numeric(y_PA_any))
  
  # Comparison 1: what's the slope of 
  fit1 <- glm(y_PA_any ~ r_PO, family = binomial, unified_df)
  fit_null <- glm(y_PA_any ~ 1, family = binomial, unified_df)
  s <- summary(fit1)
  s_null <- summary(fit_null)
  
  metrics_df <- data.frame(
    glm_slope = s$coefficients[2, 1],
    glm_slope_pval = s$coefficients[2, 4],
    pct_deviance = (s_null$deviance - s$deviance) / s_null$deviance,
    cor = cor(unified_df$y_PA_any, unified_df$r_PO),
    iter = iter,
    scenario = scenario
  )

  metrics_df
}

