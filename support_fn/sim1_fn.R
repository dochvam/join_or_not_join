library(spOccupancy)
library(tidyverse)

sim1_params <- c("n1", "n2", "J", "xi", "beta0", "beta1", "alpha0_1", 
                 "alpha1_1", "alpha0_2", "alpha1_2", "zeta", "sigma")
sim1_defaults <- c(n1 = 50,
                   n2 = 50,
                   J = 3,
                   xi = 1,
                   beta0 = 0.5,
                   beta1 = 1,
                   alpha0_1 = 0,
                   alpha1_1 = 1,
                   alpha0_2 = -0.5,
                   alpha1_2 = 0.5,
                   zeta = 0,
                   sigma = 0)

simulate_data_sim1 <- function(n1, 
                               n2, 
                               J, 
                               xi,
                               beta0,
                               beta1, 
                               alpha0_1, 
                               alpha1_1,
                               alpha0_2, 
                               alpha1_2,
                               zeta, 
                               sigma) {
  stopifnot(xi <= 1 && xi >= 0)
  
  x_1 <- runif(n1, -xi, xi)
  x_2 <- runif(n2, -1, 1)
  
  w_1 <- matrix(runif(n1 * J, -1, 1), nrow = n1)
  w_2 <- matrix(runif(n2 * J, -1, 1), nrow = n2)
  
  epsilon <- rnorm(n2, mean = x_2 * zeta, sd = sigma)
  
  psi_1 <- plogis(beta0 + beta1 * x_1)
  psi_2 <- plogis(beta0 + beta1 * x_2 + epsilon)
  
  z_1 <- rbinom(n1, 1, psi_1)
  z_2 <- rbinom(n2, 1, psi_2)
  
  p_1 <- plogis(alpha0_1 + alpha1_1 * w_1)
  p_2 <- plogis(alpha0_2 + alpha1_2 * w_2)
  
  y_1 <- p_1
  for (i in 1:nrow(y_1)) {
    y_1[i, ] <- rbinom(J, 1, p_1[i, ] * z_1[i])
  }
  
  y_2 <- p_2
  for (i in 1:nrow(y_2)) {
    y_2[i, ] <- rbinom(J, 1, p_2[i, ] * z_2[i])
  }
  
  # Return formatted for spOccupancy
  return(list(
    y = list(
      y_1 = y_1, y_2 = y_2
    ),
    occ.covs = rbind(data.frame(site = 1:n1, x = x_1), 
                     data.frame(site = n1 + 1:n2, x = x_2)),
    det.covs = list(
      list(w_1 = w_1), list(w_2 = w_2)
    ),
    sites = list(
      1:n1, n1 + 1:n2
    ),
    true_values = data.frame(
      param = sim1_params,
      truth = c(n1, n2, J, xi, beta0, beta1, alpha0_1, 
                alpha1_1, alpha0_2, alpha1_2, zeta, sigma)
    )
  ))
}


run_one_sim1 <- function(iter, 
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
                         sigma = 0) {
  set.seed(scenario + 1000 * iter + 14908)
  
  dat <- simulate_data_sim1(n1, n2, J, xi, beta0, beta1, alpha0_1, alpha1_1, 
                            alpha0_2, alpha1_2, zeta, sigma)
  dat_OOS <- simulate_data_sim1(n1, n2, J, xi = 1, beta0, beta1, alpha0_1, alpha1_1, 
                                alpha0_2, alpha1_2, zeta, sigma)
  
  dat_one <- list(
    y = dat$y[[1]],
    occ.covs = dat$occ.covs[1:n1, ],
    det.covs = dat$det.covs[[1]]
  )
  
  fit_joint <- spOccupancy::intPGOcc(~x, list(~w_1, ~w_2), data = dat, 
                                     n.burn = 2500, n.samples = 5000,
                                     n.chains = 2,
                                     verbose = F)
  # browser()
  
  fit_one <- spOccupancy::PGOcc(~x, ~w_1, data = dat_one, 
                                  n.burn = 2500, n.samples = 5000,
                                  n.chains = 2, verbose = F)
  
  estimation_result <- bind_rows(
    data.frame(
      param = c("beta0", "beta1", "alpha0_1", "alpha1_1"),
      mean = colMeans(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2])),
      sd   = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2]), 2, sd),
      Q025 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2]), 2, quantile, probs = 0.025),
      Q975 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2]), 2, quantile, probs = 0.975)
    ) %>% 
      left_join(dat$true_values, by = "param") %>% 
      mutate(type = "joint"),
    data.frame(
      param = c("beta0", "beta1", "alpha0_1", "alpha1_1"),
      mean = colMeans(cbind(fit_one$beta.samples, fit_one$alpha.samples)),
      sd   = apply(cbind(fit_one$beta.samples, fit_one$alpha.samples), 2, sd),
      Q025 = apply(cbind(fit_one$beta.samples, fit_one$alpha.samples), 2, quantile, probs = 0.025),
      Q975 = apply(cbind(fit_one$beta.samples, fit_one$alpha.samples), 2, quantile, probs = 0.975)
    ) %>% 
      left_join(dat$true_values, by = "param") %>% 
      mutate(type = "one")
  )
  
  # TODO next: calculate CV using PA-1 for both joint and single-dataset model
  ncv_sites <- 30
  stopifnot(ncv_sites <= n1)
  
  cv_result <- calculate_CV_sim1(ncv_sites, ests = estimation_result, 
                                    dat_OOS = dat_OOS)
  
  runtime_result <- data.frame(
    type = c("joint", "one"),
    vaule = unname(c(fit_joint$run.time[3], fit_one$run.time[3]))
  )
  
  return(list(estimation_result = estimation_result %>% mutate(iter = iter),
              cv_result = cv_result %>% mutate(iter = iter), 
              runtime_result = runtime_result %>% mutate(iter = iter)))
}


run_many_sim1 <- function(specs_df_onerow, nsim, cl = NULL) {
  
  sim1_values <- sim1_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim1_params) {
      sim1_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  if (!is.null(cl)) {
    result_list <- parLapply(cl, X = 1:nsim, fun = run_one_sim1,
                             scenario = specs_df_onerow$scenario,
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
      result_list[[i]] <- run_one_sim1(iter = i,
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
  
  estimation_result_list <- cv_result_list <- runtime_result_list <- list()
  for (i in 1:nsim) {
    estimation_result_list[[i]] <- result_list[[i]]$estimation_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    cv_result_list[[i]] <- result_list[[i]]$cv_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    runtime_result_list[[i]]<- result_list[[i]]$runtime_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
  }
  
  return(list(
    estimation_result = bind_rows(estimation_result_list),
    cv_result = bind_rows(cv_result_list),
    runtime_result = bind_rows(runtime_result_list)
  ))
}


calculate_CV_sim1 <- function(ncv_sites, ests, dat_OOS) {
  res <- data.frame()
  for (this_type in unique(ests$type)) {
    this_ests <- ests %>% filter(type == this_type)
    pred_psi <- plogis(this_ests$mean[this_ests$param == "beta0"] +
     this_ests$mean[this_ests$param == "beta1"] * dat_OOS$occ.covs[1:ncv_sites, "x"])
    pred_p <- plogis(this_ests$mean[this_ests$param == "alpha0_1"] +
     this_ests$mean[this_ests$param == "alpha1_1"] * dat_OOS$det.covs[[1]]$w_1[1:ncv_sites, ])
    
    predicted_value <- matrix(pred_psi, nrow = ncv_sites, ncol = ncol(pred_p)) * 
      pred_p
    
    brier_score <- mean((predicted_value - dat_OOS$y[[1]][1:ncv_sites, ])^2)
    
    res <- bind_rows(res, data.frame(brier_score = brier_score, type = this_type))
  }
  res
}


