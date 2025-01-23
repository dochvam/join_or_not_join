library(spOccupancy)
library(tidyverse)
library(nimbleEcology)
library(MCMCvis)

sim2_params <- c("S", "nPA", "J", "xi", "eta", "q",
                 "beta0", "beta1", "alpha0", "alpha1",
                 "theta0", "theta1", "phi", "zeta", "sigma")
sim2_defaults <- c(S = 100,
                   nPA = 30,
                   J = 3,
                   xi = 1,
                   eta = 30,
                   q = 0.7,
                   beta0 = 0.5,
                   beta1 = 1,
                   alpha0 = 0,
                   alpha1 = 1,
                   theta0 = -3,
                   theta1 = 1,
                   phi = 0.25,
                   zeta = 0,
                   sigma = 0)

simulate_data_sim2 <- function(S,
                               nPA,
                               J,
                               xi,
                               eta,
                               q,
                               beta0,
                               beta1,
                               alpha0,
                               alpha1,
                               theta0,
                               theta1,
                               phi,
                               zeta,
                               sigma) {
  stopifnot(xi <= 1)
  
  # Latent process
  x <- -1 + 2/S * (1:S)
  lambda <- exp(beta0 + beta1 * x)
  
  # PA data
  cells_to_sample <- 1:(S * xi)
  grid_cells <- sample(cells_to_sample, size = nPA, replace = TRUE)
  w <- matrix(runif(nPA * J, -1, 1), nrow = nPA)
  
  x_PA <- x[grid_cells]
  psi <- icloglog(log(lambda[grid_cells]))
  z <- rbinom(nPA, 1, prob = psi)
  p <- icloglog(alpha0 + alpha1 * w)
  
  y_PA <- p
  for (i in 1:nrow(y_PA)) {
    y_PA[i, ] <- rbinom(J, 1, p[i, ] * z[i])
  }
  
  # PO data
  epsilon <- rnorm(S, mean = x * zeta, sd = sigma)
  mu <- exp(theta0 + theta1 * log(lambda) + epsilon)
  
  E_PO <- rpois(S, eta) * rbinom(S, 1, q)
  
  y_PO <- rnbinom(S, mu = mu * E_PO, size = 1/phi)
  
  # Return formatted for NIMBLE
  return(list(
    data = list(
      y_PA = y_PA,
      x_PA = x_PA,
      w = w,
      E_PO = E_PO,
      y_PO = y_PO,
      x = x
    ),
    constants = list(
      J = J, nPA = nPA, S = S
    ),
    inits = list(
      beta0 = sim2_defaults["beta0"],
      beta1 = sim2_defaults["beta1"],
      alpha0 = sim2_defaults["alpha0"],
      alpha1 = sim2_defaults["alpha1"],
      theta0 = sim2_defaults["theta0"],
      theta1 = sim2_defaults["theta1"],
      phi = sim2_defaults["phi"]
    ),
    true_values = data.frame(
      param = sim2_params,
      truth = c(S, nPA, J, xi, eta, q, beta0, beta1, alpha0, alpha1, 
                theta0, theta1, phi, zeta, sigma)
    )
  ))
}


# run_one_sim2 <- function(iter,
#                          S = 100,
#                          nPA = 30,
#                          J = 3,
#                          xi = 1,
#                          eta = 30,
#                          q = 0.7,
#                          beta0 = 0.5,
#                          beta1 = 1,
#                          alpha0 = 0,
#                          alpha1 = 1,
#                          theta0 = -3,
#                          theta1 = 1,
#                          phi = 0.25,
#                          zeta = 0,
#                          sigma = 0) {
# 
#   dat <- simulate_data_sim2(n1, n2, J, xi, beta0, beta1, alpha0_1, alpha1_1,
#                             alpha0_2, alpha1_2, zeta, sigma)
#   dat_OOS <- simulate_data_sim2(n1, n2, J, xi, beta0, beta1, alpha0_1, alpha1_1,
#                                 alpha0_2, alpha1_2, zeta, sigma)
# 
#   dat_one <- list(
#     y = dat$y[[1]],
#     occ.covs = dat$occ.covs[1:n1, ],
#     det.covs = dat$det.covs[[1]]
#   )
# 
#   fit_joint <- spOccupancy::intPGOcc(~x, list(~w_1, ~w_2), data = dat,
#                                      n.burn = 2500, n.samples = 5000,
#                                      n.chains = 2,
#                                      verbose = F)
#   browser()
# 
#   fit_one <- spOccupancy::PGOcc(~x, ~w_1, data = dat_one,
#                                 n.burn = 2500, n.samples = 5000,
#                                 n.chains = 2, verbose = F)
# 
#   estimation_result <- bind_rows(
#     data.frame(
#       param = c("beta0", "beta1", "alpha0_1", "alpha1_1"),
#       mean = colMeans(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2])),
#       Q025 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2]), 2, quantile, probs = 0.025),
#       Q975 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:2]), 2, quantile, probs = 0.975)
#     ) %>%
#       left_join(dat$true_values, by = "param") %>%
#       mutate(type = "joint"),
#     data.frame(
#       param = c("beta0", "beta1", "alpha0_1", "alpha1_1"),
#       mean = colMeans(cbind(fit_one$beta.samples, fit_one$alpha.samples)),
#       Q025 = apply(cbind(fit_one$beta.samples, fit_one$alpha.samples), 2, quantile, probs = 0.025),
#       Q975 = apply(cbind(fit_one$beta.samples, fit_one$alpha.samples), 2, quantile, probs = 0.975)
#     ) %>%
#       left_join(dat$true_values, by = "param") %>%
#       mutate(type = "one")
#   )
# 
#   # TODO next: calculate CV using PA-1 for both joint and single-dataset model
#   ncv_sites <- 30
#   stopifnot(ncv_sites <= n1)
# 
#   cv_result <- calculate_CV_sim2(ncv_sites, ests = estimation_result,
#                                     dat_OOS = dat_OOS)
# 
#   runtime_result <- data.frame(
#     type = c("one", "joint"),
#     vaule = unname(c(fit_joint$run.time[3], fit_one$run.time[3]))
#   )
# 
#   return(list(estimation_result = estimation_result %>% mutate(iter = iter),
#               cv_result = cv_result %>% mutate(iter = iter),
#               runtime_result = runtime_result %>% mutate(iter = iter)))
# }


run_many_sim2 <- function(specs_df_onerow, nsim) {
  set.seed(specs_df_onerow$scenario + 14908)
  
  sim2_values <- sim2_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim2_params) {
      sim2_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  ncv_sites <- 30
  stopifnot(ncv_sites <= sim2_values["nPA"])
  
  # SETUP: simulate one dummy dataset to build a NIMBLE model
  dat_dummy <- simulate_data_sim2(S = sim2_values["S"],
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
  
  
  model_code <- nimbleCode({
    # PA data
    for (i in 1:nPA) {
      cloglog(psi[i]) <- beta0 + beta1 * x_PA[i]
      for (j in 1:J) {
        cloglog(p[i, j]) <- alpha0 + alpha1 * w[i, j]
      }
      y_PA[i, 1:J] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J], len = J)
    }
    # Priors
    beta0 ~ dnorm(0, sd = 2.72)
    beta1 ~ dnorm(0, sd = 2.72)
    alpha0 ~ dnorm(0, sd = 2.72)
    alpha1 ~ dnorm(0, sd = 2.72)
    
    # PO data
    if (type == "joint") {
      for (s in 1:S) {
        log(lambda[s]) <- beta0 + beta1 * x[s]
        log(mu[s]) <- theta0 + theta1 * log(lambda[s])
        
        y_PO[s] ~ dnbinom(size = 1/phi, prob = 1 / (1 + phi*mu[s]*E_PO[s]))
      }
      # Priors
      theta0 ~ dnorm(0, sd = 2.72)
      theta1 ~ dnorm(0, sd = 2.72)
      phi ~ dgamma(1, 1)
    }
  })
  dat_dummy$constants$type <- "joint"
  joint_mod <- nimbleModel(model_code, 
                     inits = dat_dummy$inits, 
                     constants = dat_dummy$constants,
                     data = dat_dummy$data)
  joint_mcmc <- buildMCMC(joint_mod)
  joint_complist <- compileNimble(joint_mod, joint_mcmc)
  
  dat_dummy$constants$type <- "one"
  single_mod <- nimbleModel(model_code, 
                     inits = dat_dummy$inits, 
                     constants = dat_dummy$constants,
                     data = dat_dummy$data)
  single_mcmc <- buildMCMC(single_mod)
  single_complist <- compileNimble(single_mod, single_mcmc)
  
  estimation_result_list <- cv_result_list <- runtime_list <- list()
  # Loop over iters, simulating data and reusing the NIMBLE model to 
  for (sim in 1:nsim) {
    # Simulate a dataset
    dat <- simulate_data_sim2(S = sim2_values["S"],
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
    dat_OOS <- simulate_data_sim2(S = sim2_values["S"],
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
    
    # Update the data in the model
    for (i in 1:length(dat$data)) {
      joint_complist$joint_mod[[names(dat$data)[i]]] <- dat$data[[i]]
      if (names(dat$data)[i] %in% c("y_PA", "x_PA", "w")) {
        single_complist$single_mod[[names(dat$data)[i]]] <- dat$data[[i]]
      }
    }
    
    joint_time <- system.time(
      samples_joint <- runMCMC(joint_complist$joint_mcmc, niter = 20000, nburnin = 5000,
                             nchains = 2, thin = 1, progressBar = F))
    single_time <- system.time(
      samples_single <- runMCMC(single_complist$single_mcmc, niter = 20000, nburnin = 5000,
                             nchains = 2, thin = 1, progressBar = F))
    
    summary_joint <- MCMCsummary(samples_joint) %>% mutate(type = "joint")
    summary_joint$param <- rownames(summary_joint)
    rownames(summary_joint) <- NULL
    summary_one <- MCMCsummary(samples_single) %>% mutate(type = "one")
    summary_one$param <- rownames(summary_one)
    rownames(summary_one) <- NULL
    
    estimation_result_list[[sim]] <- bind_rows(
      summary_joint, summary_one
    ) %>% 
      mutate(iter = sim, scenario = specs_df_onerow$scenario) %>% 
      left_join(dat$true_values, by= "param")
    
    cv_result_list[[sim]] <- 
       calculate_CV_sim2(ncv_sites, ests = estimation_result_list[[sim]], 
                         dat_OOS = dat_OOS) %>% 
      mutate(iter = sim, scenario = specs_df_onerow$scenario)
    
    runtime_list[[sim]] <- data.frame(
      type = c("joint", "one"),
      vaule = unname(c(joint_time[3], single_time[3]))
    ) %>% 
      mutate(iter = sim, scenario = specs_df_onerow$scenario)
  }
  
  return(list(
    estimation_result = bind_rows(estimation_result_list),
    cv_result = bind_rows(cv_result_list),
    runtime_result = bind_rows(runtime_list)
  ))
}

run_many_sim2_parallel_wrapper <- function(scenario, nsim, specs_df) {
  run_many_sim2(specs_df_onerow = specs_df[scenario, ], nsim)
}


calculate_CV_sim2 <- function(ncv_sites, ests, dat_OOS) {
  res <- data.frame()
  for (this_type in unique(ests$type)) {
    this_ests <- ests %>% filter(type == this_type)
    pred_psi <- icloglog(this_ests$mean[this_ests$param == "beta0"] +
                         this_ests$mean[this_ests$param == "beta1"] * dat_OOS$data$x_PA[1:ncv_sites])
    pred_p <- icloglog(this_ests$mean[this_ests$param == "alpha0"] +
                       this_ests$mean[this_ests$param == "alpha1"] * dat_OOS$data$w[1:ncv_sites, ])
    
    predicted_value <- matrix(pred_psi, nrow = ncv_sites, ncol = ncol(pred_p)) * 
      pred_p
    
    brier_score <- mean((predicted_value - dat_OOS$data$y_PA[1:ncv_sites, ])^2)
    
    res <- bind_rows(res, data.frame(brier_score = brier_score, type = this_type))
  }
  res
}


