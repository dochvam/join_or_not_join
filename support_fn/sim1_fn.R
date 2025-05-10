library(spOccupancy)
library(tidyverse)

sim1_params <- c("n1", "n2", "J", "xi", "beta0", "beta1", "alpha0_1", 
                 "alpha1_1", "alpha0_2", "alpha1_2", "zeta", "sigma")
sim1_defaults <- c(n1 = 50,
                   n2 = 100,
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

split_mtx <- function(mtx, n) {
  split(seq_len(nrow(mtx)), rep(1:n, each = nrow(mtx)/n)) |>
    lapply(function(i) mtx[i, , drop = FALSE])
}

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
                         sigma = 0,
                         ppc = F,
                         seed) {
  set.seed(seed + (iter * 17))

  ncv_sites <- 50
  
  dat <- simulate_data_sim1(n1, n2, J, xi, beta0, beta1, alpha0_1, alpha1_1, 
                            alpha0_2, alpha1_2, zeta, sigma)
  dat_OOS <- simulate_data_sim1(n1 = ncv_sites, n2 = ncv_sites, J, xi = 1, beta0,
                                beta1, alpha0_1, alpha1_1, 
                                alpha0_2, alpha1_2, zeta, sigma)
  
  dat_one <- list(
    y = dat$y[[1]],
    occ.covs = dat$occ.covs[1:n1, ],
    det.covs = dat$det.covs[[1]]
  )
  dat_two <- list(
    y = dat$y[[2]],
    occ.covs = dat$occ.covs[(n1+1):(n1 + n2), ],
    det.covs = dat$det.covs[[2]]
  )
  
  fit_joint <- spOccupancy::intPGOcc(~x, list(~w_1, ~w_2), data = dat, 
                                     n.burn = 2500, n.samples = 5000,
                                     n.chains = 2,
                                     verbose = F)
  # browser()
  
  fit_one <- spOccupancy::PGOcc(~x, ~w_1, data = dat_one, 
                                  n.burn = 2500, n.samples = 5000,
                                  n.chains = 2, verbose = F)
  fit_two <- spOccupancy::PGOcc(~x, ~w_2, data = dat_two, 
                                  n.burn = 2500, n.samples = 5000,
                                  n.chains = 2, verbose = F)
  
  if (ppc) {
    ppc_res_ft <- ppcOcc(fit_joint, "freeman-tukey", 1)
    ppc_res_cs <- ppcOcc(fit_joint, "chi-squared", 1)
    ppc_res_dev <- ppcOcc_Deviance(object = fit_joint)
    
    ppc_result <- data.frame(
      type = c("D1", "D2"),
      ppc_pval_ft = c(
        mean(ppc_res_ft$fit.y.rep[[1]] > ppc_res_ft$fit.y[[1]]),
        mean(ppc_res_ft$fit.y.rep[[2]] > ppc_res_ft$fit.y[[2]])
      ),
      ppc_pval_cs = c(
        mean(ppc_res_cs$fit.y.rep[[1]] > ppc_res_cs$fit.y[[1]]),
        mean(ppc_res_cs$fit.y.rep[[2]] > ppc_res_cs$fit.y[[2]])
      ),
      ppc_pval_dev = c(
        mean(ppc_res_dev$fit.y.rep[[1]] > ppc_res_dev$fit.y[[1]]),
        mean(ppc_res_dev$fit.y.rep[[2]] > ppc_res_dev$fit.y[[2]])
      )
    )
    
  } else {
    ppc_result <- data.frame(
      type = c("D1", "D2"),
      ppc_pval_ft = c(NA, NA),
      ppc_pval_cs = c(NA, NA),
      ppc_pval_dev = c(NA, NA)
    )
  }
  
  estimation_result <- bind_rows(
    data.frame(
      param = c("beta0", "beta1", "alpha0_1", "alpha1_1", "alpha0_2", "alpha1_2"),
      mean = colMeans(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:4])),
      sd   = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:4]), 2, sd),
      Q025 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:4]), 2, quantile, probs = 0.025),
      Q975 = apply(cbind(fit_joint$beta.samples, fit_joint$alpha.samples[, 1:4]), 2, quantile, probs = 0.975)
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
      mutate(type = "one"),
    data.frame(
      param = c("beta0", "beta1", "alpha0_2", "alpha1_2"),
      mean = colMeans(cbind(fit_two$beta.samples, fit_two$alpha.samples)),
      sd   = apply(cbind(fit_two$beta.samples, fit_two$alpha.samples), 2, sd),
      Q025 = apply(cbind(fit_two$beta.samples, fit_two$alpha.samples), 2, quantile, probs = 0.025),
      Q975 = apply(cbind(fit_two$beta.samples, fit_two$alpha.samples), 2, quantile, probs = 0.975)
    ) %>% 
      left_join(dat$true_values, by = "param") %>% 
      mutate(type = "two")
  )
  
  cv_result_d1 <- calculate_CV_sim1(ncv_sites, ests = estimation_result, 
                                    dat_OOS = dat_OOS, d = 1)  
  cv_result_d2 <- calculate_CV_sim1(ncv_sites, ests = estimation_result, 
                                    dat_OOS = dat_OOS, d = 2)

  runtime_result <- data.frame(
    type = c("joint", "one"),
    value = unname(c(fit_joint$run.time[3], fit_one$run.time[3])),
    min_ESS = c(min(unlist(fit_joint$ESS)), min(unlist(fit_one$ESS)))
  ) %>% 
    mutate(min_ESS_persec = min_ESS / value)
  
  return(list(estimation_result = estimation_result %>% mutate(iter = iter),
              cv_result = bind_rows(cv_result_d1, cv_result_d2) %>% mutate(iter = iter), 
              runtime_result = runtime_result %>% mutate(iter = iter),
              ppc_result = ppc_result %>% mutate(iter = iter)))
}


run_many_sim1 <- function(specs_df_onerow, nsim, cl = NULL, ppc = FALSE) {
  stopifnot(all(colnames(specs_df_onerow) %in% c(names(sim1_defaults), "scenario", "seed")))
  
  sim1_values <- sim1_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim1_params) {
      sim1_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  if (!is.null(cl)) {
    result_list <- parLapply(cl, X = 1:nsim, fun = run_one_sim1,
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
              sigma = sim1_values["sigma"],
              ppc = ppc)
    
  } else {
    result_list <- list()
    for (i in 1:nsim) {
      result_list[[i]] <- run_one_sim1(iter = i,
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
                                       sigma = sim1_values["sigma"],
                                       ppc = ppc) 
    }
  }
  
  estimation_result_list <- cv_result_list <- runtime_result_list <- ppc_result_list <- list()
  for (i in 1:nsim) {
    estimation_result_list[[i]] <- result_list[[i]]$estimation_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    cv_result_list[[i]] <- result_list[[i]]$cv_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    runtime_result_list[[i]]<- result_list[[i]]$runtime_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    ppc_result_list[[i]]<- result_list[[i]]$ppc_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
  }
  
  return(list(
    estimation_result = bind_rows(estimation_result_list),
    cv_result = bind_rows(cv_result_list),
    runtime_result = bind_rows(runtime_result_list),
    ppc_result = bind_rows(ppc_result_list)
  ))
}


run_one_sim1_randomspecs <- function(iter, main_seed) {
  set.seed(main_seed + (17 * iter))
  
  n1 <- sample(20:200, 1)
  n2 <- sample(10 * 2:200, 1)
  J <- sample(2:6, 1)
  xi <- runif(1, 0.1, 1)
  beta0 <- runif(1, -2, 2)
  beta1 <- runif(1, -2, 2)
  alpha0_1 <- runif(1, -2, 2)
  alpha0_2 <- runif(1, -2, 2)
  alpha1_1 <- runif(1, -2, 2)
  alpha1_2 <- runif(1, -2, 2)
  zeta <- runif(1, -2, 2)
  sigma <- runif(1, 0, 5)
  
  result_list <- 
    run_one_sim1(iter, scenario = NA, n1, n2, J, xi, beta0, beta1, alpha0_1, 
                 alpha1_1, alpha0_2, alpha1_2, zeta, sigma, seed = main_seed)
  
  result_list$estimation_result <- result_list$estimation_result %>% 
    mutate(MSE = sd^2 + (mean - truth)^2,
           abs_error = mean - truth,
           covered = Q975 > truth & Q025 < truth)

  result_list$specs_df <- data.frame(
    n1 = n1,
    n2 = n2,
    J = J,
    xi = xi,
    beta0 = beta0,
    beta1 = beta1,
    alpha0_1 = alpha0_1,
    alpha0_2 = alpha0_2,
    alpha1_1 = alpha1_1,
    alpha1_2 = alpha1_2,
    zeta = zeta,
    sigma = sigma,
    iter = iter,
    main_seed = main_seed
  )
  result_list
}

run_many_sim_randomspecs <- function(nsim, cl = NULL) {
  main_seed <- floor(runif(1, 0, 1) * 1000000)
  
  if (!is.null(cl)) {
    result_list <- parLapply(cl, 1:nsim, fun = run_one_sim1_randomspecs, main_seed = main_seed)
  } else {
    result_list <- lapply(1:nsim, fun = run_one_sim1_randomspecs, main_seed = main_seed)
  }
  
  estimation_result_list <- cv_result_list <- runtime_result_list <- specs_df_list <- list()
  for (i in 1:nsim) {
    estimation_result_list[[i]] <- result_list[[i]]$estimation_result
    cv_result_list[[i]] <- result_list[[i]]$cv_result
    runtime_result_list[[i]] <- result_list[[i]]$runtime_result
    specs_df_list[[i]] <- result_list[[i]]$specs_df
  }
  
  return(list(
    estimation_result = bind_rows(estimation_result_list),
    cv_result = bind_rows(cv_result_list),
    runtime_result = bind_rows(runtime_result_list),
    specs_df = bind_rows(specs_df_list)
  ))
}


calculate_CV_sim1 <- function(ncv_sites, ests, dat_OOS, d) {
  res <- data.frame()
  
  if (d == 1) {
    for (this_type in c("joint", "one")) {
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
  } else if (d == 2) {
    for (this_type in c("joint", "two")) {
      this_ests <- ests %>% filter(type == this_type)
      pred_psi <- plogis(this_ests$mean[this_ests$param == "beta0"] +
                           this_ests$mean[this_ests$param == "beta1"] * dat_OOS$occ.covs[1:ncv_sites, "x"])
      pred_p <- plogis(this_ests$mean[this_ests$param == "alpha0_2"] +
                         this_ests$mean[this_ests$param == "alpha1_2"] * dat_OOS$det.covs[[2]]$w_2[1:ncv_sites, ])
      
      predicted_value <- matrix(pred_psi, nrow = ncv_sites, ncol = ncol(pred_p)) * 
        pred_p
      
      brier_score <- mean((predicted_value - dat_OOS$y[[2]][1:ncv_sites, ])^2)
      
      res <- bind_rows(res, data.frame(brier_score = brier_score, type = this_type))
    }
  } else {
    stop("d must be 1 or 2")
  }
  res %>% mutate(d=d) %>% return()
}



ppcOcc_Deviance <- function(object) {
  
  fit_stat <- "Deviance"
  call <- match.call()
  
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  
  out <- list()
  
  y <- object$y
  
  n.data <- length(y)
  sites <- object$sites
  X.p <- object$X.p
  p.det.long <- sapply(X.p, function(a) dim(a)[2])
  J.long <- sapply(y, nrow)
  occ.prob.all <- object$psi.samples
  fitted.out <- spOccupancy:::fitted.intPGOcc(object)
  y.rep.all <- fitted.out$y.rep.samples
  det.prob.all <- fitted.out$p.samples
  fit.y.list <- list()
  fit.y.rep.list <- list()
  fit.y.group.quants.list <- list()
  fit.y.rep.group.quants.list <- list()
  
  for (q in 1:n.data) {
    y.rep.samples <- y.rep.all[[q]] 
    
    z.samples <- object$z.samples[, sites[[q]], drop = FALSE]
    alpha.indx.r <- unlist(sapply(1:n.data, function(a) rep(a, p.det.long[a])))
    alpha.samples <- object$alpha.samples[, alpha.indx.r == q, drop = FALSE]
    # Get detection probability
    det.prob <- det.prob.all[[q]]
    n.samples <- object$n.post * object$n.chains
    fit.y <- rep(NA, n.samples)
    fit.y.rep <- rep(NA, n.samples)
    e <- 0.0001
    # Do the stuff 
    
    ### The stuff I modified is mostly here vvvvvvv ###
    fit.y.list[[q]] <- fit.y.rep.list[[q]] <- numeric(n.samples)
    for (i in 1:n.samples) {
      lprob_real <- lprob_rep <- numeric(length(nrow(y[[q]])))
      for (s in 1:nrow(y[[q]])) {
        
        ## Calculate log prob of the observed data given p, psi
        this_y <- y[[q]][s, ]
        this_detprob <-  det.prob[i, s, !is.na(this_y)]
        this_y <- this_y[!is.na(this_y)]
        
        if (sum(this_y, na.rm = T) > 0) {
          lprob_real[s] <- log(occ.prob.all[i, sites[[q]][s]]) +
            log(prod(this_y * this_detprob + (1 - this_y) * (1 - this_detprob)))
        } else {
          lprob_real[s] <- log(
            (1 - occ.prob.all[i, sites[[q]][s]]) +
              (occ.prob.all[i, sites[[q]][s]] * prod(1 - this_detprob))
          )
        }
        
        ## Calculate log prob of the sample data given p, psi
        this_yrep <- y.rep.all[[q]][i, s, ]
        this_yrep <- this_yrep[!is.na(this_yrep)]
        if (sum(this_yrep, na.rm = T) > 0) {
          lprob_rep[s] <- log(occ.prob.all[i, sites[[q]][s]]) +
            log(prod(this_yrep * this_detprob + (1-this_yrep) * (1-this_detprob)))
        } else {
          lprob_rep[s] <- log(
            (1 - occ.prob.all[i, sites[[q]][s]]) +
              (occ.prob.all[i, sites[[q]][s]] * prod(1 - this_detprob))
          )
        }
      }
      
      # Deviance = -2 * sum(ll), by tradition
      fit.y.list[[q]][i] <- -2 * sum(lprob_real)
      fit.y.rep.list[[q]][i] <- -2 * sum(lprob_rep)
    }
    ### The stuff I modified is mostly here ^^^^^^^ ###
  }
  
  out$fit.y <- fit.y.list
  out$fit.y.rep <- fit.y.rep.list
  
  out$fit.stat <- "Deviance"
  out$class <- class(object)
  out$call <- call
  out$n.samples <- object$n.samples
  out$n.burn <- object$n.burn
  out$n.thin <- object$n.thin
  out$n.post <- object$n.post
  out$n.chains <- object$n.chains
  
  class(out) <- 'ppcOcc'
  
  return(out)
}
