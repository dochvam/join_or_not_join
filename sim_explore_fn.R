library(tidyverse)
library(unmarked)
library(nimbleEcology)


#### Scenario codes:
#' 1 = Matches the model to estimate
#' 2 = PO effort is observed with random error
#' 3 = PO effort is observed with error confounded with an intensity covar.
#' 4 = PO data are totally unrelated to covariates
#' 5 = PO data have different relationships to covariates than camera data

#### Function to simulate data ####
simulate_joint_data <- function(sim_scenario,
                                dataset_ID,
                                holdout_frac,
                                grid_dim = 100,
                                ncamera = 50,
                                nreps_percam = 3,
                                PO_agg_factor = 4,
                                PO_avg_effort = 50,
                                PO_zi = 0.3,
                                intensity_intercept = -0.3,
                                intensity_b1 = 1,
                                intensity_b2 = 1,
                                intensity_b3 = 0,
                                camera_detInt = 0.5,
                                camera_detb1 = 0.3,
                                PO_b1 = -0.5,
                                PO_b2 = -1,
                                PO_b3 = -1,
                                PO_theta0 = -6,
                                PO_theta1 = 1,
                                PO_overdisp = 0.5,
                                PO_error_sd = 0.5) {

  true_params <- data.frame(
    param = c("intensity_int",
              "intensity_beta[1]",
              "intensity_beta[2]",
              "intensity_beta[3]",
              "cam_p_int",
              "det_beta[1]",
              "det_beta[2]",
              "theta0",
              "theta1",
              "overdisp",
              "log_overdisp",
              "PO_error_sd"),
    value = c(intensity_intercept,
              intensity_b1,
              intensity_b2,
              intensity_b3,
              camera_detInt,
              camera_detb1,
              0,
              PO_theta0,
              PO_theta1,
              PO_overdisp,
              log(PO_overdisp),
              PO_error_sd)
  )
  
  stopifnot(sim_scenario %in% 1:5)
  
  
  # True intensity process of interest:
  process_df <- expand.grid(x = 1:grid_dim, y = 1:grid_dim) %>% 
    mutate(PO_agg_x = floor((x - 1) / PO_agg_factor) + 1,
           PO_agg_y = floor((y - 1) / PO_agg_factor) + 1) %>% 
    mutate(cov1 = rnorm(grid_dim^2), 
           cov2 = rnorm(grid_dim^2),
           cov3 = as.numeric(scale(y))) %>% 
    mutate(
      log_intensity = intensity_intercept + 
        intensity_b1 * cov1 +
        intensity_b2 * cov2 +
        intensity_b3 * cov3
    )
  
  
  # Simulate camera trap locations, occupancies, and detection histories
  camera_df <- data.frame(
    cam.x = runif(ncamera, 0, grid_dim),
    cam.y = runif(ncamera, 0, grid_dim)
  ) %>% 
    mutate(x = ceiling(cam.x), y = ceiling(cam.y),
           covC = rnorm(ncamera)) %>% 
    left_join(process_df, by = c("x", "y")) %>% 
    mutate(psi = nimble::icloglog(log_intensity),
           p = nimble::icloglog(camera_detInt + covC * camera_detb1)) %>% 
    mutate(z = rbinom(ncamera, 1, psi))
  
  visit_cols <- paste0("V", 1:nreps_percam)
  for (i in 1:nreps_percam) {
    camera_df[[visit_cols[i]]] <- 
      rbinom(ncamera, 1, prob = camera_df$p * camera_df$z)
  }
  
  
  
  #### Simulate PO-like count data, with effort
  
  if (sim_scenario == 1) {
    po_df <- process_df %>% 
      group_by(PO_agg_x, PO_agg_y) %>% 
      summarize(mu = exp(PO_theta0 + PO_theta1 * log(sum(exp(log_intensity)))),
                .groups = "drop")
  } else if (sim_scenario == 2) {
    po_df <- process_df %>% 
      group_by(PO_agg_x, PO_agg_y) %>% 
      summarize(log_sum_exp_llam = log(sum(exp(log_intensity))),
                .groups = "drop") %>% 
      mutate(mu = exp(PO_theta0 + PO_theta1 * log_sum_exp_llam +
                        PO_error_sd * rnorm(n(), 0, 1)))
  } else if (sim_scenario == 3) {
    po_df <- process_df %>% 
      group_by(PO_agg_x, PO_agg_y) %>% 
      summarize(log_sum_exp_llam = log(sum(exp(log_intensity))),
                mean_cov1 = mean(cov1),
                .groups = "drop") %>% 
      mutate(mu = exp(PO_theta0 + PO_theta1 * log_sum_exp_llam +
                        PO_error_sd * rnorm(n(), mean_cov1, 0.3)))
    
  } else if (sim_scenario == 4) {
    po_df <- process_df %>% 
      distinct(PO_agg_x, PO_agg_y) %>% 
      mutate(mu = exp(PO_theta0 + PO_error_sd * rnorm(n(), 0, 1)))
  } else if (sim_scenario == 5) {
    process_df$log_intensity_2 <- intensity_intercept + 
      PO_b1 * process_df$cov1 +
      PO_b2 * process_df$cov2 +
      PO_b3 * process_df$cov3
    
    po_df <- process_df %>% 
      group_by(PO_agg_x, PO_agg_y) %>% 
      summarize(mu = exp(PO_theta0 + PO_theta1 * log(sum(exp(log_intensity_2)))),
                .groups = "drop")
  }
  po_df <- po_df %>% 
    mutate(effort = rpois(nrow(po_df), PO_avg_effort) * 
           rbinom(nrow(po_df), 1, 1 - PO_zi)) %>% 
    mutate(y = rnbinom(nrow(po_df), mu = mu * effort, size = 1 / PO_overdisp))
  
  
  #### Put together cell-level info ####
  
  cell_dat <- process_df %>% 
    arrange(PO_agg_x, PO_agg_y) %>% 
    group_by(PO_agg_x, PO_agg_y) %>% 
    mutate(agg_cell_ID = cur_group_id()) %>% 
    ungroup()
  PO_cell_start <- PO_cell_end <- numeric(length(unique(cell_dat$agg_cell_ID)))
  
  for (i in 1:length(PO_cell_start)) {
    PO_cell_start[i] <- min(which(cell_dat$agg_cell_ID == i))
    PO_cell_end[i]   <- max(which(cell_dat$agg_cell_ID == i))
  }
  
  po_df <- left_join(po_df, distinct(cell_dat, PO_agg_x, PO_agg_y, agg_cell_ID),
                     by = c("PO_agg_x", "PO_agg_y")) %>% 
    arrange(agg_cell_ID)
  
  
  
  # Separate the holdout data
  datlist_holdout <- list()
  
  holdout_inds_camera <- sample(1:nrow(camera_df), 
                                size = nrow(camera_df) * holdout_frac, 
                                replace = F)
  
  datlist_holdout <- list()
  datlist_holdout$camera_df <- camera_df[holdout_inds_camera, ]
  camera_df <- camera_df[-holdout_inds_camera, ]
  
  holdout_inds_PO <- sample(1:nrow(po_df), 
                            size = nrow(po_df) * holdout_frac, 
                            replace = F)
  datlist_holdout$po_df <- po_df[holdout_inds_PO, ]
  datlist_holdout$PO_cell_start <- PO_cell_start[holdout_inds_PO]
  datlist_holdout$PO_cell_end <- PO_cell_end[holdout_inds_PO]
  
  po_df         <- po_df[-holdout_inds_PO, ]
  PO_cell_start <- PO_cell_start[-holdout_inds_PO]
  PO_cell_end   <- PO_cell_end[-holdout_inds_PO]
  
  
  return(list(
    camera_df = camera_df,
    cell_dat = cell_dat,
    po_df = po_df,
    PO_cell_start = PO_cell_start,
    PO_cell_end = PO_cell_end,
    true_params = true_params,
    visit_cols = visit_cols,
    datlist_holdout = datlist_holdout,
    dataset_ID = dataset_ID
  ))
}




simulate_joint_data_v2 <- function(sim_scenario,
                                dataset_ID,
                                grid_dim = 100,
                                ncamera = 50,
                                ncamera_holdout = 25,
                                nreps_percam = 3,
                                PO_agg_factor = 4,
                                PO_avg_effort = 50,
                                PO_zi = 0.3,
                                intensity_intercept = -0.3,
                                intensity_b1 = 1,
                                camera_detInt = 0.5,
                                camera_detb1 = 0.3,
                                cam_pct_xcover = 1,
                                PO_bias = 0,
                                PO_theta0 = -6,
                                PO_theta1 = 1,
                                PO_overdisp = 0.5,
                                PO_error_sd = 0.5) {

  true_params <- data.frame(
    param = c("intensity_int",
              "intensity_beta",
              "cam_p_int",
              "det_beta",
              "theta0",
              "theta1",
              "overdisp",
              "log_overdisp",
              "PO_error_sd"),
    value = c(intensity_intercept,
              intensity_b1,
              camera_detInt,
              camera_detb1,
              PO_theta0,
              PO_theta1,
              PO_overdisp,
              log(PO_overdisp),
              PO_error_sd)
  )
  
  PO_b1 <- intensity_b1 + PO_bias
  
  # True intensity process of interest:
  process_df <- expand.grid(x = 1:grid_dim, y = 1:grid_dim) %>% 
    mutate(PO_agg_x = floor((x - 1) / PO_agg_factor) + 1,
           PO_agg_y = floor((y - 1) / PO_agg_factor) + 1) %>% 
    mutate(cov1 = runif(grid_dim^2, -1, 1),
           covC = rnorm(grid_dim^2, 0, 1)) %>% 
    mutate(
      log_intensity = intensity_intercept + 
        intensity_b1 * cov1
    )
  
  possible_cells <- process_df %>% 
    filter(cov1 < -1 + 2*cam_pct_xcover)
  
  # Simulate camera trap locations, occupancies, and detection histories
  camera_df <- possible_cells %>% 
    sample_n(ncamera, replace = TRUE) %>% 
    mutate(camera_id = row_number()) %>% 
    mutate(psi = nimble::icloglog(log_intensity),
           p = nimble::icloglog(camera_detInt + covC * camera_detb1)) %>% 
    mutate(z = rbinom(ncamera, 1, psi))
  
  visit_cols <- paste0("V", 1:nreps_percam)
  for (i in 1:nreps_percam) {
    camera_df[[visit_cols[i]]] <- 
      rbinom(ncamera, 1, prob = camera_df$p * camera_df$z)
  }
  
  
  
  #### Simulate PO-like count data, with effort
  process_df$log_intensity_2 <- intensity_intercept + 
    PO_b1 * process_df$cov1
  
  po_df <- process_df %>% 
    group_by(PO_agg_x, PO_agg_y) %>% 
    summarize(log_sum_exp_llam = log(sum(exp(log_intensity_2))),
              .groups = "drop") %>% 
    mutate(mu = exp(PO_theta0 + PO_theta1 * log_sum_exp_llam +
                        PO_error_sd * rnorm(n(), 0, 1)))
    
    
  po_df <- po_df %>% 
    mutate(effort = rpois(nrow(po_df), PO_avg_effort) * 
           rbinom(nrow(po_df), 1, 1 - PO_zi)) %>% 
    mutate(y = rnbinom(nrow(po_df), mu = mu * effort, size = 1 / PO_overdisp))
  
  
  #### Put together cell-level info ####
  
  cell_dat <- process_df %>% 
    arrange(PO_agg_x, PO_agg_y) %>% 
    group_by(PO_agg_x, PO_agg_y) %>% 
    mutate(agg_cell_ID = cur_group_id()) %>% 
    ungroup()
  PO_cell_start <- PO_cell_end <- numeric(length(unique(cell_dat$agg_cell_ID)))
  
  for (i in 1:length(PO_cell_start)) {
    PO_cell_start[i] <- min(which(cell_dat$agg_cell_ID == i))
    PO_cell_end[i]   <- max(which(cell_dat$agg_cell_ID == i))
  }
  
  po_df <- left_join(po_df, distinct(cell_dat, PO_agg_x, PO_agg_y, agg_cell_ID),
                     by = c("PO_agg_x", "PO_agg_y")) %>% 
    arrange(agg_cell_ID)
  
  
  
  # Separate the holdout data
  datlist_holdout <- list()
  camera_df_holdout <- process_df %>% 
    sample_n(size = ncamera_holdout) %>% 
    mutate(psi = nimble::icloglog(log_intensity),
           p = nimble::icloglog(camera_detInt + covC * camera_detb1)) %>% 
    mutate(z = rbinom(ncamera_holdout, 1, psi))
  
  visit_cols <- paste0("V", 1:nreps_percam)
  for (i in 1:nreps_percam) {
    camera_df_holdout[[visit_cols[i]]] <- 
      rbinom(ncamera_holdout, 1, prob = camera_df_holdout$p * camera_df_holdout$z)
  }
  
  datlist_holdout$camera_df <- camera_df_holdout
  
  return(list(
    camera_df = camera_df,
    cell_dat = cell_dat,
    po_df = po_df,
    PO_cell_start = PO_cell_start,
    PO_cell_end = PO_cell_end,
    true_params = true_params,
    visit_cols = visit_cols,
    datlist_holdout = datlist_holdout,
    dataset_ID = dataset_ID,
    sim_scenario = sim_scenario
  ))
}

joint_inits <- function(modtype, ncovInt, ncovP, MLE = FALSE) {
  if (modtype == "joint") {
    inits_list <- list(
      intensity_int = 0,
      intensity_beta = rnorm(ncovInt),
      cam_p_int = 0,
      det_beta = rnorm(ncovP),
      theta0 = -6, theta1 = 1, overdisp = 1
    )
  } else if (modtype == "camera_only") {
    inits_list <- list(
      intensity_int = 0,
      intensity_beta = rnorm(ncovInt),
      cam_p_int = 0,
      det_beta = rnorm(ncovP)
    )
  } else {
    inits_list <- list(
      intensity_int = 0,
      intensity_beta = rnorm(ncovInt),
      theta0 = -6, theta1 = 1, overdisp = 1
    )
  }
  
  if (MLE) {
    inits_list[["overdisp"]] <- NULL
    inits_list$log_overdisp <- 0
  }
  
  inits_list
}

# Set this up so that we can use the same compiled model with different
# simulated datasets, as long as they have the same data dimensions
# fit_nimblemodel <- function(datlist, nsim, modtype, holdout_frac = 0,
#                             ni = 5000, nb = 2000, nc = 2, nt = 1) {
#   
#   stopifnot(modtype %in% c("joint", "camera_only", "PO_only"))
#   
#   if (nsim == 1 && length(datlist) > 1) {
#     datlist <- list(datlist)
#   }
#   
#   datlist_holdout <- list()
#   
#   for (i in 1:nsim) {
#     holdout_inds <- sample(1:nrow(datlist[[i]]$camera_df), 
#                            size = nrow(datlist[[i]]$camera_df) * holdout_frac, 
#                            replace = F)
#     
#     datlist_holdout[[i]] <- list()
#     datlist_holdout[[i]]$camera_df <- datlist[[i]]$camera_df[holdout_inds, ]
#     
#     datlist[[i]]$camera_df <- datlist[[i]]$camera_df[-holdout_inds, ]
#   }
#   
#   modcode <- nimbleCode({
#     if (modtype == "joint" || modtype == "camera_only") {
#       for (i in 1:ncam) {
#         cloglog(psi[i]) <- intensity_int + inprod(intensity_beta[1:ncovInt], cam_dat_int[i, 1:ncovInt])
#         cloglog(p[i]) <- cam_p_int + inprod(det_beta[1:ncovP], cam_dat_p[i, 1:ncovP])
#         
#         cam_y[i, 1:nrep] ~ dOcc_s(probOcc = psi[i], probDetect = p[i], len = nrep)
#       }
#       
#       # Camera-specific priors 
#       cam_p_int ~ dnorm(0, sd = 5)
#       for (i in 1:ncovP) {
#         det_beta[i] ~ dnorm(0, sd = 5)
#       }
#       
#       if (ncam_holdout > 0) {
#         for (i in 1:ncam_holdout) {
#           cloglog(psi_holdout[i]) <- intensity_int + 
#             inprod(intensity_beta[1:ncovInt], cam_dat_int_holdout[i, 1:ncovInt])
#           
#           cloglog(p_holdout[i]) <- cam_p_int + 
#             inprod(det_beta[1:ncovP], cam_dat_p_holdout[i, 1:ncovP])
#           
#           pp_brier_camera_bysite[i] <- calcBrierScore_oneSite(
#             y = cam_y_holdout[i, 1:nrep],
#             psi = psi_holdout[i],
#             p = p_holdout[i],
#             len = nrep
#           )
#         }
#         pp_brier_camera <- sum(pp_brier_camera_bysite[1:ncam_holdout])
#       }
#     }
#     
#     if (modtype == "joint" || modtype == "PO_only") {
#       for (i in 1:nPOcell) {
#         log(mu[i]) <- theta0 + theta1 * 
#           log(sum(exp(
#             predicted_log_intensity[cell_start[i]:cell_end[i]]
#           )))
#         
#         PO_y[i] ~ dnbinom(size = 1 / overdisp,
#                           prob = 1 / (1 + overdisp * effort[i] * mu[i]))
#       }
#       
#       # PO-specific priors
#       theta0 ~ dnorm(-5, sd = 50)
#       theta1 ~ dnorm(1, sd = 2.75)
#       overdisp ~ dgamma(shape = 1, scale = 5)
#     }
#     
#     # Intensity priors
#     intensity_int ~ dnorm(0, sd = 20)
#     for (i in 1:ncovInt) {
#       intensity_beta[i] ~ dnorm(0, sd = 5)
#     }
#     
#     
#     # Predict and calculate RMSE
#     for (i in 1:nCell) {
#       predicted_log_intensity[i] <- intensity_int + 
#         inprod(intensity_beta[1:ncovInt], cell_dat_int[i, 1:ncovInt])
#       cell_sqerr[i] <- (predicted_log_intensity[i] - true_log_intensity[i])^2
#     }
#     RMSE <- sqrt(mean(cell_sqerr[1:nCell]))
#   })
#   
#   p_covs <- c("cov1", "covC")
#   int_covs <- c("cov1", "cov2", "cov3")
#   
#   ### Build the model
#   mod <- nimbleModel(
#     code = modcode, 
#     constants = list(
#       modtype = modtype,
#       nrep = length(datlist[[1]]$visit_cols),
#       nCell = nrow(datlist[[1]]$cell_dat),
#       nPOcell = length(datlist[[1]]$PO_cell_start),
#       ncam = nrow(datlist[[1]]$camera_df),
#       ncam_holdout = nrow(datlist_holdout[[1]]$camera_df),
#       cell_start = datlist[[1]]$PO_cell_start,
#       cell_end = datlist[[1]]$PO_cell_end,
#       ncovP = 2,
#       ncovInt = 3
#     ),
#     data = list(
#       cell_dat_int = as.matrix(datlist[[1]]$cell_dat[, int_covs]),
#       
#       cam_dat_int  = as.matrix(datlist[[1]]$camera_df[, int_covs]),
#       cam_dat_p = as.matrix(datlist[[1]]$camera_df[, p_covs]),
#       cam_y = as.matrix(datlist[[1]]$camera_df[, datlist[[1]]$visit_cols]),
#       
#       cam_dat_int_holdout  = as.matrix(datlist_holdout[[1]]$camera_df[, int_covs]),
#       cam_dat_p_holdout = as.matrix(datlist_holdout[[1]]$camera_df[, p_covs]),
#       cam_y_holdout = as.matrix(datlist_holdout[[1]]$camera_df[, datlist[[1]]$visit_cols]),
#       
#       true_log_intensity = as.numeric(datlist[[1]]$cell_dat$log_intensity),
#       
#       PO_y  = as.numeric(datlist[[1]]$po_df$y),
#       effort = as.numeric(datlist[[1]]$po_df$effort)
#     ),
#     inits = joint_inits(modtype, 3, 2)
#   )
#   
#   mcmcConf <- configureMCMC(mod)
#   mcmcConf$addMonitors("RMSE")
#   
#   if (holdout_frac > 0) {
#     if (modtype %in% c("joint", "camera_only")) {
#       mcmcConf$addMonitors("pp_brier_camera")
#     }
#     # if (modtype %in% c("joint", "PO_only")) {
#     #   mcmcConf$addMonitors("pp_brier_PO")
#     # }
#   }
#   
#   if (modtype == "joint") {
#     mcmcConf$removeSampler(c("theta0", "theta1"))
#     mcmcConf$addSampler(target = c("theta0", "theta1"), type = "AF_slice")
#   }
#   
#   mcmc <- buildMCMC(mcmcConf)
#   complist <- compileNimble(mod, mcmc)
#   
#   summary_list <- list()
#   for (i in 1:nsim) {
#   # Update data in the model
#     complist$mod$setData("effort" = datlist[[i]]$po_df$effort)
#     complist$mod$setData("PO_y" = datlist[[i]]$po_df$y)
#     complist$mod$setData("cell_dat_int" = as.matrix(datlist[[i]]$cell_dat[, int_covs]))
#     complist$mod$setData("cam_y" = datlist[[i]]$camera_df[, datlist[[i]]$visit_cols])
#     complist$mod$setData("cam_dat_int" = as.matrix(datlist[[i]]$camera_df[, int_covs]))
#     complist$mod$setData("cam_dat_p" = as.matrix(datlist[[i]]$camera_df[, p_covs]))
#     complist$mod$setData("true_log_intensity" = as.numeric(datlist[[i]]$cell_dat$log_intensity))
#     
#     if (holdout_frac > 0) {
#       complist$mod$setData("cam_y_holdout" = datlist_holdout[[i]]$camera_df[, datlist[[i]]$visit_cols])
#       complist$mod$setData("cam_dat_int_holdout" = as.matrix(datlist_holdout[[i]]$camera_df[, int_covs]))
#       complist$mod$setData("cam_dat_p_holdout" = as.matrix(datlist_holdout[[i]]$camera_df[, p_covs]))
#     }
#       
#   # run the MCMC
#     samples <- runMCMC(complist$mcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc,
#             inits = joint_inits(modtype, 3, 2), samplesAsCodaMCMC = TRUE)
#     
#     summary_list[[i]] <- MCMCvis::MCMCsummary(samples) %>% 
#       mutate(modtype = modtype, dataset = datlist[[i]]$dataset_ID)
#     summary_list[[i]]$param <- rownames(summary_list[[i]])
#     rownames(summary_list[[i]]) <- NULL
#   }
#   
#   return(bind_rows(
#     summary_list
#   ))
# }


fit_nimblemodel_MLE <- function(datlist, nsim = NULL, modtype, 
                                holdout_frac = 0, 
                                sim_scenario = NULL, progress = TRUE) {
  
  stopifnot(modtype %in% c("joint", "camera_only", "PO_only"))
  
  # if (is.null(sim_scenario)) sim_scenario <- datlist[[1]]$sim_scenario
  
  if (is.null(nsim)) {
    nsim <- length(datlist)
  }
  
  if (nsim == 1 && length(datlist) > 1) {
    datlist <- list(datlist)
  }
  
  datlist_holdout <- list()
  for (i in 1:nsim) {
    datlist_holdout[[i]] <- datlist[[i]]$datlist_holdout
  }
  
  #### Define the NIMBLE model to estimate ####
  modcode_MLE <- nimbleCode({
    if (modtype == "joint" || modtype == "camera_only") {
      for (i in 1:ncam) {
        cloglog(psi[i]) <- intensity_int + intensity_beta * cam_dat_int[i]
        cloglog(p[i]) <- cam_p_int + det_beta * cam_dat_p[i]
        
        cam_y[i, 1:nrep] ~ dOcc_s(probOcc = psi[i], probDetect = p[i], len = nrep)
      }
    }
    
    if (modtype == "joint" || modtype == "PO_only") {
      for (i in 1:nPOcell) {
        for (j in 1:repsPerPOCell) {
          predicted_log_intensity[i, j] <- intensity_int + 
            intensity_beta * cell_dat_int[(i-1)*repsPerPOCell + j]
        }
        
        log(mu[i]) <- theta0 + theta1 * 
          log(sum(exp(
            predicted_log_intensity[i, 1:repsPerPOCell]
          )))
        
        PO_y[i] ~ dnbinom(size = 1 / exp(log_overdisp),
                          prob = 1 / (1 + exp(log_overdisp) * effort[i] * mu[i]))
      }
    }
  })
  
  p_covs <- c("covC")
  int_covs <- c("cov1")
  
  # Make sure the order of cells matches PO data order
  order_vec <- datlist[[1]]$po_df$agg_cell_ID
  this_cell_dat <- datlist[[1]]$cell_dat[datlist[[1]]$cell_dat$agg_cell_ID %in% order_vec, ]
  stopifnot(all(unique(this_cell_dat$agg_cell_ID) == order_vec))

  ### Build the model
  mod <- nimbleModel(
    code = modcode_MLE, 
    constants = list(
      modtype = modtype,
      nrep = length(datlist[[1]]$visit_cols),
      # nCell = nrow(datlist[[1]]$cell_dat),
      nPOcell = length(datlist[[1]]$PO_cell_start),
      repsPerPOCell = sum(datlist[[1]]$cell_dat$PO_agg_x == 1 &
                          datlist[[1]]$cell_dat$PO_agg_y == 1),
      ncam = nrow(datlist[[1]]$camera_df)
    ),
    data = list(
      cell_dat_int =  as.numeric(unlist(this_cell_dat[, int_covs])),
      
      cam_dat_int  = as.numeric(unlist(datlist[[1]]$camera_df[, "cov1"])),
      cam_dat_p = as.numeric(unlist(datlist[[1]]$camera_df[, "covC"])),
      cam_y = as.matrix(datlist[[1]]$camera_df[, datlist[[1]]$visit_cols]),

      PO_y  = as.numeric(datlist[[1]]$po_df$y),
      effort = as.numeric(datlist[[1]]$po_df$effort)
    ),
    inits = joint_inits(modtype, 1, 1, MLE = TRUE)
  )
  cmod <- compileNimble(mod)
  
  if (modtype == "joint") {
    params <- c("intensity_int", "intensity_beta", 
                "cam_p_int", "det_beta",
                "theta0", "theta1", "log_overdisp") %>% 
      cmod$expandNodeNames(returnScalarComponents = TRUE)
    
  } else if (modtype == "camera_only") {
    params <- c("intensity_int", "intensity_beta", 
                "cam_p_int", "det_beta") %>% 
      cmod$expandNodeNames(returnScalarComponents = TRUE)
  } else if (modtype == "PO_only") {
    params <- c("intensity_int", "intensity_beta", "log_overdisp") %>% 
      cmod$expandNodeNames(returnScalarComponents = TRUE)
  }
  
  sc <- setAndCalculate(mod, params)
  csc <- compileNimble(sc)
  
  summary_list <- list()
  rmse_list <- list()
  CV_list <- list()
  
  if (progress) pb <- progress::progress_bar$new(total = nsim)
  for (i in 1:nsim) {
    if (progress) pb$tick()
    # if (i == 9) browser()
  # Update data in the model
    if (modtype %in% c("joint", "PO_only")) {
    # Make sure the order of cells matches PO data order
      order_vec <- datlist[[i]]$po_df$agg_cell_ID
      this_cell_dat <- datlist[[i]]$cell_dat[datlist[[i]]$cell_dat$agg_cell_ID %in% order_vec, ]
      stopifnot(all(unique(this_cell_dat$agg_cell_ID) == order_vec))
      cmod$setData("cell_dat_int" = as.numeric(unlist(this_cell_dat[, int_covs])))
      
      cmod$setData("effort" = datlist[[i]]$po_df$effort)
      cmod$setData("PO_y" = datlist[[i]]$po_df$y)
    }
    if (modtype %in% c("joint", "camera_only")) {
      cmod$setData("cam_y" = datlist[[i]]$camera_df[, datlist[[i]]$visit_cols])
      cmod$setData("cam_dat_int" = as.numeric(unlist(datlist[[i]]$camera_df[, int_covs])))
      cmod$setData("cam_dat_p" = as.numeric(unlist(datlist[[i]]$camera_df[, p_covs])))
    }

  # Gather true values
    true_values_raw <- datlist[[i]]$true_params %>% rename(true_value = value)

  # run the MLE algorithm
    this_fit <- get_MLE(csc, params)
    
    summary_list[[i]] <- this_fit %>% 
      mutate(dataset_ID = datlist[[i]]$dataset_ID,
             modtype = modtype,
             scenario = datlist[[i]]$sim_scenario) %>% 
      left_join(true_values_raw, by = "param") %>% 
      mutate(covered = est - 1.96*se < true_value & est + 1.96*se > true_value)
    # browser()
    rmse_list[[i]] <- data.frame(
      dataset_ID = datlist[[i]]$dataset_ID,
      modtype = modtype,
      scenario = datlist[[i]]$sim_scenario,
      goal = c("prediction", "inference_betas", "inference_all"),
      RMSE = c(
        calc_RMSE_prediction(this_fit, datlist[[i]]$cell_dat, int_covs),
        calc_RMSE_inference(this_fit, datlist[[i]]$true_params, betas_only = TRUE),
        calc_RMSE_inference(this_fit, datlist[[i]]$true_params, betas_only = FALSE)
      )
    )
    
    CV_list[[i]] <- calc_CV_error(estimate_df = this_fit, 
                                  holdout_dat = datlist_holdout[[i]], 
                                  cell_dat = datlist[[i]]$cell_dat,
                                  int_covs = int_covs, 
                                  p_covs = p_covs,
                                  visit_cols = datlist[[1]]$visit_cols,
                                  modtype = modtype) %>% 
      mutate(dataset_ID = datlist[[i]]$dataset_ID,
             modtype = modtype,
             scenario = datlist[[i]]$sim_scenario)
  }
  
  return(list(
    estimates = bind_rows(summary_list),
    rmse = bind_rows(rmse_list),
    CV = bind_rows(CV_list)
  ))
}


get_inits_MLE <- function(params) {
  lapply(params, function(x) {
    if (x == "theta0") {
      -6
    } else if (x == "theta1") {
      1
    } else if (x == "overdisp") {
      0.1
    } else {
      0
    }
  }) %>% 
    unlist()
}

get_MLE <- function(csc, params) {
  inits <- get_inits_MLE(params)
  
  re_optim_out <- NA
  tryCatch({
    optim_out <- optim(par = inits, fn = csc$run,
                       # method = "BFGS",
                       control = list(fnscale = -1, maxit = 10000),
                       hessian = TRUE)
    re_optim_out <- optim(par = optim_out$par, fn = csc$run,
                          # method = "BFGS",
                          control = list(fnscale = -1, maxit = 10000),
                          hessian = TRUE)
  }, error = function(e) { })
  
  if (!is.list(re_optim_out) && is.na(re_optim_out)) {
    return(data.frame(
      param = params,
      est = NA,
      se = NA
    ))
  }
  
  std_error <- rep(NA, length(re_optim_out$par))
  tryCatch({
    std_error <- sqrt(diag(solve(-re_optim_out$hess)))
  }, error = function(e) { })
  
  return(data.frame(
    param = params,
    est = re_optim_out$par,
    se = std_error
  ))
}



calcBrierScore_oneSite <- nimbleFunction(
  run = function(y = double(1),
                 psi = double(0),
                 p = double(0),
                 len = integer(0)){

    x <- numeric(len)
    for (i in 1:len) {
      x[i] <- y[i] * (1 - psi * p)^2 + (1 - y[i]) * (psi * p)^2
    }
    
    return(sum(x))
    returnType(double(0))
  })


calc_RMSE_prediction <- function(estimate_df, cell_dat, int_covs) {
  
  intensity_int <- estimate_df$est[estimate_df$param == "intensity_int"]
  intensity_beta <- estimate_df$est[grepl("intensity_beta", estimate_df$param)]
  nCell <- nrow(cell_dat)
  
  predicted_log_intensity <- as.numeric(intensity_int + 
    as.matrix(cell_dat[, int_covs]) * intensity_beta)
  
  cell_sqerr <- (predicted_log_intensity - cell_dat$log_intensity)^2

  RMSE <- sqrt(mean(cell_sqerr))
  
  return(RMSE)
}

calc_RMSE_inference <- function(estimate_df, true_values, betas_only = TRUE) {
  
  df <- left_join(estimate_df, true_values, by = "param")
  
  if (betas_only) {
    df <- df %>% filter(grepl("intensity_beta", param))
  }
  
  ((df$est - df$value)^2 + df$se^2) %>% 
    mean() %>% 
    sqrt() %>% 
    return()
}


calc_CV_error <- function(estimate_df, holdout_dat, cell_dat, 
                          int_covs, p_covs, visit_cols, modtype) {
  
  CV_df <- list()
  
  intensity_int <- estimate_df$est[estimate_df$param == "intensity_int"]
  intensity_beta <- estimate_df$est[grepl("intensity_beta", estimate_df$param)]
  
  if (modtype %in% c("joint", "camera_only")) {
    cam_p_int  <- estimate_df$est[grepl("cam_p_int", estimate_df$param)]
    det_beta <- estimate_df$est[grepl("det_beta", estimate_df$param)]
  }
  if (modtype == "joint") {
    log_overdisp <- estimate_df$est[grepl("log_overdisp", estimate_df$param)]
    theta0 <- estimate_df$est[grepl("theta0", estimate_df$param)]
    theta1 <- estimate_df$est[grepl("theta1", estimate_df$param)]
  }
  if (modtype == "PO_only") {
    log_overdisp <- estimate_df$est[grepl("log_overdisp", estimate_df$param)]
    theta0 <- 0
    theta1 <- 1
  }
  
  if (modtype %in% c("joint", "camera_only")) {
    predicted_occu <- 
      (intensity_int + as.matrix(holdout_dat$camera_df[, int_covs]) * intensity_beta) %>% 
      as.numeric() %>% 
      icloglog()
    
    predicted_det <- 
      (cam_p_int + as.matrix(holdout_dat$camera_df[, p_covs]) * det_beta) %>% 
      as.numeric() %>% 
      icloglog()
    
    brier_score <- numeric(nrow(holdout_dat$camera_df))
    for (i in 1:length(brier_score)) {
      y <- holdout_dat$camera_df[i, visit_cols] 
      brier_score[i] <- mean(unlist((predicted_occu[i] * predicted_det[i] - y)^2))
    }

    CV_df[[1]] <- data.frame(type = "camera", score = "brier", value = mean(brier_score))
  }
  
  # if (modtype %in% c("joint", "PO_only")) {
  #   PO_dat <- holdout_dat$po_df
  #   mu <- numeric(nrow(PO_dat))
  #   
  #   for (i in 1:nrow(PO_dat)) {
  #     predicted_log_intensity <- (intensity_int + 
  #       as.matrix(cell_dat[cell_dat$agg_cell_ID == PO_dat$agg_cell_ID[i], int_covs]) * intensity_beta) %>% 
  #       as.numeric()
  #     
  #     mu[i] <- exp(theta0 + theta1 * log(sum(exp(predicted_log_intensity))))
  #     
  #     
  #   }
  #   
  #   ll <- sum(dnbinom(x = PO_dat$y, size = 1 / exp(log_overdisp),
  #           prob = 1 / (1 + exp(log_overdisp) * PO_dat$effort * mu), log = TRUE))
  #   ll_saturated <- sum(dnbinom(x = PO_dat$y, size = 1 / exp(log_overdisp),
  #           prob = 1 / (1 + exp(log_overdisp) * PO_dat$y), log = TRUE))
  #   
  #   deviance <- 2 * (ll_saturated - ll)
  #   
  #   CV_df[[2]] <- data.frame(
  #     type = "PO",
  #     score = c("nll", "deviance"),
  #     value = c(-ll, deviance)
  #   )
  # }
  
  return(bind_rows(CV_df))
}

