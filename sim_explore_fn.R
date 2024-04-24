library(tidyverse)
library(unmarked)
library(nimbleEcology)


#### Scenario codes:
#' 1 = Matches the model to estimate
#' 2 = PO effort is observed with random error
#' 3 = PO effort is observed with error confounded with cov3
#' 4 = PO data have different relationships to covars.

#### Function to simulate data ####
simulate_joint_data <- function(sim_scenario,
                                dataset_ID,
                                grid_dim = 100,
                                ncamera = 50,
                                nreps_percam = 3,
                                PO_agg_factor = 4,
                                PO_avg_effort = 50,
                                PO_zi = 0.3,
                                intensity_intercept = -0.3,
                                intensity_b1 = 0,
                                intensity_b2 = 0.75,
                                intensity_b3 = -1,
                                camera_detInt = 0.5,
                                camera_detb1 = 0.3,
                                PO_theta0 = -6,
                                PO_theta1 = 1,
                                PO_overdisp = 0.5) {

  true_params <- data.frame(
    param = c("intensity_intercept",
              "intensity_b1",
              "intensity_b2",
              "intensity_b3",
              "camera_detInt",
              "camera_detb1",
              "PO_theta0",
              "PO_theta1",
              "PO_overdisp"),
    value = c(intensity_intercept,
              intensity_b1,
              intensity_b2,
              intensity_b3,
              camera_detInt,
              camera_detb1,
              PO_theta0,
              PO_theta1,
              PO_overdisp)
  )
  
  stopifnot(sim_scenario %in% 1)
  
  
  # True intensity process of interest:
  process_df <- expand.grid(x = 1:grid_dim, y = 1:grid_dim) %>% 
    mutate(PO_agg_x = floor((x - 1) / PO_agg_factor) + 1,
           PO_agg_y = floor((y - 1) / PO_agg_factor) + 1) %>% 
    mutate(cov1 = rnorm(grid_dim^2), 
           cov2 = rnorm(grid_dim^2),
           cov3 = scale(sqrt(2*y))) %>% 
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
  
  po_df <- process_df %>% 
    group_by(PO_agg_x, PO_agg_y) %>% 
    summarize(mu = exp(PO_theta0 + PO_theta1 * log(sum(exp(log_intensity)))),
              .groups = "drop")
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
  
  return(list(
    camera_df = camera_df,
    cell_dat = cell_dat,
    po_df = po_df,
    PO_cell_start = PO_cell_start,
    PO_cell_end = PO_cell_end,
    true_params = true_params,
    visit_cols = visit_cols,
    dataset_ID = dataset_ID
  ))
}

joint_inits <- function(modtype, ncovInt, ncovP) {
  if (modtype == "joint") {
    inits_list <- list(
      intensity_int = 0,
      intensity_beta = rnorm(ncovInt),
      cam_p_int = 0,
      det_beta = rnorm(ncovP),
      theta0 = -5, theta1 = 1, overdisp = 1
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
      theta0 = -5, theta1 = 1, overdisp = 1
    )
  }
  
  inits_list
}

# Set this up so that we can use the same compiled model with different
# simulated datasets, as long as they have the same data dimensions
fit_nimblemodel <- function(datlist, nsim, modtype,
                            ni = 5000, nb = 2000, nc = 2, nt = 1) {
  
  stopifnot(modtype %in% c("joint", "camera_only", "PO_only"))
  
  if (nsim == 1 && length(datlist) > 1) {
    datlist <- list(datlist)
  }  
  
  
  #### Define the NIMBLE model to estimate ####
  modcode <- nimbleCode({
    if (modtype == "joint" || modtype == "camera_only") {
      for (i in 1:ncam) {
        cloglog(psi[i]) <- intensity_int + inprod(intensity_beta[1:ncovInt], cam_dat_int[i, 1:ncovInt])
        cloglog(p[i]) <- cam_p_int + inprod(det_beta[1:ncovP], cam_dat_p[i, 1:ncovP])
        
        cam_y[i, 1:nrep] ~ dOcc_s(probOcc = psi[i], probDetect = p[i], len = nrep)
      }
      
      # Camera-specific priors 
      cam_p_int ~ dnorm(0, sd = 5)
      for (i in 1:ncovP) {
        det_beta[i] ~ dnorm(0, sd = 5)
      }
    }
    
    if (modtype == "joint" || modtype == "PO_only") {
      for (i in 1:nPOcell) {
        log(mu[i]) <- theta0 + theta1 * 
          log(sum(exp(
            predicted_log_intensity[cell_start[i]:cell_end[i]]
          )))
        
        PO_y[i] ~ dnbinom(size = 1 / overdisp,
                          prob = 1 / (1 + overdisp * effort[i] * mu[i]))
      }
      
      # PO-specific priors
      theta0 ~ dnorm(-5, sd = 50)
      theta1 ~ dnorm(1, sd = 2.75)
      overdisp ~ dgamma(shape = 1, scale = 5)
    }
    
    # Intensity priors
    intensity_int ~ dnorm(0, sd = 20)
    for (i in 1:ncovInt) {
      intensity_beta[i] ~ dnorm(0, sd = 5)
    }
    
    
    # Predict and calculate RMSE
    for (i in 1:nCell) {
      predicted_log_intensity[i] <- intensity_int + 
        inprod(intensity_beta[1:ncovInt], cell_dat_int[i, 1:ncovInt])
      cell_sqerr[i] <- (predicted_log_intensity[i] - true_log_intensity[i])^2
    }
    RMSE <- sqrt(mean(cell_sqerr[1:nCell]))
  })
  
  p_covs <- c("cov1", "covC")
  int_covs <- c("cov1", "cov2", "cov3")
  
  ### Build the model
  mod <- nimbleModel(
    code = modcode, 
    constants = list(
      modtype = modtype,
      nrep = length(datlist[[1]]$visit_cols),
      nCell = nrow(datlist[[1]]$cell_dat),
      nPOcell = length(datlist[[1]]$PO_cell_start),
      ncam = nrow(datlist[[1]]$camera_df),
      cell_start = datlist[[1]]$PO_cell_start,
      cell_end = datlist[[1]]$PO_cell_end,
      ncovP = 2,
      ncovInt = 3
    ),
    data = list(
      cell_dat_int = as.matrix(datlist[[1]]$cell_dat[, int_covs]),
      cam_dat_int  = as.matrix(datlist[[1]]$camera_df[, int_covs]),
      cam_dat_p = as.matrix(datlist[[1]]$camera_df[, p_covs]),
      
      true_log_intensity = as.numeric(datlist[[1]]$cell_dat$log_intensity),
      
      cam_y = as.matrix(datlist[[1]]$camera_df[, datlist[[1]]$visit_cols]),
      PO_y  = as.numeric(datlist[[1]]$po_df$y),
      effort = as.numeric(datlist[[1]]$po_df$effort)
    ),
    inits = joint_inits(modtype, 3, 2)
  )
  
  mcmcConf <- configureMCMC(mod)
  mcmcConf$addMonitors("RMSE")
  
  if (modtype == "joint") {
    mcmcConf$removeSampler(c("theta0", "theta1"))
    mcmcConf$addSampler(target = c("theta0", "theta1"), type = "AF_slice")
  }
  
  mcmc <- buildMCMC(mcmcConf)
  complist <- compileNimble(mod, mcmc)
  
  summary_list <- list()
  for (i in 1:nsim) {
    
  # Update data in the model
    complist$mod$setData("effort" = datlist[[i]]$po_df$effort)
    complist$mod$setData("PO_y" = datlist[[i]]$po_df$y)
    complist$mod$setData("cam_y" = datlist[[i]]$camera_df[, datlist[[i]]$visit_cols])
    complist$mod$setData("cell_dat_int" = as.matrix(datlist[[i]]$cell_dat[, int_covs]))
    complist$mod$setData("cam_dat_int" = as.matrix(datlist[[i]]$camera_df[, int_covs]))
    complist$mod$setData("cam_dat_p" = as.matrix(datlist[[i]]$camera_df[, p_covs]))
    complist$mod$setData("true_log_intensity" = as.numeric(datlist[[i]]$cell_dat$log_intensity))
    
  # run the MCMC
    samples <- runMCMC(complist$mcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc,
            inits = joint_inits(modtype, 3, 2), samplesAsCodaMCMC = TRUE)
    
    summary_list[[i]] <- MCMCvis::MCMCsummary(samples) %>% 
      mutate(modtype = modtype, dataset = datlist[[i]]$dataset_ID)
    summary_list[[i]]$param <- rownames(summary_list[[i]])
    rownames(summary_list[[i]]) <- NULL
  }
  
  return(bind_rows(
    summary_list
  ))
}






