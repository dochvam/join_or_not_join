library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

source("worked_example/preprocess_data.R")

estimation_result <- list()

for (species in c("coyote", "cottontail")) {
  dethist <- switch(species,
                    coyote = dethist_CL,
                    cottontail = dethist_SF)
  inat_col <- switch(species,
                     coyote = "n_inat_CL",
                     cottontail = "n_inat_SF")
  this_grid_cell_df <- grid_cell_df %>% filter(n_inat > 0)
  
  # Define NIMBLE code
  we1_code <- nimbleCode({
    # PA data
    for (i in 1:nPA) {
      cloglog(psi[i]) <- beta0 + beta1 * x_PA[i]
      for (j in 1:J[i]) {
        cloglog(p[i, j]) <- alpha0 + alpha1 * w[i] + alpha2 * x_PA[i]
      }
      y_PA[i, 1:J[i]] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J[i]], len = J[i])
    }
    # Priors
    beta0 ~  dnorm(0, sd = 2.72)
    beta1 ~  dnorm(0, sd = 2.72)
    alpha0 ~ dnorm(0, sd = 2.72)
    alpha1 ~ dnorm(0, sd = 2.72)
    alpha2 ~ dnorm(0, sd = 2.72)
    
    # PO data
    if (type == "joint") {
      for (s in 1:S) {
        log(lambda[s]) <- beta0 + beta1 * x[s]
        log(mu[s])     <- theta0 + theta1 * log(lambda[s])
        
        y_PO[s] ~ dnbinom(size = 1/phi, prob = 1 / (1 + phi*mu[s]*E_PO[s]))
      }
      # Priors
      theta0 ~ dnorm(0, sd = 10)
      theta1 <- 1
      phi ~ dgamma(1, 1)
    }
  })
  
  
  # Build and compile the joint model
  joint_mod <- nimbleModel(
    code = we1_code,
    constants = list(
      type = "joint",
      J = J_vec,
      x_PA = this_depl$ndvi,
      w = as.numeric(scale(this_depl$log_roaddist)),
      E_PO = this_grid_cell_df$n_inat,
      x = as.numeric(this_grid_cell_df$ndvi),
      nPA = nrow(dethist),
      S = nrow(this_grid_cell_df)
    ),
    data = list(y_PA = dethist,
                y_PO = this_grid_cell_df[[inat_col]]),
    inits = list(
      beta0 = 0, beta1 = 0, alpha0 = 0, alpha1 = 0, alpha2 = 0,
      theta0 = 0, theta1 = 0, phi = 0.1
    )
  )
  joint_mcmc <- buildMCMC(joint_mod)
  joint_complist <- compileNimble(joint_mod, joint_mcmc)
  
  # Build and compile the single-dataset model
  single_mod <- nimbleModel(
    code = we1_code,
    constants = list(
      type = "single",
      J = J_vec,
      x_PA = this_depl$ndvi,
      w = as.numeric(scale(this_depl$log_roaddist)),
      E_PO = grid_cell_df$n_inat,
      x = as.numeric(grid_cell_df$ndvi),
      nPA = nrow(dethist),
      S = nrow(grid_cell_df)
    ),
    data = list(y_PA = dethist,
                y_PO = grid_cell_df[[inat_col]]),
    inits = list(
      beta0 = 0, beta1 = 0, alpha0 = 0, alpha1 = 0, theta0 = 0, theta1 = 0,
      phi = 0.1
    )
  )
  single_mcmc <- buildMCMC(single_mod)
  single_complist <- compileNimble(single_mod, single_mcmc)
  
  joint_time <- system.time(
    samples_joint <- runMCMC(joint_complist$joint_mcmc, niter = 50000, nburnin = 5000,
                             nchains = 5, thin = 1))
  single_time <- system.time(
    samples_single <- runMCMC(single_complist$single_mcmc, niter = 50000, nburnin = 5000,
                              nchains = 5, thin = 1))
  
  summary_joint <- MCMCsummary(samples_joint) %>% mutate(type = "joint")
  summary_joint$param <- rownames(summary_joint)
  rownames(summary_joint) <- NULL
  summary_one <- MCMCsummary(samples_single) %>% mutate(type = "one")
  summary_one$param <- rownames(summary_one)
  rownames(summary_one) <- NULL
  
  estimation_result[[species]] <- bind_rows(
    summary_joint, summary_one
  ) %>% mutate(species = species)
}

bind_rows(estimation_result) %>% 
  write_csv("worked_example/joint_model_results_wNDVIdet.csv")


results <- read_csv("worked_example/joint_model_results_wNDVIdet.csv") 

results %>% 
  filter(param %in% c("beta0", "beta1")) %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.1)) +
  theme_bw() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~species, nrow = 2) + coord_flip() + xlab("") +
  ylab("Estimate (95%CI)")

results %>% 
  filter(param %in% c("theta0", "theta1")) %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.1)) +
  theme_bw() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~species, nrow = 2) + coord_flip() + xlab("") +
  ylab("Estimate (95%CI)")

results %>% 
  filter(param %in% c("alpha0", "alpha1", "alpha2")) %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.1)) +
  theme_bw() + 
  theme(axis.ticks = element_blank()) +
  facet_wrap(~species, nrow = 2) + coord_flip() + xlab("") +
  ylab("Estimate (95%CI)")


