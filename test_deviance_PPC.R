library(spOccupancy)

#### Define the deviance PPC function ####
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



#### Deviance PPC reproducible example ####
set.seed(325697)

J.x <- 15
J.y <- 15
J.all <- J.x * J.y

# Number of data sources.
n.data <- 2
# Sites for each data source. 
J.obs <- sample(ceiling(0.2 * J.all):ceiling(0.5 * J.all), n.data, replace = TRUE)
# Replicates for each data source.
n.rep <- list()
for (i in 1:n.data) {
  n.rep[[i]] <- sample(1:4, size = J.obs[i], replace = TRUE)
}
# Occupancy covariates
beta <- c(0.5, 1, -3)
p.occ <- length(beta)
# Detection covariates
alpha <- list()
for (i in 1:n.data) {
  alpha[[i]] <- runif(2, -1, 1)
}
p.det.long <- sapply(alpha, length)
p.det <- sum(p.det.long)
sigma.sq <- 2
phi <- 3 / .5
sp <- TRUE

# Simulate occupancy data. 
dat <- simIntOcc(n.data = n.data, J.x = J.x, J.y = J.y, J.obs = J.obs, 
                 n.rep = n.rep, beta = beta, alpha = alpha, sp = FALSE)

# Inflate D2 sites with 1s to induce some error, which hopefully shows up in PPC
for (i in 1:(nrow(dat$y[[2]]) / 2)) {
  dat$y[[2]][i, !is.na(dat$y[[2]][i, ])] <- 1
}


y <- dat$y
X <- dat$X.obs
X.p <- dat$X.p
sites <- dat$sites

# Package all data into a list
occ.covs <- X[, 2, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list()
# Add covariates one by one
det.covs[[1]] <- list(det.cov.1.1 = X.p[[1]][, , 2]) 
det.covs[[2]] <- list(det.cov.2.1 = X.p[[2]][, , 2]) 
data.list <- list(y = y, 
                  occ.covs = occ.covs,
                  det.covs = det.covs, 
                  sites = sites)


# Fit the model
fit <- intPGOcc(~occ.cov,
                list(~det.cov.1.1, ~det.cov.2.1), 
                data = data.list, 
                n.burn = 2500, 
                n.samples = 5000,
                n.chains = 2, 
                verbose = T)

# Run the deviance 
dev_ppc_result <- ppcOcc_Deviance(fit)

# Make sure these samples look sane
hist(dev_ppc_result$fit.y[[2]])
hist(dev_ppc_result$fit.y.rep[[2]])
hist(dev_ppc_result$fit.y.rep[[2]] - dev_ppc_result$fit.y[[2]])

summary(dev_ppc_result)
# p-value isn't actually far from 0.5




