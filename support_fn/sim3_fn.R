# Load necessary packages
library(PointedSDMs)
library(sp)
library(spatstat)
library(tidyverse)
library(INLA)
library(inlabru)
library(raster)
library(terra)
library(gstat)

sim3_params <- c("beta0_1", "beta0_2", "beta1", "zeta", "sigma")
sim3_defaults <- c(beta0_1 = -12,
                   beta0_2 = -12,
                   beta1 = 1,
                   zeta = 0,
                   sigma = 0)

sim3_modelOptions <- list(control.inla = 
                       list(int.strategy = 'ccd',
                            cmin = 0), 
                     verbose = FALSE,
                     safe = TRUE)


sim_crs <- "EPSG:32119"

# Define study area extent
extent_x <- c(0, 5000)
extent_y <- c(0, 5000)
res <- 50  # Resolution of covar raster
x_coords <- seq(extent_x[1] + res/2, extent_x[2], by = res)
y_coords <- seq(extent_y[1] + res/2, extent_y[2], by = res)


# Create spatially varying environmental layers
x <- outer(x_coords, y_coords, function(x, y) pmax(x/3000, y / 3000))
x_r <- raster(x)
extent(x_r) <- c(extent_x, extent_y)
crs(x_r) <- sim_crs

bbox <- Polygon(extent(x_r)) # bbox

x_r <- rast(x_r)
ext(x_r) <- c(extent_x, extent_y)
crs(x_r) <- sim_crs
names(x_r) <- "x1"

mesh <- fm_mesh_2d_inla(boundary = inla.sp2segment(bbox),
                        max.edge = max(extent_x) / 25,
                        crs = crs(x_r))
bbox_vect <- st_as_sfc(st_bbox(x_r))


simulate_grf <- function(input_raster, range = 500, sill = 1) {
  # Convert raster to spatial points
  coords <- as.data.frame(input_raster, xy = T)[, 1:2]
  colnames(coords) <- c("x", "y")
  
  # Define the variogram model
  vgm_model <- vgm(psill = sill, model = "Gau", range = range, nugget = 0)
  
  # Simulate Gaussian Random Field
  grf_sim <- gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE, 
                   beta = 0, model = vgm_model, nmax = 30)
  
  grf_values <- predict(grf_sim, coords, nsim = 1)
  
  # Create output raster
  output_raster <- input_raster
  values(output_raster) <- grf_values$sim1
  
  return(output_raster)
}

simulate_data_sim3 <- function(beta0_1,
                               beta0_2,
                               beta1,
                               zeta,
                               sigma) {
  
  ### Generate covariate data
  x <- outer(x_coords, y_coords, function(x, y) pmax(x/3000, y / 3000))
  x_r <- rast(x)
  ext(x_r) <- c(extent_x, extent_y)
  crs(x_r) <- sim_crs
  names(x_r) <- "x1"
  
  v_r <- simulate_grf(x_r, sill = 0.25)
  x_r <- simulate_grf(x_r)
  
  names(v_r) <- "v"
  
  covar_brick <- c(x_r)
  names(covar_brick) <- c("x")
  
  #### Simulate for dataset 1
  # Calculate intensity
  intensity1 <- exp(beta0_1 + beta1 * x_r + v_r)
  ext(intensity1) <- c(extent_x, extent_y)
  crs(intensity1) <- sim_crs
  names(intensity1) <- "intensity1"
  intensity1_mtx <- as.matrix(intensity1, wide = T)
  intensity1_mtx <- intensity1_mtx[nrow(intensity1):1, ] # This needs to be flipped for some reason
  
  # Simulate presence points using an inhomogeneous Poisson process with `spatstat`
  win <- as.owin(c(xmin = extent_x[1], xmax = extent_x[2], 
                   ymin = extent_y[1], ymax = extent_y[2]))
  lambda <- as.im(intensity1_mtx, W = win)  # Intensity function for Poisson process
  
  pp1 <- rpoispp(lambda)
  
  # Convert points to a SpatialPointsDF
  y_PO1 <- st_as_sf(data.frame(x = pp1$x, y = pp1$y),
                    coords = c("x", "y"), crs = crs(x_r))
  
  #### Dataset 2
  # Dataset 2: Calculate intensity
  epsilon_r <- x_r
  values(epsilon_r) <- rnorm(ncell(epsilon_r), values(x_r) * zeta, sigma)
  
  intensity2 <- exp(beta0_1 + (beta1) * (x_r) + (v_r) + epsilon_r)
  ext(intensity2) <- c(extent_x, extent_y)
  crs(intensity2) <- sim_crs
  names(intensity2) <- "intensity2"
  intensity2_mtx <- as.matrix(intensity2, wide = T)
  intensity2_mtx <- intensity2_mtx[nrow(intensity2):1, ] # This needs to be flipped for some reason
  
  # Simulate presence points using an inhomogeneous Poisson process with `spatstat`
  win <- as.owin(c(xmin = extent_x[1], xmax = extent_x[2], 
                   ymin = extent_y[1], ymax = extent_y[2]))
  lambda <- as.im(intensity2_mtx, W = win)  # Intensity function for Poisson process
  
  pp2 <- rpoispp(lambda)
  
  # Convert points to a SpatialPointsDF
  y_PO2 <- st_as_sf(data.frame(x = pp2$x, y = pp2$y),
                    coords = c("x", "y"), crs = crs(x_r))
  
  input_data <- list(d1 = y_PO1, d2 = y_PO2)
  
  return(list(input_data = input_data, covar_brick = c(x_r, v_r),
              true_intensity = as.data.frame(intensity1, xy = T)))
}




run_one_sim3 <- function(iter, 
                         scenario = 0,
                         beta0_1,
                         beta0_2,
                         beta1,
                         zeta,
                         sigma,
                         seed) {
  set.seed(seed + (iter * 17))
  
  capture <- capture.output(
    sim_output <- simulate_data_sim3(beta0_1, beta0_2, beta1, zeta, sigma)
  )
  
  #### Joint model
  joint_runtime <- system.time({
    
    model <- startISDM(sim_output$input_data, 
                       Boundary = bbox_vect,
                       Projection = crs(sim_output$covar_brick), 
                       Mesh = mesh,
                       spatialCovariates = sim_output$covar_brick,
                       Formulas = list(covariateFormula = ~x1))
    
    model$specifySpatial(sharedSpatial = TRUE,
                         constr = TRUE,
                         prior.sigma = c(0.1, 0.05),
                         prior.range = c(200, 0.1))
    

    estimate <- fitISDM(data = model, 
                        options = sim3_modelOptions)
    
    capture <- capture.output(summary_mtx <- summary(estimate)$fixed)
    
    estimation_result_joint <- summary_mtx %>% 
      as.data.frame() %>% 
      dplyr::select(mean, sd, Q025 = `0.025quant`, Q975 = `0.975quant`) %>% 
      mutate(iter = iter, scenario = scenario, type = "joint",
             param = c("beta0_1", "beta0_2", "beta1"),
             truth = c(beta0_1, beta0_2, beta1))
  })
  
  
  ### Single-dataset model
  single_runtime <- system.time({
    
    model <- startISDM(list(d1 = sim_output$input_data$d1), 
                       Boundary = bbox_vect,
                       Projection = crs(sim_output$covar_brick), 
                       Mesh = mesh,
                       spatialCovariates = sim_output$covar_brick,
                       Formulas = list(covariateFormula = ~x1))
    
    model$specifySpatial(sharedSpatial = TRUE,
                         constr = TRUE,
                         prior.sigma = c(0.1, 0.05),
                         prior.range = c(200, 0.1))
    
    estimate <- fitISDM(data = model, 
                        options = sim3_modelOptions)
    
    capture <- capture.output(summary_mtx <- summary(estimate)$fixed)
    
    estimation_result_single <- summary_mtx %>% 
      as.data.frame() %>% 
      dplyr::select(mean, sd, Q025 = `0.025quant`, Q975 = `0.975quant`) %>% 
      mutate(iter = iter, scenario = scenario, type = "one",
             param = c("beta0_1", "beta1"),
             truth = c(beta0_1, beta1))
  })
  
  estimation_result <- bind_rows(estimation_result_single, estimation_result_joint)
  rownames(estimation_result) <- NULL
  
  runtime_result <- data.frame(
    type = c("joint", "one"),
    vaule = unname(c(joint_runtime[3], single_runtime[3]))
  )
  

  return(list(estimation_result = estimation_result %>% mutate(iter = iter),
              runtime_result = runtime_result %>% mutate(iter = iter)))
}



run_many_sim3 <- function(specs_df_onerow, nsim, cl = NULL) {
  stopifnot(all(colnames(specs_df_onerow) %in% c(names(sim3_defaults), "scenario", "seed")))
  
  sim3_values <- sim3_defaults
  for (i in 1:ncol(specs_df_onerow)) {
    if (colnames(specs_df_onerow)[i] %in% sim3_params) {
      sim3_values[colnames(specs_df_onerow)[i]] <- specs_df_onerow[, i]
    }
  }
  
  if (!is.null(cl)) {
    result_list <- parLapply(cl, X = 1:nsim, fun = run_one_sim3,
                             scenario = specs_df_onerow$scenario,
                             seed = specs_df_onerow$seed,
                             beta0_1 = sim3_values["beta0_1"],
                             beta0_2 = sim3_values["beta0_2"],
                             beta1 = sim3_values["beta1"], 
                             zeta = sim3_values["zeta"], 
                             sigma = sim3_values["sigma"])
    
  } else {
    result_list <- list()
    for (i in 1:nsim) {
      result_list[[i]] <- run_one_sim3(iter = i,
                                       scenario = specs_df_onerow$scenario,
                                       seed = specs_df_onerow$seed,
                                       beta0_1 = sim3_values["beta0_1"],
                                       beta0_2 = sim3_values["beta0_2"],
                                       beta1 = sim3_values["beta1"], 
                                       zeta = sim3_values["zeta"], 
                                       sigma = sim3_values["sigma"])
    }
  }
  
  estimation_result_list <- runtime_result_list <- list()
  for (i in 1:nsim) {
    estimation_result_list[[i]] <- result_list[[i]]$estimation_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
    runtime_result_list[[i]]<- result_list[[i]]$runtime_result %>% 
      mutate(scenario = specs_df_onerow$scenario)
  }
  
  return(list(
    estimation_result = bind_rows(estimation_result_list),
    runtime_result = bind_rows(runtime_result_list)
  ))
}

