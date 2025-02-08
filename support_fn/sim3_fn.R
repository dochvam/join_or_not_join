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

scenario <- 1
beta0_1 <- -12
# beta0_2 <- -12
beta1 <- 1.5
nsim <- 20
crs <- "EPSG:32119"


# Define study area extent
extent_x <- c(0, 5000)
extent_y <- c(0, 5000)
res <- 50  # Resolution of covar raster
# Generate environmental covariates
x_coords <- seq(extent_x[1] + res/2, extent_x[2], by = res)
y_coords <- seq(extent_y[1] + res/2, extent_y[2], by = res)

# Create spatially varying environmental layers
x <- outer(x_coords, y_coords, function(x, y) pmax(x/3000, y / 3000))
x_r <- raster(x)
extent(x_r) <- c(extent_x, extent_y)
crs(x_r) <- crs

bbox <- Polygon(extent(x_r)) # bbox

x_r <- rast(x_r)
ext(x_r) <- c(extent_x, extent_y)
crs(x_r) <- crs
names(x_r) <- "x1"

mesh <- fm_mesh_2d_inla(boundary = inla.sp2segment(bbox),
                        max.edge = max(extent_x) / 25,
                        crs = crs(x_r))
bbox_vect <- st_as_sfc(st_bbox(x_r))

estimates_list <- list()
for (i in 1:nsim) {
  x_r <- simulate_grf(x_r)
  v_r <- simulate_grf(x_r, sill = 0.25)
  names(v_r) <- "v"
  
  # Calculate intensity
  intensity_r <- exp(beta0_1 + beta1 * (x_r) + (v_r))
  ext(intensity_r) <- c(extent_x, extent_y)
  crs(intensity_r) <- crs
  names(intensity_r) <- "intensity"
  intensity_mtx <- as.matrix(intensity_r, wide = T)
  intensity_mtx <- intensity_mtx[nrow(intensity_mtx):1, ] # This needs to be flipped for some reason
  
  # Simulate presence points using an inhomogeneous Poisson process with `spatstat`
  win <- as.owin(c(xmin = extent_x[1], xmax = extent_x[2], 
                   ymin = extent_y[1], ymax = extent_y[2]))
  lambda <- as.im(intensity_mtx, W = win)  # Intensity function for Poisson process
  
  pp1 <- rpoispp(lambda)
  pp2 <- rpoispp(lambda)
  
  # Convert points to a SpatialPointsDF
  y_PO1 <- st_as_sf(data.frame(x = pp1$x, y = pp1$y),
                    coords = c("x", "y"), crs = crs(x_r))
  y_PO2 <- st_as_sf(data.frame(x = pp2$x, y = pp2$y), 
                    coords = c("x", "y"), crs = crs(x_r))
  
  input_data <- list(d1 = y_PO1, d2 = y_PO2)
  
  model <- startISDM(input_data, 
                     Boundary = bbox_vect,
                     Projection = crs(x_r), 
                     Mesh = mesh,
                     spatialCovariates = x_r)
  
  model$specifySpatial(sharedSpatial = TRUE,
                       constr = TRUE,
                       prior.sigma = c(0.1, 0.05),
                       prior.range = c(200, 0.1))
  
  
  modelOptions <- list(control.inla = 
                         list(int.strategy = 'ccd',
                              cmin = 0), 
                       verbose = FALSE,
                       safe = TRUE)
  
  estimate <- fitISDM(data = model, 
                      options = modelOptions)
  
  estimates_list[[i]] <- estimate$summary.fixed %>% 
    dplyr::select(mean, sd, Q025 = `0.025quant`, Q975 = `0.975quant`) %>% 
    mutate(iter = i, scenario = scenario,
           param = rownames(estimate$summary.fixed))
  rownames(estimates_list[[i]]) <- NULL
}


