library(spOccupancy)
library(tidyverse)
library(parallel)
library(gridExtra)
library(PointedSDMs)
library(sp)
library(spatstat)
library(INLA)
library(inlabru)
library(raster)
library(terra)
library(gstat)

source("support_fn/sim3_fn.R")


# Check that the fn works:
run_one_sim3(1, scenario = 1, beta0_1 = -13, beta0_2 = -13, beta1 = 1,
             zeta = 0, sigma = 0, seed = 7654)

# Expected output:
# $estimation_result
# mean    sd    Q025    Q975 iter scenario  type   param truth
# 1 -12.880 0.141 -13.157 -12.602    1        1   one beta0_1   -13
# 2   0.875 0.105   0.669   1.082    1        1   one   beta1     1
# 3 -12.846 0.119 -13.079 -12.614    1        1 joint beta0_1   -13
# 4 -12.722 0.114 -12.946 -12.499    1        1 joint beta0_2   -13
# 5   0.822 0.078   0.669   0.976    1        1 joint   beta1     1
# 
# $runtime_result
# type num_points vaule iter
# 1 joint        233 11.84    1
# 2   one        109  5.35    1



# Checking that we get an NA result on error:
run_one_sim3(iter = 1, scenario = 1, beta0_1 = -13, beta0_2 = NA, beta1 = 1,
             zeta = 0, sigma = 0, seed = 7654)

