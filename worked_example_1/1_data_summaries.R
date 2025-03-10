##############################################################################
#' Step 1: Summarize the observational data
#' 
#' The goal here is to characterize the dimensions of the input data as in the
#' simulation study. This should enable us to anticipate the marginal benefits
#' to our parameter estimates that we expect.
##############################################################################

#### Step 0: Load the data ####

source("worked_example_1/preprocess_data.R")

#### Step 1: Summarize the Snapshot USA PA dataset ####

# Date range:
min(cam_pts$start_date)
max(cam_pts$end_date)

median(cam_pts$end_date - cam_pts$start_date)
range(cam_pts$end_date - cam_pts$start_date)

# Number of sites:
nrow(dethist_CL)

# Number of replicates per site
quantile(rowSums(!is.na(dethist_CL)), prob = c(0, 0.25, 0.5, 0.75, 1))

# Naive occupancy rate: what pct of cites have at least one det?
# For coyote:
has_det_CL <- rowSums(dethist_CL, na.rm = T) > 0
mean(has_det_CL)

# For cottontail:
has_det_SF <- rowSums(dethist_SF, na.rm = T) > 0
mean(has_det_SF)

# Naive detection rate: What's the det. rate at sites with at least one det?
# For coyote:
mean(dethist_CL[has_det_CL, ], na.rm = T)

# For cottontail:
mean(dethist_SF[has_det_SF, ], na.rm = T)


# What's the coverage of covariate space?
diff(range(this_depl$ndvi)) / diff(range(grid_cell_df$ndvi))


#### Step 2: Summarize the iNaturalist PO dataset ####

# Number of grid cells with effort
sum(grid_cell_df$n_inat > 0)

# Total num. observations
sum(grid_cell_df$n_inat[grid_cell_df$n_inat > 0])

# iNat effort per grid cell with any
quantile(grid_cell_df$n_inat[grid_cell_df$n_inat > 0], c(0:10/10))
mean(grid_cell_df$n_inat[grid_cell_df$n_inat > 0])

sum(grid_cell_df$n_inat == 1)

# Naive det rate per unit effort:
# Coyote:
sum(grid_cell_df$n_inat_CL) / sum(grid_cell_df$n_inat)

# Cottontail:
sum(grid_cell_df$n_inat_SF) / sum(grid_cell_df$n_inat)

# Coverage of covariate space:
diff(range(grid_cell_df$ndvi[grid_cell_df$n_inat > 0])) / diff(range(grid_cell_df$ndvi))

