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


#### Run a quick simulation using these values ####
source('support_fn/sim2_fn.R')


# Based on coyote
specs_df <- data.frame(
  S = sum(grid_cell_df$n_inat > 0),
  nPA = nrow(dethist_CL), 
  J = 3,
  xi = 0.517,
  eta = 10.8,
  q = 1,
  beta0 = logit(0.6),
  beta1 = 1,
  alpha0 = logit(0.46),
  alpha1 = 0,
  theta0 = log(0.006) - logit(0.6),
  theta1 = 1
)
CL_sim <- run_many_sim2(specs_df_onerow = specs_df, nsim = 100)

saveRDS(CL_sim, file = "worked_example_1/CL_sim.RDS")

# Based on cottontail
specs_df <- data.frame(
  S = sum(grid_cell_df$n_inat > 0),
  nPA = nrow(dethist_CL), 
  J = 3,
  xi = 0.517,
  eta = 10.8,
  q = 1,
  beta0 = logit(0.6),
  beta1 = 1,
  alpha0 = logit(0.674),
  alpha1 = 0,
  theta0 = log(0.055) - logit(0.6),
  theta1 = 1
)
SF_sim <- run_many_sim2(specs_df_onerow = specs_df, nsim = 100)

saveRDS(SF_sim, file = "worked_example_1/SF_sim.RDS")


bind_rows(
  # Summarize CL
  CL_sim$estimation_result %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = (mean - truth)^2 + sd^2) %>% 
    select(MSE, param, iter, type) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint, species = "Coyote"),
  
  # Summarize SF
  SF_sim$estimation_result %>% 
    filter(param %in% c("beta0", "beta1")) %>% 
    mutate(MSE = (mean - truth)^2 + sd^2) %>% 
    select(MSE, param, iter, type) %>% 
    pivot_wider(names_from = type, values_from = MSE) %>% 
    mutate(improvement = one - joint, species = "Cottontail")
) %>% 
  group_by(species, param) %>% 
  summarize(median_improvement = median(improvement),
            improvement_rate = mean(improvement > 0))

# Make a result table