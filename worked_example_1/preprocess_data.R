library(tidyverse)
library(terra)
library(unmarked)

#### Load all observation data ####

inat_dat  <- read_csv("worked_example_1/data/WE1_iNat.csv")
this_seq  <- read_csv("worked_example_1/data/WE1_sequences.csv")
this_depl <- read_csv("worked_example_1/data/WE1_deployments.csv")

## Load the spatial layer of interest
ndvi_nyc_raw <- rast("worked_example_1/data/Avg_NDVI_Landsat_NYC.tif")


#### Spatialize and grid all data to  0.05 decimal degrees ####

main_grid <- rast(ext(ndvi_nyc_raw), resolution = 0.01)
terra::values(main_grid) <- 1:ncell(main_grid)

ndvi_nyc <- resample(ndvi_nyc_raw, main_grid, method = "bilinear")
terra::values(ndvi_nyc) <- scale(terra::values(ndvi_nyc))

inat_pts <- vect(inat_dat, geom = c("longitude", "latitude"),
                 crs = "+proj=longlat")
cam_pts  <- vect(this_depl, geom = c("longitude", "latitude"),
                 crs = "+proj=longlat")


inat_dat$grid_cell  <- as.numeric(extract(main_grid, inat_pts)[, 2])
this_depl$grid_cell <- as.numeric(extract(main_grid, cam_pts)[, 2])

inat_dat$ndvi  <- as.numeric(extract(ndvi_nyc, inat_pts)[, 2])
this_depl$ndvi <- as.numeric(extract(ndvi_nyc, cam_pts)[, 2])


grid_cell_df <- data.frame(grid_cell = 1:ncell(main_grid)) %>% 
  left_join(
    inat_dat %>% count(grid_cell) %>% rename(n_inat = n)
  ) %>% 
  left_join(
    this_depl %>% count(grid_cell) %>% rename(n_ct = n)
  )
grid_cell_df[is.na(grid_cell_df)] <- 0
grid_cell_df$ndvi <- terra::values(ndvi_nyc)

inat_effort_rast <- ct_effort_rast <- main_grid
terra::values(inat_effort_rast) <- grid_cell_df$n_inat
terra::values(ct_effort_rast) <- grid_cell_df$n_ct



#### Construct PA detection histories ####

window <- 10

dethist_SF <- dethist_CL <- matrix(nrow = nrow(this_depl), 
                     ncol = floor(max(this_depl$end_date - this_depl$start_date)/window))
J_vec <- numeric(nrow(dethist_SF))

for (i in 1:nrow(dethist_SF)) {
  J_vec[i] <- floor((this_depl$end_date[i] - this_depl$start_date[i])/window)
  
  for (j in 1:J_vec[i]) {
    dethist_SF[i, j] <- as.numeric(sum(
      this_seq$deployment_id == this_depl$deployment_id[i] & 
        this_seq$common_name == "Eastern Cottontail" &
        date(this_seq$start_time) >= this_depl$start_date[i] + window*(j - 1) &
        date(this_seq$start_time) <= this_depl$start_date[i] + window*j
    ) > 0)
    
    dethist_CL[i, j] <- as.numeric(sum(
      this_seq$deployment_id == this_depl$deployment_id[i] & 
        this_seq$common_name == "Coyote" &
        date(this_seq$start_time) >= this_depl$start_date[i] + window*(j - 1) &
        date(this_seq$start_time) <= this_depl$start_date[i] + window*j
    ) > 0)
  }
}

#### Calculate iNat species counts per cell ####

grid_cell_df <- grid_cell_df %>% 
  left_join(
    inat_dat %>% filter(grepl("Canis latrans", inat_dat$scientificName)) %>% 
      count(grid_cell) %>% rename(n_inat_CL = n)
  ) %>% 
  left_join(
    inat_dat %>% filter(grepl("Sylvilagus floridanus", inat_dat$scientificName)) %>% 
      count(grid_cell) %>% rename(n_inat_SF = n)
  )

grid_cell_df$n_inat_CL[is.na(grid_cell_df$n_inat_CL)] <- 0
grid_cell_df$n_inat_SF[is.na(grid_cell_df$n_inat_SF)] <- 0


#### Make a map of the data ####

(grid_cell_df %>% 
  left_join(as.data.frame(main_grid, xy = TRUE), by = c("grid_cell" = "lyr.1")) %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = n_inat)) +
  geom_point(aes(x, y, col = n_ct > 0, alpha = n_ct > 0)) +
  scale_color_manual("Cell contains CT", values = c(NA, "#dd3344")) +
  scale_alpha_manual("Cell contains CT", values = c(0.05, 1)) +
  scale_fill_viridis_c("Num. iNat obs.",
                       trans = "log", na.value = "#dddddd", breaks = c(1, 8, 64, 400)) +
  theme_minimal() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Map of sampling effort")) %>% 
  ggsave(filename = "worked_example_1/effort_plot.jpg", width = 5, height = 4)
