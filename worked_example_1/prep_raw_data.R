library(tidyverse)
library(rgbif)

# Get the relevant Snapshot data

snapshot_all_depl <- read_csv("../iNat_CT_joint/all_deployments.csv")
snapshot_all_seq  <- read_csv("../iNat_CT_joint/all_sequences.csv")

this_depl <- snapshot_all_depl %>% 
  filter(year(start_date) == 2023, subproject_name == "NYC") %>% 
  select(start_date, end_date, site_name, deployment_id, longitude, latitude, 
         log_roaddist, Canopy_height)

stopifnot(length(unique(this_depl$deployment_id)) == nrow(this_depl))

this_seq <- snapshot_all_seq %>% 
  filter(deployment_id %in% this_depl$deployment_id) %>% 
  filter(common_name %in% c("Eastern Cottontail", "Coyote"))

write_csv(this_depl, "worked_example_1/data/WE1_deployments.csv")
write_csv(this_seq,  "worked_example_1/data/WE1_sequences.csv")

bbox <- c(min(this_depl$longitude), max(this_depl$longitude),
          min(this_depl$latitude),  max(this_depl$latitude))
bbox[c(1, 3)] <- bbox[c(1, 3)] - 0.05
bbox[c(2, 4)] <- bbox[c(2, 4)] + 0.05


# Download all iNat mammal data in this bounding box

mammal_key <- name_backbone(name = "Mammalia", rank = "class")$usageKey

# Query GBIF for iNaturalist research-grade mammal observations within the bounding box
mammal_data <- occ_search(
  taxonKey = mammal_key,
  basisOfRecord = "HUMAN_OBSERVATION",
  datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7", # iNaturalist dataset key in GBIF
  hasCoordinate = TRUE,
  eventDate = "2022,2024", # Optional: filter by date range
  decimalLatitude = paste(bbox[3], bbox[4], sep = ","),
  decimalLongitude = paste(bbox[1], bbox[2], sep = ","),
  limit = 10000
)

# Extract the data frame
mammal_df <- mammal_data$data %>% 
  select(longitude = decimalLongitude, latitude = decimalLatitude, 
         scientificName, year, month, day, eventDate, datasetName, recordedBy)

# Display the first few rows
head(mammal_df)

# Save to CSV
write_csv(mammal_df, "worked_example_1/data/WE1_iNat.csv")


