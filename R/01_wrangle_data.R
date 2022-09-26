# Create tibbles in format suitable for analysis from raw csv files.

library(tidyverse)
library(forcats)

# Set to true to create a data set that does not include the K. Dee. 
# This is a river with few data points, but is an outlier for some covariates,
# so for robustness it is useful to make sure it does not have too strong an
# influence on the results.
OMIT_K_DEE <- FALSE

MODEL_DATA_FILE_PATH <- "data"

# #### Stock Recruitment data ####
# # Includes covariates for between stock relationships

sr_data <- readr::read_csv(paste0(MODEL_DATA_FILE_PATH,'/stock_recruitment.csv'))

if(OMIT_K_DEE) {
  # Scotland + Lune but remove the K. Dee
  sr_data <- sr_data %>% filter(site != 'River Dee (Kirkcudbrightshire)')
  file_name <- "/model_data_no_dee.rds"
}else{
  # Scotland + Lune
  file_name <- "/model_data.rds"
}

sr_data <- sr_data %>% 
  mutate(site = as.factor(as.character(site))) %>% 
  # Convert to egg density
  mutate(S = s/area,
         R = r/area,
         # Turn categoricals to factors
         spring_fish_raw = spring_fish,
         spring_fish = factor(spring_fish),
         coast_raw = coast,
         coast = factor(coast))

readr::write_rds(file = paste0(MODEL_DATA_FILE_PATH,file_name),
                 sr_data)



# #### All river covariates ####
# Between stock covariate data for all rivers

all_river_data <- readr::read_csv(paste0(MODEL_DATA_FILE_PATH,'/river_covariates.csv'))

all_river_data <- all_river_data %>% 
  # All the assessable are necessarily whole catchment, set spring_fish to 0
  mutate(spring_fish = factor('n',levels = c('n','y'))) %>% 
  # rename river to site for compatibility
  mutate(site = river) %>% 
  # West coast if gis order greater than 55, set to 1
  mutate(coast = ifelse(gis.order > 55,'w','e')) %>% 
  arrange(site)

readr::write_rds(all_river_data,file = paste0(MODEL_DATA_FILE_PATH,'/covariate_data.rds'))

