# Load fits for the best fit model and
# those which are less parsimonious and calculate and compare
# PSIS-LOO criterion

fit_file_list <- list.files(path = "results/",pattern = "fit_sf.*sp\\.rds",
           full.names = TRUE)

# For all fits, uncomment
# fit_file_list <- list.files(path = "results/",pattern = "fit_.*.rds",
#                             full.names = TRUE)

# Iterate over fits and extract psis-loo criterion
loo_list <- list()
for(f in fit_file_list){
  print(paste("Loading fits",f))
  fit <- readRDS(f)
  loo_list[[f]] <- loo::loo(fit$stan_fit)
}

# Check output, a couple of models have 1 ok (not good) sample
# should be fine
loo_list

# Loo comparison
loo_comp_tbl <- loo::loo_compare(loo_list)
loo_comp_tbl

# Save
saveRDS(loo_comp_tbl,file = 'results/loo_comparison.rds')
