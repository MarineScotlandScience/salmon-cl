# Master file to run a model and save output
#
# Runs stan and saves the output. Can be called from the command line with
# arguments specifying the covariates to include in the model
# E.g
#
# For linux
# > cd path/to/repo
# > Rscript R/02_run_model.R -c "lat pc1"
#
# For Windows, from cmd prompt
# > cd path\to\repo
# > path\to\R\Rscript.exe R\02_run_model.R -c "lat pc1"
#

#### SETUP

rm(list=ls(all=TRUE))
gc()

library(optparse)   # For parsing command line arguments
# Parse before loading big libraries, for brevity
opt <- OptionParser()
option_list = list(
  make_option(c("-c", "--covariates"), action="store", type='character',
              default = "lat pc1", 
              help='Covariates to include in the model, must be a space delimited string 
              containing the following covariates: 
                lat  - latitude
                pc1  - land usage principal component
                lcpa - log catch per area
                sf   - spring fish population
                e.g.               "lat pc1" [default %default]'),
  make_option(c("-s", "--spatial_position_smoother"), action="store_true", 
              default=FALSE, help="Include the spatial position smoother as an additional term [default %default]"),
  make_option(c("-n", "--n_chains"), action="store", type='integer',
              default=3, help="Number of chains to run [default %default]"),
  make_option(c("-i", "--iterations"), action="store", type='integer',
              default=50, help="Number of post warmup iterations to run per chain [default %default]"),
  make_option(c("-w", "--warmup_factor"), action="store", type='integer',
              default=1, help="Number of warmup iterations = warmup_factor * iterations [default %default]"),
  make_option(c("-D", "--omit_k_dee"), action="store_true", 
              default=FALSE, help="Omit the river K. Dee from analysis [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))


library(rstan)        # For running Bayesian models
library(dplyr)        # Wrangling

source('R/utilities.R') # Shared functions for model construction

# Options for parallelisation
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Select the covariate to add to the model,
# only used in covariate models
COVARIATE_LIST <- c('lat','pc1','lcpa','sf')
# Stan model file name
STAN_MODEL <- 'stan/bhsr_ricker_p_spline.stan'
# Data location
MODEL_DATA_FILE_PATH <- "data/"


n_chains = opt$n_chains
n_iter = opt$iterations + floor(opt$warmup*opt$iterations)

covariates <- sort(stringr::str_split(opt$covariates," ")[[1]])
for(c in covariates){
  if(!c %in% COVARIATE_LIST & c!="")
    stop(paste("Unknown covariate ",c,"in arguments"))
}

print("STARTING MODEL")
print(opt)

####


#### DATA ####
# Load data
if(opt$omit_k_dee){
  sr_data <- readr::read_rds(paste0(MODEL_DATA_FILE_PATH,'model_data_no_dee.rds'))
}else{
  sr_data <- readr::read_rds(paste0(MODEL_DATA_FILE_PATH,'model_data.rds'))
}
all_river_data <- readr::read_rds(paste0(MODEL_DATA_FILE_PATH,'covariate_data.rds'))

# Create model data list
model_data <- create_model_data(sr_data,all_river_data,covariates,opt$spatial_position_smoother)

################################## Run model ###################################

# Parameters for rstan call
inits = function(){list(sigma=runif(1,0,1),
                        sigma_h=runif(1,0,1),
                        sigma_S=runif(1,0,1))}
inits_all = lapply(1:n_chains,function(x) {inits()})



# Run stan
stan_fit = stan(STAN_MODEL,data=model_data,init = inits_all,
                chains = n_chains,iter=n_iter,
                control=list(adapt_delta=0.99,
                             max_treedepth=15))

# Create results object
results = list(stan_fit = stan_fit,
               model_data=model_data,
               opt = opt,
               covariates = covariates,
               sp = opt$spatial_position_smoother,
               code = readLines(STAN_MODEL))

# Save results
results_dir <- ifelse(opt$omit_k_dee,'results/omit_k_dee/','results/')
if(!dir.exists(results_dir)) dir.create(results_dir)
f_name <- paste0(results_dir,
                 'fit_',
                 stringr::str_replace(opt$covariates,' ','_'),
                 ifelse(opt$spatial_position_smoother,'_sp',''),
                 '.rds')
saveRDS(results,file = f_name)
print("DONE")
print(opt)
