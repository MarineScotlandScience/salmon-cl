rm(list=ls(all=TRUE))
gc()

library(optparse)   # For parsing command line arguments
# Parse before loading big libraries, for brevity
opt <- OptionParser()
option_list = list(
  make_option(c("-m", "--model"), action="store",
              default = "./results/fit_sf_sp.rds",
              help="Filename for previously run model for which to run CV"),
  make_option(c("-n", "--n_chains"), action="store", type='integer',
              default=3, help="Number of chains to run [default %default]"),
  make_option(c("-i", "--iterations"), action="store", type='integer',
              default=100, help="Number of post warmup iterations to run per chain [default %default]"),
  make_option(c("-w", "--warmup_factor"), action="store", type='integer',
              default=1, help="Number of warmup iterations = warmup_factor * iterations [default %default]"),
  make_option(c("-e", "--n_eff_diag"), action="store", type='integer',
              default=1000, help="Effective sample size for diagnostic checks [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))


library(rstan)        # For running Bayesian models
library(bayesplot)    # Additional diagnostics
library(dplyr)        # Wrangling

source('R/utilities.R') # Shared functions for model construction

# Options for parallelisation
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

MODEL_DATA_FILE_PATH <- "data/"

n_chains = opt$n_chains
n_iter = opt$iterations
if(!file.exists(opt$model)){
  stop(paste("Could not find model file",opt$model)) 
}

print("STARTING CV")
print(opt)

#### DATA ####
# Load data
sr_data <- readr::read_rds(paste0(MODEL_DATA_FILE_PATH,'model_data.rds'))
all_river_data <- readr::read_rds(paste0(MODEL_DATA_FILE_PATH,'covariate_data.rds'))

# Load model
results <- readRDS(opt$model)

# Initialise
to_run <- unique(sr_data$site)

# Is there existing progress to load?
cv_filename <- stringr::str_replace(opt$model,"fit","cv")
if(file.exists(cv_filename)){
  print(paste('Loading existing CV file',cv_filename))
  cv_results <- readRDS(cv_filename)
  to_run <- to_run[!to_run %in% names(cv_results$folds)]
}else{
  # New results
  cv_results <- list()
  cv_results$covariates <- results$covariates
  cv_results$sp <- results$sp
  cv_results$model <- opt$model
  cv_results$folds <- list()
}

# Iterate through required folds
for(r in to_run){
  # Create model data list
  print(paste("Running CV holding out",r,
              " (",
              which(to_run==r),"/",
              length(to_run),
              "remaining)"))
  model_data <- create_model_data(sr_data,all_river_data,
                                  results$covariates,
                                  results$sp,
                                  holdout_river = r)
  
  # Parameters for rstan call
  inits = function(){list(sigma=runif(1,0,1),
                          sigma_h=runif(1,0,1),
                          sigma_S=runif(1,0,1))}
  inits_all = lapply(1:n_chains,function(x) {inits()})
  
  # Run stan
  stan_fit = stan(model_code = stringr::str_flatten(results$code,collapse="\n"),
                  data=model_data,
                  init = inits_all,
                  chains = n_chains,
                  iter=n_iter,
                  control=list(adapt_delta=0.99,
                               max_treedepth=15))
  
  # Were there divergences?
  np <- nuts_params(stan_fit)
  n_divergent = sum(np[np$Parameter=='divergent__',]$Value)
  
  # Check rhat and neff
  diag_pars <- c('sigma','alpha_h','alpha_S','beta_h','beta_S','tau_h','tau_S','mu_h_tilde','log_S_tilde')
  s <- summary(stan_fit,pars=diag_pars)
  n_eff <- s$summary[,'n_eff']
  r_hat <- s$summary[,'Rhat']
  
  if(n_divergent > 0) {
    print(paste0(n_divergent,' divergent transitions in model fit, moving to next fold.'))
    #TODO: Check for convergence issues etc?
  }else if(sum(abs(r_hat - 1) > 0.05) > 0){
    print(paste(sum(abs(r_hat - 1) > 0.05) > 0,'parameters have rhat > 1.05 / < 0.95, moving to next fold.'))
    print(r_hat) 
  }else if(sum(n_eff < opt$n_eff_diag)){
    print(paste(sum(n_eff < opt$n_eff_diag),'parameters have effective sample size <',opt$n_eff_diag,', moving to next fold'))
    print(n_eff)
  }else{
    # Good fit, 
    # save log_lik of heldout fold and continue
    log_lik <- extract(stan_fit,'log_lik_new')
    cv_results$folds[[r]] <- log_lik
    saveRDS(cv_results,file = cv_filename)
  }
}

if(length(cv_results$folds)!=length(unique(sr_data$site))){
  print(paste("Only",length(cv_results$folds),"of",length(unique(sr_data$site)),
              "folds complete. Run again to calculate missing folds"))
}else{
  print("CV complete") 
  print(opt)
}


