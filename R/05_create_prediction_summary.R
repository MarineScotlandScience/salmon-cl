# Construct summary tibbles for in sample and out of
# sample sites, including hypothetical sites for
# positions around the coast

## CONFIG ####
rm(list = ls())
library(tidyverse)
theme_set(ggpubr::theme_pubclean())

source('R/utilities.R')

## LOAD ####
# Load the model
results <- readRDS('results/fit_sf_sp.rds')

# Load data and models ----------------------------------------------------

# Original data used to fit the model
sr_data <- readRDS('data/model_data.rds')

# Data for all rivers in assessment
all_river_data = readRDS('./data/covariate_data.rds')

num_basis = results$model_data$num_basis
stan_fit_df = as.data.frame(results$stan_fit)

# Extract posterior samples for parameters
S_star  <- stan_fit_df %>% select(starts_with('S_star['))
beta_S  <- stan_fit_df %>% select(starts_with('beta_S['))
alpha_S <- stan_fit_df %>% select(starts_with('alpha_S['))
sigma_S <- stan_fit_df %>% select(starts_with('sigma_S'))

H_star  <- stan_fit_df %>% select(starts_with('h_star['))
beta_H  <- stan_fit_df %>% select(starts_with('beta_h['))
alpha_H <- stan_fit_df %>% select(starts_with('alpha_h['))
sigma_H <- stan_fit_df %>% select(starts_with('sigma_h'))

sigma   <- stan_fit_df %>% select(sigma)


# Construct predictions ---------------------------------------------------

# FUNCTIONS

# Predict S using global samples and specified
# data X,Y
predict_S <- function(X,Y,linpred=FALSE){
  N <- nrow(X)
  
  sigma_S_pred <- t(matrix(rep(as.matrix(sigma_S),N),ncol = N))
  
  S_pred <- X %*% t(alpha_S) + Y %*% t(beta_S)
  #S_pred <- X %*% t(alpha_S) 
  
  # If we just want linear predictor
  if(linpred) return(exp(S_pred))
  
  S_pred <- rnorm(length(S_pred),S_pred,sigma_S_pred)
  S_pred <- matrix(exp(S_pred),nrow = N)
}

# Predict H using global samples and specified
# data X,Y
predict_H <- function(X,Y,linpred=FALSE){
  N <- nrow(X)
  
  sigma_H_pred <- t(matrix(rep(as.matrix(sigma_H),N),ncol = N))
  
  H_pred <- X %*% t(alpha_H) + Y %*% t(beta_H)
  #H_pred <- X %*% t(alpha_H) 
  
  # If we just want linear predictor
  if(linpred) return(InvLogit(H_pred))
    
  H_pred <- rnorm(length(H_pred),H_pred,sigma_H_pred)
  H_pred <- matrix(InvLogit(H_pred),nrow = N)
}

# Calculate posterior predictions of R
# given S_star,h_star and predicted S
predict_ricker = function(S_star,h_star,S,linpred=FALSE){
  
  S <- matrix(S,ncol = 1)
  N <- nrow(S)
  i <- matrix(rep(1,N),ncol=1)
  
  sigma_pred <- t(matrix(rep(as.matrix(sigma),N),ncol = N))
  
  # mu = S/(1-h_star) * exp(h_star*(1-S/S_star))
  # R ~ lognormal(log(mu),sigma)
  tmp_1 <- S %*% (1/(1-h_star))
  tmp_2 <- (1-S %*% (1/S_star))
  mu <- tmp_1 * exp(i%*%h_star * tmp_2)
  
  if(linpred) mu
  
  R <- rlnorm(length(mu),log(mu),sigma_pred)
  R <- matrix(R,nrow = N)
  
  return(R)
}

# Samples to quantiles
to_quantiles_df <- function(X,var_name='S',...){
  X_quants <- t(apply(X,1,quantile,probs=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),...))
  colnames(X_quants) <- paste(var_name,c('q05','q10','q25','q50','q75','q90','q95'),sep='_')
  return(as_tibble(X_quants))
  
}

# Bayesian R2
# See R-squared for Bayesian regression models
# Gelman et al. 2018

# Using fitted variance parameter for residual variance
bayes_R2 <- function(pred,sigma,y) {
  var_fit <- apply(pred, 1, var)
  var_res <- sigma^2
  r2 <- var_fit / (var_fit + var_res)
  colnames(r2) <- "R2"
  r2
}


# Using residuals for residual variance
bayes_R2_2 <- function(y_hat,y) {
  var_fit <- apply(y_hat, 2, var)
  e <- -1*sweep(y_hat,1,y)
  var_e <- apply(e,2,var)
  r2 <- var_fit / (var_fit + var_e)
  r2 <- data.frame(R2 = r2)
  r2
}

# Prediction data for hypothetical rivers
unscaled_range = range(all_river_data$round.coast,sr_data$round.coast)
unscaled_lat = range(all_river_data$lat,sr_data$lat)

# Values of SP which are closest to 55 degree on both coasts
# for comparison with previous analyses
fifty_five <- c(1.264,7.74)

P <- 52
pred_rc <- c(seq(unscaled_range[1],
                 unscaled_range[2],
                 length.out = P-2),fifty_five)
pred_lat <- seq(unscaled_lat[1],unscaled_lat[2],length.out=P-2)

pred_data <- tibble(round.coast = rep(pred_rc,2),
                    lat = 0,
                    coast = factor('e',levels=c('e','w')),
                    spring_fish = factor(rep(c('n','y'),each=P)),
                    land.cover.pc1=0,
                    cpa=rep(c(1,2),each=P)) %>% 
  mutate(site = paste0('site_',row_number())) %>% 
  arrange(site)



# Posterior predictions, hypothetical new rivers ------------------------

# Construct design matrices for linear covariates and P-splines
X_pred <- create_model_matrix(pred_data,covariates = results$covariates)
X_pred[,'pc1'] <- 0
X_pred[,'lat'] <- 0

Y_pred <- create_coast_spline_matrix(pred_data$round.coast,s_range=unscaled_range,num_basis)

S_pred <- predict_S(X_pred,Y_pred)
H_pred <- predict_H(X_pred,Y_pred)

S_quants <- to_quantiles_df(S_pred,'S')
H_quants <- to_quantiles_df(H_pred,'H')


# Bayesian R2 -------------------------------------------------------------
# Linear predictors required for bayesian R squared values on 
# derived parameters

# R2 values
# A few issues with R2 values
# 1. R2 values for S* are slightly weird, as the S* values themselves are 
#   parameter estiamtes with uncertainty. S* (or h*) are not dependent variables 
#   that have been measured.
# 2. R2 values are slightly problematic for Bayesian models, as it is possible 
#   (with strong prior information) to have an R2 greater than 1, using the 
#   classic definition.
# 3. There exisits a 'Bayesian R2' value which can be calculated, but this does 
#   not have the same interpretation as 'proportion variance explained'. It is 
#   given as ratio of variance in predictive values, to the variance of 
#   predictive value plus the variance of the residuals. This is conditional on 
#   the model and hence cannot always be used to show improvement for model fits.

S_pp <- predict_S(results$model_data$X,results$model_data$Y,linpred = TRUE)
H_pp <- predict_H(results$model_data$X,results$model_data$Y,linpred = TRUE)

R2_S <- bayes_R2(S_pp,sigma_S)
R2_h <- bayes_R2(H_pp,sigma_H)

# PLot and summarise
R2_S %>% ggplot(aes(x=R2)) + geom_histogram()
R2_h %>% ggplot(aes(x=R2)) + geom_histogram()
summary(R2_S)
summary(R2_h)

# Posterior predictions actual new rivers -------------------------------
X_all_rivers <- create_model_matrix(all_river_data,results$covariates)
Y_all_rivers <- create_coast_spline_matrix(all_river_data$round.coast,
                                           num_basis = num_basis,
                                           s_range = unscaled_range)
S_all_rivers <- predict_S(X_all_rivers,Y_all_rivers)
H_all_rivers <- predict_H(X_all_rivers,Y_all_rivers)

S_quants_all_rivers <- to_quantiles_df(S_all_rivers,'S',na.rm=TRUE)
H_quants_all_rivers <- to_quantiles_df(H_all_rivers,'H',na.rm=TRUE)

# Posterior predictions for fit sites -----------------------------------
S_star_quants <- to_quantiles_df(t(S_star),'S')
H_star_quants <- to_quantiles_df(t(H_star),'H')

# SR predictions for fit sites
pred_sr <- sr_data %>% group_by(site) %>% 
  summarise(S_max = max(S),
            S_pred = seq(0,S_max,length.out = P)) %>% 
  ungroup()


summary_sr_data <- list()
y_hat <- matrix(nrow = nrow(sr_data),ncol = nrow(sigma))
for(i in 1:nlevels(pred_sr$site)){
  
  # Predictions and summary for new data
  sub_data <- pred_sr %>% filter(site==levels(site)[i])
  R <- predict_ricker(S_star[,i],H_star[,i],sub_data$S_pred)
  summary_sr_data[[i]] <- sub_data %>% bind_cols(to_quantiles_df(R,'R'))
  
  # y_hat estimates for orig data (for use in Bayes R2)
  y_hat_data_i <- which(sr_data$site==levels(sr_data$site)[i])
  y_hat[y_hat_data_i,] <- predict_ricker(S_star[,i],H_star[,i],sr_data[y_hat_data_i,]$S,linpred = TRUE)
  
}
summary_sr_data <- bind_rows(summary_sr_data)

R2_all <- bayes_R2_2(y_hat,sr_data$R)
R2_all %>% ggplot(aes(x=R2)) + geom_histogram()

resid <- sweep(log(y_hat),1,log(sr_data$R))
resid_tbl <- to_quantiles_df(resid,var_name = 'resid') %>% 
  mutate(y = 1:n(),
         site = sr_data$site,
         R = sr_data$R)
resid_tbl %>% ggplot(aes(x=R,y=resid_q50,ymin=resid_q05,ymax=resid_q95,col=site)) + 
  geom_errorbar(stat = 'identity',alpha=0.3) + geom_point()

resid_tbl %>% ggplot(aes(x=R,y=resid_q50,col=site)) + 
   geom_point()

resid_tbl %>% ggplot(aes(sample=resid_q50,col=site)) + 
  stat_qq() + stat_qq_line()



# CONSOLIDATE 

pred_data <- pred_data %>% 
  bind_cols(S_quants) %>% 
  bind_cols(H_quants) %>% 
  mutate(spring_fish = factor(spring_fish))

all_river_data <- all_river_data %>% 
  select(site,lat,coast,cpa,spring_fish,round.coast) %>% 
  bind_cols(S_quants_all_rivers) %>% 
  bind_cols(H_quants_all_rivers)

fit_rivers <- sr_data %>% 
  group_by(site) %>% 
  summarise_all(first) %>% 
  select(site,lat,coast,cpa,spring_fish,round.coast) %>% 
  bind_cols(S_star_quants) %>% 
  bind_cols(H_star_quants) %>% 
  ungroup()

all_river_data <- all_river_data %>% 
  filter(!(site %in% fit_rivers$site)) 

summary_data <- bind_rows(new_site = all_river_data,
                          fit_site = fit_rivers,
                          pred_site = pred_data,.id='prediction_type')

# SAVE 
saveRDS(summary_data,'data/final_model_predictions_summary.rds')
saveRDS(summary_sr_data,'data/final_model_predictions_sr_summary.rds')

# Plots -------------------------------------------------------------------


## Plots to check ####
summary_data %>% filter(prediction_type=='pred_site') %>% 
  ggplot(aes(x=round.coast,ymin=S_q05,ymax=S_q95,y=S_q50,fill=spring_fish)) + 
  geom_ribbon(alpha=0.1) +
  geom_line(aes(col=spring_fish)) + 
  geom_errorbar(data=summary_data %>% filter(prediction_type=='fit_site'),
                aes(col=spring_fish)) +
  geom_point(data=summary_data %>% filter(prediction_type=='fit_site'),
             aes(col=spring_fish)) +
  scale_y_log10()

summary_data %>% filter(prediction_type=='pred_site') %>% 
  ggplot(aes(x=round.coast,ymin=H_q05,ymax=H_q95,y=H_q50,fill=spring_fish)) + 
  geom_ribbon(alpha=0.1) +
  geom_line(aes(col=spring_fish)) + 
  geom_errorbar(data=summary_data %>% filter(prediction_type=='fit_site'),
                aes(col=spring_fish)) +
  geom_point(data=summary_data %>% filter(prediction_type=='fit_site'),
             aes(col=spring_fish)) 

summary_data %>% 
  filter(prediction_type %in% c('fit_site','new_site'),
         spring_fish=='n',
         site!='River Lune') %>% 
  arrange(round.coast) %>% 
  mutate(site = fct_reorder(site,round.coast)) %>% 
  ggplot(aes(x=site,ymin=S_q05,lower=S_q25,middle=S_q50,
             upper=S_q75,ymax=S_q95))+
  geom_boxplot(stat = 'identity')


summary_sr_data %>% 
  ggplot(aes(x=S_pred,y=R_q50,ymin=R_q05,ymax=R_q95)) + 
  geom_ribbon(alpha=0.2) + 
  geom_line() + 
  geom_vline(data = summary_data %>% filter(prediction_type=='fit_site'),aes(xintercept=S_q50))+
  geom_point(data = sr_data,aes(x=S,y=R),inherit.aes = FALSE) +
  facet_wrap(~site, scales = 'free')


# Results for ms text -----------------------------------------------------

# 80 intervals for 55 degree equivalent sites
summary_data %>% filter(round.coast %in% fifty_five, spring_fish=="n")


# Ranges of egg targets for new sites
summary_data %>% filter(prediction_type == 'new_site') %>% 
  summarise(min = min(S_q50),
            max = max(S_q50),
            min_05 = min(S_q05),
            max_05 = max(S_q05),
            min_95 = min(S_q95),
            max_95 = max(S_q95))

# Summary of covariates
s <- rstan::summary(results$stan_fit,pars=c("alpha_S",'alpha_h'),probs=c(0.05,0.95))
as.data.frame(s$summary) %>% 
  mutate(param = rep(colnames(results$model_data$X),2)) %>% 
  select(param,mean,`5%`,`95%`)

# Bayesian R2 summaries
summary(R2_S)
summary(R2_h)
