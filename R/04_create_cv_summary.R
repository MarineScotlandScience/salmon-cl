
library(dplyr)
library(ggplot2)
library(stringr)

OMIT_K_DEE <- FALSE

CV_PATH <- 'results/'
N_FOLDS <- 12
if(OMIT_K_DEE){
  CV_PATH <- 'results/omit_k_dee/'
  N_FOLDS <- 11
}

# Functions ####
log_mean_exp <- function(x){log(mean(exp(x)))}

get_spatial_type <- function(x,y){
  if(y) return("SP")
  if(str_detect(x,'coast_lat')) return("Coast + Lat")
  if(str_detect(x,'coast')) return("Coast")
  if(str_detect(x,'lat')) return("Lat")
  if(str_detect(x,'lsp')) return("log(SP)")
  return("None")
}


cv_file_list <- list.files(path = CV_PATH,pattern = "cv.*\\.rds",full.names = TRUE)

cv_list <- list()
loo_list <- list()
for(f in cv_file_list){
  
  print(paste("Loading cv results",f))
  cv_results <- readRDS(f)
  full_results <- readRDS(stringr::str_replace(f,"cv_","fit_"))
  cv_results$folds <- cv_results$folds[order(names(cv_results$folds))]
  if(length(cv_results$folds)<N_FOLDS){
    print(paste("cv results",f,"incomplete...skipping"))
    next()
  }
  
  lppds <- list()
  
  # Iterate over folds
  for(r in names(cv_results$folds)){
    # Get the fit for this fold
    # S X N (held out) log likelihoods
    cv_fit <- as.data.frame(cv_results$folds[[r]])
    
    # Joint log posterior predictive densities
    # Sum over all predictive values for this group
    # lppds (S samples for this river)
    lppds[[r]] <- rowSums(cv_fit)
  }
  
  # lppds from brms method
  # S x G df
  lppds <- do.call(cbind,lppds)
  # Expected log predictive densities
  # (average over samples)
  elpds <- apply(lppds,2,log_mean_exp)
  
  # compute effective number of parameters
  # (requires deriving likelihood for model fit on all data)
  log_lik <- rstan::extract(full_results$stan_fit,'log_lik')[[1]]
  ll_full_marg <- matrix(nrow = nrow(log_lik), ncol = length(names(cv_results$folds)))
  
  # Derive joint log likelihood over rivers
  for(r_idx in 1:N_FOLDS){
    fold_idx = full_results$model_data$river == r_idx
    ll_full_marg[,r_idx] <- rowSums(log_lik[,fold_idx,drop=FALSE])
  }
  # Average over samples
  lpds <- apply(ll_full_marg, 2, log_mean_exp)
  # Estimated number of parameters
  ps <- lpds - elpds
  
  # Construct a loo object manually
  make_loo <- function(elpds,ps){
    pointwise <- cbind(elpd_kfold = elpds, p_kfold = ps, kfoldic = -2 * elpds)
    est <- colSums(pointwise)
    se_est <- sqrt(nrow(pointwise) * apply(pointwise, 2, var))
    estimates <- cbind(Estimate = est, SE = se_est)
    
    rownames(estimates) <- colnames(pointwise)
    out <- list(estimates=estimates, pointwise=pointwise)
    out <- structure(out, class = c("kfold", "loo"))
  }
  
  loo_list[[f]] <- make_loo(elpds,ps)
  
  cv_list[[f]] <- data.frame(filename=f,
                             model=cv_results$model,
                             covariates=stringr::str_c(cv_results$covariates,collapse='_'),
                             sp=cv_results$sp)
  
  
}




# Sanitise model detailes table and add comparison
cv_table <- dplyr::bind_rows(cv_list) %>% 
  rowwise() |>
  mutate(spatial_type = get_spatial_type(covariates,sp)) |>
  ungroup() |>
  mutate(covariates = stringr::str_replace(covariates,"_lat_*","_"),
         covariates = stringr::str_replace(covariates,"lat_*",""),
         covariates = stringr::str_replace(covariates,"coast_*",""),
         covariates = stringr::str_replace(covariates,"lsp_*",""),
         covariates = gsub('_',', ', covariates),
         covariates = sub('lcpa','CPA',covariates), 
         covariates = sub('sf','CT',covariates), 
         covariates = sub('pc1','LU',covariates),
         covariates = replace(covariates,is.na(covariates),''))  

#### Yates paper approach
#
# Yates, Luke A., Zach Aandahl, Shane A. Richards, and Barry W. Brook. 2023. 
# “Cross Validation for Model Selection: A Review with Examples from Ecology.” 
# Ecological Monographs 93(1): e1557. https://doi. org/10.1002/ecm.1557

cv_plot <- function(logo_df){
  spatial_levels <-  c('Coast + Lat','Coast','None','Lat','SP','log(SP)')
  covariate_levels <- c('CT', 'LU, CT','CPA, CT','CPA, LU, CT')
  
  logo_df <- logo_df |> mutate(spatial_type = factor(spatial_type,levels = spatial_levels),
                               covariates = factor(covariates,levels=covariate_levels)) |>
    mutate(ose_model = ifelse(metric + se >= max(metric),1,0))
  
  
  logo_df |> mutate(spatial_type = factor(spatial_type,levels = spatial_levels),
                    covariates = factor(covariates,levels=covariate_levels))|>
    ggplot(aes(x = covariates)) +
    geom_point(aes(y=metric,size=ose_model),col='black', shape = 1, show.legend = F) +
    geom_point(aes(y = metric, col = spatial_type), size = 2, show.legend = F) +
    geom_linerange(aes(ymin = metric - se, ymax = metric + se, col = spatial_type), show.legend = F) +
    theme_classic() +
    theme(strip.placement = "outside", panel.border = element_blank(), 
          strip.background = element_rect(), axis.ticks.x = element_blank(),
          strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90"),
          axis.title.y = element_text(size = 8),axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~spatial_type,nrow = 1, strip.position = "top") +
    scale_color_brewer(type = "div", palette = "Dark2") + 
    scale_size(range=c(0,4)) +
    labs(y=expression(hat(S)[CV]),x=NULL)
}

make_plot_data <- function(metric_data, levels = names(metric_data)){
  best_model <- metric_data %>% purrr::map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% purrr::map_dbl(mean),
         metric_diff = metric - max(metric),
         se = metric_data %>% purrr::map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_diff = metric_data %>% purrr::map(~ .x - metric_data[[best_model]]) %>% purrr::map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model]) 
}

add_cv_details <- function(x,cv_table){ 
  x |> mutate(covariates = factor(cv_table$covariates),
              spatial_type = cv_table$spatial_type,
              filename = cv_table$filename,
              se_ose = se,
              se = se_mod)
}

# Pointwise (stock wise) elpd list to data frame
logo_pointwise <- loo_list |> 
  purrr::map("pointwise") |> 
  purrr::map_dfc(~.[,"elpd_kfold"]) |>
  relocate(all_of(cv_table$filename))

# Refine candidate list
#
# Remove models with coast 
# Enforce catchment type
ct_only_no_coast_cv_table <- cv_table |> filter(str_detect(covariates,"CT")) |> filter(!str_detect(spatial_type,"Coast"))
logo_pointwise |> select(ct_only_no_coast_cv_table$filename) |> make_plot_data() |>  add_cv_details(ct_only_no_coast_cv_table) |> cv_plot()

# Save plot
if(!dir.exists('paper')) dir.create('paper')
ggsave('paper/model_select.png',width=4,height=4)



