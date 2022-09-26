
library(dplyr)

OMIT_K_DEE <- FALSE

CV_PATH <- 'results/'
N_FOLDS <- 12
if(OMIT_K_DEE){
 CV_PATH <- 'results/omit_k_dee/'
 N_FOLDS <- 11
}

cv_file_list <- list.files(path = CV_PATH,pattern = "cv.*\\.rds",
                           full.names = TRUE)

cv_list <- list()
lpd_mat <- list()
lpd_by_river_mat <- list()
for(f in cv_file_list){
  
  print(paste("Loading cv results",f))
  cv_results <- readRDS(f)
  cv_results$folds <- cv_results$folds[order(names(cv_results$folds))]
  if(length(cv_results$folds)<N_FOLDS){
    print(paste("cv results",f,"incomplete...skipping"))
    next()
  }
  
  
  lpd_list <- list()
  lpd_by_river_list <- list()
  # Iterate over folds
  for(r in names(cv_results$folds)){
    # Get the fit for this fold
    cv_fit <- as.data.frame(cv_results$folds[[r]])
    
    # Calculate the lpd for each point
    lpd_list[[r]] <- log(colMeans(exp(as.matrix(cv_fit))))
    
    # Calculate the lpd for all data as a "single data point"
    lpd_by_river_list[[r]] <- log(mean(exp(as.matrix(cv_fit))))
  }
  
  # Pointwise lpd
  # Weight data points equally in cv
  lpd <- unlist(lpd_list)
  lpd_mat[[f]] <- lpd
  logo_cv = sum(unlist(lpd))
  logo_cv_se = sqrt(length(lpd)) * sd(lpd)
  
  # River wise lpd 
  # Weight rivers equally in cv
  lpd_by_river <- unlist(lpd_by_river_list)
  lpd_by_river_mat[[f]] <- lpd_by_river
  logo_cv_2 = sum(lpd_by_river)
  logo_cv_2_se = sqrt(length(lpd_by_river)) * sd(lpd_by_river)
  
  cv_list[[f]] <- data.frame(filename=f,
                             model=cv_results$model,
                             covariates=stringr::str_c(cv_results$covariates,collapse='_'),
                             sp=cv_results$sp,
                             logo_cv,
                             logo_cv_se,
                             logo_cv_2,
                             logo_cv_2_se)
  
  
}


elpd_diff_function <- function(elpd_mat){

  # Which is the best model  
  max_elpd <- which(apply(elpd_mat,1,sum) == max(apply(elpd_mat,1,sum)))
  # Calculate diff relative to best model
  elpd_diff <- apply(t(elpd_mat) - elpd_mat[max_elpd,],2,sum)
  # Calculate sd of difference relative to best model
  elpd_diff_se <- sqrt(ncol(elpd_mat))*apply(t(elpd_mat) - elpd_mat[max_elpd,],2,sd)
  
  return(data.frame(elpd_diff,elpd_diff_se))
}

# Calculate difference on all models
lpd_mat <- as.matrix(bind_rows(lpd_mat))
lpd_by_river_mat <- as.matrix(bind_rows(lpd_by_river_mat))

elpd_diff <- elpd_diff_function(lpd_mat)
elpd_by_river_diff <- elpd_diff_function(lpd_by_river_mat)
names(elpd_by_river_diff) <- paste0(names(elpd_by_river_diff),"_by_river")


spatial_names <- c('None','SP','Lat','Coast + Lat') 

cv_summary <- dplyr::bind_rows(cv_list) %>% 
  bind_cols(elpd_diff,elpd_by_river_diff) %>% 
  arrange(desc(elpd_diff_by_river)) %>% 
  mutate(spatial_type = replace(sp,str_detect(covariates,'lat'),2),
         spatial_type = spatial_names[spatial_type+1]) %>% 
  mutate(covariates = str_replace(covariates,"_lat_*",""),
         covariates = str_replace(covariates,"lat_*",""))

cv_summary


saveRDS(cv_summary,file = paste0(CV_PATH,'cross_validation_summary.rds'))

# Plot by river comparison with standard error
cv_table <- cv_summary %>% 
  mutate(covariates = gsub('_',', ', covariates),
         covariates = sub('lat','Lat',covariates),
         covariates = sub('lcpa','CPA',covariates), 
         covariates = sub('coast','Coast',covariates), 
         covariates = sub('sf','CT',covariates), 
         covariates = sub('pc1','LU',covariates),
         covariates = replace(covariates,is.na(covariates),'')) %>% 
  filter(!grepl('Coast \\+ Lat',spatial_type)) %>% 
  mutate(rank_2 = nrow(.) - rank(elpd_diff_by_river) + 1) %>%  
  mutate(covariates = fct_inorder(covariates)) %>% 
  mutate(elpd_diff_min = elpd_diff - elpd_diff_se,
         elpd_diff_max = elpd_diff + elpd_diff_se) %>% 
  mutate(elpd_diff_2_min = elpd_diff_by_river - elpd_diff_se_by_river,
         elpd_diff_2_max = elpd_diff_by_river + elpd_diff_se_by_river) %>% 
  arrange(-1*elpd_diff_by_river) %>% 
  mutate(covariates = fct_inorder(covariates))

# Plot river level cv score 
pl2 <- ggplot(cv_table,
              aes(covariates,elpd_diff_by_river,
                  ymin = elpd_diff_2_min,
                  ymax = elpd_diff_2_max,
                  col=spatial_type,
                  shape=spatial_type)) +
  geom_point(size=2,position = position_dodge(0.5)) + 
  geom_errorbar(position = position_dodge(0.5),width=0.5,size=.7) +
  scale_color_brewer(palette = 'Set2') + 
  scale_color_discrete(name='Spatial structure') + 
  scale_shape_discrete(name='Spatial structure') +
  labs(x='Among-stock explanatory variables',
       y=expression(Delta*elpd_[CV*"*"])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pl2

# ggplot(cv_table,
#        aes(covariates,elpd_diff,
#            ymin = elpd_diff_min,
#            ymax = elpd_diff_max,  
#            col=spatial_type,
#            shape=spatial_type)) +
#   geom_point(size=2,position = position_dodge(0.5)) + 
#   geom_errorbar(position = position_dodge(0.5),width=0.5,size=.7) +
#   scale_color_brewer(palette = 'Set2') + 
#   scale_color_discrete(name='Spatial structure') + 
#   scale_shape_discrete(name='Spatial structure') +
#   labs(x='Among-stock explanatory variables',
#        y=expression(Delta*elpd_[CV*"*"])) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


