# Construct summary plots and statistics using best model from 
# cross validation

# for in sample and out of
# sample sites

# Libraries and config ----------------------------------------------------
rm(list = ls())
library(rgdal)
library(tmap)
library(tidyverse)
theme_set(tidybayes::theme_tidybayes())

source('R/utilities.R')
source('config.R')

# Functions ---------------------------------------------------------------

to_brp_tbl <- function(x,S_mat,h_mat){
  R_mat <- S_mat/(1-h_mat)
  rr_mat <- R_mat/S_mat
  
  x |> bind_cols(to_quantiles_df(S_mat,'S'),
                 to_quantiles_df(h_mat,'h'),
                 to_quantiles_df(R_mat,'R'),
                 to_quantiles_df(rr_mat,'rr'))
}

site_fix <- function(data){
  data %>% 
    mutate(site = gsub(pattern = "River ",replacement = "",site),
           site = gsub(pattern = "SAC",replacement = "",site),
           site = replace(site,site == "Dee (Kirkcudbrightshire)","K. Dee"))
}

# Load data and models ----------------------------------------------------

## Stock data ----
# Original data used to fit the model
sr_data <- readRDS('data/modelling/model_data.rds')
# Data for all rivers in assessment
# 173 areas
all_river_data = readRDS('./data/modelling/covariate_data.rds') |> filter(!is.na(land.cover.pc1))



## Models ----
# Best model from cv
results <- readRDS('results/fit_lat_sf pc1 lcpa.rds')
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

## Prediction data for hypothetical rivers ----

pred_data <- bind_rows(
  data.frame(lat=seq(min(all_river_data$lat),max(all_river_data$lat),length.out=100),
             cpa = 1,round.coast=0,land.cover.pc1=0,coast='e',spring_fish='n') |> mutate(site=paste0('pred_lat_',lat)),
  data.frame(cpa=seq(min(all_river_data$cpa),max(all_river_data$cpa),length.out=100),
             lat = 55,round.coast=0,land.cover.pc1=0,coast='e',spring_fish='n') |> mutate(site=paste0('pred_cpa_',cpa)),
  data.frame(land.cover.pc1=seq(min(all_river_data$land.cover.pc1),max(all_river_data$land.cover.pc1),length.out=100),
             lat = 55,round.coast=0,cpa=1,coast='e',spring_fish='n') |> mutate(site=paste0('pred_pc1_',land.cover.pc1))
)

#pred_data <- pred_data |> bind_rows(pred_data |> mutate(spring_fish = 'y',site=paste0(site,'_sf')))

# Construct design matrices for linear covariates and P-splines
X_pred <- create_new_model_matrix(old_data = sr_data,
                                  new_data = pred_data,
                                  covariates = results$covariates)
row.names(X_pred) <- pred_data$site
X_pred_sf <- X_pred
X_pred_sf[,'sf'] <- 1
row.names(X_pred_sf) <- paste0(pred_data$site,'_sf')

X_pred <- rbind(X_pred,X_pred_sf)

X_pred[!grepl('pred_cpa',rownames(X_pred)),'lcpa'] <- 0
X_pred[!grepl('pred_lat',rownames(X_pred)),'lat'] <- 0
X_pred[!grepl('pred_pc1',rownames(X_pred)),'pc1'] <- 0

S_pred <- predict_S(X_pred, alpha_S = alpha_S ,sigma_S = sigma_S)
H_pred <- predict_H(X_pred, alpha_H = alpha_H ,sigma_H = sigma_H)

pred_site_tbl <- as_tibble(X_pred,rownames = 'site') |> to_brp_tbl(S_pred,H_pred) |> mutate(sf=factor(sf))

## S* predictions ####
s_plot <- function(x,var,xlab=""){
  x |> 
    ggplot(aes(x={{var}},y=S_q50,ymin=S_q05,ymax=S_q95)) + 
    geom_ribbon(aes(fill=sf),alpha=0.1) +
    geom_line(aes(col=sf)) +
    labs(y="S* eggs/m^2",x=xlab) + 
    scale_y_log10(limits=c(0.1,50)) 
  
}

ggpubr::ggarrange(pred_site_tbl %>% filter(grepl('pred_cpa',site)) |> 
                    mutate(lcpa = sd(unique(log(sr_data$cpa)))*lcpa + mean(unique(log(sr_data$cpa)))) |> 
                    s_plot(lcpa,"CPA"),
                  pred_site_tbl |> filter(grepl('pred_lat',site)) |> 
                    mutate(lat = sd(unique(sr_data$lat))*lat + mean(unique(sr_data$lat))) |> 
                    s_plot(lat,"Lat"),
                  pred_site_tbl |> filter(grepl('pred_pc1',site)) |> 
                    mutate(pc1 = sd(unique(sr_data$land.cover.pc1))*pc1 + mean(unique(sr_data$land.cover.pc1))) |> 
                    s_plot(pc1,"LU"),
                  nrow=1,common.legend = TRUE)


## Posterior predictions for fit sites -----------------------------------

# SR predictions for fit sites
pred_sr <- sr_data %>% group_by(site) %>% 
  summarise(S_max = max(S),
            S_pred = seq(0,S_max,length.out = 100)) %>% 
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

fit_site_tbl <- as_tibble(results$model_data$X,rownames="site") %>% 
  left_join(sr_data |> group_by(site) |> 
              summarise(round.coast=first(round.coast),riverSalmon=first(area),cpa=first(cpa))) |>
  to_brp_tbl(t(S_star),t(H_star)) 

## Posterior predictions actual new rivers -------------------------------
X_all_rivers <- create_new_model_matrix(old_data = sr_data,
                                        new_data = all_river_data,
                                        covariates = results$covariates)

S_all_rivers <- predict_S(X_all_rivers, alpha_S = alpha_S ,sigma_S = sigma_S)
H_all_rivers <- predict_H(X_all_rivers, alpha_H = alpha_H ,sigma_H = sigma_H)

new_site_tbl <- as_tibble(X_all_rivers,rownames = 'site') |> 
  left_join(all_river_data |> mutate(site=river) |> select(site,round.coast,riverSalmon,cpa)) |>
  to_brp_tbl(S_all_rivers,H_all_rivers) |> 
  filter(!(site %in% fit_site_tbl$site)) 


# Plot 1. Ricker fits and prediction intervals ----

S_star_median <- fit_site_tbl %>% mutate(S_star = S_q50) %>% select(site,S_star) %>% 
  site_fix()

summary_sr_data %>% site_fix() |>
  ggplot(aes(x=S_pred,y=R_q50,ymin=R_q05,ymax=R_q95)) + 
  geom_ribbon(alpha=0.2) + 
  geom_line(size=.8) + 
  geom_abline(slope = 1,size=0.2,linetype=2) +
  geom_point(data = sr_data %>% site_fix(),aes(x=S,y=R),inherit.aes = FALSE) +
  geom_vline(data = S_star_median,aes(xintercept = S_star),size=.8,linetype=2) +
  facet_wrap(~site, scales = 'free') + 
  labs(y=expression(paste("Recruitment (eggs/",m^2,")")),
       x=expression(paste('Spawing stock (eggs/',m^2,")"))) 
ggsave('paper/ricker_fit.png',width=8,height=6)


# Plot 2. S* assessment areas ####

all_sites <- bind_rows(new_site = new_site_tbl,
                       fit_site = fit_site_tbl,.id='prediction_type') |>
  mutate(site=factor(site),sf=factor(sf))

all_sites %>% 
  filter(prediction_type %in% c('fit_site','new_site'),
         sf=='0',
         site!='River Lune') %>% site_fix() |>
  mutate(site = fct_reorder(site,S_q50)) %>% 
  ggplot(aes(x=site,y=S_q50,ymin=S_q05,ymax=S_q95,col=prediction_type))+
  geom_linerange() + 
  geom_linerange(aes(ymin=S_q25,ymax=S_q75),size=1) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y=NULL) + 
  scale_y_log10()

# Plot 3. HP Fit vs predictions ----

# Out of sample predictions for all rivers in cv process

fit_dir <- 'results/cv_lat_sf pc1 lcpa_fits/'
fit_fnames <- list.files(fit_dir)
mu_new_list <- list()
hp_new_list <- list()
for(fname in fit_fnames){
  fit <- readRDS(paste0(fit_dir,fname))
  r <- str_remove(fname,'.rds')
  mu_new_list[[r]] <- fit |> tidybayes::tidy_draws() |> tidybayes::spread_draws(mu_new[n],sigma,S_star_new,h_star_new) 
  hp_new_list[[r]] <- fit |> tidybayes::spread_draws(S_star_new,h_star_new) |> 
    mutate(rr_star_new = 1/(1-h_star_new)) |> 
    tidybayes::gather_draws(S_star_new,h_star_new,rr_star_new)
}


hp_samples <- bind_rows(hp_new_list,.id='site') |>  
  bind_rows(results$stan_fit |> tidybayes::tidy_draws() |> 
              tidybayes::spread_draws(S_star[site],h_star[site]) |>
              mutate(rr_star = 1/(1-h_star)) |>
              tidybayes::unspread_draws(S_star[site],h_star[site],rr_star[site]) |>
              tidybayes::gather_draws(S_star[site],h_star[site],rr_star[site]) |> 
              mutate(site=names(hp_new_list)[site])
  )

plt_data <-  hp_samples |>
  group_by(site,.variable) |>
  #mutate(.value = ifelse(grepl('S_star',.variable),log(.value),.value)) |>
  tidybayes::median_qi(.width = c(.5,.9)) |>
  tidyr::separate_wider_delim(.variable,delim = "_",names = c('param','temp','type'),too_few = 'align_start') |>
  tidyr::replace_na(list(type='Fit')) |>
  left_join(fit_site_tbl |> select(site,lat)) |>
  site_fix() |>
  mutate(param=paste0(param,'*'),
         type=replace(type,type=="new","Prediction")) |>
  mutate(site = forcats::fct_reorder(site,lat))


ggpubr::ggarrange(plt_data |> filter(param=="S*") |>
                    ggplot(aes(x=site,y=.value)) + 
                    labs(x=NULL) +tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=type),position=position_dodge(0.6)) + 
                    #facet_wrap(~param) +
                    labs(col='',y=expression(S*"* (eggs/"*m^2*")"))+
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    scale_y_log10(),
                  plt_data |> filter(param=="h*") |>
                    ggplot(aes(x=site,y=.value)) + 
                    #facet_wrap(~param) +
                    labs(x=NULL,col='',y='h*') +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    tidybayes::geom_pointinterval(aes(ymin = .lower, ymax = .upper,col=type),position=position_dodge(0.6)),
                  common.legend = TRUE,legend = "top",nrow=1)
ggsave('paper/hyperparam_predictions.png',width=6,height=3.5,scale = 1.2)

# Plot 4. Effect sizes ----
cov_map <- list(int = 'Int',sf='CT',lat='Lat',lcpa='CPA',pc1='LU')

cov_summary <- rstan::summary(results$stan_fit,pars=c("alpha_S",'alpha_h'),probs=c(0.05,.5,0.95))
cov_summary <- as.data.frame(cov_summary$summary) %>% 
  mutate(param = rep(colnames(results$model_data$X),2)) %>% 
  select(param,mean,`5%`,`50%`,`95%`)

cov_summary

# Plot covariate coefficients
as_tibble(cov_summary,rownames = "hyperparam") |>
  mutate(hyperparam = str_match(hyperparam,"alpha_(S|h)")[,2]) |>
  mutate(hyperparam = paste0(hyperparam,'*'),
         hyperparam = factor(hyperparam,levels=c('S*','h*'))) |>
  filter(!param %in% c('coast')) |> 
  mutate(param = unlist(cov_map[param]),
         param = factor(param,levels=cov_map)) |>
  ggplot(aes(x=param,y=`50%`,ymin=`5%`,ymax=`95%`)) + 
  geom_hline(yintercept = 0,col='grey') +
  geom_point() + 
  geom_linerange() + 
  labs(y='Effect size',x=NULL) +
  facet_wrap(~hyperparam)

ggsave(filename = 'paper/effect_size.png',width = 4,height = 4)

# Maps ----
## Map data ----

# Base map of Scotland for background of plot
scotland_map <- readOGR(dsn=paste0(GIS_DATA_PATH,"pan50.shp"))

# Polygons for all assessed areas
river_shapes <- readOGR(dsn=paste0(GIS_DATA_PATH,"AR2017_Boundaries_170613.shp"))

# Map a variable
map_value <- function(x,v_name, fit_point=TRUE, fit_text = FALSE, pal='Blues',style='pretty'){
  
  # Update Moriston to have Ness values
  x[x$site=="River Moriston SAC",c(v_name,"prediction_type")] <-  x[x$site == "River Ness",c(v_name,"prediction_type")]
  
  # Map labels
  label_data <- sr_data |> group_by(site) |> 
    summarise(lat = mean(lat),lon = mean(long),CT=first(spring_fish)) |> 
    site_fix() |> 
    arrange(site) |>
    mutate(label = paste0("(",letters[1:12],")"),
           CT = case_match(CT,
                           'y'~'upper',
                           'n'~'whole'
           ),
           ymod = ifelse(site%in%c("North Esk"),-1.5,0),
           xmod = ifelse(site=="Baddoch",-1,0))
  
  label_sf <- sf::st_as_sf(label_data,coords = c('lon','lat'))
  
  v <- unlist(x[,v_name])
  pt <- x$prediction_type
  names(v) <- names(pt) <- x$site
  
  river_shapes$v <- signif(v[as.character(river_shapes$river)],2)
  river_shapes$pt <- pt[as.character(river_shapes$river)]
  
  shetland_lim <- c(1060000,1220289)
  
  sc_map <-  tm_shape(scotland_map,ylim = c(530223,1060000)) + 
    tm_graticules(labels.inside.frame=TRUE,n.y=5,n.x=4) +
    tm_fill('grey80') +
    tm_layout(frame = FALSE, outer.margins = c(0.01, 0.01, 0.01, 0.01), inner.margins = c(0, 0, 0, 0),legend.bg.color = "white") +
    tm_shape(river_shapes) +
    tm_fill("v",title=v_name,palette=pal,style = style) +
    tm_borders(col='grey60')  #
  if(fit_point){
    sc_map <- sc_map +
      tm_shape(label_sf) + 
      tm_symbols(col = 'CT',size = .5,border.col='black')
  }
  if(fit_text){
    sc_map <- sc_map + 
      tm_text("site",just=c(0.5,-.5),ymod="ymod",xmod="xmod",shadow=TRUE)
  }
  inset_map <- tm_shape(scotland_map,
                        ylim = c(1060000,1220289),
                        xlim = c(380000, 470315)) + 
    tm_graticules(labels.inside.frame=FALSE,lines=FALSE,n.y=4,n.x=2) +
    tm_fill('grey80')
  
  sc_grob <- tmap_grob(sc_map)
  inset_grob <- tmap_grob(inset_map)
  
  cowplot::ggdraw() +
    cowplot::draw_plot(sc_grob,width=1,height=1) +
    cowplot::draw_plot(inset_grob,
                       width = 0.2, height = 0.2,
                       x = 0.75, y = 0.725) 
}


## Plot 5. Covariate map ----

map_data <- all_sites %>% mutate(coast_no = as.numeric(factor(coast)))
# ar .809
map_lcpa <- map_data |> mutate(`CPA` = cpa) |> map_value("CPA",pal='Oranges',style='log10')
map_lu <- map_data |> mutate(LU = pc1) |> map_value("LU",pal='Oranges',style='cont')
map_sp <- map_data |> mutate(SP = round.coast) |> map_value("SP",pal='Oranges',style='cont')
map_wa <- map_data |> mutate(`WA (mil)` = riverSalmon/1000000) |> map_value("WA (mil)",pal='Oranges',style='log10',fit_text = TRUE)

pl_save <- ggpubr::ggarrange(map_wa,map_lcpa,map_lu,map_sp,labels = c('a)','b)','c)','d)'))
ggsave("./paper/map_model_cov.png",plot = pl_save,width=10,height=12)

## Plot 6. Prediction map ----

map_s <- map_data |> mutate(`S*` = S_q50) |> map_value("S*",pal='Blues',fit_point=FALSE)
map_h <- map_data |> mutate(`h*` = h_q50) |> map_value("h*",pal='Blues',fit_point=FALSE)
#map_rr <- map_data |> mutate("R*/S*" = rr_q50) |> plot_value("R*/S*")

pl_save <- ggpubr::ggarrange(map_s,map_h,labels = c('a)','b)'),nrow = 1)
ggsave("./paper/map_model_pred.png",plot = pl_save,width=10,height=6)

# Results for ms text -----------------------------------------------------


# Ranges of egg targets for new sites
new_site_tbl %>% 
  summarise(min = min(S_q50),
            max = max(S_q50),
            min_05 = min(S_q05),
            max_05 = max(S_q05),
            min_95 = min(S_q95),
            max_95 = max(S_q95))

