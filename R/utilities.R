# Utility functions shared across multiple scripts


create_model_data <- function(sr_data,all_river_data,covariates,sp,holdout_river=NULL){ 
  
  # Define the model matrix, 
  # X[,1] = intercept
  # X[,2] = log catch per area
  # X[,3] = land pc1
  # X[,4] = spring fish
  # X[,5] = lat
  # X[,6] = coast
  X = create_model_matrix(sr_data,covariates)

  # Define the round coast bspline matrix here
  # B-Spline for distance round the coast, non-linear spatial relationship between 
  # rivers
  # Unscaled range defines the extent of the knots. In order to extrapolate
  # to the entire data set, need to specify the range for all rivers, including
  # those without data
  num_basis <- 9
  unscaled_range = range(all_river_data$round.coast,sr_data$round.coast)
  sr_summary <- sr_data %>% group_by(site) %>% summarise(round.coast = first(round.coast))
  Y = create_coast_spline_matrix(sr_summary$round.coast,
                                 s_range=unscaled_range,
                                 num_basis = num_basis)
  # If we are not using the smoother, set coefficents to 0
  if(!sp)  Y[,] = 0
  
  if(is.null(holdout_river)){
    # Not running CV
    
    # Need dummy data for CV, to allow model to run
    # Dummy data, normally used for cross validation phase
    holdout_river <- 'Baddoch'
    holdout_index <- which(rownames(X) == holdout_river)
    
    sr_data_holdout <- sr_data %>% filter(site == holdout_river)
    X_holdout <- X[holdout_index,]
    Y_holdout <- Y[holdout_index,]

  }else{
    # River being heldout for cv
    holdout_index <- which(rownames(X) == holdout_river)
    
    # Holdout data
    sr_data_holdout <- sr_data %>% filter(site == holdout_river)
    X_holdout <- X[holdout_index,]
    Y_holdout <- Y[holdout_index,]
    
    # Remove holdout data from model
    sr_data <- sr_data %>% filter(site != holdout_river)
    X <- X[-holdout_index,]
    Y <- Y[-holdout_index,]
  }
  
  # Define data for the stan model
  model_data = list('N'=length(sr_data$S),
                    'R'=sr_data$R,
                    'S'=sr_data$S,
                    'N_rivers'=length(unique(sr_data$site)),
                    'N_X' = ncol(X),
                    'river'=as.numeric(factor(sr_data$site)),
                    'num_basis'= num_basis,
                    'X'=X,
                    'Y'=Y,
                    'S_holdout'=sr_data_holdout$S,
                    'R_holdout'=sr_data_holdout$R,
                    'X_holdout'=X_holdout,
                    'Y_holdout'=Y_holdout,
                    'N_holdout'=length(sr_data_holdout$S))
  
  return(model_data)
}


# Define the model matrix, these are now hardcoded to the only
# covariates used in final models
# All covariates are scaled
create_model_matrix = function(sr_data,covariates){
  int = rep(1,length(unique(sr_data$site)))
  lcpa = scale(tapply(log(sr_data$cpa),sr_data$site,min))[,1]
  pc1 = scale(tapply(sr_data$land.cover.pc1,sr_data$site,min))[,1]
  sf = (tapply(as.numeric(sr_data$spring_fish),sr_data$site,min) - 1) 
  lat = scale(tapply(sr_data$lat,sr_data$site,min))[,1]
  coast = (tapply(as.numeric(sr_data$coast),sr_data$site,min) - 1) 
  
  X <- cbind(int,lcpa,pc1,sf,lat,coast)
  
  # Set unused covariate data to 0
  X[,!(c('int','lcpa','pc1','sf','lat','coast') %in% c('int',covariates))] = 0
  return(X)
}

# Creates the design matrix X for a bspline given x
create_coast_spline_matrix = function(round_coast,s_range,num_basis=12,order=3){
  
  
  # Spline configurations
  # For cubic spline
  # Order = m+1 = 3
  m = order - 1
  # Parameters (number of basis functions and coefficients)
  k = num_basis
  n_knots = m + k + 2
  
  # Create knots
  xk <- create_knots(num_basis,s_range,m)
  
  # Smooth matrix
  # X[N,k] 
  X = matrix(0,nrow = length(round_coast),ncol=k)
  
  # Construct bsplines for the knots for this river
  X = sapply(1:k,bspline,x=round_coast,xk=xk,m=m)
  
  return(X)
}

# Create equally spaced knots given range and number
# of basis required
create_knots <- function(num_basis,s_range,m=2){
  n_knots = m + num_basis + 2
  
  # Define the knots 
  # Need first and last m+1 knots to exceed range of x
  xk = seq(min(s_range),max(s_range),length.out = n_knots - 2*(m+1))
  
  # Adding m+1 knots to either end with same distance between
  # knots
  xk = sort(c(xk[1] - (1:(m+1))*diff(xk)[1],xk[length(xk)] + (1:(m+1))*diff(xk)[1],xk))
  
  return(xk)
}

bspline = function(x,xk,i,m=2){
  if(m==-1){ # Base of recursion   
    res = as.numeric(x<xk[i+1] & x>=xk[i])
  } else {
    # Weighting of lower order splines
    z0 = (x - xk[i]) / (xk[i+m+1] - xk[i])
    z1 = (xk[i+m+2] - x)/ (xk[i+m+2] - xk[i+1])
    
    # Recursion to lower order (m-1) bspline
    res = z0*bspline(x,xk,i,m-1) + z1*bspline(x,xk,i+1,m-1)
    
  }
  res
}


InvLogit = function(x){ 
  # Return the inverse logit e^x/(1+e^x) of a numeric object.
  #
  # Args:
  #   x: The numeric object tofor which the inverse logit is to be calculated
  #
  # Returns:
  #   The inverse logit of x
  return(exp(x)/(1+exp(x))) 
}
