data{
  int<lower=0> N;
  int<lower=0> N_rivers;
  int<lower=0> N_X;
  vector[N] R; // Raw recruitment counts
  vector[N] S; // Raw stock counts
  int river[N];
  matrix[N_rivers,N_X] X; // Design matrix - fixed to 3 covariates, but set to 0 where unused
  int num_basis;
  matrix[N_rivers, num_basis] Y;  //  Model matrix for b-splines on river
  
   // Data for cross validation
  int<lower=0> N_holdout;
  vector[N_holdout] R_holdout;
  vector[N_holdout] S_holdout;
  vector[num_basis] Y_holdout;
  vector[N_X] X_holdout;
}

parameters{
  real<lower=0> sigma; // Error in Ricker curve
  
  vector[N_X] alpha_S;                        // Coefficients for design matrix on S_star
  real<lower=0> sigma_S;                    // Hyper param S.D. of lognormal dist for S_star
  
  vector[N_X] alpha_h;                        // Coefficients for design matrix on h_star
  real<lower=0>sigma_h;                       // Hyper param spread of dist for h_Star
  
  vector[N_rivers] mu_h_tilde;              // Parameter for non-centered implementation
  vector[N_rivers] log_S_tilde;             // Parameter for non-centered implementation
    
  vector[num_basis] beta_raw_S;      // B spline coefficients
  vector[num_basis] beta_raw_h;      // B spline coefficients
  real<lower=0> tau_S;            // Penalisation parameter for smooth
  real<lower=0> tau_h;            // Penalisation parameter for smooth

}

transformed parameters{
  // For constraint on wiggliness, 
  // use prior tau to reduce difference
  // in consecutive coefficients -> closer to
  // straight line
  vector[num_basis] beta_S;
  vector[num_basis] beta_h;
  
  vector<lower=0>[N_rivers] S_star;         // The stock at MSY point
  vector<lower=0,upper=1>[N_rivers] h_star; // Exploitation rate at MSY point
  real mu_h[N_rivers];
  real log_S[N_rivers];
  
  // Second order random walk prior
  beta_S[1] = beta_raw_S[1];
  beta_h[1] = beta_raw_h[1];
  beta_S[2] = beta_raw_S[2];
  beta_h[2] = beta_raw_h[2];
  
  for (i in 3:num_basis){
    beta_S[i] = 2*beta_S[i-1] - beta_S[i-2] + beta_raw_S[i]*tau_S;
    beta_h[i] = 2*beta_h[i-1] - beta_h[i-2] + beta_raw_h[i]*tau_h;
  }
  
  // Non-centered parameterisation of hyper parameters
  for(r in 1:N_rivers){
    mu_h[r] = X[r,]*alpha_h + Y[r,]*beta_h + sigma_h*mu_h_tilde[r];
    h_star[r] = inv_logit(mu_h[r]);
    
    log_S[r] = X[r,]*alpha_S + Y[r,]*beta_S + sigma_S*log_S_tilde[r];
    S_star[r] = exp(log_S[r]);
  }
}


model{
  // Priors
  sigma ~ cauchy(0,1);
  sigma_h ~ cauchy(0,0.5);
  sigma_S ~ cauchy(0,0.5);
  alpha_h ~ normal(0,1);
  alpha_S ~ normal(0,1);
  
  beta_raw_S ~ normal(0, 1);
  beta_raw_h ~ normal(0, 1);
  //tau_S ~ cauchy(0,1);
  //tau_h ~ cauchy(0,1);
  tau_S ~ normal(0,1);
  tau_h ~ normal(0,1);
 
  log_S_tilde ~ normal(0,1);
  mu_h_tilde ~ normal(0,1);

  // Model
  for(n in 1:N){
      R[n] ~ lognormal(h_star[river[n]] - log(1-h_star[river[n]]) + log(S[n]) - 
      h_star[river[n]]/S_star[river[n]] * S[n],sigma);
  }
}



generated quantities{
  vector[N] log_lik;

  // CV output
  vector[N_holdout] log_lik_new;
  vector[N_holdout] mu_new;
  real S_star_new;
  real h_star_new;
  real mu_h_new;
  real log_S_new;

  // Model store log_lik
  for(n in 1:N){
      log_lik[n] = lognormal_lpdf(R[n] | h_star[river[n]] - log(1-h_star[river[n]]) + log(S[n]) - h_star[river[n]]/S_star[river[n]] * S[n],sigma);
  }

  // Predictions from heirarchical model
  log_S_new = normal_rng(dot_product(X_holdout,alpha_S) + dot_product(Y_holdout,beta_S),sigma_S);
  S_star_new = exp(log_S_new);
  mu_h_new = normal_rng(dot_product(X_holdout,alpha_h) + dot_product(Y_holdout,beta_h),sigma_h);
  //mu_h_new = normal_rng(dot_product(X_holdout,alpha_h),sigma_h);
  h_star_new = inv_logit(mu_h_new);

  // Calculate CV score for out of sample data
  for(n in 1:N_holdout){
          mu_new[n] = h_star_new - log(1-h_star_new) + log(S_holdout[n]) - h_star_new/S_star_new * S_holdout[n];
          log_lik_new[n] = lognormal_lpdf(R_holdout[n] |mu_new[n] ,sigma);
  }

}


