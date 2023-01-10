//
// This Stan program defines the Underhill-Zucchini Type 5 repeated measures model
//  covariates need to be on the individual level not the observation model for logical consistncy

data {
  //responses
  int<lower=0,upper=1> flat_prior;//translated logical use flat prior on start and duration?
  int<lower=0> N_ind;//number of individuals
  int<lower=0> N_ind_rep;//number of recaptured individuals
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  int<lower=0> Nobs_replicated;//number of observations from individuals with repeat measures
  int<lower=0> Nobs_not_replicated_old;//number of old obs from unreplicated individuals
  int<lower=0> Nobs_not_replicated_moult;//number of moult obs from unreplicated individuals
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0>individual[N_moult+N_old]; //individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
  int<lower=0>replicated[Nobs_replicated];//indices of obs from individuals with repeat measures
  int<lower=0>not_replicated[(N_moult + N_old) - Nobs_replicated];//indices of obs from individuals without repeat measures
  int<lower=0>not_replicated_old[Nobs_not_replicated_old];//indices of old obs from individuals without repeat measures
  int<lower=0>not_replicated_moult[Nobs_not_replicated_moult];//indices of moult obs from individuals without repeat measures. NB these are relative to N_moult observations, NOT the full dataset
  int<lower=0>is_replicated[N_ind];
  int<lower=0>replicated_individuals[N_ind_rep];//individual id's that are replicated
  int<lower=0>replicated_id_obs[N_moult+N_old];//new running id from 1:N_ind_rep to index the individual start date estimates. Unreplicated individuals receive a zero value
  //int<lower=0> N_new;//K
  //vector[N_new] new_dates;//v_k
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_old+N_moult,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
  int year_factor[N_old+N_moult];//year factor for each observation
  int N_years;//numebr of unique year_factor values
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  //real<lower=0> tau;//moult duration
  //real<lower=0> sigma;//population varianmoult start date
  vector[N_pred_mu] beta_mu;//regression coefficients for start date
  vector[N_pred_tau] beta_tau;//regression coefficients for duration
  vector[N_pred_sigma] beta_sigma;//regression coefficients for sigma start date
  vector[N_ind_rep] mu_ind;//individual effect on sigma start date
  real<lower=0> sigma_mu_ind;//?residual variance in regression of score on date within individuals
  vector[N_years] u_year_mean;//raneff intercept on start date
  vector[N_years] u_year_duration;//raneff intercetpt on duration
  real<lower = 0> sd_year_mean;//raneff std dev
  real<lower = 0> sd_year_duration;
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
  //post-sweep random effects
  real beta_star = beta_mu[1] + mean(mu_ind) + mean(u_year_mean);
  real tau_star = beta_tau[1] + mean(u_year_duration);
  vector[N_ind_rep] mu_ind_star = mu_ind - mean(mu_ind);
  vector[N_years] u_year_mean_star = u_year_mean - mean(u_year_mean);
  vector[N_years] u_year_duration_star = u_year_duration - mean(u_year_duration);
  real finite_sd = sd(mu_ind_star);
}

// The model to be estimated.
model {
  vector[N_old] Rt;
  vector[N_old] P;
  vector[N_moult] Ru;
  vector[N_moult] q;
  vector[N_old+N_moult] mu;//start date lin pred
  vector[N_old+N_moult] tau;//duration lin pred
  vector[N_old+N_moult] sigma;//duration lin pred

  mu = X_mu * beta_mu + u_year_mean[year_factor];
//  print(mu);
  tau = X_tau * beta_tau + u_year_duration[year_factor];
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

for (i in 1:N_old) {
	if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  P[i] = 1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	} else {//standard likelihood for Type 5 model
    P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
    Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
	}
}
for (i in 1:N_moult){
  if (is_replicated[individual[i + N_old]] == 1) {
   q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old] + mu_ind[replicated_id_obs[i + N_old]], sigma_mu_ind);//replicated individuals. NB - indexing looks messy because i runs from 1:N_moult, but the function uses both vectors of the total dataset (1:(N_old+N_moult) and the moult dataset (1:N_moult))
  } else {
   Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
   q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
  }
}
mu_ind ~ normal(0, sigma[individual_first_index][replicated_individuals]);//

//raneffs
u_year_mean ~ normal(0, sd_year_mean);
u_year_duration ~ normal(0, sd_year_duration);

target += sum(log(P)) - sum(log1m(Rt[not_replicated_old])) + sum(q) - sum(log1m(Ru[not_replicated_moult]));

//priors
if (flat_prior == 1) {
 beta_mu[1] ~ uniform(0,366);
 beta_tau[1] ~ uniform(0,366);
} else {
 beta_mu[1] ~ normal(150,50)T[0,366];
 beta_tau[1] ~ normal(100,30)T[0,366];
}
beta_sigma[1] ~ normal(0,5);
sigma_mu_ind ~ normal(0,1);

sd_year_mean ~ normal(0,2);
sd_year_duration ~ normal(0,2);
}


generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
  vector[N_old] Rt;
  vector[N_old] P;
  vector[N_moult] Ru;
  vector[N_moult] q;
  vector[N_old+N_moult] mu;//start date lin pred
  vector[N_old+N_moult] tau;//duration lin pred
  vector[N_old+N_moult] sigma;//duration lin pred
  vector[N_old+N_moult] log_lik;
  vector[N_ind_rep] mu_ind_out;//individual intercepts for output
  vector[N_pred_mu] beta_mu_out;//post-swept regression coefficients for start date
  vector[N_pred_tau] beta_tau_out;//post-swept regression coefficients for duration

  if (N_pred_mu > 1){
    beta_mu_out = append_row(beta_star,beta_mu[2:N_pred_mu]);// collate post-swept intercept with remaining
  } else {
    beta_mu_out[1] = beta_star;// intercept-only model
  }

  if (N_pred_tau > 1){
    beta_tau_out = append_row(tau_star,beta_tau[2:N_pred_tau]);// collate post-swept intercept with remaining
  } else {
    beta_tau_out[1] = tau_star;// intercept-only model
  }


  mu_ind_out = mu_ind_star + beta_star;

    mu = X_mu * beta_mu + u_year_mean[year_factor];
//  print(mu);
  tau = X_tau * beta_tau + u_year_duration[year_factor];
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

for (i in 1:N_old) {
	if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  P[i] = 1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  log_lik[i] = log(P[i]);
	} else {//standard likelihood for Type 5 model
    P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
    Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
    log_lik[i] = log(P[i])- log1m(Rt[i]);
	}

}
for (i in 1:N_moult){
  if (is_replicated[individual[i + N_old]] == 1) {
   q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old] + mu_ind[replicated_id_obs[i + N_old]], sigma_mu_ind);//replicated individuals. NB - indexing looks messy because i runs from 1:N_moult, but the function uses both vectors of the total dataset (1:(N_old+N_moult) and the moult dataset (1:N_moult))
   log_lik[N_old+i]=q[i];
  } else {
   Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
   q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);
   log_lik[N_old+i]=q[i] - log1m(Ru[i]);//N.B. unlike P and R this returns a log density
  }
}
//target += sum(log(P)) - sum(log1m(Rt[not_replicated_old])) + sum(q) - sum(log1m(Ru[not_replicated_moult]));
//log_lik = append_row((log(P) - log1m(Rt)), (q - log1m(Ru)));//NAs arise because Rt and P / q and Ru are not the same length!? - slower but less error prone to assemble log_lik via loop ? - possibly a more elegant solution via multiplication with a 0-1 indicator - but currently replication is tracked by individual not by obs?!
}
