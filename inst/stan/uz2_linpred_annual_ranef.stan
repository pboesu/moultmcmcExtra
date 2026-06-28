//
// Underhill-Zucchini Type 2 model with annual random effects on start date and duration.
// No individual random effects (use when within-season recaptures are absent).
// Type 2: samples are representative of old, moulting, and new plumage birds.

data {
  // control flags
  int<lower=0,upper=1> flat_prior;  // use flat priors on start date and duration intercepts?
  real beta_sd;                      // sd for non-flat prior on regression coefficients
  int<lower=0,upper=1> lumped;       // use lumped (Type 2L) likelihood?
  int<lower=0,upper=1> llik;         // calculate and return pointwise log-likelihood?
  // sample sizes
  int<lower=0> N_old;
  int<lower=0> N_moult;
  int<lower=0> N_new;
  // responses
  vector[N_old] old_dates;
  vector[N_moult] moult_dates;
  vector<lower=0,upper=1>[N_moult] moult_indices;
  vector[N_new] new_dates;
  // annual structure
  int N_years;                           // number of unique years
  int<lower=1, upper=N_years> year_factor[N_old+N_moult+N_new]; // year assignment per observation
  int<lower=0,upper=1> ranef_mu;        // include annual random effect on start date?
  int<lower=0,upper=1> ranef_tau;       // include annual random effect on duration?
  // predictors (design matrices ordered: old, moult, new)
  int N_pred_mu;
  matrix[N_old+N_moult+N_new, N_pred_mu] X_mu;
  int N_pred_tau;
  matrix[N_old+N_moult+N_new, N_pred_tau] X_tau;
  int N_pred_sigma;
  matrix[N_old+N_moult+N_new, N_pred_sigma] X_sigma;
}

parameters {
  vector[N_pred_mu] beta_mu;       // regression coefficients for start date
  vector[N_pred_tau] beta_tau;     // regression coefficients for duration
  vector[N_pred_sigma] beta_sigma; // regression coefficients for log(sigma)
  vector[N_years] u_year_mean;     // annual random intercepts on start date
  vector[N_years] u_year_duration; // annual random intercepts on duration
  real<lower=0> sd_year_mean;
  real<lower=0> sd_year_duration;
}

transformed parameters {
  real sigma_intercept = exp(beta_sigma[1]);
  // post-sweep: marginalise intercepts over annual random effects
  real beta_star = beta_mu[1] + ranef_mu  * mean(u_year_mean);
  real tau_star  = beta_tau[1] + ranef_tau * mean(u_year_duration);
  vector[N_years] u_year_mean_star     = ranef_mu  * (u_year_mean - mean(u_year_mean));
  vector[N_years] u_year_duration_star = ranef_tau * (u_year_duration - mean(u_year_duration));
}

model {
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_new] R;
  vector[N_old+N_moult+N_new] mu;
  vector[N_old+N_moult+N_new] tau;
  vector[N_old+N_moult+N_new] sigma;

  mu    = X_mu  * beta_mu  + ranef_mu  * u_year_mean[year_factor];
  tau   = X_tau * beta_tau + ranef_tau * u_year_duration[year_factor];
  sigma = exp(X_sigma * beta_sigma);

  if (lumped == 0) {
    for (i in 1:N_old)  P[i] = 1 - Phi((old_dates[i] - mu[i]) / sigma[i]);
    for (i in 1:N_new)  R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult]) /
                                     sigma[i + N_old + N_moult]);
  } else {
    for (i in 1:N_old)  P[i] = (1 - Phi((old_dates[i] - mu[i]) / sigma[i])) +
                                      Phi((old_dates[i] - tau[i] - mu[i]) / sigma[i]);
    for (i in 1:N_new)  R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult]) /
                                     sigma[i + N_old + N_moult]) +
                                (1 - Phi((new_dates[i] - mu[i + N_old + N_moult]) /
                                          sigma[i + N_old + N_moult]));
  }

  for (i in 1:N_moult)
    q[i] = log(tau[i + N_old]) +
           normal_lpdf((moult_dates[i] - moult_indices[i] * tau[i + N_old]) |
                        mu[i + N_old], sigma[i + N_old]);

  target += sum(log(P)) + sum(q) + sum(log(R));

  // annual random effect priors
  u_year_mean     ~ normal(0, sd_year_mean);
  u_year_duration ~ normal(0, sd_year_duration);

  // fixed effect priors
  if (flat_prior == 1) {
    beta_mu[1]  ~ uniform(-366, 366);
    beta_tau[1] ~ uniform(0, 366);
  } else {
    beta_mu[1]  ~ normal(150, 50) T[-366, 366];
    beta_tau[1] ~ normal(100, 30) T[0, 366];
  }
  if (beta_sd > 0) {
    if (N_pred_mu > 1) {
      for (i in 2:N_pred_mu)    beta_mu[i]    ~ normal(0, beta_sd);
    }
    if (N_pred_tau > 1) {
      for (i in 2:N_pred_tau)   beta_tau[i]   ~ normal(0, beta_sd);
    }
    if (N_pred_sigma > 1) {
      for (i in 2:N_pred_sigma) beta_sigma[i] ~ normal(0, beta_sd);
    }
  }
  beta_sigma[1]    ~ normal(0, 2);
  sd_year_mean     ~ normal(0, 2);
  sd_year_duration ~ normal(0, 2);
}

generated quantities {
  // NB: code duplication for likelihood is less than ideal — refactor to Stan functions?
  vector[N_pred_mu] beta_mu_out;
  vector[N_pred_tau] beta_tau_out;
  vector[(N_old + N_moult + N_new) * llik] log_lik;

  if (N_pred_mu > 1) {
    beta_mu_out = append_row(beta_star, beta_mu[2:N_pred_mu]);
  } else {
    beta_mu_out[1] = beta_star;
  }
  if (N_pred_tau > 1) {
    beta_tau_out = append_row(tau_star, beta_tau[2:N_pred_tau]);
  } else {
    beta_tau_out[1] = tau_star;
  }

  if (llik == 1) {
    vector[N_old] P;
    vector[N_moult] q;
    vector[N_new] R;
    vector[N_old+N_moult+N_new] mu;
    vector[N_old+N_moult+N_new] tau;
    vector[N_old+N_moult+N_new] sigma;

    mu    = X_mu  * beta_mu  + ranef_mu  * u_year_mean[year_factor];
    tau   = X_tau * beta_tau + ranef_tau * u_year_duration[year_factor];
    sigma = exp(X_sigma * beta_sigma);

    if (lumped == 0) {
      for (i in 1:N_old)  P[i] = 1 - Phi((old_dates[i] - mu[i]) / sigma[i]);
      for (i in 1:N_new)  R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult]) /
                                       sigma[i + N_old + N_moult]);
    } else {
      for (i in 1:N_old)  P[i] = (1 - Phi((old_dates[i] - mu[i]) / sigma[i])) +
                                        Phi((old_dates[i] - tau[i] - mu[i]) / sigma[i]);
      for (i in 1:N_new)  R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult]) /
                                       sigma[i + N_old + N_moult]) +
                                  (1 - Phi((new_dates[i] - mu[i + N_old + N_moult]) /
                                            sigma[i + N_old + N_moult]));
    }

    for (i in 1:N_moult)
      q[i] = log(tau[i + N_old]) +
             normal_lpdf((moult_dates[i] - moult_indices[i] * tau[i + N_old]) |
                          mu[i + N_old], sigma[i + N_old]);

    log_lik = append_row(append_row(log(P), q), log(R));
  }
}
