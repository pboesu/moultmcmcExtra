#' Bayesian inference for the Type 2 moult model with repeat measures and annual random intercepts
#'
#' Fits a Type 2 (Underhill-Zucchini) moult model that includes individual-level random
#' intercepts for within-season recaptures and annual random intercepts on start date and
#' duration. Type 2 data include birds in old plumage, active moult, and new plumage.
#'
#' @export
#' @param date_column the name of the column in `data` containing sampling dates, encoded
#'   as days since an arbitrary reference date (numeric vector)
#' @param moult_index_column the name of the column in `data` containing moult indices,
#'   i.e. a numeric vector of (linearised) moult scores: 0 = old plumage, 1 = new plumage
#' @param id_column factor identifier for within-season individual recaptures
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param year_factor_column the name of a factor column in `data` containing year identifiers
#' @param flat_prior logical; use flat (uniform) priors on start date and duration
#'   intercepts? Defaults to FALSE (weakly informative normal priors)
#' @param lumped logical; use the lumped (Type 2L) likelihood where pre- and post-moult
#'   birds are treated as indistinguishable? Defaults to FALSE
#' @param use_phi_approx logical; use `Phi_approx` instead of `Phi` for computational
#'   efficiency? Defaults to FALSE
#' @param raneff_components character vector specifying which linear predictors receive
#'   annual random intercepts. Any combination of `"start"` (start date) and `"duration"`.
#'   Defaults to `c("start", "duration")` (both).
#' @param beta_sd use zero-centred normal priors for regression coefficients other than
#'   intercepts? If <= 0 Stan's default improper flat priors are used
#' @param data input data frame
#' @param init specification of initial values. Can be `"auto"` for an automatic guess
#'   based on the data, or any permitted `rstan::sampling` option
#' @param log_lik logical; retain pointwise log-likelihood in output? Enables model
#'   comparison via the `loo` package. Defaults to TRUE
#' @param standata_only logical; if TRUE, return the Stan data list without fitting
#' @param ... arguments passed to `rstan::sampling` (e.g. iter, chains)
#' @return When `standata_only = FALSE`, an object of class `moultmcmc`; otherwise the
#'   Stan data list
#' @importFrom stats model.matrix sd
#'
uz2_recap_annual_raneff <- function(moult_index_column, date_column, id_column,
                                    start_formula = ~1, duration_formula = ~1, sigma_formula = ~1,
                                    year_factor_column,
                                    raneff_components = c("start", "duration"),
                                    flat_prior = FALSE, lumped = FALSE, use_phi_approx = FALSE,
                                    beta_sd = 0,
                                    data, init = "auto", log_lik = TRUE, standata_only = FALSE, ...) {
  stopifnot(all(data[[moult_index_column]] >= 0 & data[[moult_index_column]] <= 1))
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.factor(data[[id_column]]))
  stopifnot(is.factor(data[[year_factor_column]]))
  stopifnot(is.data.frame(data))
  data <- droplevels(data)

  # order data by moult category: old (0), active moult (0,1), new (1)
  data <- data[order(data[[moult_index_column]]), ]
  # model matrices (over all N_old + N_moult + N_new rows)
  X_mu    <- model.matrix(start_formula,    data)
  X_tau   <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula,    data)

  N_old   <- sum(data[[moult_index_column]] == 0)
  N_moult <- sum(data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1)
  N_new   <- sum(data[[moult_index_column]] == 1)

  # first occurrence of each individual in the model frame
  id_first <- match(unique(data[[id_column]]), data[[id_column]])

  # replicated/non-replicated individual structure
  tab <- table(data[[id_column]])
  replicated     <- which(data[[id_column]] %in% names(tab[tab > 1]))
  is_replicated  <- as.integer(ifelse(tab > 1, 1, 0))

  raneff_mu  <- as.integer("start"    %in% raneff_components)
  raneff_tau <- as.integer("duration" %in% raneff_components)

  standata <- list(
    flat_prior       = as.integer(flat_prior),
    beta_sd          = beta_sd,
    lumped           = as.integer(lumped),
    llik             = as.integer(log_lik),
    use_phi_approx   = as.integer(use_phi_approx),
    raneff_mu        = raneff_mu,
    raneff_tau       = raneff_tau,
    N_ind            = length(unique(data[[id_column]])),
    N_ind_rep        = length(unique(as.numeric(data[[id_column]])[replicated])),
    N_old            = N_old,
    N_moult          = N_moult,
    N_new            = N_new,
    Nobs_replicated  = length(replicated),
    old_dates        = as.array(data[[date_column]][data[[moult_index_column]] == 0]),
    moult_dates      = as.array(data[[date_column]][data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1]),
    moult_indices    = as.array(data[[moult_index_column]][data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1]),
    new_dates        = as.array(data[[date_column]][data[[moult_index_column]] == 1]),
    individual       = as.numeric(data[[id_column]]),
    individual_first_index = as.array(id_first),
    is_replicated    = as.array(is_replicated),
    replicated_individuals = unique(as.numeric(data[[id_column]])[replicated]),
    year_factor      = as.integer(data[[year_factor_column]]),
    N_years          = length(unique(data[[year_factor_column]])),
    X_mu             = X_mu,
    N_pred_mu        = ncol(X_mu),
    X_tau            = X_tau,
    N_pred_tau       = ncol(X_tau),
    X_sigma          = X_sigma,
    N_pred_sigma     = ncol(X_sigma)
  )

  if (standata_only) return(standata)

  outpars <- c('beta_mu_out', 'beta_tau_out', 'beta_tau', 'beta_sigma', 'sigma_intercept',
               'sigma_mu_ind', 'beta_star', 'finite_sd', 'mu_ind_star')
  if (raneff_mu)  outpars <- c(outpars, 'u_year_mean', 'u_year_mean_star', 'sd_year_mean')
  if (raneff_tau) outpars <- c(outpars, 'u_year_duration', 'u_year_duration_star', 'sd_year_duration')
  if (log_lik)    outpars <- c(outpars, 'log_lik')

  if (init == "auto") {
    mu_start    <- mean(c(min(standata$moult_dates), max(standata$old_dates)))
    tau_start   <- max(10, max(standata$moult_dates) - max(standata$old_dates))
    sigma_start <- min(10, sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      list(
        beta_mu    = as.array(c(mu_start,           rep(0, standata$N_pred_mu    - 1))),
        beta_tau   = as.array(c(tau_start,           rep(0, standata$N_pred_tau   - 1))),
        beta_sigma = as.array(c(log(sigma_start),    rep(0, standata$N_pred_sigma - 1))),
        mu_ind     = as.array(rep(0, standata$N_ind))
      )
    }
    out <- rstan::sampling(stanmodels$uz2_recap_annual_raneff,
                           data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz2_recap_annual_raneff,
                           data = standata, init = init, pars = outpars, ...)
  }

  # rename parameters for interpretability
  names(out)[grep('beta_mu_out',    names(out))] <- paste('mean',    colnames(X_mu),    sep = '_')
  names(out)[grep('beta_tau',       names(out))] <- paste('duration',colnames(X_tau),   sep = '_')
  names(out)[grep('beta_sigma',     names(out))] <- paste('log_sd',  colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept',names(out))] <- 'sd_(Intercept)'

  out_struc <- list()
  out_struc$stanfit                  <- out
  out_struc$terms$date_column        <- date_column
  out_struc$terms$moult_index_column <- moult_index_column
  out_struc$terms$moult_cat_column   <- NA
  out_struc$terms$id_column          <- id_column
  out_struc$terms$year_column        <- year_factor_column
  out_struc$terms$start_formula      <- start_formula
  out_struc$terms$duration_formula   <- duration_formula
  out_struc$terms$sigma_formula      <- sigma_formula
  out_struc$data <- data
  out_struc$type <- "2R"
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}


#' Bayesian inference for the Type 2 moult model with annual random intercepts (no individual recaptures)
#'
#' Fits a Type 2 (Underhill-Zucchini) moult model with annual random intercepts on start
#' date and duration. No individual-level random effects; use when within-season recaptures
#' are absent. Type 2 data include birds in old plumage, active moult, and new plumage.
#'
#' @export
#' @param date_column the name of the column in `data` containing sampling dates, encoded
#'   as days since an arbitrary reference date (numeric vector)
#' @param moult_index_column the name of the column in `data` containing moult indices,
#'   i.e. a numeric vector of (linearised) moult scores: 0 = old plumage, 1 = new plumage
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param year_factor_column the name of a factor column in `data` containing year identifiers
#' @param flat_prior logical; use flat (uniform) priors on start date and duration
#'   intercepts? Defaults to FALSE (weakly informative normal priors)
#' @param lumped logical; use the lumped (Type 2L) likelihood? Defaults to FALSE
#' @param raneff_components character vector specifying which linear predictors receive
#'   annual random intercepts. Any combination of `"start"` (start date) and `"duration"`.
#'   Defaults to `c("start", "duration")` (both).
#' @param beta_sd use zero-centred normal priors for regression coefficients other than
#'   intercepts? If <= 0 Stan's default improper flat priors are used
#' @param data input data frame
#' @param init specification of initial values. Can be `"auto"` or any permitted
#'   `rstan::sampling` option
#' @param log_lik logical; retain pointwise log-likelihood in output? Defaults to TRUE
#' @param standata_only logical; if TRUE, return the Stan data list without fitting
#' @param ... arguments passed to `rstan::sampling` (e.g. iter, chains)
#' @return When `standata_only = FALSE`, an object of class `moultmcmc`; otherwise the
#'   Stan data list
#' @importFrom stats model.matrix sd
#'
uz2_linpred_annual_raneff <- function(moult_index_column, date_column,
                                      start_formula = ~1, duration_formula = ~1, sigma_formula = ~1,
                                      year_factor_column,
                                      raneff_components = c("start", "duration"),
                                      flat_prior = FALSE, lumped = FALSE,
                                      beta_sd = 0,
                                      data, init = "auto", log_lik = TRUE, standata_only = FALSE, ...) {
  stopifnot(all(data[[moult_index_column]] >= 0 & data[[moult_index_column]] <= 1))
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.factor(data[[year_factor_column]]))
  stopifnot(is.data.frame(data))
  data <- droplevels(data)

  # order data by moult category: old (0), active moult (0,1), new (1)
  data <- data[order(data[[moult_index_column]]), ]
  X_mu    <- model.matrix(start_formula,    data)
  X_tau   <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula,    data)

  N_old   <- sum(data[[moult_index_column]] == 0)
  N_moult <- sum(data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1)
  N_new   <- sum(data[[moult_index_column]] == 1)

  raneff_mu  <- as.integer("start"    %in% raneff_components)
  raneff_tau <- as.integer("duration" %in% raneff_components)

  standata <- list(
    flat_prior    = as.integer(flat_prior),
    beta_sd       = beta_sd,
    lumped        = as.integer(lumped),
    llik          = as.integer(log_lik),
    raneff_mu     = raneff_mu,
    raneff_tau    = raneff_tau,
    N_old         = N_old,
    N_moult       = N_moult,
    N_new         = N_new,
    old_dates     = as.array(data[[date_column]][data[[moult_index_column]] == 0]),
    moult_dates   = as.array(data[[date_column]][data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1]),
    moult_indices = as.array(data[[moult_index_column]][data[[moult_index_column]] >  0 & data[[moult_index_column]] < 1]),
    new_dates     = as.array(data[[date_column]][data[[moult_index_column]] == 1]),
    year_factor   = as.integer(data[[year_factor_column]]),
    N_years       = length(unique(data[[year_factor_column]])),
    X_mu          = X_mu,
    N_pred_mu     = ncol(X_mu),
    X_tau         = X_tau,
    N_pred_tau    = ncol(X_tau),
    X_sigma       = X_sigma,
    N_pred_sigma  = ncol(X_sigma)
  )

  if (standata_only) return(standata)

  outpars <- c('beta_mu_out', 'beta_tau_out', 'beta_tau', 'beta_sigma', 'sigma_intercept',
               'beta_star')
  if (raneff_mu)  outpars <- c(outpars, 'u_year_mean', 'u_year_mean_star', 'sd_year_mean')
  if (raneff_tau) outpars <- c(outpars, 'u_year_duration', 'u_year_duration_star', 'sd_year_duration')
  if (log_lik)    outpars <- c(outpars, 'log_lik')

  if (init == "auto") {
    mu_start    <- mean(c(min(standata$moult_dates), max(standata$old_dates)))
    tau_start   <- max(10, max(standata$moult_dates) - max(standata$old_dates))
    sigma_start <- min(10, sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      list(
        beta_mu    = as.array(c(mu_start,        rep(0, standata$N_pred_mu    - 1))),
        beta_tau   = as.array(c(tau_start,        rep(0, standata$N_pred_tau   - 1))),
        beta_sigma = as.array(c(log(sigma_start), rep(0, standata$N_pred_sigma - 1)))
      )
    }
    out <- rstan::sampling(stanmodels$uz2_linpred_annual_raneff,
                           data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz2_linpred_annual_raneff,
                           data = standata, init = init, pars = outpars, ...)
  }

  # rename parameters for interpretability
  names(out)[grep('beta_mu_out',    names(out))] <- paste('mean',    colnames(X_mu),    sep = '_')
  names(out)[grep('beta_tau',       names(out))] <- paste('duration',colnames(X_tau),   sep = '_')
  names(out)[grep('beta_sigma',     names(out))] <- paste('log_sd',  colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept',names(out))] <- 'sd_(Intercept)'

  out_struc <- list()
  out_struc$stanfit                  <- out
  out_struc$terms$date_column        <- date_column
  out_struc$terms$moult_index_column <- moult_index_column
  out_struc$terms$moult_cat_column   <- NA
  out_struc$terms$id_column          <- NA
  out_struc$terms$year_column        <- year_factor_column
  out_struc$terms$start_formula      <- start_formula
  out_struc$terms$duration_formula   <- duration_formula
  out_struc$terms$sigma_formula      <- sigma_formula
  out_struc$data <- data
  out_struc$type <- "2"
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
