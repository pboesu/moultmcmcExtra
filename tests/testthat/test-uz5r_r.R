library(moultmcmc)
library(moultmcmcExtra)
library(dplyr)
n_distinct(recaptures$id)
set.seed(1234)
recaps_yr <- left_join(recaptures, data.frame(id = unique(recaptures$id), year_fac = factor(sample(8, size = n_distinct(recaptures$id), replace = TRUE))))

recaps_yr %>% group_by(year_fac, id) %>% summarize(n = n()) %>% ungroup() -> summs
table(summs$year_fac, summs$n)

test_that("uz5r_ranef works", {
  uz5r_r = moultmcmcExtra::uz5_linpred_recap_annual_raneff("pfmg_sampled",
                                           date_column = "date_sampled",
                                           id_column = "id",
                                           year_factor_column = 'year_fac',
                                           data = subset(recaps_yr, pfmg_sampled != 1),
                                           log_lik = FALSE,
                                           chains = 2,cores=2,
                                           iter = 1000)
  expect_s3_class(uz5r_r, "moultmcmc")

  compare_plot_annual_raneff(
    uz5r_r
  )
  #predict(uz5r_r)
  #predict_ranef(uz5r_r)
})

test_that("moultmcmc_ranef works", {
  muz5r_r = moultmcmcExtra::moultmcmc_ranef("pfmg_sampled",
                                            date_column = "date_sampled",
                                            id_column = "id",
                                            year_factor_column = 'year_fac',
                                            type = 5,
                                            data = subset(recaps_yr, pfmg_sampled != 1),
                                            log_lik = FALSE,
                                            chains = 2,cores=2,
                                            iter = 1000)
  expect_s3_class(muz5r_r, "moultmcmc")

  # compare_plot_annual_raneff(
  #   uz5r_r, muz5r_r
  # )
  #predict(uz5r_r)
  #predict_ranef(uz5r_r)
})
