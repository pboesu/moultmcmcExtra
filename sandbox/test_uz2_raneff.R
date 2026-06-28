#install from r-universe
install.packages("moultmcmcExtra", repos = "https://pboesu.r-universe.dev")

library(moult)
library(moultmcmc)
library(moultmcmcExtra)
library(ggplot2)
data(weavers)

mscores <- substr(weavers$Moult, 1, 9)

## convert moult scores to proportion of feather mass grown
feather.mass <- c(10.4, 10.8, 11.5, 12.8, 14.4, 15.6, 16.3, 15.7, 15.7)
weavers$pfmg <- ms2pfmg(mscores, feather.mass)

## days since 1 August
weavers$day <- date2days(weavers$RDate, dateformat = "yyyy-mm-dd", startmonth = 8)
ssex <- ifelse(weavers$Sex == 1 | weavers$Sex == 3, "male",
               ifelse(weavers$Sex == 2 | weavers$Sex == 4, "female", NA))
weavers$ssex <- as.factor(ssex)
weavers$year_fac <- as.factor(weavers$Year)
ggplot(weavers, aes(x=day, y=pfmg)) + geom_point()

table(weavers$Year)
table(is.na(weavers$pfmg))

weavers <- dplyr::filter(weavers, !is.na(pfmg))


## moult model with duration and mean start date depending on sex
mmf <- moult(pfmg ~ day | ssex | ssex, data = weavers, type = 3)
summary(mmf)

uz2l = moultmcmc("pfmg",
                 date_column = "day",
                 id_column = NULL,
                 start_formula = ~ssex,
                 duration_formula = ~ssex,
                 type=2,
                 lump_non_moult = TRUE,
                 log_lik = FALSE,
                 data = weavers,
                 chains = 2,cores=2,
                 iter = 1000)
summary_table(uz2l)
moult_plot(uz2l, newdata = data.frame(ssex = c('male','female')))
moult_plot(mmf, newdata = data.frame(ssex = c('male','female')))
uz2l_ranef = moultmcmc_ranef(
  "pfmg",
  date_column = "day",
  ranef_factor_column = 'year_fac',
  start_formula = ~ssex,
  duration_formula = ~ssex,
  type = 2,
  lump_non_moult = TRUE,
  data = weavers,
  log_lik = FALSE,
  chains = 2,
  cores = 2,
  iter = 1000
)

uz2l_ranef_start = moultmcmc_ranef(
  "pfmg",
  date_column = "day",
  ranef_factor_column = 'year_fac',
  start_formula = ~ssex,
  duration_formula = ~ssex,
  type = 2,
  lump_non_moult = TRUE,
  ranef_components = c("start"),
  data = weavers,
  log_lik = FALSE,
  chains = 2,
  cores = 2,
  iter = 1000
)

uz2l_ranef_duration = moultmcmc_ranef(
  "pfmg",
  date_column = "day",
  ranef_factor_column = 'year_fac',
  start_formula = ~ssex,
  duration_formula = ~ssex,
  type = 2,
  lump_non_moult = TRUE,
  ranef_components = c("duration"),
  data = weavers,
  log_lik = FALSE,
  chains = 2,
  cores = 2,
  iter = 1000
)
moult_plot(uz2l_ranef, newdata = data.frame(ssex = c('male','female')))
summary_table(uz2l_ranef)
compare_plot_annual_ranef(uz2l_ranef, uz2l_ranef_start,uz2l_ranef_duration, mmf, uz2l)
