# moultmcmc: Bayesian inference for moult phenology models

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/pboesu/moultmcmcExtra/workflows/R-CMD-check/badge.svg)](https://github.com/pboesu/moultmcmcExtra/actions)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/moultmcmc)]()
[![moultmcmcExtra status badge](https://pboesu.r-universe.dev/badges/moultmcmcExtra)](https://pboesu.r-universe.dev)
<!-- badges: end -->

# Introduction
The [`moultmcmc` package](https://github.com/pboesu/moultmcmc) implements a Bayesian inference framework for models of avian primary moult data.

`moultmcmcExtra` implements further model extensions, which are currently outside the scope of `moultmcmc`.

# Installation
The package `moultmcmcExtra` is built around pre-compiled [Stan](https://mc-stan.org/) models, the easiest and quickest way of installing it is to install the package from [R-universe](https://pboesu.r-universe.dev/) use the following code:

```r
install.packages("moultmcmcExtra", repos = "https://pboesu.r-universe.dev")
```
On Mac and Windows systems this will make use of pre-compiled binaries, which means the models can be run  without having to install a C++ compiler. On Linux this will install the package from a source tarball. Because of the way the Stan models are currently structured, compilation from source can be a lengthy process (several minutes), depending on system setup and compiler toolchain.

To install `moultmcmcExtra` from the github source (not generally recommended for Windows users) use the following code. This requires a working C++ compiler and a working installation of [rstan](https://mc-stan.org/rstan):

```r
#not generally recommended for Windows/MacOS users
install.packages("remotes")
remotes::install_github("pboesu/moultmcmcExtra")
```


