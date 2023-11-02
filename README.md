
# bayesurv

<!-- badges: start -->
<!-- badges: end -->

The goal of bayesurv is to implement the Gibbs sampler method to group-tested data with misspecification based on the semiparametric probit model.

## Installation

You can install the development version of bayesurv like so:

``` r
devtools::install_github("jihyunk1114/bayesurv")
```

## Arguments
Y: individual status by group test result (length: sample size)

Yi: group status (length: number of groups)

X: covariates in matrix form

c: test time

grid: can define grid and find corresponding baseline survival function

theta0, Sigma0, m0, v0, a0, b0: For priors

alpha: Sensitivity

beta: Specificity

order: order for knots (usually 2 or 3)

knots: If NULL default, can set up own knots

maxiter: maximum interation number

burn.in: burn-in number

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bayesurv)
Gibbs(Y=Y,X=X,c=c,groupID=groupID,alpha=0.9,beta=0.8)
```

