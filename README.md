
<!-- README.md is generated from README.Rmd. Please edit that file -->

# subsampling

<!-- badges: start -->

[![R-CMD-check](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The R package `subsampling` provides optimal subsampling methods for
common statistical models such as glm, softmax(multinomial) regression,
rare event logistic regression and quantile regression model.

## Installation

You can install the development version of subsampling from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dqksnow/Subsampling")
```

## Getting Started

The [Online document](https://dqksnow.github.io/Subsampling/) provides a
guidance for quick start.

- [Generalized Linear
  Model](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html).

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(subsampling)
set.seed(1)
N <- 1e4
beta0 <- rep(-0.5, 7)
d <- length(beta0) - 1
corr <- 0.5
sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
X <- MASS::mvrnorm(N, rep(0, d), sigmax)
colnames(X) <- paste("V", 1:ncol(X), sep = "")
P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
Y <- rbinom(N, 1, P)
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ .
n.plt <- 200
n.ssp <- 600
ssp.results <- ssp.glm(formula = formula,
                       data = data,
                       n.plt = n.plt,
                       n.ssp = n.ssp,
                       family = "quasibinomial",
                       criterion = "optL",
                       sampling.method = "poisson",
                       likelihood = "weighted"
                       )
summary(ssp.results)
#> Model Summary
#> 
#> Call:
#> 
#> ssp.glm(formula = formula, data = data, n.plt = n.plt, n.ssp = n.ssp, 
#>     family = "quasibinomial", criterion = "optL", sampling.method = "poisson", 
#>     likelihood = "weighted")
#> 
#> Subsample Size:
#>                                
#> 1       Total Sample Size 10000
#> 2 Expected Subsample Size   600
#> 3   Actual Subsample Size   651
#> 4   Unique Subsample Size   651
#> 5  Expected Subample Rate    6%
#> 6    Actual Subample Rate 6.51%
#> 7    Unique Subample Rate 6.51%
#> 
#> Coefficients:
#> 
#>           Estimate Std. Error z value Pr(>|z|)
#> Intercept  -0.4092     0.0795 -5.1486  <0.0001
#> V1         -0.5861     0.0949 -6.1791  <0.0001
#> V2         -0.4514     0.1066 -4.2343  <0.0001
#> V3         -0.5557     0.1005 -5.5283  <0.0001
#> V4         -0.3915     0.1006 -3.8898   0.0001
#> V5         -0.3732     0.1046 -3.5697   0.0004
#> V6         -0.6454     0.0969 -6.6589  <0.0001
```
