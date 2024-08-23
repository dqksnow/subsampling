
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(subsampling)
set.seed(1)
N <- 5e3
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
n.plt <- 100
n.ssp <- 500
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
#> 
#> Call:
#> 
#> ssp.glm(formula = formula, data = data, n.plt = n.plt, n.ssp = n.ssp, 
#>     family = "quasibinomial", criterion = "optL", sampling.method = "poisson", 
#>     likelihood = "weighted")
#> 
#> Subsample Size:
#>                                
#> 1       Total Sample Size  5000
#> 2 Expected Subsample Size   500
#> 3   Actual Subsample Size   499
#> 4   Unique Subsample Size   499
#> 5  Expected Subample Rate   10%
#> 6    Actual Subample Rate 9.98%
#> 7    Unique Subample Rate 9.98%
#> 
#> Coefficients:
#> 
#>           Estimate Std. Error z value Pr(>|z|)
#> Intercept  -0.6193     0.3142 -1.9712   0.0487
#> V1         -0.4909     0.3331 -1.4739   0.1405
#> V2         -0.6077     0.4277 -1.4209   0.1554
#> V3         -0.5586     0.3963 -1.4095   0.1587
#> V4         -0.4988     0.3496 -1.4266   0.1537
#> V5         -0.4414     0.3743 -1.1792   0.2383
#> V6         -0.4866     0.3683 -1.3212   0.1864
```
