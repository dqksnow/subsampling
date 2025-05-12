
# subsampling

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/subsampling)](https://CRAN.R-project.org/package=subsampling)
[![Total_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/subsampling)](https://CRAN.R-project.org/package=subsampling)
[![Monthly_Downloads](https://cranlogs.r-pkg.org/badges/subsampling)](https://CRAN.R-project.org/package=subsampling)
[![R-CMD-check](https://github.com/dqksnow/subsampling/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dqksnow/subsampling/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A major challenge in big data statistical analysis is the demand for
computing resources. For example, when fitting a logistic regression
model to binary response variable with $N \times d$ dimensional
covariates, the computational complexity of estimating the coefficients
using the IRLS algorithm is $O(\zeta N d^2)$, where $\zeta$ is the
number of iteriation. When $N$ is large, the cost can be prohibitive,
especially if high performance computing resources are unavailable.
Subsampling has become a widely used technique to balance the trade-off
between computational efficiency and statistical efficiency.

The R package `subsampling` provides optimal subsampling methods for
various statistical models such as generalized linear models (GLM),
softmax (multinomial) regression, rare event logistic regression and
quantile regression model. Specialized subsampling techniques are
provided to address specific challenges across different models and
datasets. With specified model assumptions and subsampling techniques,
it draws subsample from the full data, fits model on the subsample and
perform statistical inferences.

## Installation

You can install the package by

``` r
# Install from CRAN
install.packages("subsampling")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("dqksnow/subsampling")
```

## Getting Started

The [Online document](https://dqksnow.github.io/subsampling/) provides a
guidance for quick start.

- [Generalized linear
  model](https://dqksnow.github.io/subsampling/articles/ssp-logit.html).
- [Rare event logistic
  regression](https://dqksnow.github.io/subsampling/articles/ssp-relogit.html).
- [Softmax (multinomial)
  regression](https://dqksnow.github.io/subsampling/articles/ssp-softmax.html).
- [Quantile
  regression](https://dqksnow.github.io/subsampling/articles/ssp-quantreg.html).

## Example

This is an example of subsampling method on logistic regression:

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
#> 3   Actual Subsample Size   635
#> 4   Unique Subsample Size   635
#> 5  Expected Subample Rate    6%
#> 6    Actual Subample Rate 6.35%
#> 7    Unique Subample Rate 6.35%
#> 
#> Coefficients:
#> 
#>           Estimate Std. Error z value Pr(>|z|)
#> Intercept  -0.4149     0.0803 -5.1694  <0.0001
#> V1         -0.5874     0.0958 -6.1286  <0.0001
#> V2         -0.4723     0.1086 -4.3499  <0.0001
#> V3         -0.5492     0.1014 -5.4164  <0.0001
#> V4         -0.4044     0.1012 -3.9950  <0.0001
#> V5         -0.3725     0.1045 -3.5649   0.0004
#> V6         -0.6703     0.0973 -6.8859  <0.0001
```

## Acknowledgments

The development of this package was supported by the National Eye
Institute of the National Institutes of Health under Award Number
R21EY035710.
