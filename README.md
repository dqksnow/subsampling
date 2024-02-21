
# subsampling

<!-- badges: start -->
[![R-CMD-check](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of subsampling is to ...

## Meeting summary

### Feb 21

- **[1]** When calculating $Var(\hat{\beta}_{plt} - \beta_{true})$, I forget to add a term to correct its difference with $Var(\hat{\beta}_{plt} - \beta_{full})$. As a result, this term is missed in the calculation of $Var(\hat{\beta}_{cmb} - \beta_{true})$. Check this problem in softmax code and previous code.

- **[2]** Implement MSCLE with the assistance of the Julia code. First implement OptL.

- **[3]** LUC


## Installation

You can install the development version of subsampling from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dqksnow/Subsampling")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(subsampling)
## basic example code
```

