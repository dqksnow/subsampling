
# subsampling

<!-- badges: start -->
[![R-CMD-check](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dqksnow/Subsampling/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of subsampling is to ...

### Mar 20

- **[1]** Fix (n.plt+n.ssp) and then enlarge n.plt so that $Var(\hat{\beta}_{cmb} - \beta_{true})$ should be close to $Var(\hat{\beta}_{plt} - \beta_{true})$. The simulation results show that it meets expectation except when n.plt and n.ssp are both large. Try to use degree of freedom to correct.

- **[2]** Implement MSCLE. First implement OptL.

## Meeting summary

### Mar 6

- **[1]** Since $Var(\hat{\beta}_{plt} - \beta_{true})$ and $Var(\hat{\beta}_{ssp} - \beta_{true})$ work well, we should double check the calculation of $Var(\hat{\beta}_{cmb} - \beta_{true})$. Simulation strategy: fix (n.plt+n.ssp) and then enlarge n.plt so that $Var(\hat{\beta}_{cmb} - \beta_{true})$ should be close to $Var(\hat{\beta}_{plt} - \beta_{true})$.

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

