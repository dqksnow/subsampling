#' Optimal Subsampling Methods for Statistical Models
#'
#' Subsampling methods are utilized in statistical modeling for 
#' massive datasets. These methods aim to draw representative subsamples from the 
#' full dataset based on specific sampling probabilities, with the goal of 
#' maintaining inference efficiency. The sampling probabilities are tailored to 
#' particular objectives, such as minimizing the variance of the estimated 
#' coefficients or reducing prediction error. By using subsampling techniques,
#' the package balances the trade-off between computational efficiency and statistical
#' efficiency, making it a practical tool for massive data 
#' analysis.
#'
#' @name subsampling
#' @useDynLib subsampling
#'
#' @importFrom stats as.formula runif model.frame model.matrix model.response pnorm quantile coef quasibinomial
#' @importFrom quantreg rq
#' @importFrom Rcpp evalCpp
#' 
#' @docType package
#'
#'
#' @section Models Supported:
#' \itemize{
#'   \item Generalized Linear Models (GLMs)
#'   \item Softmax (Multinomial) Regression
#'   \item Rare Event Logistic Regression
#'   \item Quantile Regression
#' }
#'
#'
"_PACKAGE"
NULL
