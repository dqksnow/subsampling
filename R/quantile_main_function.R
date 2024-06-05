#' Optimal Subsampling Methods for Quantile Models
#' @details
#' Additional details... briefly introduce the idea.
#'
#' @param formula An object of class "formula" which describes the model to be
#'  fitted.
#' @param data A data frame containing the variables in the model. Usually it
#' contains a response vector and a design matrix.
#'     The binary response vector that takes the value of 0 or 1, where 1 means
#'      the event occurred.
#'     The design matrix contains predictor variables. A column representing
#'     the intercept term with all 1's will be automatically added.
#' @param tau The quantile.
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator.
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size).
#' @param B TBD
#' @param criterion The criterion of optimal subsampling probabilities.
#' @param sampling.method The sampling method for drawing the optimal subsample.
#' @param estimate.method The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator.
#'
#' @return
#' \describe{
#'   \item{model.call}{the model}
#'   \item{beta.plt}{pilot estimator}
#'   \item{beta.ssp}{optimal subsample estimator}
#'   \item{est.cov.ssp}{covariance matrix of \code{beta.ssp}}
#'   \item{index.plt}{index of pilot subsample}
#'   \item{index.ssp}{index of optimal subsample}
#' }
#'
#' @examples
#' #quantile regression
#' set.seed(1)
#' N <- 1e6
#' B <- 10
#' tau <- 0.5
#' beta.true <- rep(1, 7)
#' d <- length(beta.true) - 1
#' corr  <- 0.5
#' sigmax  <- matrix(0, d, d)
#' for (i in 1:d) for (j in 1:d) sigmax[i, j] <- corr^(abs(i-j))
#' X <- MASS::mvrnorm(N, rep(0, d), sigmax)
#' err <- rnorm(N, 0, 1) - qnorm(tau)
#' Y <- beta.true[1] + X %*% beta.true[-1] + 
#' err * rowMeans(abs(X))
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ X
#' n.plt <- 1000
#' n.ssp <- 1000
#' optL.results <- quantile.subsampling(formula,data,tau = tau,n.plt = n.plt,
#' n.ssp = n.ssp,B,criterion = 'OptL',sampling.method = 'WithReplacement',
#' estimate.method = 'Weighted',parallel = F)
#' uni.results <- quantile.subsampling(formula,data,tau = tau,n.plt = n.plt,
#' n.ssp = n.ssp,B,sampling.method = 'WithReplacement',
#' estimate.method = 'Uniform',parallel = F)

quantile.subsampling <- function(formula,
                                 data,
                                 tau,
                                 n.plt,
                                 n.ssp,
                                 B,
                                 criterion = c('OptL'),
                                 sampling.method = c('WithReplacement'),
                                 estimate.method = c('Weighted', 'Uniform'),
                                 parallel = parallel) {
  
  model.call <- match.call()
  mf <- model.frame(formula, data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
  N <- nrow(X)
  d <- ncol(X)
  if (estimate.method %in% c("Weighted")) {
    ### pilot step ###
    plt.results <- quantile.plt.estimation(X, Y, tau, N, n.plt)
    beta.plt <- plt.results$beta.plt
    Ie.full <- plt.results$Ie.full
    index.plt <- plt.results$index.plt
    ### subsampling and boot step ###
    ssp.results <- quantile.ssp.estimation(X = X,
                                           Y = Y,
                                           n.ssp = n.ssp,
                                           B = B,
                                           tau = tau,
                                           Ie.full = Ie.full,
                                           estimate.method = estimate.method,
                                           parallel = parallel)
    Betas.ssp <- ssp.results$Betas.ssp
    beta.ssp <- ssp.results$beta.ssp
    est.cov.ssp <- ssp.results$est.cov.ssp
    index.ssp <- ssp.results$index.ssp
    return(list(model.call = model.call,
                beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                est.cov.ssp = est.cov.ssp,
                index.plt = index.plt,
                index.ssp = index.ssp,
                uniduqe.index.ssp <- unique(as.vector(index.ssp))
    )
    )
  } else if (estimate.method == "Uniform"){
    ### subsampling and boot step ###
    uni.results <- quantile.ssp.estimation(X = X,
                                           Y = Y,
                                           n.ssp = n.ssp,
                                           B = B,
                                           tau = tau,
                                           Ie.full = NA,
                                           estimate.method = estimate.method,
                                           parallel = parallel)
    Betas.uni <- uni.results$Betas.ssp
    beta.uni <- uni.results$beta.ssp
    est.cov.uni <- uni.results$est.cov.ssp
    index.uni <- uni.results$index.ssp
    return(list(model.call = model.call,
                beta.plt = NA,
                beta.uni = beta.uni,
                est.cov.uni = est.cov.uni,
                index.plt = NA,
                index.uni = index.uni,
                uniduqe.index.uni <- unique(as.vector(index.uni))
                )
           )
    
    n.uni <- n.plt + n.ssp
    index.uni <- sample(N, n.uni)
    x.uni <- X[index.uni, ]
    y.uni = Y[index.uni]
    results <- quantreg::rq(y.uni ~ x.uni - 1, tau=tau)
    beta.uni <- results$coefficients
    return(list(model.call = model.call,
                index.uni = index.uni,
                beta.uni = beta.uni
                )
           )
  }
}
