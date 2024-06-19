#' Optimal Subsampling Methods for Generalized Linear Models
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
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator as well as to
#' estimate the optimal sampling probability.
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size).
#' @param family defalut = 'binomial'.
#' @param criterion The criterion of optimal subsampling probabilities,
#' currently there are three choices \code{OptA}, \code{OptL}, and \code{LCC}.
#' @param sampling.method The sampling method for drawing the optimal subsample,
#'  currently there are two choices \code{WithReplacement} and \code{Poisson}.
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator, currently there are two choices
#'  \code{Weighted} and \code{LogOddsCorrection}.
#' @param alpha Mixture proportions of optimal subsampling probability and
#' uniform subsampling probability. Default = 0.1.
#' @param b This parameter controls the upper threshold for optimal subsampling
#'  probabilities.
#'
#' @return
#' \describe{
#'   \item{beta.plt}{pilot estimator}
#'   \item{beta.ssp}{optimal subsample estimator}
#'   \item{beta.cmb}{combined estimator of \code{beta.plt} and \code{beta.ssp}}
#'   \item{var.ssp}{covariance matrix of \code{beta.ssp}}
#'   \item{var.cmb}{covariance matrix of \code{beta.cmb}}
#'   \item{index.plt}{index of pilot subsample}
#'   \item{index.ssp}{index of optimal subsample}
#' }
#'
#' @examples
#' # logistic regression
#' set.seed(1)
#' N <- 2e4
#' beta0 <- rep(-0.5, 7)
#' d <- length(beta0) - 1
#' X <- matrix(0, N, d)
#' generate_rexp <- function(x) x <- rexp(N, rate = 2)
#' X <- apply(X, 2, generate_rexp)
#' Y <- as.integer(rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))))
#' print(paste('N: ', N))
#' print(paste('sum(Y): ', sum(Y)))
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' n.plt <- 500
#' n.ssp <- 1000
#' subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'binomial', criterion = "OptL", sampling.method = 'Poisson',
#' likelihood = "LogOddsCorrection")
#' subsampling.summary(subsampling.results)
#' subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'binomial', criterion = "OptL",
#' sampling.method = 'WithReplacement', likelihood = "Weighted")
#' subsampling.summary(subsampling.results)
#' Uni.subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'binomial', criterion = 'Uniform')
#' subsampling.summary(Uni.subsampling.results)
#' ############################################################################
#' # poisson regression
#' set.seed(1)
#' N <-  1e4
#' beta0 <- rep(0.5, 7)
#' d <- length(beta0) - 1
#' X <- matrix(runif(N * d), N, d)
#' epsilon <- runif(N)
#' lambda <- exp(beta0[1] + X %*% beta0[-1])
#' Y <- rpois(N, lambda)
#' hist(Y)
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' n.plt <- 200
#' n.ssp <- 600
#' subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'poisson', criterion = "OptL", sampling.method = 'Poisson',
#' likelihood = "Weighted")
#' subsampling.summary(subsampling.results)
#' subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'poisson', criterion = "OptL", sampling.method = 'WithReplacement',
#' likelihood = "Weighted")
#' subsampling.summary(subsampling.results)
#' Uni.subsampling.results <- glm.subsampling(formula, data, n.plt, n.ssp,
#' family = 'poisson', criterion = 'Uniform')
#' subsampling.summary(Uni.subsampling.results)
#' @export

glm.subsampling <- function(formula,
                            data,
                            n.plt,
                            n.ssp,
                            family = c('binomial', 'poisson', 'gamma'),
                            criterion = c('OptL', 'OptA', 'LCC', 'Uniform'),
                            sampling.method = c('Poisson', 'WithReplacement'),
                            likelihood = c('LogOddsCorrection',
                                                'Weighted'),
                            alpha = 0.1,
                            b = 2) {

  model.call <- match.call()
  if(is.function(family)) family <- family$family
  if(missing(family)) {
    stop("Specify a valid 'family' from c('binomial','poisson',gamma')")
  }
  family <- switch(family,
                  "binomial" = binomial.expand(),
                  "poisson" = poisson.expand(),
                  "gamma" = gamma.expand())
  mf <- model.frame(formula, data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (criterion %in% c('OptL', 'OptA', 'LCC')) {

    ### pilot step ###
    plt.estimate.results <- pilot.estimate(X = X,
                                           Y = Y,
                                           n.plt = n.plt,
                                           family = family
                                          )
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    Lambda.plt <- plt.estimate.results$Lambda.plt
    d.psi <- plt.estimate.results$d.psi
    index.plt <- plt.estimate.results$index.plt

    ### subsampling step ###
    ssp.results <- subsampling(X = X,
                               Y = Y,
                               n.ssp = n.ssp,
                               alpha = alpha,
                               b = b,
                               criterion = criterion,
                               likelihood = likelihood,
                               sampling.method = sampling.method,
                               p.plt = p.plt,
                               ddL.plt.correction = ddL.plt.correction,
                               d.psi = d.psi,
                               index.plt = index.plt
                               )
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    ### subsample estimating step ###
    ssp.estimate.results <- subsample.estimate(x.ssp = X[index.ssp, ],
                                               y.ssp = Y[index.ssp],
                                               n.ssp = n.ssp,
                                               N = N,
                                               w.ssp = w.ssp,
                                               offset = offset,
                                               beta.plt = beta.plt,
                                             sampling.method = sampling.method,
                                             likelihood = likelihood,
                                               family = family)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    var.ssp <- ssp.estimate.results$var.ssp

    ### combining step ###
    combining.results <- combining(ddL.plt = ddL.plt,
                                   ddL.ssp = ddL.ssp,
                                   dL.sq.plt = dL.sq.plt,
                                   dL.sq.ssp = dL.sq.ssp,
                                   Lambda.plt = Lambda.plt,
                                   Lambda.ssp = Lambda.ssp,
                                   n.plt = n.plt,
                                   n.ssp = n.ssp,
                                   beta.plt = beta.plt,
                                   beta.ssp = beta.ssp)
    beta.cmb <- combining.results$beta.cmb
    var.cmb <- combining.results$var.cmb

    return(list(model.call = model.call,
                beta.plt = beta.plt,
                beta.ssp = beta.ssp,
                beta = beta.cmb,
                var.ssp = var.ssp,
                var = var.cmb,
                index.plt = index.plt,
                index = index.ssp,
                N = N,
                subsample.size.expect = n.ssp
                )
           )
  } else if (criterion == "Uniform"){
    n.uni <- n.plt + n.ssp
    index.uni <- random.index(N, n.uni)
    x.uni <- X[index.uni, ]
    y.uni = Y[index.uni]
    results.uni <- glm.coef.estimate(X = x.uni, Y = y.uni, family = family)
    beta.uni <- results.uni$beta
    linear.predictor.uni <- as.vector(x.uni %*% beta.uni)
    ddL.uni <- ddL(linear.predictor.uni,
                   x.uni,
                   weights = 1 / n.uni,
                   family = family)
    dL.sq.plt <- dL.sq(linear.predictor.uni,
                       x.uni,
                       y.uni,
                       weights = 1 / n.uni ^ 2,
                       family = family)
    var.uni <- solve(ddL.uni) %*% 
                    (dL.sq.plt * (1 + n.uni / N)) %*% solve(ddL.uni)

    return(list(model.call = model.call,
                index = index.uni,
                beta = beta.uni,
                var = var.uni,
                N = N,
                subsample.size.expect = n.uni
                )
           )
  }
}
