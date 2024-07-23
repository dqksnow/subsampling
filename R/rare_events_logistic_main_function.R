#' Optimal Subsampling for Logistic Regression Model with Rare Events Data
#' @details
#' Additional details... briefly introduce the idea of this method.
#'
#' @param formula An object of class "formula" which describes the model to be
#' fitted.
#' @param data A data frame containing the variables in the model. Usually it
#' contains a response vector and a design matrix. The binary response vector
#' that takes the value of 0 or 1, where 1 means the event occurred. The
#' design matrix contains predictor variables. A column representing the
#' intercept term with all 1's will be automatically added.
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size) drawn from those samples with \code{Y=0}.
#' @param criterion The criterion of optimal subsampling probabilities,
#' currently there are three choices \code{optA}, \code{optL}, and \code{LCC}.
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator, currently there are two choices
#'  \code{weighted} and \code{logOddsCorrection}.
#' @param alpha Mixture proportions of optimal subsampling probability and
#' uniform subsampling probability. Default = 0.1.
#' @param b This parameter controls the upper threshold for optimal subsampling
#' probabilities.
#'
#' @return
#' \describe{
#'   \item{beta.plt}{pilot estimator}
#'   \item{beta.ssp}{optimal subsample estimator}
#'   \item{beta.cmb}{combined estimator of \code{beta.plt} and \code{beta.ssp}}
#'   \item{var.ssp}{covariance matrix of \code{beta.ssp}}
#'   \item{var.cmb}{covariance matrix of \code{beta.cmb}}
#'   \item{index.plt}{index of pilot subsample}
#'   \item{index.ssp}{index of optimal subsample. The expectation of
#'   \code{length(index.ssp)} is n.ssp plus the number of rare event data.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' N <- 2 * 1e4
#' beta0 <- c(-1.8, -rep(1, 6))
#' d <- length(beta0) - 1
#' X <- matrix(0, N, d)
#' generate_rexp <- function(x) x <- rexp(N)
#' X <- apply(X, 2, generate_rexp)
#' Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
#' print(paste('N: ', N))
#' print(paste('sum(Y): ', sum(Y)))
#' n.plt <- 200
#' n.ssp <- 1000
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' subsampling.results <- relogit.ssp(formula,
#'                                      data,
#'                                      n.plt,
#'                                      n.ssp,
#'                                      criterion = 'optA',
#'                                      likelihood = 'logOddsCorrection')
#' summary(subsampling.results)


relogit.ssp <-  function(formula,
                         data,
                         n.plt,
                         n.ssp,
                         criterion = c('optL', 'optA', 'LCC', 
                                       'uniform'),
                         likelihood = c('logOddsCorrection',
                                        'weighted'),
                         alpha = 0.1,
                         b = 2) {
  model.call <- match.call()
  mf <- model.frame(formula, data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1

  if (criterion %in% c('optL', 'optA', 'LCC')){

    ## pilot step
    plt.estimate.results <- rare.pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt

    ## subsampling step
    ssp.results <- rare.subsampling(X = X,
                                    Y = Y,
                                    n.ssp = n.ssp,
                                    alpha = alpha,
                                    b = b,
                                    criterion = criterion,
                                    likelihood = likelihood,
                                    p.plt = p.plt,
                                    ddL.plt.correction = ddL.plt.correction,
                                    P.plt = P.plt,
                                    index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    ## subsample estimation step
    ssp.estimate.results <- rare.subsample.estimate(X[index.ssp, ],
                                             Y[index.ssp],
                                             n.ssp = n.ssp,
                                             N = N,
                                             w.ssp = w.ssp,
                                             offset = offset,
                                             beta.plt = beta.plt,
                                             likelihood = likelihood)
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    var.ssp <- ssp.estimate.results$var.ssp

    ## combine step
    combining.results <- rare.combining(ddL.plt = ddL.plt,
                                        ddL.ssp = ddL.ssp,
                                        dL.sq.plt = dL.sq.plt,
                                        dL.sq.ssp = dL.sq.ssp,
                                        n.plt = n.plt,
                                        n.ssp = n.ssp,
                                        beta.plt = beta.plt,
                                        beta.ssp = beta.ssp)
    beta.cmb <- combining.results$beta.cmb
    var.cmb <- combining.results$var.cmb
    results <- list(model.call = model.call,
                    beta.plt = beta.plt,
                    beta.ssp = beta.ssp,
                    beta = beta.cmb,
                    var.ssp = var.ssp,
                    var = var.cmb,
                    index.plt = index.plt,
                    index = index.ssp,
                    N = N,
                    subsample.size.expect = N1 + n.ssp
                    )
    
    class(results) <- c("relogit.ssp", "list")
    return(results)
  } else if (criterion == "uniform"){
    ## Poisson sampling
    n.uni <- N1 + n.ssp
    pi.uni <- rep(1, N)
    pi.uni[Y == 0] <- n.ssp / N0
    index.uni <- poisson.index(N, pi.uni)
    results.uni <- rare.coef.estimate(X = X[index.uni, ], Y = Y[index.uni])
    beta.uni <- results.uni$beta
    var.uni <- results.uni$cov
    beta.uni[1] <- beta.uni[1] + log(n.ssp / N0) # correct intercept
    results <- list(model.call = model.call,
                    index = index.uni,
                    beta = beta.uni,
                    var = var.uni,
                    N = N,
                    subsample.size.expect = n.uni
                    )
    class(results) <- c("relogit.ssp", "list")
    return(results)
  }
}


