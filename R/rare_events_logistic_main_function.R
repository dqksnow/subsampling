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
#' @param subset An optional vector specifying a subset of rows to be used
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size) drawn from those samples with \code{Y=0}.
#' @param criterion The criterion of optimal subsampling probabilities,
#' currently there are three choices \code{optA}, \code{optL}, and \code{LCC}.
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator, currently there are two choices
#'  \code{weighted} and \code{logOddsCorrection}.
#' @param contrasts an optional list. It specifies how categorical variables are represented in the design matrix. For example, contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')
#' @param control a list of parameters for controlling the fitting process. 
#' @param ... a list of parameters for controlling the fitting process. 
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
#' beta0 <- c(-5, -rep(0.7, 6))
#' d <- length(beta0) - 1
#' X <- matrix(0, N, d)
#' corr <- 0.5
#' sigmax <- corr ^ abs(outer(1:d, 1:d, "-"))
#' sigmax <- sigmax / 4
#' X <- MASS::mvrnorm(n = N, mu = rep(0, d), Sigma = sigmax)
#' Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
#' print(paste('N: ', N))
#' print(paste('sum(Y): ', sum(Y)))
#' n.plt <- 200
#' n.ssp <- 1000
#' data <- as.data.frame(cbind(Y, X))
#' colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""))
#' formula <- Y ~ .
#' subsampling.results <- ssp.relogit(formula = formula,
#'                                      data = data,
#'                                      n.plt = n.plt,
#'                                      n.ssp = n.ssp,
#'                                      criterion = 'optA',
#'                                      likelihood = 'logOddsCorrection')
#' summary(subsampling.results)

ssp.relogit <-  function(formula,
                         data,
                         subset = NULL,
                         n.plt,
                         n.ssp,
                         criterion = 'optL',
                         likelihood = 'logOddsCorrection',
                         control = list(...),
                         contrasts = NULL,
                         ...
                         ) {
  model.call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"),
             names(mf),
             0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  X <- model.matrix(mt, mf, contrasts)
  colnames(X)[1] <- "Intercept"
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  criterion <- match.arg(criterion, c('optL', 'optA', 'LCC', 'uniform'))
  likelihood <- match.arg(likelihood, c('logOddsCorrection', 'weighted'))
  control <- do.call("relogit.control", control)
  
  
  inputs <- list(X = X, Y = Y, N = N, N1 = N1, N0 = N0, d = d,
                 n.plt = n.plt, n.ssp = n.ssp,
                 criterion = criterion, 
                 likelihood = likelihood, 
                 control = control
                 )
  
  
  if (criterion %in% c('optL', 'optA', 'LCC')){

    ## pilot step
    plt.estimate.results <- rare.pilot.estimate(inputs, ...)
    # plt.estimate.results <- rare.pilot.estimate(X = X, Y = Y, n.plt = n.plt)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    P.plt <- plt.estimate.results$P.plt
    index.plt <- plt.estimate.results$index.plt

    ## subsampling step
    ssp.results <- rare.subsampling(inputs,
                                    p.plt = p.plt,
                                    ddL.plt.correction = ddL.plt.correction,
                                    P.plt = P.plt,
                                    index.plt = index.plt)
    # ssp.results <- rare.subsampling(X = X,
    #                                 Y = Y,
    #                                 n.ssp = n.ssp,
    #                                 alpha = alpha,
    #                                 b = b,
    #                                 criterion = criterion,
    #                                 likelihood = likelihood,
    #                                 p.plt = p.plt,
    #                                 ddL.plt.correction = ddL.plt.correction,
    #                                 P.plt = P.plt,
    #                                 index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    ## subsample estimation step
    ssp.estimate.results <- rare.subsample.estimate(inputs, 
                                                    w.ssp = w.ssp,
                                                    offset = offset,
                                                    beta.plt = beta.plt,
                                                    index.ssp = index.ssp,
                                                    ...)
    # ssp.estimate.results <- rare.subsample.estimate(X[index.ssp, ],
    #                                          Y[index.ssp],
    #                                          n.ssp = n.ssp,
    #                                          N = N,
    #                                          w.ssp = w.ssp,
    #                                          offset = offset,
    #                                          beta.plt = beta.plt,
    #                                          likelihood = likelihood)
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
    names(beta.plt) <- names(beta.ssp) <- names(beta.cmb) <- colnames(X)
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
    
    class(results) <- c("ssp.relogit", "list")
    return(results)
  } else if (criterion == "uniform"){
    ## Poisson sampling
    n.uni <- N1 + n.ssp
    pi.uni <- rep(1, N)
    pi.uni[Y == 0] <- n.ssp / N0
    index.uni <- poisson.index(N, pi.uni)
    results.uni <- rare.coef.estimate(X = X[index.uni, ],
                                      Y = Y[index.uni],
                                      ...)
    beta.uni <- results.uni$beta
    var.uni <- results.uni$cov
    beta.uni[1] <- beta.uni[1] + log(n.ssp / N0) # correct intercept
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    index = index.uni,
                    beta = beta.uni,
                    var = var.uni,
                    N = N,
                    subsample.size.expect = n.uni
                    )
    class(results) <- c("ssp.relogit", "list")
    return(results)
  }
}


