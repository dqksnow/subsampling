#' Optimal Subsampling for Logistic Regression Model with Rare Events Data
#' @description
#' Draw subsample from full dataset and fit logistic regression model on subsample.
#'
#' @param formula An object of class "formula" which describes the model to be
#'  fitted.
#' @param data A data frame containing the variables in the model.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator as well as to
#' estimate the optimal subsampling probability.
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size) drawn from those samples with \code{Y=0}. All rare events (\code{Y=1}) are included in the optimal subsample automatically.
#' @param criterion The criterion of optimal subsampling probabilities.
#' Choices include \code{optA}, \code{optL}(default), \code{LCC} and \code{uniform}. 
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator. Choices include 
#'  \code{weighted} and \code{logOddsCorrection}(default). 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control A list of parameters for controlling the sampling process. Default is \code{list(alpha=0, b=2)}.
#' @param ... A list of parameters which will be passed to \code{svyglm()}. 
#'
#' @return
#' ssp.glm returns an object of class "ssp.glm" containing the following components (some are optional):
#' \describe{
#'   \item{model.call}{model call}
#'   \item{coef.plt}{pilot estimator}
#'   \item{coef.ssp}{optimal subsample estimator.}
#'   \item{coef}{weighted combination of \code{coef.plt} and \code{coef.ssp}.}
#'   \item{cov.ssp}{covariance matrix of \code{coef.ssp}}
#'   \item{cov}{covariance matrix of \code{beta.cmb}}
#'   \item{index.plt}{index of pilot subsample in the full sample}
#'   \item{index.ssp}{index of optimal subsample in the full sample}
#'   \item{N}{number of observations in the full sample}
#'   \item{subsample.size.expect}{expected subsample size}
#'   \item{terms}{model terms}
#' }
#'
#' @details
#' Rare event stands for the number of Y=1 is rare compare to the number of Y=0 in the full sample. When \code{criterion = uniform}, it draws (n.plt+n.ssp) subsmples from the full sample with equal sampling probability. When \code{criterion = optA, optL or LCC}, observations with Y=1 are preserved and it draw n.ssp subsmples from observations with Y=0.
#' 
#' Most of the arguments and returned variables have the same meaning with \link{ssp.glm}. Also refer to [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html)
#' 
#' @references
#' Wang, H., Zhang, A., & Wang, C. (2021). Nonuniform negative sampling and log odds correction with rare events data. \emph{Advances in Neural Information Processing Systems}, \strong{34}, 19847-19859.
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
#' @export
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
  if (attr(mt, "intercept") == 1) {
    colnames(X)[1] <- "Intercept"
  }
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
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset

    ## subsample estimation step
    ssp.estimate.results <- rare.subsample.estimate(inputs, 
                                                    w.ssp = w.ssp,
                                                    offset = offset,
                                                    beta.plt = beta.plt,
                                                    index.ssp = index.ssp,
                                                    ...
                                                    )
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    cov.ssp <- ssp.estimate.results$cov.ssp

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
    cov.cmb <- combining.results$cov.cmb
    names(beta.plt) <- names(beta.ssp) <- names(beta.cmb) <- colnames(X)
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef = beta.cmb,
                    cov.ssp = cov.ssp,
                    cov = cov.cmb,
                    index.plt = index.plt,
                    index = index.ssp,
                    N = N,
                    subsample.size.expect = N1 + n.ssp,
                    terms = mt
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
                                      ...
                                      )
    beta.uni <- results.uni$beta
    cov.uni <- results.uni$cov
    beta.uni[1] <- beta.uni[1] + log(n.ssp / N0) # correct intercept
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    index = index.uni,
                    coef = beta.uni,
                    cov = cov.uni,
                    N = N,
                    subsample.size.expect = n.uni,
                    terms = mt
                    )
    class(results) <- c("ssp.relogit", "list")
    return(results)
  }
}


