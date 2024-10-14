#' Optimal Subsampling for Logistic Regression Model with Rare Events Data
#' @description
#' Draw subsample from full dataset and fit logistic regression model on subsample. For a quick start, refer to the [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-relogit.html).
#'
#'
#' @param formula A model formula object of class "formula" that describes the model to be fitted.
#' @param data A data frame containing the variables in the model. Denote \eqn{N} as the number of observations in `data`.
#' @param subset An optional vector specifying a subset of observations from `data` to use for the analysis. This subset will be viewed as the full data.
#' @param n.plt The pilot subsample size (first-step subsample size).
#' This subsample is used to compute the pilot estimator and estimate the optimal subsampling probabilities.
#' @param n.ssp The expected subsample size (the second-step subsample
#' size) drawn from those samples with \code{Y=0}. All rare events (\code{Y=1}) are included in the optimal subsample automatically.
#' @param criterion The choices include `optA`, `optL`(default), `LCC` and `uniform.`
#' 
#' - `optA` Minimizes the trace of the asymptotic covariance matrix of the subsample estimator. 
#' 
#' - `optL` Minimizes the trace of a transformation of the asymptotic covariance matrix. The computational complexity of
#' optA is \eqn{O(N d^2)} while that of optL is \eqn{O(N d)}.
#' 
#' - `LCC` Local Case-Control sampling probability, used as a baseline subsampling strategy.
#' 
#' - `uniform` Assigns equal subsampling probability
#' \eqn{\frac{1}{N}} to each observation, serving as a baseline subsampling strategy.
#' 
#' @param likelihood The likelihood function to use. Options include `weighted` and
#' `logOddsCorrection` (default). A bias-correction likelihood function is required for subsample since unequal subsampling probabilities introduce bias.
#' 
#' - `weighted` Applies a weighted likelihood function where each observation is weighted by the inverse of its subsampling probability.
#' 
#' - `logOddsCorrection` This lieklihood is available only for logistic regression model (i.e., when family is binomial or quasibinomial). It uses a conditional likelihood, where each element of the likelihood represents the probability of \eqn{Y=1}, given that this subsample was drawn.
#' 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control The argument `control` contains two tuning parameters `alpha` and `b`. 
#' 
#' - `alpha` \eqn{\in [0,1]} is the mixture weight of the user-assigned subsampling
#' probability and uniform subsampling probability. The actual subsample
#' probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. This protects the estimator from extreme small
#' subsampling probability. The default value is 0.
#' 
#' - `b` is a positive number which is used to constaint the poisson subsampling probability. `b` close to 0 results in subsampling probabilities closer to uniform probability \eqn{\frac{1}{N}}. `b=2` is the default value. See relevant references for further details.
#' 
#' @param ... A list of parameters which will be passed to \code{svyglm()}. 
#'
#' @return
#' ssp.relogit returns an object of class "ssp.relogit" containing the following components (some are optional):
#' 
#' \describe{
#'   \item{model.call}{The original function call.}
#'   \item{coef.plt}{The pilot estimator. See Details for more information.}
#'   \item{coef.ssp}{The estimator obtained from the optimal subsample.}
#'   \item{coef}{The weighted linear combination of `coef.plt` and `coef.ssp.` The combination weights depend on the relative size of `n.plt` and `n.ssp` and the estimated covariance matrices of `coef.plt` and `coef.ssp.` We blend the pilot subsample information into optimal subsample estimator since the pilot subsample has already been drawn. The coefficients and standard errors reported by summary are `coef` and the square root of `diag(cov)`.}
#'   \item{cov.ssp}{The covariance matrix of \code{coef.ssp}}
#'   \item{cov}{The covariance matrix of \code{beta.cmb}}
#'   \item{index.plt}{Row indices of pilot subsample in the full dataset.}
#'   \item{index.ssp}{Row indices of of optimal subsample in the full dataset.}
#'   \item{N}{The number of observations in the full dataset.}
#'   \item{subsample.size.expect}{The expected subsample size.}
#'   \item{terms}{The terms object for the fitted model.}
#' }
#'
#' @details
#' 'Rare event' stands for the number of observations where \eqn{Y=1} is rare compare to the number of \eqn{Y=0} in the full data. In the face of logistic regression with rare events, @wang2021nonuniform shows that the available information ties to the number of positive instances instead of the full data size. Based on this insight, one can keep all the rare instances and perform subsampling on the non-rare instances to reduce the computational cost. When \code{criterion = optA, optL or LCC}, all observations with \eqn{Y=1} are preserved and it draw `n.ssp` subsmples from observations with Y=0. When \code{criterion = uniform}, it draws (n.plt+n.ssp) subsmples from the full sample with equal sampling probability. 
#' 
#' A pilot estimator for the unknown parameter  \eqn{\beta} is required because both optA and
#' optL subsampling probabilities depend on \eqn{\beta}. This
#' is achieved by drawing half size subsample from rare observations and half from non-rare observations.
#' 
#' Most of the arguments and returned variables have similar meaning with \link{ssp.glm}. Refer to [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html)
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


