#' Optimal Subsampling Methods for Quantile Regression Model
#' @description
#' Draw subsample from full dataset and fit quantile regression model. For a quick start, refer to the [vignette](https://dqksnow.github.io/subsampling/articles/ssp-quantreg.html).
#'
#'
#' @param formula A model formula object of class "formula" that describes the model to be fitted.
#' @param data A data frame containing the variables in the model. Denote \eqn{N} as the number of observations in `data`.
#' @param subset An optional vector specifying a subset of observations from `data` to use for the analysis. This subset will be viewed as the full data.
#' @param tau The interested quantile.
#' @param n.plt The pilot subsample size (first-step subsample size).
#' This subsample is used to compute the pilot estimator and estimate the optimal subsampling probabilities.
#' @param n.ssp The expected size of the optimal subsample (second-step subsample). For `sampling.method = 'withReplacement'`, The exact subsample size is `n.ssp`. For `sampling.method = 'poisson'`, `n.ssp` is the expected subsample size. 
#' @param B The number of subsamples for the iterative sampling algorithm. Each subsample contains \code{n.ssp} observations. This allows us to estimate the covariance matrix.
#' @param boot If TRUE then perform iterative sampling algorithm and estimate the covariance matrix. If FALSE then only one subsample with size `B*n.ssp` is returned.
#' @param criterion It determines how subsampling probabilities are computed.
#' Choices include \code{optL}(default) and \code{uniform}. 
#' 
#' - `optL` Minimizes the trace of a transformation of the asymptotic covariance matrix of the subsample estimator.
#' 
#' - `uniform` Assigns equal subsampling probability
#' \eqn{\frac{1}{N}} to each observation, serving as a baseline subsampling strategy.
#' 
#' @param sampling.method The sampling method for drawing the optimal subsample.
#' Choices include \code{withReplacement} and \code{poisson}(default). `withReplacement` draws exactly `n.ssp`
#'   subsamples from size \eqn{N} full dataset with replacement, using the specified
#' subsampling probabilities. `poisson` draws observations independently by
#' comparing each subsampling probability with a realization of uniform random
#' variable  \eqn{U(0,1)}.
#' 
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator. Currently \code{weighted} is implemented which applies a weighted likelihood function where each observation is weighted by the inverse of its subsampling probability.
#' 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control The argument `control` contains two tuning parameters `alpha` and `b`. 
#' 
#' - `alpha` \eqn{\in [0,1]} is the mixture weight of the user-assigned subsampling
#' probability and uniform subsampling probability. The actual subsample
#' probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. This protects the estimator from extreme small
#' subsampling probability. The default value is 0.
#' 
#' - `b` is a positive number which is used to constraint the poisson subsampling probability. `b` close to 0 results in subsampling probabilities closer to uniform probability \eqn{\frac{1}{N}}. `b=2` is the default value.
#' See relevant references for further details.
#' 
#' @param ... A list of parameters which will be passed to \code{quantreg::rq()}. 
#' 
#' @return
#' `ssp.quantreg` returns an object of class "ssp.quantreg" containing the following components (some are optional):
#' 
#' \describe{
#'   \item{model.call}{The original function call.}
#'   \item{coef.plt}{The pilot estimator. See Details for more information.}
#'   \item{coef}{The estimator obtained from the optimal subsample.}
#'   \item{cov}{The covariance matrix of \code{coef}}
#'   \item{index.plt}{Row indices of pilot subsample in the full dataset.}
#'   \item{index.ssp}{Row indices of of optimal subsample in the full dataset.}
#'   \item{N}{The number of observations in the full dataset.}
#'   \item{subsample.size.expect}{The expected subsample size}
#'   \item{terms}{The terms object for the fitted model.}
#' }
#' 
#' @details
#' Most of the arguments and returned variables have the same meaning with \link{ssp.glm}. Refer to [vignette](https://dqksnow.github.io/subsampling/articles/ssp-logit.html)
#' 
#' A pilot estimator for the unknown parameter \eqn{\beta} is required because 
#' optL subsampling probabilities depend on \eqn{\beta}. There is no "free lunch" when determining optimal subsampling probabilities. For quantile regression, this
#' is achieved by drawing a size `n.plt` subsample with replacement from full
#' dataset, using uniform sampling probability.
#' 
#' If `boot`=TRUE, the returned value `subsample.size.expect` equals to `B*n.ssp`, and the covariance matrix for `coef` would be calculated. 
#' If `boot`=FALSE, the returned value `subsample.size.expect` equals to `B*n.ssp`, but the covariance matrix won't be estimated. 
#' 
#' @references
#' Wang, H., & Ma, Y. (2021). Optimal subsampling for quantile regression in big data. \emph{Biometrika}, \strong{108}(1), 99-112.
#' 
#' @examples
#' #quantile regression
#' set.seed(1)
#' N <- 1e4
#' B <- 5
#' tau <- 0.75
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
#' colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""))
#' formula <- Y ~ .
#' n.plt <- 200
#' n.ssp <- 100
#' optL.results <- ssp.quantreg(formula,data,tau = tau,n.plt = n.plt,
#' n.ssp = n.ssp,B = B,boot = TRUE,criterion = 'optL',
#' sampling.method = 'withReplacement',likelihood = 'weighted')
#' summary(optL.results)
#' uni.results <- ssp.quantreg(formula,data,tau = tau,n.plt = n.plt,
#' n.ssp = n.ssp,B = B,boot = TRUE,criterion = 'uniform',
#' sampling.method = 'withReplacement', likelihood = 'weighted')
#' summary(uni.results)
#' @export
ssp.quantreg <- function(formula,
                         data,
                         subset = NULL,
                         tau = 0.5,
                         n.plt,
                         n.ssp,
                         B = 5,
                         boot = TRUE,
                         criterion = 'optL',
                         sampling.method = 'withReplacement',
                         likelihood = c('weighted'),
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
  ## avoid problems with 1D arrays, but keep names
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
  criterion <- match.arg(criterion, c('optL', 'uniform'))
  sampling.method <- match.arg(sampling.method, c('poisson', 'withReplacement'))
  likelihood <- match.arg(likelihood, c('weighted'))

  ## check subsample size
  if (n.ssp * B > 0.1 * N) {
    warning("The total subsample size n.ssp*B exceeds the recommended maximum
    value (10% of full sample size).")
  }
  
  if (boot == FALSE | B == 1){
    n.ssp <- n.ssp * B
    B <- 1
    boot <- FALSE
  }
  
  control <- do.call("quantreg.control", control)
  
  ## create a list to store variables
  inputs <- list(X = X, Y = Y, tau = tau, N = N, d = d,
                 n.plt = n.plt, n.ssp = n.ssp, B = B, boot = boot, 
                 criterion = criterion, sampling.method = sampling.method,
                 likelihood = likelihood,
                 control = control
                 )
  
  if (criterion %in% c("optL")) {
    
    ## pilot step
    plt.results <- plt.estimation.quantreg(inputs, ...)
    beta.plt <- plt.results$beta.plt
    Ie.full <- plt.results$Ie.full
    index.plt <- plt.results$index.plt

    ## subsampling and boot step
    ssp.results <- ssp.estimation.quantreg(inputs,
                                           Ie.full = Ie.full,
                                           index.plt = index.plt,
                                           ...
                                           )
    Betas.ssp <- ssp.results$Betas.ssp
    beta.ssp <- ssp.results$beta.ssp
    est.cov.ssp <- ssp.results$cov.ssp
    index.ssp <- ssp.results$index.ssp
    
    names(beta.ssp) <- names(beta.plt) <- colnames(X)

    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef = beta.ssp,
                    cov = est.cov.ssp,
                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    N = N,
                    subsample.size.expect = c(n.ssp, B),
                    terms = mt
                    )
    class(results) <- c("ssp.quantreg", "list")
    return(results)
  } else if (criterion == "uniform"){
    ## subsampling and boot step
    uni.results <- ssp.estimation.quantreg(inputs, ...)
    Betas.uni <- uni.results$Betas.ssp
    beta.uni <- uni.results$beta.ssp
    est.cov.uni <- uni.results$cov.ssp
    index.uni <- uni.results$index.ssp
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    coef.plt = NA,
                    coef = beta.uni,
                    cov = est.cov.uni,
                    index = index.uni,
                    N = N,
                    subsample.size.expect = c(n.ssp, B),
                    terms = mt
                    )
    class(results) <- c("ssp.quantreg", "list")
    return(results)
  }
}
###############################################################################