#' Optimal Subsampling Methods for Quantile Regression Model
#' @description
#' This function fits generalized linear models using ....
#'
#'
#' @param formula An object of class "formula" which describes the model to be
#'  fitted.
#' @param data A data frame containing the variables in the model.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param tau The quantile to be estimated,
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator as well as to
#' estimate the optimal subsampling probability.
#' @param n.ssp The expectation optimal subsample size (the second-step subsample
#' size). For \code{sampling.method = 'withReplacement'}, \code{n.ssp} is exactly the subsample size. For \code{sampling.method = 'poisson'}, \code{n.ssp} is the expectation of subsample size. 
#' @param B The number of subsamples for the iterative sampling algorithm. Each subsample contains \code{n.ssp} observations. This allows us to estimate the covariance matrix.
#' @param boot If TRUE then perform iterative sampling algorithm and estimate the covariance matrix. If FALSE then only one subsample is returned which contains n.ssp observations.
#' @param criterion The criterion of optimal subsampling probabilities.
#' Choices include \code{optL}(default) and \code{uniform}. 
#' @param sampling.method The sampling method for drawing the optimal subsample. 
#' Choices include \code{withReplacement} and \code{poisson}(default).
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator. Currently it uses \code{weighted}. 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control A list of parameters for controlling the sampling process. Default is \code{list(alpha=0, b=2)}.
#' @param ... A list of parameters which will be passed to \code{quantreg::rq()}. 
#' 
#' @return
#' ssp.quantreg returns an object of class "ssp.quantreg" containing the following components (some are optional):
#' 
#' \describe{
#'   \item{model.call}{model call}
#'   \item{beta.plt}{pilot estimator}
#'   \item{coefficients}{optimal subsample estimator.}
#'   \item{cov}{covariance matrix of \code{coefficients}}
#'   \item{index.plt}{index of pilot subsample in the full sample}
#'   \item{index.ssp}{index of optimal subsample in the full sample}
#'   \item{N}{number of observations in the full sample}
#'   \item{subsample.size.expect}{expected subsample size}
#'   \item{terms}{model terms}
#' }
#' @details
#' Additional details... briefly introduce the idea.
#' 
#' @references
#' Wang, H., & Ma, Y. (2021). Optimal subsampling for quantile regression in big data. \emph{Biometrika}, \strong{108}(1), 99-112. \doi{https://doi.org/10.1093/biomet/asaa043}
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
  colnames(X)[1] <- "Intercept"
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
                    beta.plt = beta.plt,
                    coefficients = beta.ssp,
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
                    beta.plt = NA,
                    coefficients = beta.uni,
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