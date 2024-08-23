#' Optimal Subsampling Methods for Generalized Linear Models
#' @description
#' Draw subsample from full dataset and fit glm on subsample. Refer to [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html) for a quick start.
#'
#' @param formula An object of class "formula" which describes the model to be
#'  fitted.
#' @param data A data frame containing the variables in the model.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator as well as to
#' estimate the optimal subsampling probability.
#' @param n.ssp The expectation optimal subsample size (the second-step subsample
#' size). For \code{sampling.method = 'withReplacement'}, \code{n.ssp} is exactly the subsample size. For \code{sampling.method = 'poisson'}, \code{n.ssp} is the expectation of subsample size. 
#' @param family \code{family} can be a character string naming a family function, a family function or the result of a call to a family function. Currently \code{'binomial'}, \code{'poisson'} and \code{'Gamma'} are implemented. 
#' @param criterion The criterion of optimal subsampling probabilities.
#' Choices include \code{optA}, \code{optL}(default), \code{LCC} and \code{uniform}. 
#' @param sampling.method The sampling method for drawing the optimal subsample. 
#' Choices include \code{withReplacement} and \code{poisson}(default).
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator. Choices include 
#'  \code{weighted} and \code{logOddsCorrection}. 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control A list of parameters for controlling the sampling process. Default is \code{list(alpha=0, b=2)}.
#' @param ... A list of parameters which will be passed to \code{svyglm()}. 
#'
#' @return
#' ssp.glm returns an object of class "ssp.glm" containing the following components (some are optional):
#' 
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
#' As suggested by \code{survey::svyglm()}, for binomial and poisson families use \code{family=quasibinomial()} and \code{family=quasipoisson()} to avoid a warning "In eval(family$initialize) : non-integer #successes in a binomial glm!". The warning is due to the non-integer survey weights. The ‘quasi’ versions of the family objects give the same point estimates and do not give the warning. Subsampling methods only use point estimates from \code{svyglm()} for further computation so that would not bring problems. 
#' 
#' For Gamma family, it will only return the estimation of coefficients, not dispersion parameter.
#' 
#' \code{likelihood = logOddsCorrection} is implemented only for logistic regression (family = binomial or quasibonomial).
#'
#' In \code{control}, `alpha` is the mixture proportions of optimal subsampling probability and 
#' uniform sampling probability. `b` is the parameter controls the upper 
#' threshold for optimal subsampling probability. 
#'
#' @references
#' Wang, H. (2019). More efficient estimation for logistic regression with optimal subsamples. \emph{Journal of machine learning research}, \strong{20}(132), 1-59.
#' 
#' Ai, M., Yu, J., Zhang, H., & Wang, H. (2021). Optimal subsampling algorithms for big data regressions. \emph{Statistica Sinica}, \strong{31}(2), 749-772.
#' 
#' Wang, H., & Kim, J. K. (2022). Maximum sampled conditional likelihood for informative subsampling. \emph{Journal of machine learning research}, \strong{23}(332), 1-50.
#' 
#' @examples
#' # logistic regression
#' set.seed(2)
#' N <- 1e4
#' beta0 <- rep(-0.5, 7)
#' d <- length(beta0) - 1
#' X <- matrix(0, N, d)
#' generate_rexp <- function(x) x <- rexp(N, rate = 2)
#' X <- apply(X, 2, generate_rexp)
#' Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' n.plt <- 500
#' n.ssp <- 1000
#' subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'quasibinomial',
#' criterion = "optL",
#' sampling.method = 'poisson',
#' likelihood = "logOddsCorrection")
#' summary(subsampling.results)
#' subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'binomial', 
#' criterion = "optL",
#' sampling.method = 'withReplacement', 
#' likelihood = "weighted")
#' summary(subsampling.results)
#' Uni.subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'binomial', 
#' criterion = 'uniform')
#' summary(Uni.subsampling.results)
#' ####################
#' # poisson regression
#' set.seed(1)
#' N <-  1e4
#' beta0 <- rep(0.5, 7)
#' d <- length(beta0) - 1
#' X <- matrix(runif(N * d), N, d)
#' epsilon <- runif(N)
#' lambda <- exp(beta0[1] + X %*% beta0[-1])
#' Y <- rpois(N, lambda)
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' n.plt <- 200
#' n.ssp <- 600
#' subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'poisson',
#' criterion = "optL", 
#' sampling.method = 'poisson',
#' likelihood = "weighted")
#' summary(subsampling.results)
#' subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'poisson', 
#' criterion = "optL", 
#' sampling.method = 'withReplacement',
#' likelihood = "weighted")
#' summary(subsampling.results)
#' Uni.subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'poisson', 
#' criterion = 'uniform')
#' summary(Uni.subsampling.results)
#' ##################
#' # gamma regression
#' set.seed(1)
#' N <- 1e4
#' p <- 3
#' beta0 <- rep(0.5, p + 1)
#' d <- length(beta0) - 1
#' shape <- 2
#' X <- matrix(runif(N * d), N, d)
#' link_function <- function(X, beta0) 1 / (beta0[1] + X %*% beta0[-1])
#' scale <- link_function(X, beta0) / shape
#' Y <- rgamma(N, shape = shape, scale = scale)
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' n.plt <- 200
#' n.ssp <- 1000
#' subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'Gamma',
#' criterion = "optL", 
#' sampling.method = 'poisson',
#' likelihood = "weighted")
#' summary(subsampling.results)
#' @export
ssp.glm <- function(formula,
                    data,
                    subset = NULL,
                    n.plt,
                    n.ssp,
                    family = 'binomial',
                    criterion = 'optL',
                    sampling.method = 'poisson',
                    likelihood = 'weighted',
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
  family <- match.arg(family, c('binomial', 'poisson', 'Gamma',
                                'quasibinomial', 'quasipoisson'))
  criterion <- match.arg(criterion, c('optL', 'optA', 'LCC', 'uniform'))
  sampling.method <- match.arg(sampling.method, c('poisson', 'withReplacement'))
  likelihood <- match.arg(likelihood, c('logOddsCorrection', 'weighted'))
  control <- do.call("glm.control", control)
  
  ## family 
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  N <- nrow(X)
  d <- ncol(X)
  inputs <- list(X = X, Y = Y, N = N, d = d,
                 n.plt = n.plt, n.ssp = n.ssp,
                 criterion = criterion, sampling.method = sampling.method,
                 likelihood = likelihood, family = family,
                 control = control
                 )
  
  if (criterion %in% c('optL', 'optA', 'LCC')) {
    ## pilot step
    plt.estimate.results <- pilot.estimate(inputs, ...)
    p.plt <- plt.estimate.results$p.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    ddL.plt.correction <- plt.estimate.results$ddL.plt.correction
    Lambda.plt <- plt.estimate.results$Lambda.plt
    d.psi <- plt.estimate.results$d.psi
    index.plt <- plt.estimate.results$index.plt
    
    ## subsampling step
    ssp.results <- subsampling(inputs,
                               p.plt = p.plt,
                               ddL.plt.correction = ddL.plt.correction,
                               d.psi = d.psi,
                               index.plt = index.plt
                               )
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    offset <- ssp.results$offset
    
  
    ## subsample estimating step
    ssp.estimate.results <- subsample.estimate(inputs,
                                               w.ssp = w.ssp,
                                               offset = offset,
                                               beta.plt = beta.plt,
                                               index.ssp = index.ssp,
                                               ...
                                               )
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp <- ssp.estimate.results$cov.ssp

    ## combining step
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
    cov.cmb <- combining.results$cov.cmb
    
    names(beta.cmb) <- names(beta.ssp) <- names(beta.plt) <- colnames(X)
    
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef = beta.cmb,
                    cov.ssp = cov.ssp,
                    cov = cov.cmb,
                    index.plt = index.plt,
                    index = index.ssp,
                    N = N,
                    subsample.size.expect = n.ssp,
                    terms = mt
                    )
    class(results) <- c("ssp.glm", "list")
    return(results)
  } else if (criterion == "uniform"){
    n.uni <- n.plt + n.ssp
    if (sampling.method == 'withReplacement') {
      index.uni <- random.index(N, n.uni)
    } else if (sampling.method == 'poisson') {
      index.uni <- poisson.index(N, n.uni / N)
    }
    x.uni <- X[index.uni, ]
    y.uni = Y[index.uni]
    results.uni <- glm.coef.estimate(X = x.uni, Y = y.uni, 
                                     family = family,
                                     ...)
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
    if (sampling.method == 'withReplacement') {
      cov.uni <- solve(ddL.uni) %*% 
        (dL.sq.plt * (1 + n.uni / N)) %*% solve(ddL.uni)
    } else if (sampling.method == 'poisson') {
      cov.uni <- solve(ddL.uni) %*% dL.sq.plt %*% solve(ddL.uni)
    }
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    index = index.uni,
                    coef = beta.uni,
                    cov = cov.uni,
                    N = N,
                    subsample.size.expect = n.uni,
                    terms = mt
                    )
    class(results) <- c("ssp.glm", "list")
    return(results)
  }
}
