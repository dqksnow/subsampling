#' Optimal Subsampling Methods for Generalized Linear Models
#' @description
#' Draw subsample from full dataset and fit a generalized linear model (GLM) on the subsample. For a quick start, refer to the [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html).
#'
#' @param formula A model formula object of class "formula" that describes the model to be fitted.
#' @param data A data frame containing the variables in the model. Denote \eqn{N} as the number of observations in `data`.
#' @param subset An optional vector specifying a subset of observations from `data` to use for the analysis. This subset will be viewed as the full data.
#' @param n.plt The pilot subsample size (first-step subsample size).
#' This subsample is used to compute the pilot estimator and estimate the optimal subsampling probabilities.
#' @param n.ssp The expected size of the optimal subsample (second-step subsample). For `sampling.method = 'withReplacement'`, The exact subsample size is `n.ssp`. For `sampling.method = 'poisson'`, `n.ssp` is the expected subsample size. 
#' @param family `family` can be a character string naming a family function, a family function or the result of a call to a family function. 
#' 
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
#' @param sampling.method The sampling method to use. Options include `withReplacement`
#' and `poisson` (default). `withReplacement` draws exactly `n.ssp`
#'   subsamples from size \eqn{N} full dataset with replacement, using the specified
#' subsampling probabilities. `poisson` draws observations independently by
#' comparing each subsampling probability with a realization of uniform random
#' variable  \eqn{U(0,1)}.
#' Differences between methods:
#' 
#' - Sample size: `withReplacement` draws exactly  `n.ssp` subsamples while `poisson` draws
#' subsamples with expected size `n.ssp`, meaning the actual size may vary.
#' 
#' - Memory usage: `withReplacement` requires the entire dataset to be loaded at once, while `poisson`
#' allows for processing observations sequentially (will be implemented in future version). 
#' 
#' - Estimator performance: Theoretical results show that the `poisson` tends to get a
#' subsample estimator with lower asymptotic variance compared to the
#' `withReplacement`
#' 
#' @param likelihood The likelihood function to use. Options include `weighted` (default) and
#' `logOddsCorrection`. A bias-correction likelihood function is required for subsample since unequal subsampling probabilities introduce bias.
#' 
#' - `weighted` Applies a weighted likelihood function where each observation is weighted by the inverse of its subsampling probability.
#' 
#' - `logOddsCorrection` This lieklihood is available only for logistic regression model (i.e., when family is binomial or quasibinomial). It uses a conditional likelihood, where each element of the likelihood represents the probability of \eqn{Y=1}, given that this subsample was drawn.
#'   
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, `contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')`.
#' 
#' @param control The argument `control` contains two tuning parameters `alpha` and `b`. 
#' 
#' - `alpha` \eqn{\in [0,1]} is the mixture weight of the user-assigned subsampling
#' probability and uniform subsampling probability. The actual subsample
#' probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. This protects the estimator from extreme small
#' subsampling probability. The default value is 0.
#' 
#' - `b` is a positive number which is used to constaint the poisson subsampling probability. `b` close to 0 results in subsampling probabilities closer to uniform probability \eqn{\frac{1}{N}}. `b=2` is the default value. See relevant references for further details.
#' 
#' @param ... A list of parameters which will be passed to `svyglm()`. 
#'
#' @return
#' `ssp.glm` returns an object of class "ssp.glm" containing the following components (some are optional):
#' 
#' \describe{
#'   \item{model.call}{The original function call.}
#'   \item{coef.plt}{The pilot estimator. See Details for more information.}
#'   \item{coef.ssp}{The estimator obtained from the optimal subsample.}
#'   \item{coef}{The weighted linear combination of `coef.plt` and `coef.ssp`. The combination weights depend on the relative size of `n.plt` and `n.ssp` and the estimated covariance matrices of `coef.plt` and `coef.ssp.` We blend the pilot subsample information into optimal subsample estimator since the pilot subsample has already been drawn. The coefficients and standard errors reported by summary are `coef` and the square root of `diag(cov)`.}
#'   \item{cov.ssp}{The covariance matrix of `coef.ssp`.}
#'   \item{cov}{The covariance matrix of `coef`.}
#'   \item{index.plt}{Row indices of pilot subsample in the full dataset.}
#'   \item{index.ssp}{Row indices of of optimal subsample in the full dataset.}
#'   \item{N}{The number of observations in the full dataset.}
#'   \item{subsample.size.expect}{The expected subsample size, equals to `n.ssp` for `ssp.glm.` Note that for other functions, such as \link{ssp.relogit}, this value may differ.}
#'   \item{terms}{The terms object for the fitted model.}
#' }
#'
#' @details
#' 
#' A pilot estimator for the unknown parameter  \eqn{\beta} is required because both optA and
#' optL subsampling probabilities depend on \eqn{\beta}. There is no "free lunch" when determining optimal subsampling probabilities. Fortunately the
#' pilot estimator only needs to satisfy mild conditions. For logistic regression, this
#' is achieved by drawing a size `n.plt` subsample with replacement from full
#' dataset. The case-control subsample probability is applied, that is, \eqn{\pi_i =
#'   \frac{1}{2N_1}} for  \eqn{Y_i=1} and  \eqn{\pi_i = \frac{1}{2N_0}} for  \eqn{Y_i=0},
#'  \eqn{i=1,...,N}, where\eqn{N_0} and \eqn{N_1} are the counts of observations with \eqn{Y = 0} and \eqn{Y = 1}, respectively. For other
#' families, uniform subsampling probabilities are applied. Typically, `n.plt` is
#' relatively small compared to `n.ssp`.
#' 
#' When `criterion = 'uniform'`, there is no need to compute the pilot estimator. In this case, a size `n.plt + n.ssp` subsample will be drawn with uniform sampling probability and `coef` is the corresponding  estimator.
#' 
#' As suggested by `survey::svyglm()`, for binomial and poisson families, use `family=quasibinomial()` and `family=quasipoisson()` to avoid a warning "In eval(family$initialize) : non-integer #successes in a binomial glm!". The quasi versions of the family objects give the same point estimates and suppress the warning. Since subsampling methods only rely on point estimates from svyglm() for further computation, using the quasi families does not introduce any issues.
#' 
#' For Gamma family, `ssp.glm` returns only the estimation of coefficients, as the dispersion parameter is not estimated.
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
#' corr <- 0.5
#' sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
#' X <- MASS::mvrnorm(N, rep(0, d), sigmax)
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
#' family = 'quasibinomial', 
#' criterion = "optL",
#' sampling.method = 'withReplacement', 
#' likelihood = "weighted")
#' summary(subsampling.results)
#' Uni.subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'quasibinomial', 
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
