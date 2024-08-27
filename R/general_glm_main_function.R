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
#' size). For `sampling.method = 'withReplacement'`, `n.ssp` is exactly the subsample size. For `sampling.method = 'poisson'`, `n.ssp` is the expectation of subsample size. 
#' @param family `family` can be a character string naming a family function, a family function or the result of a call to a family function. 
#' 
#' @param criterion The choices of `criterion` include `optA`, `optL`(default), `LCC` and `uniform.`
#' 
#' - `optA` subsampling probabilities are derived by minimizing the
#' trace of the asymptotic covariance of subsample estimator. `optL` subsampling
#' probabilities are derived by minimizing the trace of a transportation of the
#' asymptotic covariance of subsample estimator. The computational complexity of
#' optA subsampling probabilities is \eqn{O(N d^2)} while that of optL is \eqn{O(N d)}.
#' 
#' - `LCC` stands for the Local Case-Control subsampling probability, serving as a baseline criterion.
#' 
#' - `uniform` assigns each observation with equal subsampling probability
#' \eqn{\frac{1}{N}}, serving as a baseline criterion.
#' 
#' @param sampling.method The options for the `sampling.method` argument include `withReplacement`
#' and `poisson` (default). `withReplacement` stands for drawing `n.ssp`
#'   subsamples from full dataset of size \eqn{N} with replacement, using the specified
#' subsampling probability. `poisson` stands for drawing subsamples one by one by
#' comparing the subsampling probability with a realization of uniform random
#' variable  \eqn{U(0,1)}. The expected number of drawed samples are \eqn{n.ssp}.
#' The main differences are:
#' 
#' - `withReplacement` draws exactly  `n.ssp` subsamples while `poisson` draws
#' subsamples with expectation `n.ssp` meaning the actual number may vary
#' slightly.
#' 
#' - `withReplacement` requires loading the full dataset at once while `poisson`
#' allows for scanning the dataset one observation at a time.
#' 
#' - Theoretical results showed that the `poisson` method tends to get a
#' subsample estimator with smaller asymptotic variance compared to the
#' `withReplacement` method.
#' 
#' @param likelihood The available choices for `likelihood` include `weighted` (default) and
#' `logOddsCorrection`. The reason we can not use an equally weighted likelihood
#' function for the subsample is that it introduces bias due to the different
#' subsampling probabilities. Therefore, we need to apply methods to correct the
#' bias.
#' 
#' - `weighted` refers to the weighted likelihood function for subsample, where
#'   each observation is weighted by the inverse of its subsampling probability.
#' 
#' - `logOddsCorrection` stands for the conditional likelihood function for the
#'   subsample. "conditional" means that each element in the likelihood function
#'   is the probability of \eqn{Y=1} given that this subsample was drawn. `likelihood = logOddsCorrection` is implemented only for logistic regression (family = binomial or quasibonomial).
#'   
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, `contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')`.
#' @param control The argument `control` contains two tuning parameters `alpha` and `b`. 
#' 
#' - `alpha` \eqn{\in [0,1]} is the mixture weights of the user assigned subsampling
#' probability and uniform subsampling probability. That is, the actual subsample
#' probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. The aim is to
#' protect the subsample estimator from those subsamples with extreme small
#' subsampling probability. The default value of `alpha` is 0.
#' 
#' - `b` is also used to constaint the subsample probability. It can be viewed as
#' the threshold for too large subsample probability. It take values
#' between \eqn{(0,\frac{N}{n})}. `b` close to 0 means subsample probabilities are
#' compressed to uniform probability \eqn{\frac{1}{N}}. `b=2` is the default value
#' and it works well for many cases.
#' 
#' @param ... A list of parameters which will be passed to `svyglm()`. 
#'
#' @return
#' `ssp.glm` returns an object of class "ssp.glm" containing the following components (some are optional):
#' 
#' \describe{
#'   \item{model.call}{model call}
#'   \item{coef.plt}{pilot estimator}
#'   \item{coef.ssp}{optimal subsample estimator}
#'   \item{coef}{the linear combination of `coef.plt` and `coef.ssp.` The combine weights depend on the relative size of `n.plt` and `n.ssp` as well as the estimated covariance matrix of `coef.plt` and `coef.ssp.` We blend the pilot subsample information into optimal subsample estimator since the pilot subsample has been drawn. The coefficients and standard errors printed by summary are `coef` and the square root of `diag(cov)`}
#'   \item{cov.ssp}{covariance matrix of `coef.ssp`}
#'   \item{cov}{covariance matrix of `coef`}
#'   \item{index.plt}{the row index of drawn pilot subsamples in the full data}
#'   \item{index.ssp}{the row index of drawn optimal subsamples in the full data}
#'   \item{N}{number of observations in the full sample}
#'   \item{subsample.size.expect}{`subsample.size.expect` is the expected subsample size which is equals to `n.ssp` when we use `ssp.glm.` In some other models like `ssp.relogit` it might be different.}
#'   \item{terms}{model terms}
#' }
#'
#' @details
#' 
#' A pilot estimator for the unknown parameter  \eqn{\beta} is required because optA and
#' optL subsampling probabilities depend on  \eqn{\beta}. Yes, no free lunch when it
#' comes to the optimal subsampling probabilities. Fortunately we only need the
#' pilot estimator to satisfy some mild conditions. For logistic regression, this
#' is achieved by drawing a size `n.plt` subsample with replacement from full
#' dataset. The case-control subsample probability is applied, that is, \eqn{\pi_i =
#'   \frac{1}{2N_1}} for  \eqn{Y_i=1} and  \eqn{\pi_i = \frac{1}{2N_0}} for  \eqn{Y_i=0},
#'  \eqn{i=1,...,N}, where  \eqn{N_0} is the count of  \eqn{Y=0} and  \eqn{N_1 = N - N_0}. For other
#' families in glm, uniform subsampling probability is used. Typically, `n.plt` is
#' relatively small compared to `n.ssp`.
#' 
#' As suggested by `survey::svyglm()`, for binomial and poisson families use `family=quasibinomial()` and `family=quasipoisson()` to avoid a warning "In eval(family$initialize) : non-integer #successes in a binomial glm!". The warning is due to the non-integer survey weights. The ‘quasi’ versions of the family objects give the same point estimates and do not give the warning. Subsampling methods only use point estimates from `svyglm()` for further computation so that would not bring problems. 
#' 
#' For Gamma family, it will only return the estimation of coefficients, not dispersion parameter.
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
