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
#' @param subset An optional vector specifying a subset of rows to be used
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator as well as to
#' estimate the optimal sampling probability.
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size).
#' @param family defalut = 'binomial'.
#' @param criterion The criterion of optimal subsampling probabilities,
#' currently there are three choices \code{optA}, \code{optL}, and \code{LCC}.
#' @param sampling.method The sampling method for drawing the optimal subsample,
#'  currently there are two choices \code{withReplacement} and \code{poisson}.
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator, currently there are two choices
#'  \code{weighted} and \code{logOddsCorrection}.
#' @param alpha Mixture proportions of optimal subsampling probability and
#' uniform subsampling probability. Default = 0.1.
#' @param b This parameter controls the upper threshold for optimal subsampling probabilities.
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
#'   \item{index.ssp}{index of optimal subsample}
#' }
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
#' family = 'binomial',
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
#' likelihood = "weighted",
#' contrasts = list(F1 = 'contr.sum'))
#' summary(subsampling.results)
#' 
#' Uni.subsampling.results <- ssp.glm(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'binomial', 
#' criterion = 'uniform')
#' summary(Uni.subsampling.results)
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
#' ############################################################################
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
#' family = 'gamma',
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
                    alpha = 0.1,
                    b = 2,
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
  family <- match.arg(family, c('binomial', 'poisson', 'gamma'))
  criterion <- match.arg(criterion, c('optL', 'optA', 'LCC', 'uniform'))
  sampling.method <- match.arg(sampling.method, c('poisson', 'withReplacement'))
  likelihood <- match.arg(likelihood, c('logOddsCorrection', 'weighted'))
  
  if(is.function(family)) family <- family$family
  if(missing(family)) {
    stop("Specify a valid 'family' from c('binomial','poisson',gamma')")
  }
  family <- switch(family,
                  "binomial" = binomial.expand(),
                  "poisson" = poisson.expand(),
                  "gamma" = gamma.expand())

  N <- nrow(X)
  d <- ncol(X)
  if (criterion %in% c('optL', 'optA', 'LCC')) {

    ## pilot step
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
    
    ## subsampling step
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

    ## subsample estimating step
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
    var.cmb <- combining.results$var.cmb
    
    names(beta.cmb) <- names(beta.ssp) <- names(beta.plt) <- colnames(X)
    
    results <- list(model.call = model.call,
                    beta.plt = beta.plt,
                    beta.ssp = beta.ssp,
                    beta = beta.cmb,
                    var.ssp = var.ssp,
                    var = var.cmb,
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
    if (sampling.method == 'withReplacement') {
      var.uni <- solve(ddL.uni) %*% 
        (dL.sq.plt * (1 + n.uni / N)) %*% solve(ddL.uni)
    } else if (sampling.method == 'poisson') {
      var.uni <- solve(ddL.uni) %*% dL.sq.plt %*% solve(ddL.uni)
    }
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    index = index.uni,
                    beta = beta.uni,
                    var = var.uni,
                    N = N,
                    subsample.size.expect = n.uni,
                    terms = mt
                    )
    class(results) <- c("ssp.glm", "list")
    return(results)
  }
}
