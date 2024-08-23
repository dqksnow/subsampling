#' Optimal Subsampling Method for Softmax(multinomial logistic) Regression Model
#' @description
#' Draw subsample from full dataset and fit softmax(multinomial logistic) regression model on subsample.
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
#' @param criterion The criterion of optimal subsampling probabilities.
#' Choices include \code{optA}, \code{optL}, \code{MSPE}(default), \code{LUC} and \code{uniform}. 
#' @param sampling.method The sampling method for drawing the optimal subsample. 
#' Choices include \code{withReplacement} and \code{poisson}(default).
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator. Choices include 
#'  \code{weighted} and \code{MSCLE}(default). 
#' @param constraint The constraint for identifiability of softmax model. Choices include 
#'  \code{baseline} and \code{summation}(default). 
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control A list of parameters for controlling the sampling process. Default is \code{list(alpha=0, b=2)}.
#' @param ... A list of parameters which will be passed to \code{nnet::multinom()}. 
#'
#' @return
#' ssp.softmax returns an object of class "ssp.softmax" containing the following components (some are optional):
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
#' @details
#' The returned coefficients are under baseline constraints, that is, the first category is baseline with zero coefficient.
#' 
#' \code{criterion = MSPE} stands for optimal subsampling probabilities by minimizing the Mean Squared Prediction Error which is immune to the choice of model constraint. \code{criterion = LUC} stands for local uncertainty sampling. Refer to Yao, Zou and Wang (2023B).
#'
#' In \code{control}, alpha is the mixture proportions of optimal subsampling probability and uniform sampling probability. b is the parameter controls the upper threshold for optimal subsampling probability. 
#'
#' Most of the arguments and returned variables have the same meaning with \link{ssp.glm}. Also refer to [vignette](https://dqksnow.github.io/Subsampling/articles/ssp-logit.html)
#'
#' @references
#' Yao, Y., & Wang, H. (2019). Optimal subsampling for softmax regression. \emph{Statistical Papers}, \strong{60}, 585-599.
#' 
#' Han, L., Tan, K. M., Yang, T., & Zhang, T. (2020). Local uncertainty sampling for large-scale multiclass logistic regression. \emph{Annals of Statistics}, \strong{48}(3), 1770-1788.
#' 
#' Yao, Y., Zou, J., & Wang, H. (2023). Optimal poisson subsampling for softmax regression. \emph{Journal of Systems Science and Complexity}, \strong{36}(4), 1609-1625.
#' 
#' Yao, Y., Zou, J., & Wang, H. (2023). Model constraints independent optimal subsampling probabilities for softmax regression. \emph{Journal of Statistical Planning and Inference}, \strong{225}, 188-201.
#' 
#' @examples
#' # softmax regression
#' d <- 3 # dim of covariates
#' K <- 2 # K + 1 classes
#' G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
#' N <- 1e4
#' beta.true.baseline <- cbind(rep(0, d), matrix(-1.5, d, K))
#' beta.true.summation <- cbind(rep(1, d), 0.5 * matrix(-1, d, K))
#' set.seed(1)
#' mu <- rep(0, d)
#' sigma <- matrix(0.5, nrow = d, ncol = d)
#' diag(sigma) <- rep(1, d)
#' X <- MASS::mvrnorm(N, mu, sigma)
#' prob <- exp( X %*% beta.true.summation)
#' prob <- prob / rowSums(prob)
#' Y <- apply(prob, 1, function(row) sample(0:K, size = 1, prob = row))
#' n.plt <- 500
#' n.ssp <- 1000
#' data <- as.data.frame(cbind(Y, X))
#' colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""))
#' head(data)
#' formula <- Y ~ . -1
#' WithRep.MSPE <- ssp.softmax(formula = formula,
#'  data = data, 
#'  n.plt = n.plt,
#'  n.ssp = n.ssp,
#'  criterion = 'MSPE', 
#'  sampling.method = 'withReplacement',
#'  likelihood = 'weighted',
#'  constraint = 'baseline')
#' summary(WithRep.MSPE)
#' @export
ssp.softmax <- function(formula, 
                        data,
                        subset,
                        n.plt,
                        n.ssp,
                        criterion = 'MSPE',
                        sampling.method = 'poisson',
                        likelihood = 'MSCLE',
                        constraint = 'summation',
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

  dimension <- dim(X)
  N <- dimension[1]
  d <- dimension[2]
  K <- length(unique(Y)) - 1
  G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d) 
  ## G: transformation matrix
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(c(1:length(Y)), Y)] <- 1
  
  criterion <- match.arg(criterion,c('optL', 'optA', 'MSPE', 'LUC', 'uniform'))
  sampling.method <- match.arg(sampling.method, c('poisson', 'withReplacement'))
  likelihood <- match.arg(likelihood, c('weighted', 'MSCLE'))
  constraint <- match.arg(constraint, c('baseline', 'summation'))
  control <- do.call("softmax.control", control)
  
  ## create a list to store variables
  inputs <- list(X = X, Y = Y, Y.matrix = Y.matrix,
                 N = N, d = d, K = K, G = G, 
                 n.plt = n.plt, n.ssp = n.ssp, 
                 criterion = criterion, sampling.method = sampling.method,
                 likelihood = likelihood, constraint = constraint,
                 control = control
                 )
  
  if (criterion %in% c('optL', 'optA', 'MSPE', 'LUC')) {
    plt.estimate.results <- softmax.plt.estimate(inputs, ...)
    p.plt <- plt.estimate.results$p.plt
    beta.plt.b <- plt.estimate.results$beta.plt
    P1.plt <- plt.estimate.results$P1.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    index.plt <- plt.estimate.results$index.plt
    cov.plt.b <- plt.estimate.results$cov.plt
    Omega.plt <- plt.estimate.results$Omega.plt
    Lambda.plt <- plt.estimate.results$Lambda.plt

    ## subsampling step
    ssp.results <- softmax.subsampling(inputs,
                                       p.plt = p.plt,
                                       ddL.plt = ddL.plt,
                                       P1.plt = P1.plt,
                                       Omega.plt = Omega.plt,
                                       index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    p.ssp <- ssp.results$p.ssp
    offsets <- ssp.results$offsets

    ## subsample estimating step
    ssp.estimate.results <-
      softmax.subsample.estimate(inputs,
                                 index.ssp = index.ssp,
                                 p.ssp = p.ssp[index.ssp],
                                 offsets = offsets,
                                 beta.plt = beta.plt,
                                 ...)
    beta.ssp.b <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp.b <- ssp.estimate.results$cov.ssp

    ## combining step
    combining.results <- softmax.combining(inputs,
                                           n.ssp = length(index.ssp),
                                           ddL.plt = ddL.plt,
                                           ddL.ssp = ddL.ssp,
                                           dL.sq.plt = dL.sq.plt,
                                           dL.sq.ssp = dL.sq.ssp,
                                           Lambda.plt = Lambda.plt,
                                           Lambda.ssp = Lambda.ssp,
                                           beta.plt = beta.plt.b,
                                           beta.ssp = beta.ssp.b)
    beta.cmb.b <- combining.results$beta.cmb
    cov.cmb.b <- combining.results$cov.cmb
    P.cmb <- combining.results$P.cmb

    # beta.plt.b <- G %*% as.vector(beta.plt.b)
    # beta.ssp.b <- G %*% as.vector(beta.ssp.b)
    # beta.cmb.b <- G %*% as.vector(beta.cmb.b)
    # cov.plt.b <- G %*% cov.plt.b %*% t(G)
    # cov.ssp.b <- G %*% cov.ssp.b %*% t(G)
    # cov.cmb.b <- G %*% cov.cmb.b %*% t(G)
    
    beta.plt <- matrix(beta.plt.b, nrow = d)
    beta.ssp = matrix(beta.ssp.b, nrow = d)
    beta <- matrix(beta.cmb.b, nrow = d)
    rownames(beta.plt) <- rownames(beta.ssp) <- rownames(beta) <- colnames(X)
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef = beta,
                    # cov.plt = cov.plt.b,
                    cov.ssp = cov.ssp.b,
                    cov = cov.cmb.b,
                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    N = N,
                    subsample.size.expect = n.ssp,
                    terms = mt
                    )
    class(results) <- c("ssp.softmax", "list")
    return(results)
  } else if (criterion == "uniform"){
    n.uni <- n.plt + n.ssp
    if (sampling.method == 'withReplacement') {
      index.uni <- random.index(N, n.uni)
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni, ...)
      beta.uni.b <- results$beta
      P.uni <-results$P1
      ddL.uni <- softmax_ddL_cpp(X = x.uni, P = P.uni[, -1], p = rep(1, n.uni),
                                 K, d, scale = N*n.uni)
      dL.sq.uni <- softmax_dL_sq_cpp(X = x.uni, 
                                     Y_matrix = Y.matrix[index.uni, ],
                                     P = P.uni[, -1], p = rep(1, n.uni), K = K, 
                                     d = d, scale = N^2*n.uni^2)
      c <- n.uni / N
      cov.uni.b <- solve(ddL.uni) %*% (dL.sq.uni * (1+c)) %*% solve(ddL.uni)
    } else if (sampling.method == 'poisson') {
      index.uni <- poisson.index(N, n.uni/N)
      p.uni <- rep(n.uni / N, length(index.uni))
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni, ...)
      beta.uni.b <- results$beta
      P.uni <- results$P1
      ddL.uni <- softmax_ddL_cpp(X = x.uni, P = P.uni[, -1], p = p.uni, K, d,
                             scale = N)
      dL.sq.uni <- softmax_dL_sq_cpp(X = x.uni, 
                                     Y_matrix = Y.matrix[index.uni, ],
                                     P = P.uni[, -1], p = p.uni, K = K, d = d,
                                     scale=(N^2))
      cov.uni.b <- solve(ddL.uni) %*% dL.sq.uni %*% solve(ddL.uni)
    }
    beta <- matrix(beta.uni.b, nrow = d)
    rownames(beta) <- colnames(X)
    results <- list(model.call = model.call,
                    index = index.uni,
                    coef = beta,
                    cov = cov.uni.b,
                    # P = P.uni,
                    N = N,
                    subsample.size.expect = n.uni,
                    terms = mt
                    )
    class(results) <- c("ssp.softmax", "list")
    return(results)
  }
}
