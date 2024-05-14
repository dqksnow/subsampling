#' Optimal Subsampling Method for Softmax Regression Models
#'
#' @details details TBD
#' @param X An object of class "formula" which describes the model to be fitted.
#' @param Y A data frame containing the variables in the model. Usually it 
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' @param n.ssp The expected optimal subsample size (the second-step subsample.
#' @param criterion The criterion of optimal subsampling probabilities
#' @param sampling.method The sampling method for drawing the optimal subsample
#' @param estimate.method The type of the maximum likelihood function used to 
#' @param constraint 1
#' @param alpha Mixture proportions of optimal subsampling probability and 
#' @param b This parameter controls the upper threshold for optimal subsampling
#'
#' @return
#' \describe{
#'   \item{beta.plt}{pilot estimator}
#'   \item{beta.ssp}{optimal subsample estimator}
#'   \item{beta.cmb}{combined estimator of \code{beta.plt} and \code{beta.ssp}}
#'   \item{var.plt}{covariance matrix of \code{beta.ssp}}
#'   \item{var.ssp}{covariance matrix of \code{beta.cmb}}
#'   \item{var.cmb}{covariance matrix of \code{beta.cmb}}
#'   \item{P.cmb}{predicted probability matrix of full observations}
#'   \item{index.plt}{index of pilot subsample}
#'   \item{index.ssp}{index of optimal subsample}
#' }

#' @export
#'
#' @examples
#' # softmax regression
#' d <- 3 # dim of covariates
#' K <- 5 # K + 1 classes
#' G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
#' N <- 5e4
#' beta.true <- 0.2 * matrix(-1, d, K)
#' beta.true.sum <- cbind(rep(1, d), beta.true)
#' set.seed(1)
#' mu <- rep(0, d)
#' sigma <- matrix(0.5, nrow = d, ncol = d)
#' diag(sigma) <- rep(1, d)
#' X <- MASS::mvrnorm(N, mu, sigma)
#' prob <- exp( X %*% beta.true.sum)
#' prob <- prob / rowSums(prob)
#' Y <- apply(prob, 1, function(row) sample(0:K, size = 1, prob = row))
#' n.plt <- 500
#' n.ssp <- 1000
#' WithRep.MSPE <- Softmax.subsampling(X, Y, n.plt, n.ssp, criterion = 'MSPE', 
#' sampling.method = 'WithReplacement', estimate.method = 'Weighted',
#' constraint = 'baseline')
#' Poi.MSPE <- Softmax.subsampling(X, Y, n.plt, n.ssp, criterion = 'MSPE',
#' sampling.method = 'Poisson', estimate.method = 'Weighted', 
#' constraint = 'baseline')
#' Poi.LUC <- Softmax.subsampling(X, Y, n.plt, n.ssp, criterion = 'LUC',
#' sampling.method = 'Poisson', estimate.method = 'MSCLE', 
#' constraint = 'baseline')
#' Poi.MSCLE <- Softmax.subsampling(X, Y, n.plt, n.ssp, criterion = 'MSPE',
#' sampling.method = 'Poisson', estimate.method = 'MSCLE', 
#' constraint = 'baseline')
#' softmax.summary(WithRep.MSPE)
#' softmax.summary(Poi.MSPE)
#' softmax.summary(Poi.LUC)
#' softmax.summary(Poi.MSCLE)

Softmax.subsampling <-
  function(X, Y, n.plt, n.ssp,
           criterion = c('OptL', 'OptA', 'MSPE', 'LUC'),
           sampling.method = c('Poisson', 'WithReplacement'),
           estimate.method = c('Weighted', 'MSCLE', 'Uniform'),
           constraint = c('baseline', 'summation'),
           alpha = 0,
           b = 2) {
  model.call <- match.call()
  # family <- match.arg(family)

  # mf <- model.frame(formula, data)
  # Y <- model.response(mf, "any")
  # if (is.character(Y) && length(unique(Y)) == 2) {
  #   levels <- unique(Y)
  #   Y <- as.integer(Y == levels[2])  
  # Assuming levels[2] is the 'success' category
  # }
  # X <- model.matrix(formula, mf)
  # colnames(X)[1] <- "intercept"
  # print(X)
  N <- nrow(X)
  d <- ncol(X)
  K <- length(unique(Y)) - 1
  G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1

  if (estimate.method %in% c("Weighted", "MSCLE")) {
    # pilot step
    plt.estimate.results <- softmax.plt.estimate(X = X, Y = Y, Y.matrix,
                                                 n.plt = n.plt, N, K, d,
                                                 criterion)
    p.plt <- plt.estimate.results$p.plt
    beta.plt.b <- plt.estimate.results$beta.plt
    P1.plt <- plt.estimate.results$P1.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    index.plt <- plt.estimate.results$index.plt
    cov.plt.b <- plt.estimate.results$cov.plt
    Omega.plt <- plt.estimate.results$Omega.plt
    Lambda.plt <- plt.estimate.results$Lambda.plt

    # subsampling step
    ssp.results <- softmax.subsampling(X = X,
                                       Y.matrix = Y.matrix,
                                       G = G,
                                       n.ssp = n.ssp,
                                       N = N, K = K, d = d,
                                       alpha = alpha,
                                       b = b,
                                       criterion = criterion,
                                       estimate.method = estimate.method,
                                       sampling.method = sampling.method,
                                       constraint = constraint,
                                       p.plt = p.plt,
                                       ddL.plt = ddL.plt,
                                       P1.plt = P1.plt,
                                       Omega.plt = Omega.plt,
                                       index.plt = index.plt)
    index.ssp <- ssp.results$index.ssp
    p.ssp <- ssp.results$p.ssp
    offset <- ssp.results$offset

    # subsample estimating step
    ssp.estimate.results <-
      softmax.subsample.estimate(X[index.ssp, ],
                                 Y[index.ssp],
                                 Y.matrix[index.ssp, ],
                                 n.ssp = length(index.ssp),
                                 index.ssp = index.ssp,
                                 p.ssp = p.ssp[index.ssp],
                                 offset = offset,
                                 beta.plt = beta.plt,
                                 sampling.method = sampling.method,
                                 estimate.method = estimate.method,
                                 N=N, K=K, d=d)
    beta.ssp.b <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp.b <- ssp.estimate.results$cov.ssp

    # combining step
    combining.results <- softmax.combining(ddL.plt = ddL.plt,
                                           ddL.ssp = ddL.ssp,
                                           dL.sq.plt = dL.sq.plt,
                                           dL.sq.ssp = dL.sq.ssp,
                                           Lambda.plt = Lambda.plt,
                                           Lambda.ssp = Lambda.ssp,
                                           n.plt = n.plt,
                                           n.ssp = length(index.ssp),
                                           beta.plt = beta.plt.b,
                                           beta.ssp = beta.ssp.b,
                                           X = X, N = N, K = K, d = d)
    beta.cmb.b <- combining.results$beta.cmb
    cov.cmb.b <- combining.results$cov.cmb
    P.cmb <- combining.results$P.cmb

    return(list(
      model.call = model.call,
      beta.plt = matrix(beta.plt.b, nrow = d),
      beta.ssp = matrix(beta.ssp.b, nrow = d),
      beta.cmb = matrix(beta.cmb.b, nrow = d),
      cov.plt = cov.plt.b,
      cov.ssp = cov.ssp.b,
      cov.cmb = cov.cmb.b,
      P.cmb = P.cmb,
      index.plt = index.plt,
      index.ssp = index.ssp,
      N = N,
      subsample.size.expect = n.ssp))
  } else if (estimate.method == "Uniform"){
    n.uni <- n.plt + n.ssp
    if (sampling.method == 'WithReplacement') {
      index.uni <- random.index(N, n.uni)
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni)
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
    } else if (sampling.method == 'Poisson') {
      index.uni <- poisson.index(N, n.uni/N)
      p.uni <- rep(n.uni / N, length(index.uni))
      x.uni <- X[index.uni, ]
      y.uni <- Y[index.uni]
      results <- softmax.coef.estimate(x.uni, y.uni)
      beta.uni.b <- results$beta
      P.uni <- results$P1
      ddL.uni <- softmax_ddL_cpp(X = x.uni, P = P.uni[, -1], p = p.uni, K, d,
                             scale = N)
      dL.sq.uni <- softmax_dL_sq_cpp(X = x.uni, 
                                     Y.matrix = Y.matrix[index.uni, ],
                                     P = P.uni[, -1], p = p.uni, K = K, d = d,
                                     scale=(N^2))
      cov.uni.b <- solve(ddL.uni) %*% dL.sq.uni %*% solve(ddL.uni)
    }
    return(list(index = index.uni,
                beta = matrix(beta.uni.b, nrow = d),
                cov = cov.uni.b,
                P = P.uni))
  }
}
