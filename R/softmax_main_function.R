#' Optimal Subsampling Method for Softmax Regression Models
#'
#' @details details TBD
#' @param formula An object of class "formula" which describes the model to be fitted.
#' @param data A data frame containing the variables in the model. Usually it 
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' @param n.ssp The expected optimal subsample size (the second-step subsample.
#' @param criterion The criterion of optimal subsampling probabilities
#' @param sampling.method The sampling method for drawing the optimal subsample
#' @param likelihood The type of the maximum likelihood function used to 
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
#' N <- 1e4
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
#' data <- as.data.frame(cbind(Y, X))
#' formula <- Y ~ .
#' WithRep.MSPE <- ssp.softmax(formula, data, n.plt, n.ssp, criterion = 'MSPE', 
#' sampling.method = 'withReplacement', likelihood = 'weighted',
#' constraint = 'baseline')
#' summary(WithRep.MSPE)

ssp.softmax <-
  function(formula, data, n.plt, n.ssp,
           criterion = c('optL', 'optA', 'MSPE', 'LUC', 'uniform'),
           sampling.method = c('poisson', 'withReplacement'),
           likelihood = c('weighted', 'MSCLE'),
           constraint = c('baseline', 'summation'),
           alpha = 0,
           b = 2) {
    
    # model.call <- match.call()
    # m <- match.call(expand.dots = FALSE)
    # m[[1L]] <- quote(stats::model.frame)
    # m <- eval.parent(m)
    # Terms <- attr(m, "terms")
    # print(Terms)
    # X <- model.matrix(Terms, m, contrasts)
    # print(head(X))
    # cons <- attr(X, "contrasts")
    # print(cons)
    # Y <- model.response(m)

  model.call <- match.call()
  mf <- model.frame(formula, data)
  Y <- model.response(mf, "any")
  X <- model.matrix(formula, mf)
  colnames(X)[1] <- "intercept"
  
  criterion <- match.arg(criterion)
  sampling.method <- match.arg(sampling.method)
  likelihood <- match.arg(likelihood)
  constraint <- match.arg(constraint)
  
  dimension <- dim(X)
  N <- dimension[1]
  d <- dimension[2]
  K <- length(unique(Y)) - 1
  G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d) 
  ## G: transformation matrix
  Y.matrix <- matrix(0, nrow = N, ncol = K)
  Y.matrix[cbind(seq_along(Y), Y)] <- 1
  
  ## create a list to store variables
  inputs <- list(X = X, Y = Y, Y.matrix = Y.matrix,
                 N = N, d = d, K = K, G = G, 
                 n.plt = n.plt, n.ssp = n.ssp, 
                 criterion = criterion, sampling.method = sampling.method,
                 likelihood = likelihood, constraint = constraint
                 )
  
  if (criterion %in% c('optL', 'optA', 'MSPE', 'LUC')) {
    ## pilot step
    # plt.estimate.results <- softmax.plt.estimate(X = X, Y = Y, Y.matrix,
    #                                              n.plt = n.plt, N, K, d,
    #                                              criterion)
    plt.estimate.results <- softmax.plt.estimate(inputs)
    
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
    ssp.results <- softmax.subsampling(X = X,
                                       Y.matrix = Y.matrix,
                                       G = G,
                                       n.ssp = n.ssp,
                                       N = N, K = K, d = d,
                                       alpha = alpha,
                                       b = b,
                                       criterion = criterion,
                                       likelihood = likelihood,
                                       sampling.method = sampling.method,
                                       constraint = constraint,
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
      softmax.subsample.estimate(X[index.ssp, ],
                                 Y[index.ssp],
                                 Y.matrix[index.ssp, ],
                                 n.ssp = length(index.ssp),
                                 index.ssp = index.ssp,
                                 p.ssp = p.ssp[index.ssp],
                                 offsets = offsets,
                                 beta.plt = beta.plt,
                                 sampling.method = sampling.method,
                                 likelihood = likelihood,
                                 N=N, K=K, d=d)
    beta.ssp.b <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp.b <- ssp.estimate.results$cov.ssp

    ## combining step
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

    results <- list(model.call = model.call,
                    beta.plt = matrix(beta.plt.b, nrow = d),
                    beta.ssp = matrix(beta.ssp.b, nrow = d),
                    beta = matrix(beta.cmb.b, nrow = d),
                    cov.plt = cov.plt.b,
                    cov.ssp = cov.ssp.b,
                    cov = cov.cmb.b,
                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    N = N,
                    subsample.size.expect = n.ssp
                    )
    class(results) <- c("ssp.softmax", "list")
    return(results)
  } else if (criterion == "uniform"){
    n.uni <- n.plt + n.ssp
    if (sampling.method == 'withReplacement') {
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
    } else if (sampling.method == 'poisson') {
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
                                     Y_matrix = Y.matrix[index.uni, ],
                                     P = P.uni[, -1], p = p.uni, K = K, d = d,
                                     scale=(N^2))
      cov.uni.b <- solve(ddL.uni) %*% dL.sq.uni %*% solve(ddL.uni)
    }
    
    results <- list(model.call = model.call,
                    index.ssp = index.uni,
                    beta = matrix(beta.uni.b, nrow = d),
                    cov = cov.uni.b,
                    # P = P.uni,
                    N = N,
                    subsample.size.expect = n.uni
                    )
    class(results) <- c("ssp.softmax", "list")
    return(results)
  }
}
