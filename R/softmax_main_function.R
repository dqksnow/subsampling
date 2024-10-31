#' Optimal Subsampling Method for Softmax (multinomial logistic) Regression Model
#' @description
#' Draw subsample from full dataset and fit softmax(multinomial logistic) regression model on the subsample. Refer to [vignette](https://dqksnow.github.io/subsampling/articles/ssp-softmax.html) for a quick start.
#' 
#' @param formula A model formula object of class "formula" that describes the model to be fitted.
#' @param data A data frame containing the variables in the model. Denote \eqn{N} as the number of observations in `data`.
#' @param subset An optional vector specifying a subset of observations from `data` to use for the analysis. This subset will be viewed as the full data.
#' @param n.plt The pilot subsample size (first-step subsample size).
#' This subsample is used to compute the pilot estimator and estimate the optimal subsampling probabilities.
#' @param n.ssp The expected size of the optimal subsample (second-step subsample). For `sampling.method = 'withReplacement'`, The exact subsample size is `n.ssp`. For `sampling.method = 'poisson'`, `n.ssp` is the expected subsample size. 
#' @param criterion The criterion of optimal subsampling probabilities.
#' Choices include \code{optA}, \code{optL}, \code{MSPE}(default), \code{LUC} and \code{uniform}.
#' 
#' - `MSPE` Minimizes the mean squared prediction error between subsample estimator and full data estimator.
#' 
#' - `optA` Minimizes the trace of the asymptotic covariance matrix of the subsample estimator. 
#' 
#' - `optL` Minimizes the trace of a transformation of the asymptotic covariance matrix, which reduces computational costs than `optA`.
#' 
#' - `LUC` Local uncertainty sampling method, serving as a baseline subsampling strategy. See Wang and Kim (2022).
#' 
#' - `uniform` Assigns equal subsampling probability
#' \eqn{\frac{1}{N}} to each observation, serving as a baseline subsampling strategy.
#' 
#' @param sampling.method The sampling method to use. 
#' Choices include \code{withReplacement} and \code{poisson}(default). `withReplacement` draws exactly `n.ssp`
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
#' @param likelihood A bias-correction likelihood function is required for subsample since unequal subsampling probabilities introduce bias. Choices include 
#'  \code{weighted} and \code{MSCLE}(default). 
#'  
#' - `weighted` Applies a weighted likelihood function where each observation is weighted by the inverse of its subsampling probability.
#' 
#' - `MSCLE` It uses a conditional likelihood, where each element of the likelihood represents the density of \eqn{Y_i} given that this observation was drawn.
#' 
#' @param constraint The constraint for identifiability of softmax model. Choices include 
#'  \code{baseline} and \code{summation}(default). The baseline constraint assumes the coefficient for the baseline category are \eqn{0}. Without loss of generality, we set the category \eqn{Y=0} as the baseline category so that \eqn{\boldsymbol{\beta}_0=0}. The summation constraint \eqn{\sum_{k=0}^{K} \boldsymbol{\beta}_k} is also used in the subsampling method for the purpose of calculating subsampling probability. These two constraints lead to different interpretation of coefficients but are equal for computing \eqn{P(Y_{i,k} = 1 \mid \mathbf{x}_i)}. The estimation of coefficients returned by `ssp.softmax()` is under baseline constraint.
#'  
#' @param contrasts An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.
#' @param control A list of parameters for controlling the sampling process. There are two tuning parameters `alpha` and `b`. Default is \code{list(alpha=0, b=2)}.
#' 
#' - `alpha` \eqn{\in [0,1]} is the mixture weight of the user-assigned subsampling
#' probability and uniform subsampling probability. The actual subsample
#' probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. This protects the estimator from extreme small
#' subsampling probability. The default value is 0.
#' 
#' - `b` is a positive number which is used to constaint the poisson subsampling probability. `b` close to 0 results in subsampling probabilities closer to uniform probability \eqn{\frac{1}{N}}. `b=2` is the default value. See relevant references for further details.
#' 
#' @param ... A list of parameters which will be passed to \code{nnet::multinom()}. 
#'
#' @return
#' ssp.softmax returns an object of class "ssp.softmax" containing the following components (some are optional):
#' \describe{
#'   \item{model.call}{The original function call.}
#'   \item{coef.plt}{The pilot estimator. See Details for more information.}
#'   \item{coef.ssp}{The estimator obtained from the optimal subsample.}
#'   \item{coef}{The weighted linear combination of `coef.plt` and `coef.ssp`, under baseline constraint. The combination weights depend on the relative size of `n.plt` and `n.ssp` and the estimated covariance matrices of `coef.plt` and `coef.ssp.` We blend the pilot subsample information into optimal subsample estimator since the pilot subsample has already been drawn. The coefficients and standard errors reported by summary are `coef` and the square root of `diag(cov)`.}
#'   \item{coef.plt.sum}{The pilot estimator under summation constrraint. `coef.plt.sum = G %*% as.vector(coef.plt)`.}
#'   \item{coef.ssp.sum}{The estimator obtained from the optimal subsample under summation constrraint. `coef.ssp.sum = G %*% as.vector(coef.ssp)`.}
#'   \item{coef.sum}{The weighted linear combination of `coef.plt` and `coef.ssp`, under summation constrraint. `coef.sum = G %*% as.vector(coef)`.}
#'   \item{cov.plt}{The covariance matrix of \code{coef.plt}.}
#'   \item{cov.ssp}{The covariance matrix of \code{coef.ssp}.}
#'   \item{cov}{The covariance matrix of \code{coef.cmb}.}
#'   \item{cov.plt.sum}{The covariance matrix of \code{coef.plt.sum}.}
#'   \item{cov.ssp.sum}{The covariance matrix of \code{coef.ssp.sum}.}
#'   \item{cov.sum}{The covariance matrix of \code{coef.sum}.}
#'   \item{index.plt}{Row indices of pilot subsample in the full dataset.}
#'   \item{index.ssp}{Row indices of of optimal subsample in the full dataset.}
#'   \item{N}{The number of observations in the full dataset.}
#'   \item{subsample.size.expect}{The expected subsample size.}
#'   \item{terms}{The terms object for the fitted model.}
#' }
#' 
#' @details
#' 
#' A pilot estimator for the unknown parameter  \eqn{\beta} is required because MSPE, optA and
#' optL subsampling probabilities depend on \eqn{\beta}. There is no "free lunch" when determining optimal subsampling probabilities. For softmax regression, this
#' is achieved by drawing a size `n.plt` subsample with replacement from full
#' dataset with uniform sampling probability.
#' 
#'
#' @references
#' Yao, Y., & Wang, H. (2019). Optimal subsampling for softmax regression. \emph{Statistical Papers}, \strong{60}, 585-599.
#' 
#' Han, L., Tan, K. M., Yang, T., & Zhang, T. (2020). Local uncertainty sampling for large-scale multiclass logistic regression. \emph{Annals of Statistics}, \strong{48}(3), 1770-1788.
#' 
#' Wang, H., & Kim, J. K. (2022). Maximum sampled conditional likelihood for informative subsampling. \emph{Journal of machine learning research}, \strong{23}(332), 1-50.
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
#' prob <- exp(X %*% beta.true.summation)
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
    beta.plt.b <- plt.estimate.results$beta.plt # under baseline constraint
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
    beta.ssp.b <- ssp.estimate.results$beta.ssp # under baseline constraint
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
    beta.cmb.b <- combining.results$beta.cmb # under baseline constraint
    cov.cmb.b <- combining.results$cov.cmb
    P.cmb <- combining.results$P.cmb

    beta.plt.sum <- G %*% as.vector(beta.plt.b)
    beta.ssp.sum <- G %*% as.vector(beta.ssp.b)
    beta.cmb.sum <- G %*% as.vector(beta.cmb.b)
    cov.plt.sum <- G %*% cov.plt.b %*% t(G)
    cov.ssp.sum <- G %*% cov.ssp.b %*% t(G)
    cov.cmb.sum <- G %*% cov.cmb.b %*% t(G)
    
    beta.plt <- matrix(beta.plt.b, nrow = d)
    beta.ssp = matrix(beta.ssp.b, nrow = d)
    beta <- matrix(beta.cmb.b, nrow = d)
    rownames(beta.plt) <- rownames(beta.ssp) <- rownames(beta) <- colnames(X)
    
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef = beta,
                    
                    coef.plt.sum = beta.plt.sum,
                    coef.ssp.sum = beta.ssp.sum,
                    coef.sum = beta.cmb.sum,
                    
                    cov.plt = cov.plt.b,
                    cov.ssp = cov.ssp.b,
                    cov = cov.cmb.b,
                    
                    cov.plt.sum = cov.plt.sum,
                    cov.sum = cov.cmb.sum,
                    cov.ssp.sum = cov.ssp.sum,
                    
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
      P.uni <- results$P1
      ddL.uni <- softmax_ddL_cpp(X = x.uni, P = P.uni[, -1], p = rep(1, n.uni),
                                 K, d, scale = n.uni)
      dL.sq.uni <- softmax_dL_sq_cpp(X = x.uni, 
                                     Y_matrix = Y.matrix[index.uni, ],
                                     P = P.uni[, -1], p = rep(1, n.uni), K = K, 
                                     d = d, scale = n.uni)
      c <- n.uni / N
      cov.uni.b <- solve(ddL.uni) %*% (dL.sq.uni * (1+c)) %*% 
        solve(ddL.uni) / n.uni
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
    
    beta.sum <- G %*% as.vector(beta.uni.b)
    
    beta <- matrix(beta.uni.b, nrow = d)
    
    cov.sum <- G %*% cov.uni.b %*% t(G)
    
    rownames(beta) <- colnames(X)
    results <- list(model.call = model.call,
                    index.ssp = index.uni,
                    coef = beta,
                    cov = cov.uni.b,
                    coef.sum = beta.sum,
                    cov.sum = cov.sum,
                    N = N,
                    subsample.size.expect = n.uni,
                    terms = mt
                    )
    class(results) <- c("ssp.softmax", "list")
    return(results)
  }
}
