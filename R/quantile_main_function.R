#' Optimal Subsampling Methods for Quantile Models
#' @export
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
#' @param tau The quantile.
#' @param n.plt The pilot subsample size (the first-step subsample size).
#' These samples will be used to estimate the pilot estimator.
#' @param n.ssp The expected optimal subsample size (the second-step subsample
#' size).
#' @param B TBD
#' @param boot TBD
#' @param criterion The criterion of optimal subsampling probabilities.
#' @param sampling.method The sampling method for drawing the optimal subsample.
#' @param likelihood The type of the maximum likelihood function used to
#' calculate the optimal subsampling estimator.
#' @param contrasts an optional list. It specifies how categorical variables are represented in the design matrix. For example, contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')
#' @param control a list of parameters for controlling the computing process. 
#' @param ... an optional list of parameters which will be passed to 
#' quantreg::rq().
#' 
#' @return
#' \describe{
#'   \item{model.call}{the model}
#'   \item{beta.plt}{pilot estimator}
#'   \item{beta.ssp}{optimal subsample estimator}
#'   \item{est.cov.ssp}{covariance matrix of \code{beta.ssp}}
#'   \item{index.plt}{index of pilot subsample}
#'   \item{index.ssp}{index of optimal subsample}
#' }
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
  
  ## mf represents the model frame that contains all the necessary variables 
  ## for model fitting, including the response and predictors.
  ## Initially, mf is the match.all() object containing all the arguments 
  ## passed to the main function
  ## After match() matching the relevant arguments, mf is subsetted to include
  ## only those relevant arguments.
  ## It is then converted to a model frame using stats::model.frame
  
  model.call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"),
             names(mf),
             0L)
  mf <- mf[c(1L, m)] # mf[[1]]) contains the original function call
  mf$drop.unused.levels <- TRUE 
  # Sometimes, not all levels are present in the subset of data being analyzed.
  # drop.unused.levels drops any unused levels of factor variables.
  
  ## when mf is evaluated, it will execute stats::model.frame(...)
  ## instead of mf[[1L]].
  mf[[1L]] <- quote(stats::model.frame)
  ## mf is just an expression before running eval(). eval() evaluate the
  ## expression within a specific environment. parent.frame() is the 
  ## environment of the caller function.
  mf <- eval(mf, parent.frame())
  ## mt extracts the terms object from the model frame mf. This object is
  ## essential for understanding the structure of the model, such as which
  ## variables are used and how they are transformed.
  mt <- attr(mf, "terms") # allow model.frame to have updated it
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
    plt.results <- quantile.plt.estimation(inputs, ...)
    beta.plt <- plt.results$beta.plt
    Ie.full <- plt.results$Ie.full
    index.plt <- plt.results$index.plt

    ## subsampling and boot step
    ssp.results <- quantile.ssp.estimation(inputs,
                                           Ie.full = Ie.full,
                                           index.plt = index.plt,
                                           ...
                                           )
    Betas.ssp <- ssp.results$Betas.ssp
    beta.ssp <- ssp.results$beta.ssp
    est.cov.ssp <- ssp.results$est.cov.ssp
    index.ssp <- ssp.results$index.ssp
    
    names(beta.ssp) <- names(beta.plt) <- colnames(X)

    results <- list(model.call = model.call,
                    beta.plt = beta.plt,
                    beta = beta.ssp,
                    est.cov = est.cov.ssp,
                    index.plt = index.plt,
                    index = index.ssp,
                    N = N,
                    subsample.size.expect = c(n.ssp, B),
                    terms = mt
                    )
    class(results) <- c("ssp.quantreg", "list")
    return(results)
  } else if (criterion == "uniform"){
    ## subsampling and boot step
    uni.results <- quantile.ssp.estimation(inputs, ...)
    Betas.uni <- uni.results$Betas.ssp
    beta.uni <- uni.results$beta.ssp
    est.cov.uni <- uni.results$est.cov.ssp
    index.uni <- uni.results$index.ssp
    names(beta.uni) <- colnames(X)
    results <- list(model.call = model.call,
                    beta.plt = NA,
                    beta = beta.uni,
                    est.cov = est.cov.uni,
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