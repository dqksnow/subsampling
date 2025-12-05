#' Balanced Subsampling Methods for Generalized Linear Models with Rare Features
#'
#' @description
#' This function inherits the weighted objective function (`likelihood = "weighted"`) and the poisson sampling method (`sampling.method = "poisson"`) from `ssp.glm`, while additionally handling rare features in the model. Rare features refer to binary covariates with low prevalence of being one
#' (expressed) and mostly zero. Such features create challenges for subsampling,
#' because it is likely to miss these rare but informative
#' observations. The balanced subsampling method upweights observations that
#' contain expressed rare features, thereby preserving estimation efficiency for
#' the coefficients of rare features. 
#'
#' A quick start guide is provided in the vignette:
#' https://dqksnow.github.io/subsampling/articles/ssp-logit-rF.html.
#'
#'
#' @param formula A model formula object.
#'
#' @param data A data frame containing the variables in the model.
#'
#' @param subset An optional vector specifying a subset of observations to be
#'   used as the full dataset. 
#'
#' @param n.plt The pilot sample size for computing the pilot estimator and
#'   estimating optimal subsampling probabilities.
#'
#' @param n.ssp The expected size of the optimal subsample. For
#'   `sampling.method = "poisson"`, the actual sample size may vary, but the
#'   expected size equals `n.ssp`.
#'
#' @param family A character string naming a family. It can be a character string naming a family function, a family function or the result of a call to a family function. 
#'
#' @param criterion The subsampling criterion. Choices include:
#'   * `"BL-Uni"` (default): probabilities proportional to the balance score.
#'   * `"R-Lopt"`: rareness-aware L-optimality.
#'   * `"Aopt"`: classical A-optimality, minimizing the trace of the asymptotic
#'     covariance matrix.
#'   * `"Lopt"`: classical L-optimality, minimizing a transformed trace of the
#'     asymptotic covariance.
#'   * `"BL-Lopt"`: balance score combined with L-optimality.
#'   * `"Uni"`: uniform sampling.
#'
#' @param sampling.method The sampling method. Currently `"poisson"` is
#'   supported, which avoids drawing repeated observations.
#'
#' @param likelihood The objective function used for estimation. Currently 
#'   `"weighted"` is supported. Each sampled observation is weighted by the
#'   inverse of its subsampling probability.
#'
#' @param rareFeature.index Column indices of binary rare features in the design
#'   matrix, coded as `1` for the rare case. For example,
#'   `c(4, 9)` indicates that columns 4 and 9 of the design matrix contain rare
#'   features.
#'
#' @param balance.plt Logical. Whether to use `"BL-Uni"` to draw the pilot
#'   sample for two-step subsampling methods. Default is `TRUE`. A good pilot
#'   estimator significantly improves the performance of the second-step
#'   subsampling; a poor pilot estimator may cause failure.
#'
#' @param balance.Y Logical. Whether to balance the binary response variable in
#'   logistic regression. If `TRUE`, the pilot sampling probability combines the
#'   balance score with the case-control method, and the negative subsampling will be performed for the second-step subsampling method. That is, automatically include all Y=1 observations to the subsample, while performing subsampling for Y=0 observations.
#'   
#' @param contrasts Optional list specifying how categorical variables are
#'   encoded in the design matrix.
#'
#' @param control A list passed to `glm.control()`. It includes:
#'   * `alpha`: mixture weight between optimal and uniform probabilities, giving  
#'     \eqn{\pi = (1 - \alpha)\pi^{opt} + \alpha \pi^{uni}}.
#'
#' @param ... Additional arguments passed to `svyglm()`.
#'
#'
#' @return
#' An object of class `"ssp.glm.rF"` containing the following fields.
#'
#' \describe{
#'
#' \item{model.call}{The original function call.}
#'
#' \item{coef.plt}{Pilot estimator obtained from the pilot subsample.}
#'
#' \item{coef.ssp}{Estimator obtained from the optimal subsample.}
#'
#' \item{coef.cmb}{Combined estimator using the union of pilot and optimal subsamples.
#'   If `criterion = "BL-Uni"` or `"Uni"`, this equals the pilot and subsample
#'   estimators because only one sampling step is used.}
#'
#' \item{cov.plt}{Estimated covariance matrix of `coef.plt`.}
#'
#' \item{cov.ssp}{Estimated covariance matrix of `coef.ssp`.}
#'
#' \item{cov.cmb}{Estimated covariance matrix of `coef.cmb`.}
#'
#' \item{N}{Number of observations in the full dataset.}
#'
#' \item{subsample.size.expect}{Expected subsample size (`n.ssp`).}
#'
#' \item{subsample.size.actual}{Actual number of subsampled observations.}
#'
#' \item{full.rare.count}{For each rare feature, return the counts of ones in the full dataset.}
#'
#' \item{rare.count.plt}{Vector of length K giving, for each rare feature,
#'   the number of ones in the pilot subsample.}
#'
#' \item{rare.count.ssp}{Same as above, computed for the optimal subsample.}
#'
#' \item{rare.count.cmb}{Same as above, computed for the combined subsample.}
#'
#' \item{index.plt}{Row indices of the pilot subsample within the full dataset.}
#'
#' \item{index.ssp}{Row indices of the optimal subsample within the full
#'   dataset.}
#'
#' \item{index.cmb}{Union of pilot and optimal subsample row indices.}
#'
#' \item{rareFeature.index}{Column indices of rare features in the design
#'   matrix (same as the user input).}
#'
#' \item{comp.time}{Total computation time for computing sampling probabilities,
#'   drawing subsamples, and fitting the subsample estimator.}
#'
#' \item{terms}{The `terms` object for the fitted model.}
#'
#' }
#'
#'
#' @details
#' See the package vignette for more details and examples.
#' @examples
#' # logistic regression
#' set.seed(2)
#' N <- 1e4
#' d_rare <- 3
#' d_cont <- 2
#' p_rare <- c(0.01, 0.02, 0.05)
#' beta0 <- c(0.5, rep(0.5, d_rare), rep(0.5, d_cont)) 
#' corr <- 0.5
#' sigmax  <- matrix(corr, d_cont, d_cont) + diag(1-corr, d_cont)
#' X <- MASS::mvrnorm(N, rep(0, d_cont), sigmax)
#' Z <- do.call(cbind, lapply(seq_along(p_rare), function(i) {
#' rbinom(N, 1, p_rare[i])
#' }))
#' X <- cbind(Z, X)
#' P <- 1 / (1 + exp(-(beta0[1] + X %*% beta0[-1])))
#' Y <- as.integer(rbinom(N, 1, P))
#' colnames(X) <- paste0("X", 1:(d_rare + d_cont))
#' rareFeature.index <- c(1:d_rare)
#' data <- data.frame(Y = Y, X)
#' formula <- Y ~ .
#' n.plt <- 300
#' n.ssp <- 2000
#' BL.Uni.results <- ssp.glm.rF(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'quasibinomial',
#' criterion = 'BL-Uni',
#' sampling.method = 'poisson',
#' likelihood = 'weighted',
#' balance.plt = TRUE,
#' balance.Y = FALSE,
#' rareFeature.index = rareFeature.index
#' )
#' summary(BL.Uni.results)
#' R.Lopt.results <- ssp.glm.rF(formula = formula, 
#' data = data, 
#' n.plt = n.plt,
#' n.ssp = n.ssp,
#' family = 'quasibinomial',
#' criterion = 'R-Lopt',
#' sampling.method = 'poisson',
#' likelihood = 'weighted',
#' balance.plt = TRUE,
#' balance.Y = FALSE,
#' rareFeature.index = rareFeature.index
#' )
#' summary(R.Lopt.results)
#' @export
ssp.glm.rF <- function(formula,
                       data,
                       subset = NULL,
                       n.plt,
                       n.ssp,
                       family = 'binomial',
                       criterion = 'BL-Uni',
                       sampling.method = 'poisson',
                       likelihood = 'weighted',
                       balance.plt = TRUE,
                       balance.Y = FALSE,
                       rareFeature.index = NULL,
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
    name.Y <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(name.Y)) names(Y) <- name.Y
  }
  
  X <- model.matrix(mt, mf, contrasts)
  if (attr(mt, "intercept") == 1) {
    colnames(X)[1] <- "Intercept"
    if(!is.null(rareFeature.index)) {
      # the intercept is added to the left
      rareFeature.index <- rareFeature.index + 1
    }
  }


  criterion <- match.arg(criterion, c('Lopt', 'Aopt', 
                                      'Uni', 'BL-Uni',
                                      'R-Lopt', 'BL-Lopt'))
  sampling.method <- match.arg(sampling.method, c('poisson', 'withReplacement'))
  likelihood <- match.arg(likelihood, c('weighted'))
  control <- do.call("glm.control", control)
  
  ## handle family 
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  canonical_families <- c("binomial", "poisson", "Gamma", "gaussian",
                          "quasibinomial", "quasipoisson")
  if (!family$family %in% canonical_families) {
    stop(sprintf("Unsupported GLM family '%s'. Only (%s) are supported",
                 family$family,
                 paste(canonical_families, collapse = ", ")))
  }
  
  
  N <- nrow(X)
  d <- ncol(X)
  subsample.size.expect <- n.ssp
  
  
  ## compute balancing score
  bl.information <- balance_score(rareFeature.index, X)
  bl <- bl.information$bl
  DN <- bl.information$DN
  full.rare.count <- bl.information$rare.count
  Y.count <- bl.information$Y.count
  

  if(balance.plt && is.null(rareFeature.index)) {
    cat('Caution: argument rareFeature.index is NULL. 
    The pilot sample will be drawn by uniform sampling probability.',
        '\n')
  }
  
  inputs <- list(X = X, Y = Y, N = N, d = d,
                 n.plt = n.plt, n.ssp = n.ssp,
                 criterion = criterion, sampling.method = sampling.method,
                 likelihood = likelihood, family = family,
                 bl = bl, DN = DN,
                 balance.plt = balance.plt,
                 balance.Y = balance.Y,
                 rareFeature.index = rareFeature.index, 
                 control = control
                 )
  
  if (criterion %in% c('Lopt', 'Aopt', 'R-Lopt', 'BL-Lopt')) {
    t_start <- proc.time()[3]
    ## pilot step
    plt.estimate.results <- rF.pilot.estimate(inputs, 
                                           ...
                                           )
    pi.plt <- plt.estimate.results$pi.plt
    pi.plt.N <- plt.estimate.results$pi.plt.N
    varphi.plt = plt.estimate.results$varphi.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    dL.sq.plt <- plt.estimate.results$dL.sq.plt
    Lambda.plt <- plt.estimate.results$Lambda.plt
    linear.predictor <- plt.estimate.results$linear.predictor # dimension: N
    index.plt <- plt.estimate.results$index.plt
    cov.plt <- plt.estimate.results$cov.plt
    n.plt <- length(index.plt)

    ## subsampling step
    ssp.results <- rF.subsampling(inputs,
                               pi.plt = pi.plt,
                               varphi.plt = varphi.plt,
                               ddL.plt = ddL.plt,
                               linear.predictor = linear.predictor,
                               index.plt = index.plt
                               )
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    varphi.ssp <- ssp.results$varphi.ssp
    pi.ssp.N <- ssp.results$pi.ssp.N
    offset <- ssp.results$offset
    n.ssp <- length(index.ssp)
    
    ## subsample estimating step
    ssp.estimate.results <- rF.subsample.estimate(inputs,
                                               w.ssp = w.ssp,
                                               varphi.ssp = varphi.ssp,
                                               offset = offset,
                                               beta.plt = beta.plt,
                                               index.ssp = index.ssp,
                                               ...
                                               )
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    dL.sq.ssp <- ssp.estimate.results$dL.sq.ssp
    Lambda.ssp <- ssp.estimate.results$Lambda.ssp
    cov.ssp <- ssp.estimate.results$cov.ssp
    
    
    ## combining step
    combining.results <- rF.combining(inputs,
                                   index.plt = index.plt, # dim: n.plt
                                   index.ssp = index.ssp, # dim: n.ssp
                                   pi.plt = pi.plt.N, # dim: N
                                   pi.ssp = pi.ssp.N # dim: N
                                   )
    index.cmb <- combining.results$index.cmb
    beta.cmb <- combining.results$beta.cmb
    cov.cmb <- combining.results$cov.cmb
    
    
    ## prepare for displaying results 
    names(beta.cmb) <- names(beta.ssp) <- names(beta.plt) <- colnames(X)
    
    rare.counts.results.plt <- rare_counts(X, 
                                           index.plt, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.ssp <- rare_counts(X, 
                                           index.ssp, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.cmb <- rare_counts(X, 
                                           index.cmb, 
                                           rareFeature.index,
                                           Y)
    
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef.cmb = beta.cmb,
                    cov.plt = cov.plt,
                    cov.ssp = cov.ssp,
                    cov.cmb = cov.cmb,

                    N = N,
                    subsample.size.expect = subsample.size.expect,
                    subsample.size.actual = length(index.ssp),
                    full.rare.count = full.rare.count,
                    
                    rare.count.plt = rare.counts.results.plt$rare.count,
                    rare.count.ssp = rare.counts.results.ssp$rare.count,
                    rare.count.cmb = rare.counts.results.cmb$rare.count,

                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    index.cmb = index.cmb,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    terms = mt
                    )
    class(results) <- c("ssp.glm.rF", "list")
    return(results)
  } else if (criterion %in% c("Uni", "BL-Uni")){
    if (criterion == "Uni") inputs$balance.plt = FALSE
    # n.uni <- n.plt + n.ssp # sum up the expected size of pilot and subsample 
    inputs$n.uni <- n.ssp
    t_start <- proc.time()[3]
    uni.estimate.results <- rF.uniform.estimate(inputs, ...)
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    pi.uni <- uni.estimate.results$pi.uni
    varphi.uni = uni.estimate.results$varphi.uni
    beta.uni <- uni.estimate.results$beta.uni
    ddL.uni <- uni.estimate.results$ddL.uni
    dL.sq.uni <- uni.estimate.results$dL.sq.uni
    Lambda.uni <- uni.estimate.results$Lambda.uni
    linear.predictor <- uni.estimate.results$linear.predictor
    index.uni <- uni.estimate.results$index.uni
    cov.uni <- uni.estimate.results$cov.uni


    rare.counts.results <- rare_counts(X, index.uni, rareFeature.index, Y)
    
    
    results <- list(model.call = model.call,
                    coef.plt = beta.uni,
                    coef.ssp = beta.uni,
                    coef.cmb = beta.uni,
                    cov.plt = cov.uni,
                    cov.ssp = cov.uni,
                    cov.cmb = cov.uni,
                    N = N,
                    subsample.size.expect = inputs$n.uni,
                    subsample.size.actual = length(index.uni),
                    full.rare.count = full.rare.count,
                    
                    rare.count.plt = rare.counts.results$rare.count,
                    rare.count.ssp = rare.counts.results$rare.count,
                    rare.count.cmb = rare.counts.results$rare.count,

                    index.plt = index.uni,
                    index.ssp = index.uni,
                    index.cmb = index.uni,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    terms = mt
    )
    class(results) <- c("ssp.glm.rF", "list")
    return(results)
  }
}
