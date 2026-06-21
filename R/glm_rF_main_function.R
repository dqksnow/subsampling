#' Balanced Subsampling Methods for Generalized Linear Models with Rare Features
#'
#' @description
#' Rare features are binary covariates with low prevalence of being one. Because
#' uniform or classical optimal subsampling can miss expressed rare-feature
#' observations or produce unstable pilot estimates, this function uses
#' rarity-aware sampling probabilities to preserve information for estimating
#' rare-feature coefficients.
#'
#' The function extends `ssp.glm` by supporting 
#' rarity-aware designs, optional response balancing for binary outcomes,
#' weighted or unweighted pilot objectives, and a combined estimator based on the
#' union of the pilot and second-step subsamples.
#'
#' @param formula A model formula object.
#'
#' @param data A data frame containing the variables in the model.
#'
#' @param subset An optional vector specifying a subset of observations to be
#'   used as the full dataset.
#'
#' @param n.plt The expected pilot sample size for two-step methods. For
#'   one-step methods (`criterion = "Uni"` or `"BL-Uni"`), the expected sample
#'   size is `n.plt + n.ssp`.
#'
#' @param n.ssp The expected second-step subsample size. For Poisson
#'   subsampling, the actual sample size may vary.
#'
#' @param family A character string naming a family, a family function, or the
#'   result of a call to a family function. Supported families include
#'   `"binomial"`, `"quasibinomial"`, `"poisson"`, `"quasipoisson"`, `"gaussian"`,
#'   and `"Gamma"`.
#'
#' @param criterion The subsampling criterion. Choices include:
#'   * `"BL-Uni"` (default): probabilities proportional to the balance score.
#'   * `"Uni"`: uniform Poisson subsampling.
#'   * `"Lopt"`: classical L-optimality.
#'   * `"Aopt"`: classical A-optimality.
#'   * `"R-Lopt"`: rareness-aware L-optimality.
#'   * `"BL-Lopt"`: balance score combined with L-optimality.
#'
#' @param sampling.method The sampling method. Currently only `"poisson"` is
#'   supported.
#'
#' @param objective.weight.plt Objective weighting for the pilot fit. Use
#'   `"weighted"` for inverse-probability weighting or `"unweighted"` for an
#'   unweighted pilot objective. Unweighted pilot fitting is not allowed when the
#'   pilot sampling probability depends on the response.
#'
#' @param objective.weight Objective weighting for the one-step or second-step
#'   fit. Two-step methods currently require `"weighted"`.
#'
#' @param control A list passed to `glm.control()`. Supported entries include:
#'   * `alpha`: mixture weight between optimal and uniform probabilities.
#'   * `b`: pilot-based truncation tuning parameter.
#'   * `poi.method`: `"exact"` or `"estimated"` for Poisson probability
#'     normalization.
#'
#' @param contrasts Optional list specifying how categorical variables are
#'   encoded in the design matrix.
#'
#' @param balance.X.plt Logical. Whether to use balance-score sampling for the
#'   pilot sample in two-step methods.
#'
#' @param balance.Y.plt Logical. Whether to balance the binary response in the
#'   pilot sample. Ignored for non-binary response families.
#'
#' @param balance.Y.ssp Logical. For one-step `"Uni"` and `"BL-Uni"` methods,
#'   whether to allocate the expected sample size across `Y = 0` and `Y = 1`
#'   groups in a case-control style. Ignored for two-step optimality criteria and
#'   non-binary response families.
#'
#' @param balance.Y.all Logical. Whether to include all `Y = 1` observations and
#'   subsample from `Y = 0`. Ignored for non-binary response families.
#'
#' @param record.stage.time Logical. Whether to store timing for major internal
#'   stages in the returned object.
#'
#' @param rareFeature.index Rare-feature columns. Numeric values follow the same
#'   convention as the original data/model variables: if the model contains an
#'   intercept, the function internally shifts the indices to account for the
#'   intercept column in the design matrix. Character values are matched to
#'   design-matrix column names. If `NULL`, rare binary features are detected
#'   automatically using `rareThreshold`.
#'
#' @param rareThreshold Prevalence threshold used to automatically identify rare
#'   binary features, and to warn when user-supplied rare features have
#'   prevalence at or above the threshold.
#'
#' @param na.action Currently accepted for interface compatibility.
#'
#' @param ... Additional arguments passed to `glm.fit()` or `lm.wfit()`.
#'
#' @return
#' An object of class `"ssp.glm.rF"` containing fitted coefficients, covariance
#' estimates, selected row indices, rare-feature counts, response-composition
#' summaries, and optional stage timings.
#'
#' @details
#' Two-step criteria (`"Lopt"`, `"Aopt"`, `"R-Lopt"`, and `"BL-Lopt"`) draw a
#' pilot sample, compute second-step Poisson probabilities, fit the second-step
#' weighted GLM, and then refit on the union of the pilot and second-step
#' samples. One-step criteria (`"Uni"` and `"BL-Uni"`) draw a single Poisson
#' subsample with expected size `n.plt + n.ssp`.
#'
#' @examples
#' set.seed(2)
#' N <- 1000
#' Z1 <- rbinom(N, 1, 0.04)
#' Z2 <- rbinom(N, 1, 0.07)
#' X1 <- rnorm(N)
#' X2 <- rnorm(N)
#' eta <- 0.5 + 0.5 * Z1 + 0.5 * Z2 + 0.5 * X1 + 0.5 * X2
#' Y <- rbinom(N, 1, plogis(eta))
#' data <- data.frame(Y, Z1, Z2, X1, X2)
#'
#' fit_bl <- ssp.glm.rF(
#'   Y ~ .,
#'   data = data,
#'   n.plt = 100,
#'   n.ssp = 150,
#'   family = "quasibinomial",
#'   criterion = "BL-Uni",
#'   rareFeature.index = 1:2
#' )
#' summary(fit_bl)
#'
#' fit_rl <- ssp.glm.rF(
#'   Y ~ .,
#'   data = data,
#'   n.plt = 100,
#'   n.ssp = 150,
#'   family = "quasibinomial",
#'   criterion = "R-Lopt",
#'   balance.X.plt = TRUE,
#'   rareFeature.index = c("Z1", "Z2")
#' )
#' summary(fit_rl)
#'
#' @export
ssp.glm.rF <- function(formula,
                    data,
                    subset = NULL,
                    n.plt,
                    n.ssp,
                    family = 'binomial',
                    criterion = 'BL-Uni',
                    sampling.method = 'poisson',
                    objective.weight.plt = 'weighted',
                    objective.weight = 'weighted',
                    control = list(...),
                    contrasts = NULL,
                    balance.X.plt = FALSE,
                    balance.Y.plt = FALSE,
                    balance.Y.ssp = FALSE,
                    balance.Y.all = FALSE,
                    record.stage.time = FALSE,
                    rareFeature.index = NULL,
                    rareThreshold = 0.09,
                    na.action = getOption("na.action"),
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
    if (is.numeric(rareFeature.index)) {
      rareFeature.index <- rareFeature.index + 1L
    }
  }
  
  if (is.character(family)) {
    family_lookup <- c(
      binomial = "binomial",
      poisson = "poisson",
      Gamma = "Gamma",
      gaussian = "gaussian",
      Gaussian = "gaussian",
      quasibinomial = "quasibinomial",
      quasipoisson = "quasipoisson"
    )
    family <- match.arg(family, names(family_lookup))
    family <- family_lookup[[family]]
  }
  criterion <- match.arg(criterion, c('Lopt', 'Aopt', 
                                      'Uni', 'BL-Uni',
                                      'R-Lopt', 'BL-Lopt'))
  sampling.method <- match.arg(sampling.method, c('poisson'))
  objective.weight.plt <- match.arg(
    objective.weight.plt,
    c('weighted', 'unweighted')
  )
  objective.weight <- match.arg(
    objective.weight,
    c('weighted', 'unweighted')
  )
  control <- do.call("rF.glm.control", control)
  
  ## family 
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (!(family$family %in% c("binomial", "quasibinomial"))) {
    requested_y_balancing <- balance.Y.plt || balance.Y.ssp || balance.Y.all
    balance.Y.plt <- FALSE
    balance.Y.ssp <- FALSE
    balance.Y.all <- FALSE
    if (requested_y_balancing) {
      warning("only consider balancing Y when Y is binary")
    }
  }
  
  if (criterion %in% c("Lopt", "Aopt", "R-Lopt", "BL-Lopt") && balance.Y.ssp) {
    warning(
      paste0(
        "'balance.Y.ssp' is currently only used for 'Uni' and 'BL-Uni'; ",
        "ignoring it for criterion '", criterion, "'."
      )
    )
    balance.Y.ssp <- FALSE
  }

  is.two.stage <- criterion %in% c("Lopt", "Aopt", "R-Lopt", "BL-Lopt")
  is.one.stage <- criterion %in% c("Uni", "BL-Uni")

  if (is.one.stage &&
      objective.weight == "unweighted" &&
      (balance.Y.ssp || balance.Y.all)) {
    stop(
      paste0(
        "'objective.weight = \"unweighted\"' is only valid for one-step methods ",
        "when the sampling probability does not depend on Y."
      )
    )
  }

  if (is.two.stage &&
      objective.weight.plt == "unweighted" &&
      balance.Y.plt) {
    stop(
      paste0(
        "'objective.weight.plt = \"unweighted\"' is not allowed when ",
        "'balance.Y.plt = TRUE' because the pilot sampling depends on Y."
      )
    )
  }

  if (is.two.stage && objective.weight != "weighted") {
    stop(
      paste0(
        "'objective.weight' must be 'weighted' for two-step methods ",
        "('Lopt', 'Aopt', 'R-Lopt', 'BL-Lopt')."
      )
    )
  }
  
  N <- nrow(X)
  d <- ncol(X)
  subsample.size.expect <- n.ssp
  full.Y.sum <- sum(Y, na.rm = TRUE)
  full.Y.mean <- mean(Y, na.rm = TRUE)
  stage.time <- numeric(0)
  
  ## compute balancing score
  stage_start <- proc.time()[3]
  bl.information <- rF.balance_score(X, Y, family, 
                                  rareFeature.index = rareFeature.index,
                                  threshold = rareThreshold)
  stage.time["balance_score"] <- proc.time()[3] - stage_start
  bl <- bl.information$bl
  DN <- bl.information$DN
  full.rare.count <- bl.information$rare.count
  Y.count <- bl.information$Y.count
  rareFeature.index <- bl.information$rareFeature.index
  Y_subgroup_prev <- bl.information$Y_subgroup_prev
  if (attr(mt, "intercept") == 1) {
    colnames(X)[1] <- "Intercept"
  }

  inputs <- list(X = X, Y = Y, N = N, d = d,
                 Y.count = Y.count,
                 n.plt = n.plt, n.ssp = n.ssp,
                 criterion = criterion, sampling.method = sampling.method,
                 objective.weight.plt = objective.weight.plt,
                 objective.weight = objective.weight,
                 family = family,
                 bl = bl, DN = DN,
                 balance.X.plt = balance.X.plt,
                 balance.Y.plt = balance.Y.plt,
                 balance.Y.ssp = balance.Y.ssp,
                 balance.Y.all = balance.Y.all,
                 rareFeature.index = rareFeature.index, 
                 control = control
                 )
  
  if (criterion %in% c('Lopt', 'Aopt', 'R-Lopt', 'BL-Lopt')) {
    t_start <- proc.time()[3]
    ## pilot step
    stage_start <- proc.time()[3]
    plt.estimate.results <- rF.pilot.estimate(inputs, 
                                           ...
                                           )
    stage.time["pilot.estimate"] <- proc.time()[3] - stage_start
    pi.plt <- plt.estimate.results$pi.plt
    pi.plt.N <- plt.estimate.results$pi.plt.N
    varphi.plt = plt.estimate.results$varphi.plt
    beta.plt <- plt.estimate.results$beta.plt
    rF.ddL.plt <- plt.estimate.results$rF.ddL.plt
    rF.ddL.plt.design <- plt.estimate.results$rF.ddL.plt.design
    if (is.null(rF.ddL.plt.design)) {
      rF.ddL.plt.design <- rF.ddL.plt
    }
    linear.predictor <- plt.estimate.results$linear.predictor # dimension: N
    index.plt <- plt.estimate.results$index.plt
    cov.plt <- plt.estimate.results$cov.plt
    n.plt <- length(index.plt)

    ## subsampling step
    stage_start <- proc.time()[3]
    ssp.results <- rF.subsampling(inputs,
                               pi.plt = pi.plt,
                               varphi.plt = varphi.plt,
                               rF.ddL.plt = rF.ddL.plt.design,
                               linear.predictor = linear.predictor,
                               index.plt = index.plt
                               )
    stage.time["subsampling"] <- proc.time()[3] - stage_start
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    varphi.ssp <- ssp.results$varphi.ssp
    pi.ssp.N <- ssp.results$pi.ssp.N
    n.ssp <- length(index.ssp)
    
    ## subsample estimating step
    stage_start <- proc.time()[3]
    ssp.estimate.results <- rF.subsample.estimate(inputs,
                                               w.ssp = w.ssp,
                                               varphi.ssp = varphi.ssp,
                                               beta.plt = beta.plt,
                                               index.ssp = index.ssp,
                                               ...
                                               )
    stage.time["subsample.estimate"] <- proc.time()[3] - stage_start
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    beta.ssp <- ssp.estimate.results$beta.ssp
    cov.ssp <- ssp.estimate.results$cov.ssp
    
    
    ## combining step: refit on the union of pilot and second-step subsamples.
    stage_start <- proc.time()[3]
    rF.combining.union.results <- rF.combining.union(inputs,
                                   index.plt = index.plt, # dim: n.plt
                                   index.ssp = index.ssp, # dim: n.ssp
                                   pi.plt = pi.plt.N, # dim: N
                                   pi.ssp = pi.ssp.N # dim: N
                                   )
    stage.time["combining.union"] <- proc.time()[3] - stage_start
    index.cmb.union <- rF.combining.union.results$index.cmb
    beta.cmb.union <- rF.combining.union.results$beta.cmb
    cov.cmb.union <- rF.combining.union.results$cov.cmb
    beta.cmb <- beta.cmb.union
    cov.cmb <- cov.cmb.union
    index.cmb <- index.cmb.union

    ## prepare for displaying results 
    names(beta.cmb.union) <- names(beta.cmb) <-
      names(beta.ssp) <- names(beta.plt) <- colnames(X)
    
    
    
    stage_start <- proc.time()[3]
    rare.counts.results.plt <- rF.rare_counts(X, 
                                           index.plt, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.ssp <- rF.rare_counts(X, 
                                           index.ssp, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.cmb.union <- rF.rare_counts(X, 
                                           index.cmb.union, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.full <- rF.rare_counts(X, 
                                           seq_len(N), 
                                           rareFeature.index,
                                           Y)
    stage.time["result.assembly"] <- proc.time()[3] - stage_start
    
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef.cmb = beta.cmb,
                    coef.cmb.union = beta.cmb.union,
                    cov.plt = cov.plt,
                    cov.ssp = cov.ssp,
                    cov.cmb = cov.cmb,
                    cov.cmb.union = cov.cmb.union,
                    N = N,
                    family.name = family$family,
                    subsample.size.expect = subsample.size.expect,
                    subsample.size.actual = length(index.ssp),
                    full.Y.count = Y.count,
                    full.Y.sum = full.Y.sum,
                    full.Y.mean = full.Y.mean,
                    full.rare.count = full.rare.count,
                    rows.with.rare.full = rare.counts.results.full$rows.with.rare,
                    Y_subgroup_prev = Y_subgroup_prev,
                    
                    Y.count.plt = rare.counts.results.plt$Y.count,
                    Y.proportion.plt = rare.counts.results.plt$Y.proportion,
                    rare.count.plt = rare.counts.results.plt$rare.count,
                    rare.proportion.plt = rare.counts.results.plt$rare.proportion,
                    rare.counts.plt = rare.counts.results.plt$rare.counts,
                    rows.with.rare.plt = rare.counts.results.plt$rows.with.rare,
                    
                    Y.count.ssp = rare.counts.results.ssp$Y.count,
                    Y.proportion.ssp = rare.counts.results.ssp$Y.proportion,
                    rare.count.ssp = rare.counts.results.ssp$rare.count,
                    rare.proportion.ssp = rare.counts.results.ssp$rare.proportion,
                    rare.counts.ssp = rare.counts.results.ssp$rare.counts,
                    rows.with.rare.ssp = rare.counts.results.ssp$rows.with.rare,
                    
                    Y.count.cmb.union = rare.counts.results.cmb.union$Y.count,
                    Y.proportion.cmb.union = rare.counts.results.cmb.union$Y.proportion,
                    rare.count.cmb.union = rare.counts.results.cmb.union$rare.count,
                    rare.proportion.cmb.union = rare.counts.results.cmb.union$rare.proportion,
                    rare.counts.cmb.union = rare.counts.results.cmb.union$rare.counts,
                    rows.with.rare.cmb.union = rare.counts.results.cmb.union$rows.with.rare,
                    
                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    index.cmb = index.cmb,
                    index.cmb.union = index.cmb.union,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    stage.time = if (record.stage.time) stage.time else NULL,
                    terms = mt
                    )
    class(results) <- c("ssp.glm.rF", "list")
    return(results)
  } else if (criterion %in% c("Uni", "BL-Uni")){
    if (criterion == "Uni") inputs$balance.X.plt = FALSE
    inputs$n.uni <- n.plt + n.ssp
    t_start <- proc.time()[3]
    
    stage_start <- proc.time()[3]
    uni.estimate.results <- rF.uniform.estimate(inputs, ...)
    stage.time["uniform.estimate"] <- proc.time()[3] - stage_start
    
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    beta.uni <- uni.estimate.results$beta.uni
    index.uni <- uni.estimate.results$index.uni
    cov.uni <- uni.estimate.results$cov.uni


    stage_start <- proc.time()[3]
    rare.counts.results <- rF.rare_counts(X, index.uni, rareFeature.index, Y)
    rare.counts.results.full <- rF.rare_counts(X, seq_len(N), rareFeature.index, Y)
    stage.time["result.assembly"] <- proc.time()[3] - stage_start
    
    
    results <- list(model.call = model.call,
                    coef.plt = beta.uni,
                    coef.ssp = beta.uni,
                    coef.cmb = beta.uni,
                    coef.cmb.union = beta.uni,
                    cov.plt = cov.uni,
                    cov.ssp = cov.uni,
                    cov.cmb = cov.uni,
                    cov.cmb.union = cov.uni,
                    N = N,
                    family.name = family$family,
                    subsample.size.expect = inputs$n.uni,
                    subsample.size.actual = length(index.uni),
                    full.Y.count = Y.count,
                    full.Y.sum = full.Y.sum,
                    full.Y.mean = full.Y.mean,
                    full.rare.count = full.rare.count,
                    rows.with.rare.full = rare.counts.results.full$rows.with.rare,
                    Y_subgroup_prev = Y_subgroup_prev,
                    # plt, ssp and cmb are the same.
                    Y.count.plt = rare.counts.results$Y.count,
                    Y.proportion.plt = rare.counts.results$Y.proportion,

                    rare.count.plt = rare.counts.results$rare.count,
                    rare.proportion.plt = rare.counts.results$rare.proportion,
                    rare.counts.plt = rare.counts.results$rare.counts,
                    rows.with.rare.plt = rare.counts.results$rows.with.rare,
                    
                    Y.count.ssp = rare.counts.results$Y.count,
                    Y.proportion.ssp = rare.counts.results$Y.proportion,
                    rare.count.ssp = rare.counts.results$rare.count,
                    rare.proportion.ssp = rare.counts.results$rare.proportion,
                    rare.counts.ssp = rare.counts.results$rare.counts,
                    rows.with.rare.ssp = rare.counts.results$rows.with.rare,
                    
                    Y.count.cmb.union = rare.counts.results$Y.count,
                    Y.proportion.cmb.union = rare.counts.results$Y.proportion,
                    rare.count.cmb.union = rare.counts.results$rare.count,
                    rare.proportion.cmb.union = rare.counts.results$rare.proportion,
                    rare.counts.cmb.union = rare.counts.results$rare.counts,
                    rows.with.rare.cmb.union = rare.counts.results$rows.with.rare,
                    
                    index.plt = index.uni,
                    index.ssp = index.uni,
                    index.cmb = index.uni,
                    index.cmb.union = index.uni,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    stage.time = if (record.stage.time) stage.time else NULL,
                    terms = mt
    )
    class(results) <- c("ssp.glm.rF", "list")
    return(results)
  }
}


###############################################################################
