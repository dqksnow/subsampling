###############################################################################
rF.rare_counts <- function(X, index, rareFeature.index, Y = NULL) {
  
  if (is.null(rareFeature.index)) {
    result <- list(
      rare.counts = NA,
      rare.count = NA,
      unique.rare.count = NA,
      rare.proportion = NA,
      Y.count = if (!is.null(Y)) sum(Y[index], na.rm = TRUE) else NA,
      Y.proportion = if (!is.null(Y)) mean(Y[index], na.rm = TRUE) else NA
    )
    
  } else {
    X_subset <- X[index, rareFeature.index, drop = FALSE]
    
    if (is.null(dim(X_subset))) { 
      # only one rare feature column
      ones_count <- X_subset
      rare_count <- sum(X_subset == 1)
      unique_rare_count <- sum(X[unique(index), rareFeature.index] == 1)
      rare_proportion <- sum(X_subset > 0) / length(X_subset)
      
    } else {
      # multiple rare features
      ones_count <- rowSums(X_subset)
      rare_count <- colSums(X_subset == 1)
      unique_rare_count <- colSums(X[unique(index),
                                     rareFeature.index, drop = FALSE] == 1)
      rare_proportion <- sum(ones_count > 0) / nrow(X_subset)
    }
    
    # distribution of #rare features per observation
    counts <- sapply(0:length(rareFeature.index), 
                     function(k) sum(ones_count == k))
    names(counts) <- paste0(0:length(rareFeature.index), "_ones")
    
    
    rows_with_rare <- index[rowSums(X_subset) > 0]
    
    
    result <- list(
      rare.counts = counts,
      rare.count = rare_count,
      unique.rare.count = unique_rare_count,
      rare.proportion = rare_proportion,
      Y.count = if (!is.null(Y)) sum(Y[index], na.rm = TRUE) else NA,
      Y.proportion = if (!is.null(Y)) mean(Y[index], na.rm = TRUE) else NA,
      rows.with.rare = rows_with_rare
    )
  }
  
  return(result)
}


###############################################################################
rF.balance_score <- function(X, Y = NULL, family,
                          rareFeature.index = NULL,
                          threshold = 0.05) {
  
  # Handle Y if needed
  fam_name <- if(is.character(family)) family else family$family
  is_binomial <- fam_name %in% c("binomial", "quasibinomial")
  Y_count <- NULL
  Y_subgroup_prev <- NULL
  if (!is.null(Y) && is_binomial) {
    Y_count <- sum(Y, na.rm = TRUE)
  }
  
  
  N <- nrow(X)
  d <- ncol(X)
  
  # find binary variables
  cs <- colSums(X, na.rm = TRUE)
  css <- colSums(X^2, na.rm = TRUE)
  is_bin <- (cs > 0) & (cs < N) & (abs(cs - css) < 1e-6)
  bin_idx <- which(is_bin) # binary column index of X
  DN <- rep(1, d)
  
  if (length(bin_idx) == 0) { # no binary features in the data
    return(list(bl = rep(1, N),
                DN = DN,
                Y.count = Y_count, 
                binary.index = NULL,
                rareFeature.index = NULL))
  }
  
  bin_prev <- cs[bin_idx] / N
  
  if (is.null(rareFeature.index)) {
    rare_idx <- bin_idx[bin_prev < threshold]
  } else {
    if (is.character(rareFeature.index)) {
      if (is.null(colnames(X))) {
        stop("character 'rareFeature.index' requires column names in the model matrix")
      }
      matched_idx <- match(rareFeature.index, colnames(X))
      if (any(is.na(matched_idx))) {
        missing_names <- unique(rareFeature.index[is.na(matched_idx)])
        stop(
          paste0(
            "unknown rare feature names: ",
            paste(missing_names, collapse = ", ")
          )
        )
      }
      rare_idx <- matched_idx
    } else if (is.numeric(rareFeature.index)) {
      if (any(is.na(rareFeature.index)) ||
          any(rareFeature.index < 1) ||
          any(rareFeature.index > d) ||
          any(abs(rareFeature.index - round(rareFeature.index)) > .Machine$double.eps^0.5)) {
        stop("'rareFeature.index' must contain valid column indices of the model matrix")
      }
      rare_idx <- as.integer(rareFeature.index)
    } else {
      stop("'rareFeature.index' must be NULL, numeric indices, or column names")
    }
    
    rare_idx <- unique(rare_idx)
    non_binary_idx <- setdiff(rare_idx, bin_idx)
    if (length(non_binary_idx) > 0) {
      non_binary_names <- colnames(X)[non_binary_idx]
      if (is.null(non_binary_names)) {
        non_binary_names <- as.character(non_binary_idx)
      }
      stop(
        paste0(
          "the following supplied rare features are not binary: ",
          paste(non_binary_names, collapse = ", ")
        )
      )
    }
    
    supplied_prev <- bin_prev[match(rare_idx, bin_idx)]
    nonrare_idx <- rare_idx[supplied_prev >= threshold]
    if (length(nonrare_idx) > 0) {
      nonrare_names <- colnames(X)[nonrare_idx]
      if (is.null(nonrare_names)) {
        nonrare_names <- as.character(nonrare_idx)
      }
      warning(
        paste0(
          "the following supplied rare features have prevalence >= ",
          format(threshold, scientific = FALSE),
          ": ",
          paste0(
            nonrare_names,
            " (",
            format(round(supplied_prev[supplied_prev >= threshold], 4), nsmall = 4),
            ")",
            collapse = ", "
          )
        )
      )
    }
  }
  
  if (length(rare_idx) == 0) {
    return(list(
      bl = rep(1, N),
      DN = DN,
      binary.index = bin_idx,
      rareFeature.index = NULL,
      Y.count = Y_count,
      rare.count = numeric(0),
      Y_subgroup_prev = Y_subgroup_prev
    ))
  }
  
  X_rf <- X[, rare_idx, drop = FALSE]
  rare_prev <- cs[rare_idx] / N
  DN[rare_idx] <- sqrt(rare_prev)
  
  d_rare <- length(rare_idx)
  V_vec <- (1 - 2 * rare_prev) / rare_prev
  bl <- as.vector(d_rare + (X_rf %*% V_vec))
  
  rare_count <- cs[rare_idx]
  
  
  if (!is.null(Y) && is_binomial) {
    # rowsum is not rowSums
    counts <- rowsum(X_rf, Y) 
    Y_subgroup_prev <- counts / as.vector(table(Y))
  }

  list(
    bl = bl, 
    DN = DN, 
    binary.index = bin_idx,
    rareFeature.index = rare_idx, 
    Y.count = Y_count,
    rare.count = rare_count,
    Y_subgroup_prev = Y_subgroup_prev
  )
}


###############################################################################
rF.na_rate <- function(x) {
  sum(is.na(x)) / length(x)
}
###############################################################################
rF.glm.coef.estimate <- function(X,
                              Y,
                              start = rep(0, ncol(X)),
                              weights = 1,
                              family,
                              ...) {
  X <- as.matrix(X)
  if (length(weights) == 1L) {
    weights <- rep(weights, nrow(X))
  }
  
  if (family$family == "gaussian" && family$link == "identity") {
    results <- stats::lm.wfit(x = X, y = Y, w = weights)
    beta <- results$coefficients
  } else {
    results <- stats::glm.fit(
      x = X,
      y = Y,
      weights = weights,
      start = start,
      family = family,
      intercept = FALSE
    )
    beta <- results$coefficients
  }
  cov <- NULL
  
  return(list(beta = beta,
              cov = cov))
}

###############################################################################
rF.poisson.index <- function (N, pi) {
  return(which(runif(N) <= pi))
}

###############################################################################
rF.calculate.nm <- function(X, Y, rF.ddL.plt, linear.predictor, criterion, 
                              bl, DN = 1){
  
  if (criterion == "BL-Lopt"){
    nm <- sqrt(rowSums(X^2))
    nm <- bl * abs(Y - linear.predictor) * nm
  } else if (criterion == "Aopt"){ # original Aopt
    nm <- sqrt(rowSums((X %*% t(solve(rF.ddL.plt)))^2))
    nm <- abs(Y - linear.predictor) * nm # numerator
  } else if (criterion == "Lopt"){ # original Lopt
    nm <- sqrt(rowSums(X^2))
    nm <- abs(Y - linear.predictor) * nm
  } else if (criterion == "R-Lopt"){
    DN_inv <- 1 / (DN^2)
    X <- sweep(X, 2, DN_inv, "*")
    nm <- sqrt(rowSums(X^2))
    nm <- abs(Y - linear.predictor) * nm
  }
  
  return(nm)
}
rF.poi.prob.exact <- function(nm, n.ssp, N){
  ## compute exact H and exact denominator.
  ## then get the exact poisson sampling probabilities
  h.ssp <- n.ssp / N
  E.nm <- mean(nm) # denominator, the expectation of numerator
  H <- NULL
  g <- sum((nm / E.nm) > (1 / h.ssp))
  # g: the count of sampling probabilities larger than 1
  if((g != 0) & (!is.na(g))){
    # find the largest g elements in nm
    largest_g_sum <- sum(
      sort(nm, partial = (N-g+1):N)[(N-g+1):N]
    )
    H <- (E.nm*N - largest_g_sum) / (n.ssp - g)
    # re-compute nm for those exceed threshold H 
    nm <- pmin(nm, H)
    # print(paste('exact H: ', H))
  }

  return(list(H = H,
              E.nm = mean(nm),
              nm = nm))
}
###############################################################################
rF.poi.prob.estimated <- function(nm,
                               index.plt,
                               pi.plt,
                               n.ssp,
                               N,
                               b = 2) {
  if (length(index.plt) == 0) {
    stop("pilot-based truncation requires at least one pilot observation")
  }
  if (length(index.plt) != length(pi.plt)) {
    stop("'index.plt' and 'pi.plt' must have the same length")
  }
  if (any(!is.finite(pi.plt)) || any(pi.plt <= 0)) {
    stop("'pi.plt' must contain finite positive inclusion probabilities")
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0) {
    stop("'b' must be a positive finite scalar")
  }

  h.ssp <- n.ssp / N
  nm.plt <- nm[index.plt]
  trunc.quantile <- 1 - n.ssp / (b * N)
  trunc.quantile <- max(0, min(1, trunc.quantile))
  H <- as.numeric(
    stats::quantile(nm.plt, probs = trunc.quantile, names = FALSE, na.rm = TRUE)
  )
  nm.adj <- pmin(nm, H)
  E.nm <- sum(pmin(nm.plt, H) / pi.plt) / N
  if (!is.finite(E.nm) || E.nm <= 0) {
    stop("pilot-based estimate of mean(nm) must be positive and finite")
  }
  pi.ssp <- h.ssp * nm.adj / E.nm
  sum.pi <- sum(pmin(pi.ssp, 1))

  list(
    H = H,
    E.nm = E.nm,
    nm = nm.adj,
    iter = 1L,
    b = b,
    sum.pi = sum.pi
  )
}
###############################################################################
## second derivative of log likelihood function
rF.ddL <- function (eta, X, weights = 1, family) {
  # eta is linear predictor plus offsets(if have)
  # linkinv(eta):
  # for binomial = exp(eta)/(1+exp(eta))
  # poisson = exp(eta)
  # Gamma(link = "inverse") = 1/eta
  # variance(mu):
  # for binomial = mu(1 - mu)
  # for poisson = mu
  # for Gamma(link = "inverse") = mu^2
  variance <- family$variance
  linkinv  <- family$linkinv
  
  dlinear.predictor <- variance(linkinv(eta))
  
  rF.ddL <- t(X) %*% (X * (dlinear.predictor * weights))
  return(rF.ddL)
}

## square of the first derivative of log likelihood function
rF.cov.score <- function (eta, X, Y, weights = 1, family) {
  variance <- family$variance
  linkinv  <- family$linkinv
  
  temp <- (Y - linkinv(eta))^2
  
  rF.cov.score <- t(X) %*% (X * (temp * weights))
  return(rF.cov.score)
}
###############################################################################
rF.pilot.estimate <- function(inputs, ...){
  X <- inputs$X
  Y <- inputs$Y
  n.plt <- inputs$n.plt
  family <- inputs$family
  N <- inputs$N
  bl <- inputs$bl
  objective.weight.plt <- inputs$objective.weight.plt
  linkinv  <- family$linkinv
  balance.X.plt <- inputs$balance.X.plt
  balance.Y.plt <- inputs$balance.Y.plt
  Y.count <- inputs$Y.count
  
  h.plt <- n.plt / N
  if (family[["family"]] %in% c('binomial', 'quasibinomial')){
    ## for binary response, three strategies: 
    ## uniform; balancing X; balancing both X and Y.
    
    if (balance.X.plt && length(bl) != 1L) {
      # bl.mean <- mean(bl)
      # varphi.plt <-  bl / bl.mean
      if(balance.Y.plt == TRUE){
        # balancing X and case control Y
        # will not consider truncation H, since n.plt is supposed to be small.
        N1 <- Y.count
        N0 <- N - N1
        if (N0 < n.plt / 2 | N1 < n.plt / 2) {
          warning(
            paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.")
          )
        }
        w1 <- sum(bl[Y == 1])  # sum of bl in Y=1
        w0 <- sum(bl[Y == 0])  # sum of bl in Y=0
        den <- (1 - Y) * w0 + Y * w1
        varphi.plt <- (N / 2) * bl / den
        
        pi.plt <- pmin(varphi.plt * h.plt, 1) # dimension: N
        
        index.plt <- rF.poisson.index(N, pi.plt)
        pi.plt.N <- pi.plt  # dimension: N
        pi.plt <- pi.plt[index.plt]  # dimension: n.plt
        varphi.plt <- varphi.plt[index.plt]   # dimension: n.plt
        w.plt <- 1 / pi.plt
      } else {
        # only balancing X
        bl.mean <- mean(bl)
        varphi.plt <-  bl / bl.mean
        pi.plt <- pmin(varphi.plt * h.plt, 1) # dimension: N
        index.plt <- rF.poisson.index(N, pi.plt)
        pi.plt.N <- pi.plt
        pi.plt <- pi.plt[index.plt]
        varphi.plt <- varphi.plt[index.plt]
        w.plt <- 1 / pi.plt
      }
    } else{
      # does not consider balancing X; 
      # do uniform sampling or case control sampling
      if(balance.Y.plt == TRUE){
        N1 <- Y.count
        N0 <- N - N1
        varphi.plt <- (1 - Y) * (N/(2*N0)) + Y * (N/(2*N1))
        pi.plt <- h.plt * varphi.plt
        index.plt <- rF.poisson.index(N, pi.plt)
        pi.plt.N <- pi.plt
        pi.plt <- pi.plt[index.plt]
        varphi.plt <- varphi.plt[index.plt]
        w.plt <- 1 / pi.plt
      } else {
        # purely uniform 
        index.plt <- rF.poisson.index(N, h.plt)
        pi.plt <- rep(n.plt/N, length(index.plt)) # dimension: n.plt
        pi.plt.N <- rep(n.plt/N, N)
        varphi.plt <- rep(1, length(index.plt)) # all 1
        w.plt <- 1 / pi.plt
      }
    }
    
    n.plt <- length(index.plt) # actuall n.plt
    x.plt <- X[index.plt,] # dimension: n.plt
    Y.plt <- Y[index.plt] # dimension: n.plt
    fit.weights.plt <- if (objective.weight.plt == "weighted") {
      w.plt
    } else {
      rep(1, n.plt) # do unweighted estimation if the weights does not
      # reply on the response
    }
    info.weights.plt <- w.plt
    beta.plt <- rF.glm.coef.estimate(X = x.plt, Y = Y.plt, 
                                  weights = fit.weights.plt,
                                  family = family, ...)$beta
    linear.predictor.plt <- as.vector(x.plt %*% beta.plt) # dimension: n.plt
    
    rF.ddL.plt <- (1/N) * rF.ddL(eta = linear.predictor.plt,
                          X = x.plt,
                          weights = fit.weights.plt,
                          family = family)
    rF.cov.score.plt <- (n.plt / (N^2)) * rF.cov.score(eta = linear.predictor.plt,
                                                 x.plt,
                                                 Y.plt,
                                                 weights = fit.weights.plt^2,
                                                 family = family)
    rF.ddL.plt.design <- (1/N) * rF.ddL(eta = linear.predictor.plt,
                                  X = x.plt,
                                  weights = info.weights.plt,
                                  family = family)
    # rF.ddL.plt.design is for approximating full rF.ddL so we use weights.
    linear.predictor <- linkinv(X %*% beta.plt) # dimension: N
    rF.ddL.plt.inv <- solve(rF.ddL.plt)
    cov.plt <- (1/n.plt) * rF.ddL.plt.inv %*% rF.cov.score.plt %*% rF.ddL.plt.inv
  } else { # for non-binary-response models:
    if (balance.X.plt && length(bl) != 1L) {
      # balancing X
      bl.mean <- mean(bl)
      varphi.plt <- bl / bl.mean
      pi.plt <- pmin(varphi.plt * h.plt, 1)
      index.plt <- rF.poisson.index(N, pi.plt)
      pi.plt.N <- pi.plt
      pi.plt <- pi.plt[index.plt]
      varphi.plt <- varphi.plt[index.plt]
      w.plt <- 1 / pi.plt
    } else {
      # uniform
      index.plt <- rF.poisson.index(N, h.plt)
      pi.plt <- rep(n.plt / N, length(index.plt))
      pi.plt.N <- rep(n.plt / N, N)
      varphi.plt <- rep(1, length(index.plt))
      w.plt <- 1 / pi.plt
    }
    
    n.plt <- length(index.plt)
    x.plt <- X[index.plt, , drop = FALSE]
    Y.plt <- Y[index.plt]
    fit.weights.plt <- if (objective.weight.plt == "weighted") {
      w.plt
    } else {
      rep(1, n.plt)
    }
    info.weights.plt <- w.plt
    beta.plt <- rF.glm.coef.estimate(
      X = x.plt,
      Y = Y.plt,
      weights = fit.weights.plt,
      family = family,
      ...
    )$beta
    linear.predictor.plt <- as.vector(x.plt %*% beta.plt)
    rF.ddL.plt <- (1 / N) * rF.ddL(
      eta = linear.predictor.plt,
      X = x.plt,
      weights = fit.weights.plt,
      family = family
    )
    rF.cov.score.plt <- (n.plt / (N^2)) * rF.cov.score(
      eta = linear.predictor.plt,
      x.plt,
      Y.plt,
      weights = fit.weights.plt^2,
      family = family
    )
    # rF.ddL.plt.design is for approximating full rF.ddL so we use weights.
    rF.ddL.plt.design <- (1 / N) * rF.ddL(
      eta = linear.predictor.plt,
      X = x.plt,
      weights = info.weights.plt,
      family = family
    )
    linear.predictor <- linkinv(X %*% beta.plt)
    rF.ddL.plt.inv <- solve(rF.ddL.plt)
    cov.plt <- (1 / n.plt) * rF.ddL.plt.inv %*% rF.cov.score.plt %*% rF.ddL.plt.inv
  }
  return(list(pi.plt = pi.plt, 
              pi.plt.N = pi.plt.N, # dim:N
              varphi.plt = varphi.plt,
              beta.plt = beta.plt,
              rF.ddL.plt = rF.ddL.plt,
              rF.ddL.plt.design = rF.ddL.plt.design,
              rF.cov.score.plt = rF.cov.score.plt,
              linear.predictor = linear.predictor,
              index.plt = index.plt,
              cov.plt = cov.plt
              )
         )
}
###############################################################################
rF.subsampling <- function(inputs,
                        pi.plt,
                        varphi.plt,
                        rF.ddL.plt,
                        linear.predictor,
                        index.plt) {
  X <- inputs$X
  Y <- inputs$Y
  family <-  inputs$family
  n.ssp <- inputs$n.ssp
  N <- inputs$N
  control <- inputs$control
  alpha <- control$alpha
  b <- control$b
  poi.method <- control$poi.method
  criterion <- inputs$criterion
  sampling.method <- inputs$sampling.method
  bl <- inputs$bl
  DN <- inputs$DN
  balance.Y.all <- inputs$balance.Y.all
  h.ssp <- n.ssp / N
  
  ## When use poisson sampling method to draw pilot sample, 
  ## length(index.plt) might be different with original n.plt
  ## so here we reset n.plt.
  n.plt <- length(index.plt)
  w.ssp <- NA
  
  nm <- rF.calculate.nm(X, Y, rF.ddL.plt, linear.predictor, criterion, 
                          bl, DN) # numerator. dimension: N

  get_poi_prob_results <- function(nm_values,
                                   pilot_index,
                                   pilot_pi,
                                   target_n,
                                   population_size) {
    if (poi.method == "exact") {
      return(rF.poi.prob.exact(nm_values, target_n, population_size))
    }
    rF.poi.prob.estimated(
      nm = nm_values,
      index.plt = pilot_index,
      pi.plt = pilot_pi,
      n.ssp = target_n,
      N = population_size,
      b = b
    )
  }
  
  
  ### this is for binomial only
  if (balance.Y.all) {
    # include all Y=1 and do rF.subsampling for Y=0
    N0 <- N - sum(Y)
    control_positions <- which(Y == 0)
    pilot_control_global <- index.plt[Y[index.plt] == 0]
    pilot_control_local <- match(pilot_control_global, control_positions)
    pilot_control_pi <- pi.plt[Y[index.plt] == 0]
    poi.prob.results <- get_poi_prob_results(
      nm_values = nm[Y == 0],
      pilot_index = pilot_control_local,
      pilot_pi = pilot_control_pi,
      target_n = n.ssp,
      population_size = N0
    )
    nm.control <- poi.prob.results$nm # dim: N0
    h.ssp.control <- n.ssp / N0
    varphi.control <- nm.control / poi.prob.results$E.nm
    varphi.control <- (1 - alpha) * varphi.control + alpha
    pi.ssp <- rep(NA, N)
    pi.ssp[Y==1] <- 1  # maybe not efficient
    pi.ssp[Y==0] <- h.ssp.control * varphi.control # [Y==0]
  } else {
    poi.prob.results <- get_poi_prob_results(
      nm_values = nm,
      pilot_index = index.plt,
      pilot_pi = pi.plt,
      target_n = n.ssp,
      population_size = N
    )
    nm <- poi.prob.results$nm # dim: N
    varphi.ssp <- nm / poi.prob.results$E.nm
    varphi.ssp <- (1 - alpha) * varphi.ssp + alpha
    pi.ssp <- h.ssp * varphi.ssp
  }
  
  index.ssp <- rF.poisson.index(N, pi.ssp)
  pi.ssp.N <- pi.ssp
  pi.ssp <- pi.ssp[index.ssp] # dim: actual n.ssp
  w.ssp <- 1 / pi.ssp # dim: actual n.ssp
  varphi.ssp = pi.ssp / h.ssp
  
  
  return(list(index.ssp = index.ssp,
              w.ssp = w.ssp,
              varphi.ssp = varphi.ssp,
              pi.ssp = pi.ssp,
              pi.ssp.N = pi.ssp.N
              )
         )
}
###############################################################################
rF.subsample.estimate <- function(inputs,
                               w.ssp,
                               varphi.ssp,
                               beta.plt,
                               index.ssp,
                               ...) {
  x.ssp <- inputs$X[index.ssp, ]
  y.ssp = inputs$Y[index.ssp]
  n.ssp <- length(index.ssp) # re-calculate the actual second-step size
  N <- inputs$N
  sampling.method <- inputs$sampling.method
  family <- inputs$family

  
  results.ssp <- rF.glm.coef.estimate(x.ssp,
                                   y.ssp,
                                   weights = w.ssp,
                                   family = family,
                                   ...)
  beta.ssp <- results.ssp$beta
  linear.predictor.ssp <- as.vector(x.ssp %*% beta.ssp)
  
  if (sampling.method == 'poisson') {
    rF.ddL.ssp <- (1/N) * rF.ddL(linear.predictor.ssp,
                           x.ssp,
                           weights = w.ssp,
                           family = family)
    rF.cov.score.ssp <- (n.ssp / (N^2)) * rF.cov.score(linear.predictor.ssp,
                                                 x.ssp,
                                                 y.ssp,
                                                 weights = w.ssp^2,
                                                 family = family)
  }
  rF.ddL.ssp.inv <- solve(rF.ddL.ssp)
  cov.ssp <- (1/n.ssp) * rF.ddL.ssp.inv %*% rF.cov.score.ssp %*% rF.ddL.ssp.inv 
  return(list(beta.ssp = beta.ssp,
              rF.ddL.ssp = rF.ddL.ssp,
              rF.cov.score.ssp = rF.cov.score.ssp,
              cov.ssp = cov.ssp
              )
         )
}
###############################################################################
rF.combining.union <- function(inputs,
                      index.plt, # dim: n.plt
                      index.ssp, # dim: n.ssp
                      pi.plt, # dim: N
                      pi.ssp, # dim: N
                      ...) {
  X <- inputs$X # dim: N
  Y <- inputs$Y# dim: N
  family <- inputs$family
  N <- inputs$N

  index.cmb <- union(index.plt, index.ssp)
  n.cmb <- length(index.cmb)
  pi.plt <- pi.plt[index.cmb] # becomes dim: n.cmb
  pi.ssp <- pi.ssp[index.cmb] # becomes dim: n.cmb
  pi.cmb <- pi.plt + pi.ssp - pi.plt * pi.ssp
  X.cmb <- X[index.cmb, ]
  Y.cmb <- Y[index.cmb]
  w.cmb <- 1 / pi.cmb
  
  results.cmb <- rF.glm.coef.estimate(X.cmb,
                                   Y.cmb,
                                   weights = w.cmb,
                                   family = family,
                                   ...)
  beta.cmb <- results.cmb$beta
  linear.predictor.cmb <- as.vector(X.cmb %*% beta.cmb)
  
  rF.ddL.cmb <- (1/N) * rF.ddL(linear.predictor.cmb,
                         X.cmb,
                         weights = w.cmb,
                         family = family)
  rF.cov.score.cmb <- (n.cmb / (N^2)) * rF.cov.score(linear.predictor.cmb,
                                               X.cmb,
                                               Y.cmb,
                                               weights = w.cmb^2,
                                               family = family)
  rF.ddL.cmb.inv <- solve(rF.ddL.cmb)
  cov.cmb <- (1 / n.cmb) * rF.ddL.cmb.inv %*% rF.cov.score.cmb %*% rF.ddL.cmb.inv  

  
  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              index.cmb = index.cmb
              )
         )
}
###############################################################################
rF.uniform.estimate <- function(inputs, ...){
  X <- inputs$X
  Y <- inputs$Y
  n.uni <- inputs$n.uni
  family <- inputs$family
  N <- inputs$N
  bl <- inputs$bl
  linkinv  <- family$linkinv
  balance.Y.ssp <- inputs$balance.Y.ssp
  balance.Y.all <- inputs$balance.Y.all
  Y.count <- inputs$Y.count
  criterion <- inputs$criterion
  objective.weight <- inputs$objective.weight
  
  h.uni <- n.uni / N # rF.subsampling rate 
  
  if (family[["family"]] %in% c('binomial', 'quasibinomial')){
    
    ## for binary response, it optionally takes the rarity of Y into account.
    if (criterion == 'BL-Uni' && length(bl) != 1L) {
      bl.mean <- mean(bl)
      varphi.uni <-  bl / bl.mean # normalized to mean 1 
      pi.uni <- h.uni * varphi.uni # un-truncated sampling prob
      
      ## truncate those pi.uni > 1, while enforce the summation to be n
      if(balance.Y.all == TRUE){
        # Use BL to balance X and use all Y=1 
        # do truncation for Y=0, while keeping all Y=1
        pi.uni.Y0 <- pi.uni[Y==0]
        g <- sum(pi.uni.Y0 > 1) 
        if (g == 0) { # no truncation needed
          pi.uni <- Y + (1 - Y) * pi.uni # dim: N
          index.uni <- rF.poisson.index(N, pi.uni)
          pi.uni <- pmin(pi.uni[index.uni], 1) # dim: n.uni
          Y.uni <- Y[index.uni]
          varphi.uni <- (Y.uni + (1 - Y.uni) * pi.uni) / h.uni # dim: n.uni
          w.uni <- 1 / pi.uni
        } else {
          bl.Y0 <- bl[Y==0]
          N0 <- length(bl.Y0) 
          top_g <- sort(bl.Y0, partial = (N0-g+1):N0)[(N0-g+1):N0]
          B_g   <- sum(top_g)
          H     <- (sum(bl.Y0) - B_g) / (n.uni - g)
          bl.Y0.trunc   <- pmin(bl.Y0, H) # truncation
          varphi.uni.Y0 <-  bl.Y0.trunc / mean(bl.Y0.trunc) # re-normalization
          pi.uni[Y==0] <- h.uni * varphi.uni.Y0 # sums to 1, max <=1
          pi.uni[Y==1] <- 1
          index.uni <- rF.poisson.index(N, pi.uni)
          pi.uni <- pi.uni[index.uni]
          Y.uni <- Y[index.uni]
          varphi.uni <- (Y.uni + (1 - Y.uni) * pi.uni) / h.uni # dim: n.uni
          w.uni <- 1 / pi.uni
        }
      } else if (balance.Y.ssp == TRUE){
        # Use BL to balance X and CC to balance Y,
        # while keeping total expected sample size equal to n.uni
        
        N1 <- Y.count
        N0 <- N - N1
        
        if (N1 < n.uni / 2) {
          warning("n.uni/2 exceeds the number of Y=1 in the full data;
                  taking all Y=1 and allocating the remainder to Y=0.")
          target <- c("0" = n.uni - N1, "1" = N1)
        } else if (N0 < n.uni / 2) {
          warning("n.uni/2 exceeds the number of Y=0 in the full data;
                  taking all Y=0 and allocating the remainder to Y=1.")
          target <- c("0" = N0, "1" = n.uni - N0)
        } else {
          target <- c("0" = n.uni / 2, "1" = n.uni / 2)
        }
        
        pi.uni <- numeric(N)
        
        for (yy in 0:1) {
          idx <- (Y == yy)
          bl.yy <- bl[idx]
          N.yy <- sum(idx)
          n.yy <- target[as.character(yy)]
          
          if (n.yy >= N.yy) {
            pi.uni[idx] <- 1
          } else {
            pi.yy <- n.yy * bl.yy / sum(bl.yy)
            
            if (max(pi.yy) <= 1) {
              pi.uni[idx] <- pi.yy
            } else {
              g <- sum(pi.yy > 1)
              top_g <- sort(bl.yy, partial = (N.yy - g + 1):N.yy)[(N.yy - g + 1):N.yy]
              H <- (sum(bl.yy) - sum(top_g)) / (n.yy - g)
              bl.yy <- pmin(bl.yy, H)
              pi.uni[idx] <- n.yy * bl.yy / sum(bl.yy)
            }
          }
        }
        
        index.uni <- rF.poisson.index(N, pi.uni)
        pi.uni.N <- pi.uni
        pi.uni <- pi.uni[index.uni]
        varphi.uni <- pi.uni / h.uni
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
      } else { 
        # Use BL to balance X only 
        g <- sum(pi.uni > 1)
        if (g == 0) { # no truncation needed
          index.uni <- rF.poisson.index(N, pi.uni)
          pi.uni <- pi.uni[index.uni]
          Y.uni <- Y[index.uni]
          varphi.uni <- pi.uni / h.uni
          w.uni <- 1 / pi.uni
        } else {
          top_g <- sort(bl, partial = (N-g+1):N)[(N-g+1):N]
          B_g   <- sum(top_g)
          H     <- (N*bl.mean - B_g) / (n.uni - g)
          bl.trunc   <- pmin(bl, H) # truncation
          varphi.uni <-  bl.trunc / mean(bl.trunc) # re-normalization
          pi.uni <- varphi.uni * h.uni   # sums to 1, max <=1
          index.uni <- rF.poisson.index(N, pi.uni)
          pi.uni <- pi.uni[index.uni]
          w.uni <- 1 / pi.uni
          varphi.uni <- pi.uni / h.uni
          Y.uni <- Y[index.uni]
        }
      }
    } else {
      ########################
      # criterion == 'Uni'
      # do poisson sampling with uniform / CC sampling prob
      if (balance.Y.all == TRUE) { # include all Y=1 and do uniform for Y=0
        pi.uni <- Y + (1 - Y) * h.uni
        index.uni <- rF.poisson.index(N, pi.uni)
        pi.uni <- pi.uni[index.uni] # dimension: n.uni
        varphi.uni <- pi.uni / h.uni
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
        
      } else if (balance.Y.ssp == TRUE) {
        N1 <- Y.count
        N0 <- N - N1
        
        if (N1 < n.uni / 2) {
          warning("n.uni/2 exceeds the number of Y=1 in the full data; 
                  taking all Y=1 and sampling the remainder from Y=0.")
          pi.uni <- Y + (1 - Y) * ((n.uni - N1) / N0)
        } else if (N0 < n.uni / 2) {
          warning("n.uni/2 exceeds the number of Y=0 in the full data; 
                  taking all Y=0 and sampling the remainder from Y=1.")
          pi.uni <- (1 - Y) + Y * ((n.uni - N0) / N1)
        } else {
          varphi.uni <- (1 - Y) * (N / (2 * N0)) + Y * (N / (2 * N1))
          pi.uni <- h.uni * varphi.uni
        }
        
        pi.uni <- pmin(pi.uni, 1) # size: N
        index.uni <- rF.poisson.index(N, pi.uni)
        
        pi.uni <- pi.uni[index.uni] # size: Exp(n.uni)
        varphi.uni <- pi.uni / h.uni
        
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
      } else { # purely uniform
        index.uni <- rF.poisson.index(N, h.uni)
        pi.uni <- rep(h.uni, length(index.uni)) # dimension: n.uni
        varphi.uni <- rep(1, length(index.uni)) # all 1
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
      }
    }
    
    
    n.uni <- length(index.uni) # actuall n.uni
    X.uni <- X[index.uni,] # dimension: n.uni
    fit.weights.uni <- if (objective.weight == "weighted") {
      w.uni
    } else {
      rep(1, n.uni)
    }
    beta.uni <- rF.glm.coef.estimate(X = X.uni, Y = Y.uni, 
                                  weights = fit.weights.uni,
                                  family = family, ...)$beta
    linear.predictor.uni <- as.vector(X.uni %*% beta.uni) # dimension: n.uni
    
    # rF.ddL.uni and rF.cov.score.uni: the pilot estimation of score^2 and information
    if (objective.weight == "weighted") {
      rF.ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                     X = X.uni,
                     weights = (1 / (varphi.uni * n.uni)),
                     family = family)
      rF.cov.score.uni <- rF.cov.score(eta = linear.predictor.uni,
                         X.uni,
                         Y.uni,
                         weights = 1 / ((varphi.uni^2) * n.uni),
                         family = family)
    } else {
      rF.ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                     X = X.uni,
                     weights = rep(1 / n.uni, n.uni),
                     family = family)
      rF.cov.score.uni <- rF.cov.score(eta = linear.predictor.uni,
                         X.uni,
                         Y.uni,
                         weights = rep(1 / n.uni, n.uni),
                         family = family)
    }
    linear.predictor <- NA # dimension: N
    rF.ddL.uni.inv <- solve(rF.ddL.uni)
    cov.uni <- (1 / n.uni) * rF.ddL.uni.inv %*% rF.cov.score.uni %*% rF.ddL.uni.inv
  } else {
    ## for other families: only consider X rarity, not Y rarity
    if (criterion == 'BL-Uni' && length(bl) != 1L) {
      bl.mean <- mean(bl)
      varphi.uni <- bl / bl.mean
      pi.uni <- h.uni * varphi.uni
      g <- sum(pi.uni > 1)
      if (g == 0) {
        index.uni <- rF.poisson.index(N, pi.uni)
        pi.uni <- pi.uni[index.uni]
        varphi.uni <- pi.uni / h.uni
        w.uni <- 1 / pi.uni
      } else {
        top_g <- sort(bl, partial = (N - g + 1):N)[(N - g + 1):N]
        B_g <- sum(top_g)
        H <- (N * bl.mean - B_g) / (n.uni - g)
        bl.trunc <- pmin(bl, H)
        varphi.uni <- bl.trunc / mean(bl.trunc)
        pi.uni <- varphi.uni * h.uni
        index.uni <- rF.poisson.index(N, pi.uni)
        pi.uni <- pi.uni[index.uni]
        varphi.uni <- pi.uni / h.uni
        w.uni <- 1 / pi.uni
      }
    } else {
      index.uni <- rF.poisson.index(N, n.uni / N)
      pi.uni <- rep(n.uni/N, length(index.uni)) # dimension: n.uni
      varphi.uni <- rep(1, length(index.uni)) # all 1
      w.uni <- 1 / pi.uni
    }
    
    n.uni <- length(index.uni) # actuall n.uni
    X.uni <- X[index.uni,] # dimension: n.uni
    Y.uni <- Y[index.uni] # dimension: n.uni
    fit.weights.uni <- if (objective.weight == "weighted") {
      w.uni
    } else {
      rep(1, n.uni)
    }
    beta.uni <- rF.glm.coef.estimate(X = X.uni, Y = Y.uni, 
                                  weights = fit.weights.uni,
                                  family = family, ...)$beta
    linear.predictor.uni <- as.vector(X.uni %*% beta.uni) # dimension: n.uni
    
    # rF.ddL.uni and rF.cov.score.uni: the pilot estimation of gradient^2 and information
    if (objective.weight == "weighted") {
      rF.ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                     X = X.uni,
                     weights = (1 / (varphi.uni * n.uni)),
                     family = family)
      rF.cov.score.uni <- rF.cov.score(eta = linear.predictor.uni,
                         X.uni,
                         Y.uni,
                         weights = 1 / ((varphi.uni^2) * n.uni),
                         family = family)
    } else {
      rF.ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                     X = X.uni,
                     weights = rep(1 / n.uni, n.uni),
                     family = family)
      rF.cov.score.uni <- rF.cov.score(eta = linear.predictor.uni,
                         X.uni,
                         Y.uni,
                         weights = rep(1 / n.uni, n.uni),
                         family = family)
    }
    linear.predictor <- NA # dimension: N
    rF.ddL.uni.inv <- solve(rF.ddL.uni)
    cov.uni <- (1 / n.uni) * rF.ddL.uni.inv %*% rF.cov.score.uni %*% rF.ddL.uni.inv
  }
  return(list(pi.uni = pi.uni,
              varphi.uni = varphi.uni,
              beta.uni = beta.uni,
              rF.ddL.uni = rF.ddL.uni,
              rF.cov.score.uni = rF.cov.score.uni,
              linear.predictor = linear.predictor,
              index.uni = index.uni,
              cov.uni = cov.uni
              )
         )
}
###############################################################################
rF.glm.control <- function(alpha = 0, b = 2, poi.method = c("exact", "estimated"), ...)
{
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("sampling probability weight 'alpha' must between [0, 1]")
  if(!is.numeric(b) || b < 0)
    stop("sampling probability threshold 'b' must > 0")
  poi.method <- match.arg(poi.method)
  list(alpha = alpha, b = b, poi.method = poi.method)
}
###############################################################################
rF.format_p_values <- function(p.values, threshold = 0.0001) {
  formatted <- sapply(p.values, function(p.value) {
    if (p.value < threshold) {
      return(sprintf("<%.4f", threshold))
    } else {
      return(sprintf("%.4f", p.value))
    }
  })
  return(formatted)
}
###############################################################################
#' @export
summary.ssp.glm.rF <- function(object, ...) {
  ### collect results
  coef <- object$coef.ssp
  se <- sqrt(diag(object$cov.ssp))
  N <- object$N
  family.name <- object$family.name
  n.plt.actual <- length(object$index.plt)
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index.ssp)
  n.ssp.unique <- length(unique(object$index.ssp))
  n.cmb.actual <- length(object$index.cmb.union)
  overlap.count <- length(intersect(object$index.plt, object$index.ssp))
  subsample.rate.expect <- (n.ssp.expect / N) * 100
  subsample.rate.actual <- (n.ssp.actual / N) * 100
  subsample.rate.unique <- (n.ssp.unique / N) * 100
  
  safe_pct <- function(num, den, digits = 2) {
    if (is.null(num) || is.null(den) || length(num) == 0 || den == 0) {
      return(NA_character_)
    }
    paste0(round(100 * num / den, digits), "%")
  }
  overlap.rate.ssp <- safe_pct(overlap.count, n.ssp.actual)
  overlap.rate.plt <- safe_pct(overlap.count, n.plt.actual)
  full.rare.row.rate <- if (!is.null(object$rows.with.rare.full)) {
    safe_pct(length(object$rows.with.rare.full), N)
  } else {
    NA_character_
  }
  
  y_summary <- data.frame(
    Sample = c("Pilot", "Subsample", "Combined union", "Full data"),
    Size = c(n.plt.actual, n.ssp.actual, n.cmb.actual, N),
    Rare_row_rate = c(
      safe_pct(length(object$rows.with.rare.plt), n.plt.actual),
      safe_pct(length(object$rows.with.rare.ssp), n.ssp.actual),
      safe_pct(length(object$rows.with.rare.cmb.union), n.cmb.actual),
      full.rare.row.rate
    ),
    check.names = FALSE
  )
  if (family.name %in% c("binomial", "quasibinomial")) {
    y_summary <- data.frame(
      y_summary["Sample"],
      y_summary["Size"],
      Y_count = c(object$Y.count.plt,
                  object$Y.count.ssp,
                  object$Y.count.cmb.union,
                  object$full.Y.count),
      Y_rate = c(
        safe_pct(object$Y.count.plt, n.plt.actual),
        safe_pct(object$Y.count.ssp, n.ssp.actual),
        safe_pct(object$Y.count.cmb.union, n.cmb.actual),
        safe_pct(object$full.Y.count, N)
      ),
      y_summary["Rare_row_rate"],
      check.names = FALSE
    )
  }
  
  ### header
  cat("Model Summary\n")
  cat("\nCall:\n\n")
  print(object$model.call)
  cat("\n")
  
  ### design info
  cat("Design Summary:\n")
  size_table <- data.frame(
    Variable = c(
      'Pilot Sample Size',
      'Expected Subsample Size',
      'Actual Subsample Size',
      'Unique Subsample Size',
      'Combined Union Size',
      'Overlap Count',
      'Overlap Rate in Pilot',
      'Overlap Rate in Subsample',
      'Expected Subsample Rate',
      'Actual Subsample Rate',
      'Unique Subsample Rate'
    ),
    Value = c(
      n.plt.actual,
      n.ssp.expect,
      n.ssp.actual,
      n.ssp.unique,
      n.cmb.actual,
      overlap.count,
      overlap.rate.plt,
      overlap.rate.ssp,
      paste0(round(subsample.rate.expect, 2), "%"),
      paste0(round(subsample.rate.actual, 2), "%"),
      paste0(round(subsample.rate.unique, 2), "%")
    )
  )
  print(size_table, row.names = FALSE)
  cat("\n")
  
  if (family.name %in% c("binomial", "quasibinomial")) {
    cat("Sample Composition (Logistic Regression):\n")
  } else {
    cat("Sample Composition:\n")
  }
  print(y_summary, row.names = FALSE)
  cat("\n")
  
  ### Rare feature summary
  if (!is.null(object$rareFeature.index)) {
    cat("Subsample Rare-Feature Summary:\n")
    
    rare_names <- names(object$rare.count.plt)
    K <- length(rare_names)
    rare_table <- data.frame(
      Feature = rare_names,
      Subsample_count = as.numeric(object$rare.count.ssp),
      Subsample_rate = vapply(
        object$rare.count.ssp,
        function(x) safe_pct(x, n.ssp.actual),
        character(1)
      ),
      Full_count = as.numeric(object$full.rare.count),
      Full_rate = vapply(
        object$full.rare.count,
        function(x) safe_pct(x, N, digits = 4),
        character(1)
      ),
      Coverage_of_full = vapply(
        seq_len(K),
        function(k) safe_pct(object$rare.count.ssp[k], object$full.rare.count[k]),
        character(1)
      ),
      check.names = FALSE
    )
    print(rare_table, row.names = FALSE)
    cat("\n")
    
    if (!is.null(object$rare.counts.ssp) && length(object$rare.counts.ssp) > 0) {
      cat("Subsample Rare-Pattern Distribution:\n")
      pattern_table <- data.frame(
        Rare_features_in_row = names(object$rare.counts.ssp),
        Count = as.numeric(object$rare.counts.ssp),
        Rate = vapply(
          object$rare.counts.ssp,
          function(x) safe_pct(x, n.ssp.actual),
          character(1)
        ),
        check.names = FALSE
      )
      print(pattern_table, row.names = FALSE)
      cat("\n")
    }
  }
  
  
  ### Coefficients
  print_coef_block <- function(title, coef, cov) {
    cat(title, "\n\n", sep = "")
    se <- sqrt(diag(cov))
    coef_table <- data.frame(
      Estimate = round(coef, digits = 4),
      `Std. Error` = round(se, digits = 4),
      `z value` = round(coef / se, digits = 4),
      `Pr(>|z|)` = rF.format_p_values(
        2 * (1 - pnorm(abs(coef / se))),
        threshold = 0.0001
      ),
      check.names = FALSE
    )
    rownames(coef_table) <- names(coef)
    print(coef_table)
  }
  
  print_coef_block("Pilot Coefficients:", object$coef.plt, object$cov.plt)
  cat("\n")
  print_coef_block("Second-Step Coefficients:", object$coef.ssp, object$cov.ssp)
  cat("\n")
  print_coef_block("Combined-Union Coefficients:",
                   object$coef.cmb.union,
                   object$cov.cmb.union)
}
