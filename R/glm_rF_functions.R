
###############################################################################
rare_counts <- function(X, index, rareFeature.index, Y = NULL) {
  
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
balance_score <- function(rareFeature.index, 
                          X) {
  n <- nrow(X)
  p <- ncol(X)
  DN <- rep(1, p)
  prop_1 <- NULL
  
  
  ### return default value 
  if (is.null(rareFeature.index) || length(rareFeature.index) == 0) {
    return(list(
      bl = rep(1, n),
      DN = DN,
      rare.count = NULL
    ))
  }
  
  ### Rare features summary
  X_rf <- X[, rareFeature.index, drop = FALSE]
  Z_mean_rf <- colMeans(X_rf)
  DN[rareFeature.index] <- sqrt(Z_mean_rf)
  rare_count <- Z_mean_rf * n

  ### bl score formula:
  ### L1 norm bl_i = sum(|Z_ij - Z_bar_j| / Z_bar_j)
  ### equivalent formula: bl_i = d_rare + sum(Z_ij * V_j)
  
  d_rare <- length(rareFeature.index)
  V_vec <- (1 - 2 * Z_mean_rf) / Z_mean_rf
  bl <- d_rare + X_rf %*% V_vec
  bl <- as.vector(bl)

  ### a more direct but less efficient way is:
  # diff <- abs(sweep(X_rf, 2, Z_mean_rf, "-"))
  # weighted_diff <- sweep(diff, 2, Z_mean_rf, "/")
  # bl <- rowSums(weighted_diff)
  
  ### previous balance score:
  # diff <- abs(sweep(X_rf, 2, Z_mean_rf, "-"))
  # bl <- rowSums(diff)
  
  list(
    bl = bl,
    DN = DN,
    rare.count = rare_count
  )
}
###############################################################################
na_rate <- function(x) {
  sum(is.na(x)) / length(x)
}
###############################################################################
rF.glm.coef.estimate <- function(X,
                              Y,
                              start = rep(0, ncol(X)),
                              weights = 1,
                              family,
                              ...) {
  
  data <- as.data.frame(cbind(Y, X))
  
  ## use '-1' in the formula to avoid adding the intercept again.
  formula <- as.formula(paste(colnames(data)[1], "~",
                              paste(colnames(data)[-1], collapse = "+"), "-1"))
  
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  

  results <- survey::svyglm(formula,
                            design = design,
                            # start = start,
                            family = family,
                            ...)
  
  beta <- results$coefficients

  return(list(beta = beta))
}
###############################################################################
halfhalf.index <- function(N, Y, n.plt) {
  N1 <- sum(Y)
  N0 <- N - N1
  if (N0 < n.plt / 2 | N1 < n.plt / 2) {
    warning(paste("n.plt/2 exceeds the number of Y=1 or Y=0 in the full data.",
                  "All rare events will be drawn into the pilot sample.",
                  "Consider using rare.logistic.subsampling()."))
  }
  n.plt.0 <- min(N0, n.plt / 2)
  n.plt.1 <- min(N1, n.plt / 2)
  if (n.plt.0 == n.plt.1) {
    index.plt <- c(sample(which(Y == 0), n.plt.0),
                   sample(which(Y == 1), n.plt.1))
  } else if (n.plt.0 < n.plt.1) {
    index.plt <- c(which(Y == 0),
                   sample(which(Y == 1),
                          n.plt.1))
  } else if (n.plt.0 > n.plt.1) {
    index.plt <- c(sample(which(Y == 0),
                          n.plt.0),
                   which(Y == 1))
  }
  return(list(index.plt = index.plt, 
              n.plt.0 = n.plt.0, 
              n.plt.1 = n.plt.1
  )
  )
}
###############################################################################
random.index <- function (N, n, p = NULL) {
  ifelse(
    is.null(p),
    index <- sample(N, n, replace = TRUE), # faster than using specific prob
    index <- sample(N, n, replace = TRUE, prob = p)
  )
  return(as.vector(index))
}
poisson.index <- function (N, pi) {
  return(which(runif(N) <= pi))
}

###############################################################################
rF.calculate.nm <- function(X, Y, ddL.plt, linear.predictor, criterion, 
                              bl, DN = 1){
  
  if (criterion == "BL-Lopt"){
    nm <- sqrt(rowSums(X^2))
    nm <- bl * abs(Y - linear.predictor) * nm
  } else if (criterion == "Aopt"){ # original Aopt
    nm <- sqrt(rowSums((X %*% t(solve(ddL.plt)))^2))
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
###############################################################################
poi.prob.estimated <- function(nm, index.plt, n.ssp, N, b){
  ## compute approximated H and denominator.
  ## then get the estimated poisson sampling probabilities
  H <- quantile(nm, 1 - n.ssp / (b * N))
  nm <- pmin(nm, H)
  return(list(H = H,
              nm = nm))
}

###############################################################################
poi.prob.exact <- function(nm, n.ssp, N){
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
              nm = nm))
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
  
  ddL <- t(X) %*% (X * (dlinear.predictor * weights))
  return(ddL)
}

## square of the first derivative of log likelihood function
rF.dL.sq <- function (eta, X, Y, weights = 1, family) {
  variance <- family$variance
  linkinv  <- family$linkinv
  
  temp <- (Y - linkinv(eta))^2
  
  dL.sq <- t(X) %*% (X * (temp * weights))
  return(dL.sq)
}
###############################################################################
rF.pilot.estimate <- function(inputs, ...){
  X <- inputs$X
  Y <- inputs$Y
  n.plt <- inputs$n.plt
  family <- inputs$family
  N <- inputs$N
  bl <- inputs$bl
  DN <- inputs$DN
  linkinv  <- family$linkinv
  balance.plt <- inputs$balance.plt
  balance.Y <- inputs$balance.Y
  
  h.plt <- n.plt / N
  if (family[["family"]] %in% c('binomial', 'quasibinomial')){
    ## for binary response, three strategies: 
    ## uniform; balancing X; balancing both X and Y.
    if (balance.plt && length(bl) != 1L) {
      # bl.mean <- mean(bl)
      # varphi.plt <-  bl / bl.mean
      if(balance.Y == TRUE){
        # balancing both Y and X
        # will not consider truncation, since n.plt is supposed to be small.
        N1 <- sum(Y)
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
        index.plt <- poisson.index(N, pi.plt)
        pi.plt.N <- pi.plt  # dimension: N
        pi.plt <- pi.plt[index.plt]  # dimension: n.plt
        varphi.plt <- varphi.plt[index.plt]   # dimension: n.plt
        w.plt <- 1 / pi.plt
      } else { # only balancing X
        bl.mean <- mean(bl)
        varphi.plt <-  bl / bl.mean
        pi.plt <- pmin(varphi.plt * h.plt, 1) # dimension: N
        index.plt <- poisson.index(N, pi.plt)
        pi.plt.N <- pi.plt
        pi.plt <- pi.plt[index.plt]
        varphi.plt <- varphi.plt[index.plt]
        w.plt <- 1 / pi.plt
      }
    } else{ # do uniform sampling
      index.plt <- poisson.index(N, h.plt)
      pi.plt <- rep(n.plt/N, length(index.plt)) # dimension: n.plt
      pi.plt.N <- rep(n.plt/N, N)
      varphi.plt <- rep(1, length(index.plt)) # all 1
      w.plt <- 1 / pi.plt
    }
    
    n.plt <- length(index.plt) # actuall n.plt
    x.plt <- X[index.plt,] # dimension: n.plt
    Y.plt <- Y[index.plt] # dimension: n.plt
    ## weighted likelihood
    beta.plt <- rF.glm.coef.estimate(X = x.plt, Y = Y.plt, 
                                  weights = w.plt,
                                  family = family, ...)$beta
    linear.predictor.plt <- as.vector(x.plt %*% beta.plt) # dimension: n.plt
    
    # ddL.plt and dL.sq.plt: the pilot estimation of gradient^2 and information
    ddL.plt <- rF.ddL(eta = linear.predictor.plt,
                   X = x.plt,
                   weights = (1 / (varphi.plt*n.plt)),
                   family = family)
    dL.sq.plt <- rF.dL.sq(eta = linear.predictor.plt,
                       x.plt,
                       Y.plt,
                       weights = 1 / ((varphi.plt^2) * n.plt),
                       family = family)
    Lambda.plt <- 0 # placeholder for SWR
    linear.predictor <- linkinv(X %*% beta.plt) # dimension: N
    ddL.plt.inv <- solve(ddL.plt)
    cov.plt <- (1 / n.plt) * ddL.plt.inv %*% dL.sq.plt %*% ddL.plt.inv
  } else {
    ## for other families
  }
  return(list(pi.plt = pi.plt, 
              pi.plt.N = pi.plt.N, # dim:N
              varphi.plt = varphi.plt,
              beta.plt = beta.plt,
              ddL.plt = ddL.plt,
              dL.sq.plt = dL.sq.plt,
              Lambda.plt = Lambda.plt,
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
                        ddL.plt,
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
  criterion <- inputs$criterion
  likelihood <- inputs$likelihood
  sampling.method <- inputs$sampling.method
  rareFeature.index <- inputs$rareFeature.index
  bl <- inputs$bl
  DN <- inputs$DN
  balance.Y <- inputs$balance.Y
  
  h.ssp <- n.ssp / N
  
  ## When use poisson sampling method to draw pilot sample, 
  ## length(index.plt) might be different with original n.plt
  ## so here we reset n.plt.
  n.plt <- length(index.plt)
  w.ssp <- NA
  
  nm <- rF.calculate.nm(X, Y, ddL.plt, linear.predictor, criterion, 
                          bl, DN) # numerator. dimension: N
  
  if (balance.Y) {
    N0 <- N - sum(Y)
    poi.prob.results <- poi.prob.exact(nm[Y==0], n.ssp, N0)
    # poi.prob.results <- poi.prob.estimated(nm, index.plt, n.ssp, N, b)
    nm <- poi.prob.results$nm # dim: N0
    varphi.ssp <- nm / mean(nm)
    pi.ssp <- rep(NA, N)
    pi.ssp[Y==1] <- 1  # maybe not efficient
    pi.ssp[Y==0] <- h.ssp * varphi.ssp # [Y==0]
  } else {
    poi.prob.results <- poi.prob.exact(nm, n.ssp, N)
    # poi.prob.results <- poi.prob.estimated(nm, index.plt, n.ssp, N, b)
    nm <- poi.prob.results$nm # dim: N
    varphi.ssp <- nm / mean(nm)
    pi.ssp <- h.ssp * varphi.ssp
  }
  
  index.ssp <- poisson.index(N, pi.ssp)
  pi.ssp.N <- pi.ssp
  pi.ssp <- pi.ssp[index.ssp] # dim: actual n.ssp
  w.ssp <- 1 / pi.ssp # dim: actual n.ssp
  varphi.ssp = pi.ssp / h.ssp
  
  
  return(list(index.ssp = index.ssp,
              w.ssp = w.ssp,
              varphi.ssp = varphi.ssp,
              pi.ssp.N = pi.ssp.N
              )
         )
}
###############################################################################
rF.subsample.estimate <- function(inputs,
                               w.ssp,
                               varphi.ssp,
                               offset,
                               beta.plt,
                               index.ssp,
                               ...) {
  x.ssp <- inputs$X[index.ssp, ]
  y.ssp = inputs$Y[index.ssp]
  n.ssp <- length(index.ssp) # re-calculate
  N <- inputs$N
  sampling.method <- inputs$sampling.method
  likelihood <- inputs$likelihood
  family <- inputs$family
  DN <- inputs$DN
  
  
  results.ssp <- rF.glm.coef.estimate(x.ssp,
                                   y.ssp,
                                   weights = w.ssp,
                                   family = family,
                                   ...)
  beta.ssp <- results.ssp$beta
  linear.predictor.ssp <- as.vector(x.ssp %*% beta.ssp)
  
  if (sampling.method == 'poisson') {
    ddL.ssp <- rF.ddL(linear.predictor.ssp,
                   x.ssp,
                   weights = (1 / (varphi.ssp*n.ssp)),
                   family = family)
    dL.sq.ssp <- rF.dL.sq(linear.predictor.ssp,
                       x.ssp,
                       y.ssp,
                       weights = (1 / ((varphi.ssp^2) * n.ssp)),
                       family = family)
    Lambda.ssp <- 0 # placeholder
  } else if (sampling.method == "withReplacement") {
    ddL.ssp <- rF.ddL(linear.predictor.ssp,
                   x.ssp,
                   weights = w.ssp / (N * n.ssp),
                   family = family)
    dL.sq.ssp <- rF.dL.sq(linear.predictor.ssp,
                       x.ssp,
                       y.ssp,
                       weights = w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2),
                       family = family)
    c <- n.ssp / N
    Lambda.ssp <- c * rF.dL.sq(linear.predictor.ssp,
                            x.ssp,
                            y.ssp,
                            weights = w.ssp / (N * n.ssp ^ 2),
                            family = family)
  }
  ddL.ssp.inv <- solve(ddL.ssp)
  cov.ssp <- (1 / n.ssp) * ddL.ssp.inv %*% dL.sq.ssp %*% ddL.ssp.inv 
  return(list(beta.ssp = beta.ssp,
              ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
              cov.ssp = cov.ssp,
              Lambda.ssp = Lambda.ssp
              )
         )
}
###############################################################################
# combining <- function(ddL.plt,
#                       ddL.ssp,
#                       dL.sq.plt,
#                       dL.sq.ssp,
#                       Lambda.plt,
#                       Lambda.ssp,
#                       n.plt,
#                       n.ssp,
#                       beta.plt,
#                       beta.ssp) {
#   
#   Info.plt <- n.plt * ddL.plt
#   Info.ssp <- n.ssp * ddL.ssp
#   
#   Info.solve <- solve(Info.plt + Info.ssp)
#   beta.cmb <- c(Info.solve %*% (Info.plt %*% beta.plt + Info.ssp %*% beta.ssp))
#   
#   cov.cmb <- Info.solve %*% 
#     (n.plt * (dL.sq.plt + Lambda.plt) + 
#        n.ssp * (dL.sq.ssp + Lambda.ssp)) %*% 
#     Info.solve
#   
#   return(list(beta.cmb = beta.cmb,
#               cov.cmb = cov.cmb
#               )
#          )
# }
rF.combining <- function(inputs,
                      index.plt, # dim: n.plt
                      index.ssp, # dim: n.ssp
                      pi.plt, # dim: N
                      pi.ssp, # dim: N
                      ...) {
  X <- inputs$X # dim: N
  Y <- inputs$Y# dim: N
  combine <- inputs$combine
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
  cov.cmb.test <- results.cmb$cov
  linear.predictor.cmb <- as.vector(X.cmb %*% beta.cmb)
  
  ddL.cmb <- rF.ddL(linear.predictor.cmb,
                 X.cmb,
                 weights = 1 / ((n.cmb^2/N) * pi.cmb),
                 family = family)
  dL.sq.cmb <- rF.dL.sq(linear.predictor.cmb,
                     X.cmb,
                     Y.cmb,
                     weights = 1 / (n.cmb * (n.cmb/N)^2 * pi.cmb^2),
                     family = family)
  ddL.cmb.inv <- solve(ddL.cmb)
  cov.cmb <- (1 / n.cmb) * ddL.cmb.inv %*% dL.sq.cmb %*% ddL.cmb.inv  
  
  

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
  DN <- inputs$DN
  linkinv  <- family$linkinv
  balance.Y <- inputs$balance.Y
  criterion <- inputs$criterion
  
  h.uni <- n.uni / N # subsampling rate 
  
  if (family[["family"]] %in% c('binomial', 'quasibinomial')){
    ## for binary response, it optionally takes the rarity of Y into account.

    if (criterion == 'BL-Uni' && length(bl) != 1L) {
      bl.mean <- mean(bl)
      varphi.uni <-  bl / bl.mean # normalized to mean 1 
      pi.uni <- h.uni * varphi.uni # un-truncated sampling prob
      ## truncate those pi.uni > 1, while enforce the summation to be n
      if(balance.Y == TRUE){
        # do truncation for Y=0, while keeping all Y=1
        pi.uni.Y0 <- pi.uni[Y==0]
        g <- sum(pi.uni.Y0 > 1)
        if (g == 0) { # no truncation needed
          pi.uni <- Y + (1 - Y) * pi.uni # dim: N
          index.uni <- poisson.index(N, pi.uni)
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
          index.uni <- poisson.index(N, pi.uni)
          pi.uni <- pi.uni[index.uni]
          Y.uni <- Y[index.uni]
          varphi.uni <- (Y.uni + (1 - Y.uni) * pi.uni) / h.uni # dim: n.uni
          w.uni <- 1 / pi.uni
        }
      } else{ # balancing X only
        g <- sum(pi.uni > 1)
        if (g == 0) { # no truncation needed
          index.uni <- poisson.index(N, pi.uni)
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
          index.uni <- poisson.index(N, pi.uni)
          pi.uni <- pi.uni[index.uni]
          w.uni <- 1 / pi.uni
          varphi.uni <- pi.uni / h.uni
          Y.uni <- Y[index.uni]
        }
      }
    } else{ # do poisson sampling with uniform sampling prob 
      if(balance.Y == TRUE){
        pi.uni <- Y + (1-Y) * h.uni
        index.uni <- poisson.index(N, pi.uni)
        pi.uni <- pi.uni[index.uni] # dimension: n.uni
        varphi.uni <- pi.uni / h.uni
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
      } else {
        index.uni <- poisson.index(N, h.uni)
        pi.uni <- rep(h.uni, length(index.uni)) # dimension: n.uni
        varphi.uni <- rep(1, length(index.uni)) # all 1
        w.uni <- 1 / pi.uni
        Y.uni <- Y[index.uni]
      }
    }
    
    
    n.uni <- length(index.uni) # actuall n.uni
    X.uni <- X[index.uni,] # dimension: n.uni
    ## weighted likelihood
    beta.uni <- rF.glm.coef.estimate(X = X.uni, Y = Y.uni, 
                                  weights = w.uni,
                                  family = family, ...)$beta
    linear.predictor.uni <- as.vector(X.uni %*% beta.uni) # dimension: n.uni
    
    # ddL.uni and dL.sq.uni: the pilot estimation of gradient^2 and information
    ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                   X = X.uni,
                   weights = (1 / (varphi.uni*n.uni)),
                   family = family)
    dL.sq.uni <- rF.dL.sq(eta = linear.predictor.uni,
                       X.uni,
                       Y.uni,
                       weights = 1 / ((varphi.uni^2) * n.uni),
                       family = family)
    Lambda.uni <- 0 # placeholder for SWR
    linear.predictor <- NA # dimension: N
    ddL.uni.inv <- solve(ddL.uni)
    cov.uni <- (1 / n.uni) * ddL.uni.inv %*% dL.sq.uni %*% ddL.uni.inv
  } else {
    ## for other families
    # do poisson sampling with uniform sampling prob 
    index.uni <- poisson.index(N, n.uni / N)
    pi.uni <- rep(n.uni/N, length(index.uni)) # dimension: n.uni
    varphi.uni <- rep(1, length(index.uni)) # all 1
    w.uni <- 1 / pi.uni
    
    n.uni <- length(index.uni) # actuall n.uni
    X.uni <- X[index.uni,] # dimension: n.uni
    Y.uni <- Y[index.uni] # dimension: n.uni
    ## weighted likelihood
    beta.uni <- rF.glm.coef.estimate(X = X.uni, Y = Y.uni, 
                                  weights = w.uni,
                                  family = family, ...)$beta
    linear.predictor.uni <- as.vector(X.uni %*% beta.uni) # dimension: n.uni
    
    # ddL.uni and dL.sq.uni: the pilot estimation of gradient^2 and information
    ddL.uni <- rF.ddL(eta = linear.predictor.uni,
                   X = X.uni,
                   weights = (1 / (varphi.uni*n.uni)),
                   family = family)
    dL.sq.uni <- rF.dL.sq(eta = linear.predictor.uni,
                       X.uni,
                       Y.uni,
                       weights = 1 / ((varphi.uni^2) * n.uni),
                       family = family)
    Lambda.uni <- 0 # placeholder for SWR
    linear.predictor <- NA # dimension: N
    ddL.uni.inv <- solve(ddL.uni)
    cov.uni <- (1 / n.uni) * ddL.uni.inv %*% dL.sq.uni %*% ddL.uni.inv
  }
  return(list(pi.uni = pi.uni,
              varphi.uni = varphi.uni,
              beta.uni = beta.uni,
              ddL.uni = ddL.uni,
              dL.sq.uni = dL.sq.uni,
              Lambda.uni = Lambda.uni,
              linear.predictor = linear.predictor,
              index.uni = index.uni,
              cov.uni = cov.uni
              )
         )
}
###############################################################################
glm.control <- function(alpha = 0, b = 2, ...)
{
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("sampling probability weight 'alpha' must between [0, 1]")
  if(!is.numeric(b) || b < 0)
    stop("sampling probability threshold 'b' must > 0")
  list(alpha = alpha, b = b)
}
###############################################################################
format_p_values <- function(p.values, threshold = 0.0001) {
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
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index.ssp)
  n.ssp.unique <- length(unique(object$index.ssp))
  subsample.rate.expect <- (n.ssp.expect / N) * 100
  subsample.rate.actual <- (n.ssp.actual / N) * 100
  subsample.rate.unique <- (n.ssp.unique / N) * 100
  
  ### header
  cat("Model Summary\n")
  cat("\nCall:\n\n")
  print(object$model.call)
  cat("\n")
  
  ### Subsample info
  cat("Subsample Summary:\n")
  size_table <- data.frame(
    Variable = c(
      'Total Sample Size',
      'Expected Subsample Size',
      'Actual Subsample Size',
      'Unique Subsample Size',
      'Expected Subsample Rate',
      'Actual Subsample Rate',
      'Unique Subsample Rate'
    ),
    Value = c(
      N,
      n.ssp.expect,
      n.ssp.actual,
      n.ssp.unique,
      paste0(round(subsample.rate.expect, 2), "%"),
      paste0(round(subsample.rate.actual, 2), "%"),
      paste0(round(subsample.rate.unique, 2), "%")
    )
  )
  print(size_table, row.names = FALSE)
  cat("\n")
  
  ### Rare feature summary
  if (!is.null(object$rareFeature.index)) {
    cat("Rare Features Summary (pilot, subsample, full):\n")
    
    rare_names <- names(object$rare.count.plt)
    K <- length(rare_names)
    
    n_plt <- length(object$index.plt)
    n_ssp <- length(object$index.ssp)
    n_full <- N
    
    ### avoid 0/0 in BL-Uni or Uni
    safe_prop <- function(count, denom) {
      if (denom == 0) return("NA")
      paste0(round(count / denom * 100, 2), "%")
    }
    
    ### iilot block
    plt_block <- rbind(
      size_pilot   = rep(n_plt, K),
      count1_pilot = object$rare.count.plt,
      prop_pilot   = sapply(object$rare.count.plt,
                            function(x) safe_prop(x, n_plt)),
      coverage_pilot = sapply(seq_len(K),
                              function(k) safe_prop(object$rare.count.plt[k],
                                                    object$full.rare.count[k]))
    )
    colnames(plt_block) <- rare_names
    
    # bubsample block
    ssp_block <- rbind(
      size_subsample = rep(n_ssp, K),
      count1_ssp     = object$rare.count.ssp,
      prop_ssp       = sapply(object$rare.count.ssp,
                              function(x) safe_prop(x, n_ssp)),
      coverage_ssp   = sapply(seq_len(K),
                              function(k) safe_prop(object$rare.count.ssp[k],
                                                    object$full.rare.count[k]))
    )
    colnames(ssp_block) <- rare_names
    
    ### full-data block
    full_block <- rbind(
      size_full   = rep(n_full, K),
      count1_full = object$full.rare.count,
      rare_rate   = paste0(round(object$full.rare.count / N * 100, 4), "%")
    )
    colnames(full_block) <- rare_names
    
    ### combine with blank rows for readability
    blank_row <- matrix("", nrow = 1, ncol = K)
    colnames(blank_row) <- rare_names
    rownames(blank_row) <- " "
    
    final_table <- rbind(plt_block,
                         blank_row,
                         ssp_block,
                         blank_row,
                         full_block)
    
    print(final_table)
    cat("\n")
  }
  
  
  ### Coefficients
  cat("Coefficients:\n\n")
  coef_table <- data.frame(
    Estimate   = round(coef, digits = 4),
    `Std. Error` = round(se, digits = 4),
    `z value`  = round(coef / se, digits = 4),
    `Pr(>|z|)` = format_p_values(
      2 * (1 - pnorm(abs(coef / se))),
      threshold = 0.0001
    ),
    check.names = FALSE
  )
  rownames(coef_table) <- names(coef)
  print(coef_table)
}
