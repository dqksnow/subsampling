###############################################################################
softmax.coef.estimate <- function(X, y, weights = NULL, offsets = NULL, 
                                  ...){
  print(paste('Message from nnet::multinom: '))
    if (is.null(offsets)) {
    ifelse(is.null(weights),
           fit <- nnet::multinom(y ~ X - 1, ...),
           fit <- nnet::multinom(y ~ X - 1, weights = weights, ...))
  } else {
    fit <- nnet::multinom(y ~ X + offset(offsets) - 1, ...)
  }
  return(list(beta = t(coef(fit)), # return a d*K matrix
              P1 = fit$fitted.values)) # return a N*(K+1) matrix
}
###############################################################################
pbeta.multi <- function(X, beta, offsets = NULL){ # beta is a Kd*1 vector
  if (is.null(offsets)) {
    P1 <- exp(cbind(0, X %*% beta)) # K+1 class
    P1 <- P1 / rowSums(P1)
  } else {
    P1 <- exp(cbind(0, X %*% beta) + offsets)
    P1 <- P1 / rowSums(P1)
  }
}
###############################################################################
softmax.ddL <- function(X, P, p, K, d, scale){
  ## This is R version. cpp version is used.
  ddL <- matrix(0, nrow = K * d, ncol = K * d)
  X.t <- t(X)
  for (i in 1:K){
    for (j in i:K){
      ifelse(j==i,
             XPj <- X * ((P[, i] - P[, i]*P[, j]) / p),
             XPj <- X * ((- P[, i]*P[, j]) / p)
      )
      ddL[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        ddL[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XPj
    }
  }
  ddL <- ddL / scale
  return(ddL)
}
###############################################################################
softmax.dL.sq <- function(X, Y.matrix, P, p, K, d, scale){ 
  ## This is R version. cpp version is used.
  S <- Y.matrix - P
  dL <- matrix(0, nrow = K * d, ncol = K * d)
  X.t <- t(X)
  p.sq <- p^2
  for (i in 1:K){
    for (j in i:K){
      XSij <- X * (S[, i] * S[, j] / p.sq)
      dL[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        dL[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XSij
    }
  }
  dL <- dL / scale
  return(dL)
}
###############################################################################
softmax.Omega <- function(X, P1, p, K, d, scale){
  ## This is R version. cpp version is used.
  Omega <- matrix(0, nrow = K * d, ncol = K * d)
  P.sq <- rowSums(P1^2) # N * 1
  P0 <- P1[, -1]
  X.t <- t(X)
  for (i in 1:K){
    for (j in i:K){
      if (j == i){
        XP <- X * (((P.sq + 1 - 2 * P0[, i]) * P0[, i]^2) / p)
      } else {
        XP <- X * ((P.sq - P0[, i] - P0[, j]) * (P0[, i] * P0[, j])  / p)
      }
      Omega[(1+(i-1)*d):(i*d), (1+(j-1)*d):(j*d)] <-
        Omega[(1+(j-1)*d):(j*d), (1+(i-1)*d):(i*d)] <- X.t %*% XP
    }
  }
  Omega <- Omega / scale
  return(Omega)
}
###############################################################################
softmax.plt.estimate <- function(inputs, ...){
  
  X <- inputs$X
  Y <- inputs$Y
  Y.matrix <- inputs$Y.matrix
  n.plt <- inputs$n.plt
  N <- inputs$N
  K <- inputs$K
  d <- inputs$d
  criterion <- inputs$criterion
  
  index.plt <- random.index(N, n.plt)
  p.plt <- rep(1 / N, n.plt) # pilot sampling probability, uniform
  x.plt <- X[index.plt, ]
  y.plt <- Y[index.plt]
  results <- softmax.coef.estimate(x.plt, y.plt, ...)
  beta.plt <- results$beta
  P1.plt <-  results$P1 # n.plt*(K+1) matrix
  ddL.plt <- softmax_ddL_cpp(X = x.plt, P = P1.plt[, -1], p = rep(1, n.plt), K,
                             d, scale = n.plt)
  
  dL.sq.plt <- softmax_dL_sq_cpp(X = x.plt, Y_matrix = Y.matrix[index.plt, ],
                                 P = P1.plt[, -1], p = rep(1, n.plt), K = K,
                                 d = d, scale = n.plt)
  c <- n.plt / N
  Lambda.plt <- c * dL.sq.plt
  cov.plt <- solve(ddL.plt) %*% (dL.sq.plt + Lambda.plt) %*% 
             solve(ddL.plt) / n.plt

  if (criterion == "MSPE") {
    Omega.plt <- softmax_Omega_cpp(x.plt, P1 = P1.plt, p = p.plt, K, d,
                                   scale = N*n.plt)
  } else {
    Omega.plt <- NA
  }
  P1.plt.N <- pbeta.multi(X, beta.plt) # N*K
  return(
    list(
      p.plt = p.plt,
      beta.plt = beta.plt,
      P1.plt = P1.plt.N,
      ddL.plt = ddL.plt,
      dL.sq.plt = dL.sq.plt,
      index.plt = index.plt,
      cov.plt = cov.plt,
      Omega.plt = Omega.plt,
      Lambda.plt = Lambda.plt
    )
  )
}
###############################################################################
softmax.calculate.nm <- function(ddL.plt, Omega.plt, sixi, G,
                                 criterion, constraint){
  if (criterion == "optA") {
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi %*% solve(ddL.plt))^2))
    } else if (constraint == "summation"){
      temp <- sixi %*% solve(ddL.plt) %*% t(G)
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "optL"){
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi)^2))
    } else if (constraint == "summation"){
      temp <- sixi %*% t(G %*% solve(t(G) %*% G))
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "MSPE") {
    temp <- sixi %*% solve(ddL.plt) %*% expm::sqrtm(Omega.plt)
    nm <- sqrt(rowSums(temp^2))
  }

  return(nm)
}
###############################################################################
softmax.calculate.offsets <- function(P.ssp, X.ssp,
                                      H = NA,
                                      dm = NA,
                                      G, ddL.plt, Omega.plt,
                                      criterion, 
                                      constraint = NA,
                                      alpha, N,
                                      index.plt,
                                      n.plt
                                      ){
  # only compute offsets for subsample, not for full data.
  P.ssp.dim <- dim(P.ssp) # P.ssp = P1.plt[index.ssp, ], dim= n.ssp * (K+1)
  n.ssp <- P.ssp.dim[1]
  K <- P.ssp.dim[2] - 1
  d <- ncol(X.ssp)
  offsets <- matrix(NA, nrow = n.ssp, ncol = (K+1))

  if (criterion == 'LUC') {
    q <- pmax(apply(P.ssp, 1, max), 0.5)
    etaK <- P.ssp == q
    gamma <- N / n.ssp
    m <- gamma >= 2 * q
    piK1 <- 2 * m * (q + (1 - 2 * q) * etaK) / gamma + (1 - m) *
      (1 + (1 - gamma) / (gamma - q) * etaK)
    offsets <- log(piK1)
  } else { # for optL, optA, MSPE
    P0.ssp <- P.ssp[, -1]
    for (k in 1:(K+1)) { # when all subsample have y[k]=1
      if (k == 1) {
        P.temp <- -P0.ssp 
        sixi <- P.temp[, rep(seq(K), each = d)] * X.ssp[, rep(seq(d), K)]
      } else {
        P.temp <- - P0.ssp
        P.temp[, (k-1)] <- 1 + P.temp[, (k-1)]
        sixi <- P.temp[, rep(seq(K), each = d)] * X.ssp[, rep(seq(d), K)]
      }

      if (criterion == "optA") {
        if (constraint == "baseline"){
          temp <- sixi %*% solve(ddL.plt)
          nm <- sqrt(rowSums(temp^2))
        } else if (constraint == "summation"){
          temp <- sixi %*% solve(ddL.plt) %*% t(G)
          nm <- sqrt(rowSums(temp^2))
        }
      } else if (criterion == "optL"){
        if (constraint == "baseline"){
          nm <- sqrt(rowSums(sixi^2))
        } else if (constraint == "summation"){
          temp <- sixi %*% t(G %*% solve(t(G) %*% G))
          nm <- sqrt(rowSums(temp^2))
        }
 
      } else if (criterion == "MSPE") {
        temp <- sixi %*% solve(ddL.plt) %*% expm::sqrtm(Omega.plt)
        nm <- sqrt(rowSums(temp^2))
      }

      # threshold H and denominator dm are estimated on pilot sample
      nm[nm > H] <- H

      offsets[, k] <- log(
        pmin(n.ssp * ((1 - alpha) * nm / dm + alpha / N), 1)
        )
    }
  }
  return(offsets)
}
###############################################################################
softmax.subsampling <- function(inputs,
                                p.plt, ddL.plt, P1.plt, Omega.plt, index.plt) {
  X <- inputs$X
  Y.matrix <- inputs$Y.matrix
  G <- inputs$G
  n.plt <- inputs$n.plt
  n.ssp <- inputs$n.ssp
  N <- inputs$N
  K <- inputs$K
  d <- inputs$d
  criterion <- inputs$criterion
  likelihood <- inputs$likelihood
  sampling.method <- inputs$sampling.method
  constraint <- inputs$constraint
  control <- inputs$control
  alpha <- control$alpha
  b <- control$b
  
  if (criterion == 'LUC') {
    Y.matrix.1 <- cbind(1-rowSums(Y.matrix), Y.matrix) #add class 1 to Y.matrix
    q <- pmax(apply(P1.plt, 1, max), 0.5)
    etaK <- rowSums(Y.matrix.1 * P1.plt) == q
    gamma <- N / n.ssp
    m <- (gamma >= 2 * q)
    piK1 <- 2 * m * (q + (1 - 2 * q) * etaK) / gamma + (1 - m) *
      (1 + (1 - gamma) / (gamma - q) * etaK)
    p.ssp <- pmin(n.ssp * piK1 / sum(piK1), 1) # poisson sampling probability
    index.ssp <- poisson.index(N, p.ssp)
    if (sampling.method == "withReplacement") {
      stop("The 'LUC' criterion + 'withReplacement' sampling method
           has not been implemented yet.")
    }
    if (likelihood == 'MSCLE') {
      offsets <- softmax.calculate.offsets(P1.plt[index.ssp,],
                                           X[index.ssp,],
                                           G = G,
                                           ddL.plt = ddL.plt,
                                           Omega.plt = Omega.plt,
                                           criterion = criterion,
                                           alpha = alpha,
                                           N = N,
                                           index.plt = index.plt,
                                           n.plt = n.plt
                                           )
    } else if (likelihood == 'weighted') {
      offsets <- NA
    }
  } else { # for optA, optL, MSPE
    sixi <- (Y.matrix-P1.plt[, -1])[, rep(seq(K), each = d)] * 
      X[,rep(seq(d), K)]
    nm <- softmax.calculate.nm(ddL.plt, Omega.plt, sixi,
                               G, criterion, constraint)
    if (sampling.method == "withReplacement"){
      dm <- sum(nm) # denominator
      p.ssp <- (1 - alpha) * nm / dm + alpha / N # N*1
      index.ssp <- random.index(N, n.ssp, p.ssp)
      if (likelihood == 'MSCLE') {
        stop("The 'MSCLE' likelihood + 'withReplacement' sampling method
           has not been implemented yet.")
      }
      offsets <- NA
    } else if (sampling.method == "poisson"){
      # threshold H is estimated by pilot sample
      H <- quantile(nm[index.plt], 1-n.ssp/(b*N))
      nm[nm > H] <- H
      # denominator dm is estimated by pilot sample
      dm <- (sum(nm[index.plt] /p.plt) / n.plt) * (n.plt/(n.plt - K*d))
      p.ssp <- pmin(n.ssp * ((1 - alpha) * nm / dm + alpha / N), 1)
      index.ssp <- poisson.index(N, p.ssp)
      # calculate offsets
      if (likelihood == 'MSCLE') {
        offsets <- softmax.calculate.offsets(P1.plt[index.ssp, ],
                                           X[index.ssp, ],
                                           H = H,
                                           dm = dm,
                                           G = G,
                                           ddL.plt = ddL.plt,
                                           Omega.plt = Omega.plt,
                                           criterion = criterion,
                                           constraint = constraint,
                                           alpha = alpha,
                                           N = N,
                                           index.plt = index.plt,
                                           n.plt = n.plt
                                           )
      } else if (likelihood == 'weighted') {
        offsets <- NA
      }
    }
  }
  return(list(index.ssp = index.ssp,
              p.ssp = p.ssp,
              offsets = offsets))
}
###############################################################################
softmax.subsample.estimate <- function(inputs,
                                       index.ssp, p.ssp, offsets, beta.plt,
                                       ...) {
  x.ssp <- inputs$X[index.ssp, ]
  y.ssp <- inputs$Y[index.ssp]
  y.matrix.ssp <- inputs$Y.matrix[index.ssp, ]
  n.ssp <- length(index.ssp)
  N <- inputs$N
  K <- inputs$K
  d <- inputs$d
  likelihood <- inputs$likelihood
  sampling.method <- inputs$sampling.method

  if (likelihood == "weighted"){
    results <- softmax.coef.estimate(x.ssp, y.ssp, weights = 1 / p.ssp,
                                     ...)
    beta.ssp <- results$beta
    P.ssp <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "withReplacement") {
      ddL.ssp <- softmax_ddL_cpp(X = x.ssp, P = P.ssp, p = p.ssp,
                                 K = K, d = d, scale = N*n.ssp)
      dL.sq.ssp <- softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                     P = P.ssp, p = p.ssp,
                                     K = K, d = d, scale = (N^2)*n.ssp)
      c <- n.ssp / N
      Lambda.ssp <- c * softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                          P = P.ssp, p = sqrt(p.ssp),
                                          K = K, d = d, scale = N*n.ssp)
    } else if (sampling.method == 'poisson') {
      ddL.ssp <- softmax_ddL_cpp(X = x.ssp, P = P.ssp, p = p.ssp,
                                 K = K, d = d, scale = N)
      dL.sq.ssp <- softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                     P = P.ssp, p = p.ssp,
                                     K = K, d = d, scale = N^2 / n.ssp)
      Lambda.ssp <- 0
    }
  } else if (likelihood == 'MSCLE'){
    results <- softmax.coef.estimate(x.ssp, y.ssp, offsets = offsets,
                                     ...)
    beta.ssp <- results$beta
    # P.ssp <- pbeta.multi(x.ssp, beta.ssp, offsets)[, -1] # same
    P.ssp <- (results$P1)[, -1] # n.ssp*K matrix
    if (sampling.method == "withReplacement") {
      stop("The 'MSCLE' likelihood + 'withReplacement' sampling method
           has not been implemented yet.")
    } else if (sampling.method == 'poisson') {
      ddL.ssp <- softmax_ddL_cpp(X = x.ssp, P = P.ssp, p = rep(1, n.ssp),
                                 K = K, d = d, scale = n.ssp)
      dL.sq.ssp <- softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                     P = P.ssp, p = rep(1, n.ssp),
                                     K = K, d = d, scale = n.ssp)
      Lambda.ssp <- 0
    }
  }
  cov.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp) * 
    (1 / n.ssp)
  return(list(beta.ssp = beta.ssp,
              P.ssp = P.ssp,
              ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
              Lambda.ssp = Lambda.ssp,
              cov.ssp = cov.ssp)
  )
}
###############################################################################
softmax.combining <- function(inputs, n.ssp,
                              ddL.plt, ddL.ssp, dL.sq.plt, dL.sq.ssp,
                              Lambda.plt, Lambda.ssp, beta.plt, beta.ssp) {
  n.plt <- inputs$n.plt
  X <- inputs$X
  N <- inputs$N
  K <- inputs$K
  d <- inputs$d
  
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  dL.sq.plt <- n.plt * dL.sq.plt
  dL.sq.ssp <- n.ssp * dL.sq.ssp
  Lambda.plt <- n.plt * Lambda.plt
  Lambda.ssp <- n.ssp * Lambda.ssp

  ddL.inv <- solve(ddL.plt + ddL.ssp)
  beta.cmb <- ddL.inv %*% (ddL.plt %*% c(beta.plt) + ddL.ssp %*% c(beta.ssp))
  beta.cmb <- matrix(beta.cmb, nrow = d)

  
  cov.cmb <- ddL.inv %*% (dL.sq.plt + Lambda.plt + dL.sq.ssp + Lambda.ssp) %*% 
              ddL.inv

  P.cmb <- matrix(NA, nrow = N, ncol = (K+1))
  P.cmb <- pbeta.multi(X, beta.cmb)

  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              P.cmb = P.cmb)
  )
}
###############################################################################
softmax.control <- function(alpha = 0, b = 2, ...){
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("sampling probability weight 'alpha' must between [0, 1]")
  if(!is.numeric(b) || b < 0)
    stop("sampling probability threshold 'b' must > 0")
  list(alpha = alpha, b = b)
}
###############################################################################
#' @export
summary.ssp.softmax <- function(object, ...) {
  dimension <- dim(object$coef)
  d <- dimension[1]
  K <- dimension[2]
  coef <- object$coef
  se <- matrix(sqrt(diag(object$cov)), nrow = d, ncol=K)
  rownames(se) <- rownames(coef) 
  N <- object$N
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index.ssp)
  n.ssp.unique <- length(unique(object$index.ssp))
  subsample.rate.expect <- (n.ssp.expect / N) * 100
  subsample.rate.actual <- (n.ssp.actual / N) * 100
  subsample.rate.unique <- (n.ssp.unique / N) * 100
  cat("Model Summary\n\n")
  cat("\nCall:\n")
  cat("\n")
  print(object$model.call)
  cat("\n")
  cat("Subsample Size:\n")
  size_table <- data.frame(
    'Variable' = c(
      'Total Sample Size',
      'Expected Subsample Size',
      'Actual Subsample Size',
      'Unique Subsample Size',
      'Expected Subample Rate',
      'Actual Subample Rate',
      'Unique Subample Rate'
    ),
    'Value' = c(
      N,
      n.ssp.expect,
      n.ssp.actual,
      n.ssp.unique,
      paste0(subsample.rate.expect, "%"),
      paste0(subsample.rate.actual, "%"),
      paste0(subsample.rate.unique, "%")
    )
  )
  colnames(size_table) <- NULL
  rownames(size_table) <- NULL
  print(size_table)
  cat("\n")
  cat("Coefficients:\n")
  cat("\n")
  print(coef)
  cat("\n")
  cat("Std. Errors:\n")
  cat("\n")
  print(se)
  # Add more summary information as needed
}

