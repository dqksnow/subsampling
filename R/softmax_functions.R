###############################################################################
softmax.coef.estimate <- function(X, y, weights = NULL){
  ifelse(is.null(weights),
         fit <- nnet::multinom(y ~ X - 1),
         fit <- nnet::multinom(y ~ X - 1, weights = weights))
  return(list(beta = t(coef(fit)), # return a d*K matrix
              P1 = fit$fitted.values)) # return a N*(K+1) matrix
}
###############################################################################
pbeta.multi <- function(X, beta){ # beta is K*d vector
  exp(X %*% beta - matrixStats::rowLogSumExps(cbind(0, X %*% beta)))
}
###############################################################################
softmax.ddL <- function(X, P, p, K, d, scale){
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
softmax.plt.estimate <- function(X, Y, Y.matrix, n.plt, N, K, d, criterion){
  index.plt <- random.index(N, n.plt)
  p.plt <- rep(1 / N, n.plt) # plt sampling probability
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  results <- softmax.coef.estimate(x.plt, y.plt)
  beta.plt <- results$beta
  P1.plt <-  results$P1 # n.plt*(K+1) matrix
  ddL.plt <- softmax_ddL_cpp(X = x.plt, P = P1.plt[, -1], p = rep(1, n.plt), K,
                             d, scale = n.plt)
  dL.sq.plt <- softmax_dL_sq_cpp(X = x.plt, Y_matrix = Y.matrix[index.plt, ],
                                 P = P1.plt[, -1], p = rep(1, n.plt), K = K,
                                 d = d, scale = n.plt^2)
  c <- n.plt / N
  Lambda.plt <- c * softmax_dL_sq_cpp(X = x.plt, 
                                      Y_matrix = Y.matrix[index.plt, ],
                                      P = P1.plt[, -1], p = rep(1, n.plt),
                                      K = K, d = d, scale = n.plt^2)
  cov.plt <- solve(ddL.plt) %*% (dL.sq.plt + Lambda.plt) %*% solve(ddL.plt)
  # cov.plt <- solve(ddL.plt) %*% dL.sq.plt %*% solve(ddL.plt)
  if (criterion == "MSPE") {
    Omega.plt <- softmax_Omega_cpp(x.plt, P1 = P1.plt, p = p.plt, K, d,
                                   scale = N*n.plt)
  } else {
    Omega.plt <- NA
  }
  P.plt <- pbeta.multi(X, beta.plt) # N*K
  return(
    list(
      p.plt = p.plt,
      beta.plt = beta.plt,
      P.plt = P.plt,
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
softmax.calculate.nm <- function(X, Y, ddL.plt, Omega.plt, sixi, G,
 criterion, constraint){
  if (criterion == "OptA") {
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi %*% solve(ddL.plt))^2))
    } else if (constraint == "summation"){
      temp <- sixi %*% solve(ddL.plt) %*% t(G)
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "OptL"){
    if (constraint == "baseline"){
      nm <- sqrt(rowSums((sixi)^2))
    } else if (constraint == "summation"){
      tempG <- t(G %*% solve(t(G) %*% G))
      temp <- sixi %*% tempG
      nm <- sqrt(rowSums(temp^2))
    }
  } else if (criterion == "MSPE") {
    temp <- sixi %*% solve(ddL.plt) %*% expm::sqrtm(Omega.plt)
    nm <- sqrt(rowSums(temp^2))
  }
  return(nm)
}
###############################################################################
softmax.subsampling <- function(X, Y, Y.matrix, G, n.ssp, N, K, d, alpha,
                                b, criterion, estimate.method, sampling.method,
                                constraint, p.plt, ddL.plt, P.plt, Omega.plt,
                                index.plt) {
  sixi <- (Y.matrix-P.plt)[, rep(seq(K), each = d)] * X[,rep(seq(d), K)]
  nm <- softmax.calculate.nm(X, Y, ddL.plt, Omega.plt, sixi,
                             G, criterion, constraint)
  if (sampling.method == "WithReplacement"){
    dm <- sum(nm) # denominator
    p.ssp <- (1 - alpha) * nm / dm + alpha / N # N*1
    index.ssp <- random.index(N, n.ssp, p.ssp)
  } else if (sampling.method == "Poisson"){
    H <- quantile(nm[index.plt], 1-n.ssp/(b*N)) # threshold
    nm[nm > H] <- H
    dm <- (sum(nm[index.plt] /p.plt) / n.plt) * (n.plt/(n.plt - K*d))
    p.ssp <- n.ssp * ((1 - alpha) * nm / dm + alpha / N)
    index.ssp <- poisson.index(N, p.ssp)
  }
  return(list(index.ssp = index.ssp,
              p.ssp = p.ssp))
}
###############################################################################
softmax.subsample.estimate <- function(x.ssp, y.ssp, y.matrix.ssp, n.ssp,
                                       index.ssp, p.ssp, beta.plt,
                                       sampling.method, estimate.method,
                                       N, K, d) {
  results <- softmax.coef.estimate(x.ssp, y.ssp, weights = 1 / p.ssp)
  beta.ssp <- results$beta
  P.ssp <- (results$P1)[, -1] # n.ssp*K matrix
  if (sampling.method == "WithReplacement") {
    ddL.ssp <- softmax_ddL_cpp(X = x.ssp, P = P.ssp, p = p.ssp,
                               K = K, d = d, scale = N*n.ssp)
    dL.sq.ssp <- softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                   P = P.ssp, p = p.ssp,
                                   K = K, d = d, scale = (N^2)*(n.ssp^2))
    c <- n.ssp / N
    Lambda.ssp <- c * softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                        P = P.ssp, p = sqrt(p.ssp),
                                        K = K, d = d, scale = N*(n.ssp^2))
  } else if (sampling.method == 'Poisson') {
    ddL.ssp <- softmax_ddL_cpp(X = x.ssp, P = P.ssp, p = p.ssp,
                               K = K, d = d, scale = N)
    dL.sq.ssp <- softmax_dL_sq_cpp(X = x.ssp, Y_matrix = y.matrix.ssp,
                                   P = P.ssp, p = p.ssp,
                                   K = K, d = d, scale = N^2)
    Lambda.ssp <- 0
  }
  cov.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp)
  return(list(beta.ssp = beta.ssp,
              ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
              Lambda.ssp = Lambda.ssp,
              cov.ssp = cov.ssp)
  )
}
###############################################################################
softmax.combining <- function(ddL.plt, ddL.ssp, dL.sq.plt, dL.sq.ssp,
                              Lambda.plt, Lambda.ssp, n.plt, n.ssp, beta.plt, 
                              beta.ssp, X, N, K, d) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  dL.sq.plt <- n.plt^2 * dL.sq.plt
  dL.sq.ssp <- n.ssp^2 * dL.sq.ssp
  Lambda.plt <- n.plt^2 * Lambda.plt
  Lambda.ssp <- n.ssp^2 * Lambda.ssp

  ddL.inv <- solve(ddL.plt + ddL.ssp)
  beta.cmb <- ddL.inv %*% (ddL.plt %*% c(beta.plt) + ddL.ssp %*% c(beta.ssp))
  beta.cmb <- matrix(beta.cmb, nrow = d)
  cov.cmb <- ddL.inv %*% (dL.sq.plt + Lambda.plt + dL.sq.ssp + Lambda.ssp) %*% 
              ddL.inv

  P.cmb <- matrix(NA, nrow = N, ncol = K+1)
  P.cmb[, -1] <- pbeta.multi(X, beta.cmb)
  P.cmb[, 1] <- 1 - rowSums(P.cmb[, -1])

  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb,
              P.cmb = P.cmb)
  )
}
###############################################################################
#' Softmax Main results summary
#'
#' @param object A list object output by the main function, which contains the
#'  results of the estimation of the parameters, the estimation of the
#'  variance, subsample size, etc.
#'
#' @return A series of data.frame will be printed.
#' @export
#'
#' @examples
#' #TBD

softmax.summary <- function(object) {
  dimension <- dim(object$beta.cmb)
  d <- dimension[1]
  K <- dimension[2]
  coef <- object$beta.cmb
  se <- matrix(sqrt(diag(object$cov.cmb)), nrow = d, ncol=K)
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
  cat("Std. Errors::\n")
  cat("\n")
  print(se)
  # Add more summary information as needed
}

