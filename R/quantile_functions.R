###############################################################################
quantile.plt.estimation <- function(X, Y, tau, N, n.plt){
  index.plt <- sample(N, n.plt)
  p.plt <- rep(1 / N, n.plt) # plt sampling probability
  x.plt <- X[index.plt,]
  y.plt <- Y[index.plt]
  results <- quantreg::rq(y.plt ~ x.plt - 1, tau=tau)
  # results <- quantreg::rq(y.plt ~ x.plt - 1, tau=tau, weights = 1/p.plt)
  beta.plt <- results$coefficients
  Ie.full <- (results$residuals < 0)
  msg <- results$message
  return(
    list(
      # p.plt = p.plt,
      beta.plt = beta.plt,
      Ie.full = Ie.full,
      index.plt = index.plt
    )
  )
}
###############################################################################
quantile.iter.sampling.par <- function(j) { #this is used for parLapply
  index.ssp <- sample(1:N, n.ssp, replace = TRUE, prob = p.ssp)
  x.ssp <- X[index.ssp, ]
  y.ssp <- Y[index.ssp]
  ## currently we only consider Lopt and weighted likelihood
  fit <- quantreg::rq(y.ssp ~ x.ssp - 1, tau=tau, weights=1 / p.ssp[index.ssp],
                      method="pfn")
  return(list(beta.ssp = fit$coefficients,
              index.ssp <- index.ssp
  ))
}
###############################################################################
quantile.iter.sampling <- function(X, Y, N, n.ssp, p.ssp, tau) {
  index.ssp <- sample(1:N, n.ssp, replace = TRUE, prob = p.ssp)
  x.ssp <- X[index.ssp, ]
  y.ssp <- Y[index.ssp]
  ## currently we only consider Lopt and weighted likelihood
  fit <- quantreg::rq(y.ssp ~ x.ssp - 1, tau=tau, weights=1 / p.ssp[index.ssp],
                      method="pfn")
  return(list(beta.ssp = fit$coefficients,
              index.ssp <- index.ssp
  ))
}
###############################################################################
quantile.ssp.estimation <- function(X,
                                    Y,
                                    tau,
                                    n.ssp,
                                    B,
                                    Ie.full = NA,
                                    estimate.method,
                                    parallel) {
  N <- nrow(X)
  p <- ncol(X)
  
  if (estimate.method == "Uniform"){
    p.ssp <- rep(1, N)
  } else if (estimate.method %in% c("Weighted")) {
    # p.ssp <- sqrt((tau - Ie.full)^2 * rowSums(X^2)) # slower
    p.ssp <- abs(tau - Ie.full) * sqrt(rowSums(X^2))
    p.ssp <- p.ssp / sum(p.ssp)
  }
  
  Betas.ssp <- matrix(NA, nrow = p, ncol = B)
  Index.ssp <- matrix(NA, nrow = n.ssp, ncol = B)
  
  if (parallel == TRUE) {
    cl <- parallel::makeCluster(detectCores()) # how many CPU cores are called
    parallel::clusterExport(cl=cl,
                  varlist=c('N', 'n.ssp', 'p.ssp', 'X', 'Y', 'tau'),
                  envir=environment())
    results <- parallel::parLapply(cl, 1:B, quantile.iter.sampling.par)
    parallel::stopCluster(cl)
    for(i in 1:B){
      Betas.ssp[, i] <- results[[i]][[1]]
      Index.ssp[, i] <- results[[i]][[2]]
    }
  } else if (parallel == FALSE) {
    for(i in 1:B){
      index.ssp <- sample(1:N, n.ssp, replace = TRUE, prob = p.ssp)
      fit <- quantreg::rq(Y[index.ssp] ~ X[index.ssp, ] - 1,
                          tau=tau, weights=1 / p.ssp[index.ssp])
      Betas.ssp[, i] <- fit$coefficients
      Index.ssp[, i] <- index.ssp
    }
  }
  
  beta.ssp.mean <- rowMeans(Betas.ssp)
  beta.ssp.centered <- Betas.ssp - beta.ssp.mean
  mean.outer.prod <- beta.ssp.centered %*% t(beta.ssp.centered)
  r.ef <- ifelse(estimate.method == "Uniform",
                 1 - (n.ssp*B-1)/N/2,
                 1 - sum(p.ssp^2) * (n.ssp*B - 1) / 2) # could be negative
  est.cov.ssp <- mean.outer.prod / (r.ef * B*(B-1))
  return(list(Betas.ssp = Betas.ssp,
              beta.ssp = beta.ssp.mean,
              est.cov.ssp = est.cov.ssp,
              index.ssp = Index.ssp
  ))
}
###############################################################################
