glm.coef.estimate <- function(X,
                              Y,
                              offset = NULL,
                              start = rep(0, ncol(X)),
                              weights = 1,
                              family,
                              ...) {

  data <- as.data.frame(cbind(Y, X))
  formula <- as.formula(paste(colnames(data)[1], "~",
                   paste(colnames(data)[-1], collapse = "+"), "-1"))
  ## use '-1' to avoid adding the intercept again.
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  ifelse(all(is.null(offset)),
         results <- survey::svyglm(as.formula(formula),
                                   design = design,
                                   # start = start,
                                   family = family,
                                   ...),
         results <- survey::svyglm(as.formula(formula),
                                   design = design,
                                   # start = start,
                                   offset = offset,
                                   family = family,
                                   ...)
         )
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
random.index <- function (N, n, p = NULL) {
  ifelse(
    is.null(p),
    index <- sample(N, n, replace = TRUE),
    index <- sample(N, n, replace = TRUE, prob = p)
  )
  return(as.vector(index))
}
poisson.index <- function (N, pi) {
  return(which(runif(N) <= pi))
}
###############################################################################
calculate.offset <- function (X,
                              N,
                              d.psi,
                              alpha,
                              ddL.plt.correction,
                              H = NULL,
                              NPhi = NULL,
                              n.ssp = NULL,
                              criterion,
                              sampling.method) {
  # only compute offsets for subsample, not for full data.
  if (criterion == "optA") {
    norm <- sqrt(rowSums((X %*% t(solve(ddL.plt.correction)))^2))
    nm.1 <- abs(1 - d.psi) * norm
    nm.0 <- abs(d.psi) * norm
  } else if (criterion == "optL") {
    norm <- sqrt(rowSums(X^2))
    nm.1 <- abs(1 - d.psi) * norm
    nm.0 <- abs(d.psi) * norm
  } else if (criterion == "LCC") {
    nm.1 <- abs(1 - d.psi)
    nm.0 <- abs(d.psi)
  }
  
  # threshold H is estimated on pilot sample
  nm.1[nm.1 > H] <- H
  nm.0[nm.0 > H] <- H
  
  # denominator NPhi is estimated on pilot sample
  if (sampling.method == 'withReplacement') {
    stop("Currently only the 'logOddsCorrection' likelihood with
         'poisson' sampling method has been implemented.")
  } else if (sampling.method == 'poisson') {
    pi.1 <- pmin(n.ssp * ((1 - alpha) * nm.1 / NPhi + alpha / N), 1)
    pi.0 <- pmin(n.ssp * ((1 - alpha) * nm.0 / NPhi + alpha / N), 1)
  }
  offset <- log(pi.1 / pi.0)
  return(offset)
}
###############################################################################
calculate.nm <- function(X, Y, ddL.plt.correction, d.psi, criterion){
  if (criterion == "optA"){
    nm <- sqrt(rowSums((X %*% t(solve(ddL.plt.correction)))^2))
    nm <- abs(Y - d.psi) * nm # numerator
  } else if (criterion == "optL"){
    nm <- sqrt(rowSums(X^2))
    nm <- abs(Y - d.psi) * nm
  } else if (criterion == "LCC"){
    nm <- abs(Y - d.psi)
  }
  return(nm)
}
###############################################################################
## second derivative of log likelihood function
ddL <- function (eta, X, weights = 1, offset = NULL, family) {
  # linkinv(eta):
  # (eta is linear predictor plus offset(if have))
  # for binomial = exp(eta)/(1+exp(eta))
  # poisson() = exp(eta)
  # Gamma(link = "inverse") = 1/eta
  # variance(mu):
  # for logit = mu(1 - mu)
  # for poisson() = mu
  # for Gamma(link = "inverse") = mu^2
  variance <- family$variance
  linkinv  <- family$linkinv
  
  if (all(is.null(offset))) {
    dd.psi <- variance(linkinv(eta))
  } else {
    dd.psi <- variance(linkinv(eta + offset))
  }
  ddL <- t(X) %*% (X * (dd.psi * weights))

  return(ddL)
}
## square of the first derivative of log likelihood function
dL.sq <- function (eta, X, Y, weights = 1, offset = NULL, family) {
  variance <- family$variance
  linkinv  <- family$linkinv
  if (all(is.null(offset))) {
    temp <- (Y - linkinv(eta))^2
  } else {
    temp <- (Y - linkinv(eta + offset))^2
  }
  dL.sq <- t(X) %*% (X * (temp * weights))
  return(dL.sq)
}
###############################################################################
pilot.estimate <- function(inputs, ...){
  X <- inputs$X
  Y <- inputs$Y
  n.plt <- inputs$n.plt
  family <-  inputs$family
  N <- inputs$N
  family <- inputs$family
  variance <- family$variance
  linkinv  <- family$linkinv

  if (family[["family"]] %in% c('binomial', 'quasibinomial')){
    N1 <- sum(Y)
    N0 <- N - N1
    ## This is case control sampling with replacement for binary Y logistic
    ## regression.
    ## We can also use uniform sampling with rep, half half sampling 
    ## or poisson sampling.
    p.plt <- ifelse(Y == 1, 1/(2*N1), 1/(2*N0))
    index.plt <- random.index(N, n.plt, p = p.plt)
    x.plt <- X[index.plt,]
    y.plt <- Y[index.plt]
    p.plt <- p.plt[index.plt]
    ## weighted likelihood
    beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, weights = 1 / p.plt,
                                  family = family, ...)$beta
    linear.predictor.plt <- as.vector(x.plt %*% beta.plt)
    
    ddL.plt <- ddL.plt.correction <- ddL(eta = linear.predictor.plt,
                                         X = x.plt,
                                         weights = (1 / (p.plt*N*n.plt)),
                                         family = family)
    
    dL.sq.plt <- dL.sq(eta = linear.predictor.plt,
                       x.plt,
                       y.plt,
                       weights = (1 / (p.plt*N*n.plt)^2),
                       family = family)
    c <- n.plt / N
    Lambda.plt <- c * dL.sq(eta = linear.predictor.plt,
                            x.plt,
                            y.plt,
                            weights = (1 / (p.plt*N*n.plt^2)),
                            family = family)
    d.psi <- linkinv(X %*% beta.plt) # N dimension
  } else {
    ## This is uniform sampling with replacement.
    index.plt <- random.index(N, n.plt)
    x.plt <- X[index.plt,]
    y.plt <- Y[index.plt]
    p.plt <- rep(1 / N, n.plt)
    beta.plt <- glm.coef.estimate(X = x.plt, Y = y.plt, family = family,
                                  ...)$beta
    linear.predictor.plt <- as.vector(x.plt %*% beta.plt)
    ddL.plt <- ddL.plt.correction <- ddL(eta = linear.predictor.plt, 
                                         x.plt, weights = 1 / n.plt,
                                         family = family)
    dL.sq.plt <- dL.sq(eta = linear.predictor.plt,
                       x.plt, y.plt, weights = 1 / n.plt^2, family = family)
    c <- n.plt / N
    Lambda.plt <- c * dL.sq.plt
    d.psi <- linkinv(X %*% beta.plt)
  }
  return(list(p.plt = p.plt,
              beta.plt = beta.plt,
              ddL.plt = ddL.plt,
              dL.sq.plt = dL.sq.plt,
              ddL.plt.correction = ddL.plt.correction,
              Lambda.plt = Lambda.plt,
              d.psi = d.psi,
              index.plt = index.plt
              )
         )
}
###############################################################################
subsampling <- function(inputs,
                        p.plt,
                        ddL.plt.correction,
                        d.psi,
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
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt)
  
  ## If use half half sampling, length(index.plt) might be smaller than 
  ## n.plt, so here we reset n.plt.
  
  w.ssp <- offset <- NA
  nm <- calculate.nm(X, Y, ddL.plt.correction, d.psi, criterion) # numerator
  if (sampling.method == "withReplacement"){
    dm <- sum(nm) # denominator
    p.ssp <- (1 - alpha) * nm / dm + alpha / N
    index.ssp <- random.index(N, n.ssp, p.ssp)
    if (likelihood == 'logOddsCorrection') {
      stop("Currently only the 'logOddsCorrection' estimate method with
         'withReplacement' sampling method has not been implemented.")
    } else if (likelihood == 'weighted') {
      w.ssp <- 1 / p.ssp[index.ssp]
    }
  } else if (sampling.method == "poisson"){
    H <- quantile(nm[index.plt], 1 - n.ssp / (b * N)) 
    # threshold H is estimated on pilot sample
    nm[nm > H] <- H
    NPhi <- sum(nm[index.plt] / p.plt) / n.plt
    # denominator NPhi is estimated on pilot sample
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    index.ssp <- poisson.index(N, p.ssp)
    if (likelihood == 'logOddsCorrection') {
      offset <- calculate.offset(X = X[index.ssp, ],
                                 N = N,
                                 d.psi = d.psi[index.ssp],
                                 alpha = alpha,
                                 ddL.plt.correction = ddL.plt.correction,
                                 criterion = criterion,
                                 sampling.method = sampling.method,
                                 NPhi = NPhi,
                                 n.ssp = n.ssp,
                                 H = H)
    } else if (likelihood == 'weighted') {
      w.ssp <- 1 / pmin(p.ssp[index.ssp], 1)
    }
  }
  return(list(index.ssp = index.ssp,
              offset = offset,
              w.ssp = w.ssp
              )
         )
}
###############################################################################
subsample.estimate <- function(inputs,
                               w.ssp,
                               offset,
                               beta.plt,
                               index.ssp,
                               ...) {
  x.ssp <- inputs$X[index.ssp, ]
  y.ssp = inputs$Y[index.ssp]
  n.ssp <- inputs$n.ssp
  N <- inputs$N
  sampling.method <- inputs$sampling.method
  likelihood <- inputs$likelihood
  family <- inputs$family
  
  if (likelihood == "weighted") {
    results.ssp <- glm.coef.estimate(x.ssp,
                                     y.ssp,
                                     weights = w.ssp,
                                     family = family,
                                     ...)
    beta.ssp <- results.ssp$beta
    linear.predictor.ssp <- as.vector(x.ssp %*% beta.ssp)
    if (sampling.method == 'poisson') {
      ddL.ssp <- ddL(linear.predictor.ssp,
                     x.ssp,
                     weights = w.ssp / N,
                     family = family)
      dL.sq.ssp <- dL.sq(linear.predictor.ssp,
                         x.ssp,
                         y.ssp,
                         weights = w.ssp ^ 2 / N ^ 2,
                         family = family)

      Lambda.ssp <- 0 # placeholder
    } else if (sampling.method == "withReplacement") {
      ddL.ssp <- ddL(linear.predictor.ssp,
                     x.ssp,
                     weights = w.ssp / (N * n.ssp),
                     family = family)
      dL.sq.ssp <- dL.sq(linear.predictor.ssp,
                         x.ssp,
                         y.ssp,
                         weights = w.ssp ^ 2 / (N ^ 2 * n.ssp ^ 2),
                         family = family)
      c <- n.ssp / N
      Lambda.ssp <- c * dL.sq(linear.predictor.ssp,
                              x.ssp,
                              y.ssp,
                              weights = w.ssp / (N * n.ssp ^ 2),
                              family = family)
    }
  } else if (likelihood == 'logOddsCorrection') {
    if (!(family[["family"]]  %in% c('binomial', 'quasibinomial'))){
      stop("Currently 'logOddsCorrection' likelihood can only work for logistic
           regression.")
    }
    results.ssp <- glm.coef.estimate(X = x.ssp,
                                     Y = y.ssp,
                                     start = beta.plt,
                                     offset = offset,
                                     family = family,
                                     ...)
    beta.ssp <- results.ssp$beta
    linear.predictor.ssp <- as.vector(x.ssp %*% beta.ssp)
    ddL.ssp <- ddL(linear.predictor.ssp,
                   x.ssp,
                   weights = 1 / n.ssp,
                   offset = offset,
                   family = family)
    dL.sq.ssp <- dL.sq(linear.predictor.ssp,
                       x.ssp,
                       y.ssp,
                       weights = 1 / n.ssp ^ 2,
                       offset = offset,
                       family = family)
    Lambda.ssp <- 0 # placeholder
  }
  cov.ssp <- solve(ddL.ssp) %*% (dL.sq.ssp + Lambda.ssp) %*% solve(ddL.ssp)

  return(list(beta.ssp = beta.ssp,
              ddL.ssp = ddL.ssp,
              dL.sq.ssp = dL.sq.ssp,
              cov.ssp = cov.ssp,
              Lambda.ssp = Lambda.ssp
              )
         )
}
###############################################################################
combining <- function(ddL.plt,
                      ddL.ssp,
                      dL.sq.plt,
                      dL.sq.ssp,
                      Lambda.plt,
                      Lambda.ssp,
                      n.plt,
                      n.ssp,
                      beta.plt,
                      beta.ssp) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  dL.sq.plt <- n.plt ^ 2 * dL.sq.plt
  dL.sq.ssp <- n.ssp ^ 2 * dL.sq.ssp
  Lambda.plt <- n.plt ^ 2 * Lambda.plt
  Lambda.ssp <- n.ssp ^ 2 * Lambda.ssp
  MNsolve <- solve(ddL.plt + ddL.ssp)
  beta.cmb <- c(MNsolve %*% (ddL.plt %*% beta.plt + ddL.ssp %*% beta.ssp))
  cov.cmb <- MNsolve %*% 
    (dL.sq.plt + Lambda.plt + dL.sq.ssp + Lambda.ssp) %*% MNsolve
  return(list(beta.cmb = beta.cmb,
              cov.cmb = cov.cmb
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
summary.ssp.glm <- function(object, ...) {
  coef <- object$coef
  se <- sqrt(diag(object$cov))
  N <- object$N
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index)
  n.ssp.unique <- length(unique(object$index))
  subsample.rate.expect <- (n.ssp.expect / N) * 100
  subsample.rate.actual <- (n.ssp.actual / N) * 100
  subsample.rate.unique <- (n.ssp.unique / N) * 100
  cat("Model Summary\n")
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
  coef_table <- data.frame(
    Estimate = round(coef, digits = 4),
    `Std. Error` = round(se, digits = 4),
    `z value` = round(coef / se, digits = 4),
    `Pr(>|z|)` = format_p_values(2 * (1 - pnorm(abs(coef / se))),
                                 threshold = 0.0001),
    check.names = FALSE
  )
  rownames(coef_table) <- names(coef)
  print(coef_table)
  # Add more summary information as needed
}
