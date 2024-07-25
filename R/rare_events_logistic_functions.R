###############################################################################
rare.coef.estimate <- function(X,
                               Y,
                               offset = NULL,
                               start = rep(0, ncol(X)),
                               weights = 1) {
  data <- as.data.frame(cbind(Y, X))
  formula <- as.formula(paste(colnames(data)[1], "~",
                              paste(colnames(data)[-1], collapse = "+"), "-1"))
  ## use '-1' to avoid adding intercept column again.
  design <- survey::svydesign(ids =  ~ 1,
                              weights =  ~ weights,
                              data = data)
  fit <- ifelse(is.null(offset),
                results <- survey::svyglm(formula = formula,
                                          design = design,
                                          start = start,
                                          family = quasibinomial(link="logit")),
                results <- survey::svyglm(formula = formula,
                                          design = design,
                                          start = start,
                                          offset = offset,
                                          family = quasibinomial(link="logit")))
  beta <- results$coefficients
  cov <- results$cov.unscaled
  pbeta <- as.vector(results$fitted.values)
  return (list(beta = beta,
              cov = cov,
              pbeta = pbeta
              )
         )
}
###############################################################################
rare.ddL <- function(X, P, w = 1){
  phi <- P * (1 - P)
  return(t(X) %*% (X * (phi * w)))
}
rare.dL.sq <- function(X, Y, P, w = 1){
  dL.sq <- (Y - P)^2
  return(t(X) %*% (X * (dL.sq * w)))
}
pbeta <- function(X, beta, offset = NA){
  ifelse(is.na(offset),
         p <- 1 / (1 + exp(-as.vector(X %*% beta))),
         p <- 1 / (1 + exp(-as.vector(X %*% beta) - offset)))
  return (p)
}
rare.calculate.nm <- function(X, Y, ddL.plt.correction, P.plt, criterion){
  if (criterion == "optA"){
    nm <- sqrt(rowSums((X %*% t(solve(ddL.plt.correction)))^2)) # norm
    nm <- P.plt * sqrt(1 - P.plt) * nm # numerator
  } else if (criterion == "optL"){
    nm <- sqrt(rowSums(X^2))
    nm <- P.plt * sqrt(1 - P.plt) * nm
  } else if (criterion == "LCC"){
    nm <- abs(Y - P.plt)
  }
  return(nm)
}
###############################################################################
rare.pilot.estimate <- function(X, Y, n.plt){
  N <- nrow(X)
  d <- ncol(X)
  N1 <- sum(Y)
  N0 <- N - N1
  ## half half sampling
  halfhalf.index.results <- halfhalf.index(N, Y, n.plt)
  index.plt <- halfhalf.index.results$index.plt
  n.plt.0 <- halfhalf.index.results$n.plt.0
  n.plt.1 <- halfhalf.index.results$n.plt.1
  x.plt <- X[index.plt, ]
  y.plt <- Y[index.plt]
  p.plt <- c(rep(1 / (2 * N0), n.plt.0), rep(1 / (2 * N1), n.plt.1))
  ## p: sampling probability
  ## P: Prob(Y=1 | X)
  ## unweighted likelihood estimation:
  results.plt <- rare.coef.estimate(X = x.plt, Y = y.plt)
  beta.plt <- results.plt$beta
  pbeta.plt <- pbeta(x.plt, beta.plt)
  ddL.plt <- rare.ddL(x.plt, pbeta.plt, 1 / n.plt)
  dL.sq.plt <- rare.dL.sq(x.plt, y.plt, pbeta.plt, 1 / n.plt ^ 2)
  ## correct intercept:
  beta.plt[1] <- beta.plt[1] - log(N0 / N1)
  P.plt <- pbeta(X, beta.plt)
  ddL.plt.correction <- rare.ddL(x.plt, P.plt[index.plt], 1 / n.plt)

  return(list(p.plt = p.plt,
              beta.plt = beta.plt,
              ddL.plt = ddL.plt,
              dL.sq.plt = dL.sq.plt,
              ddL.plt.correction = ddL.plt.correction,
              P.plt = P.plt,
              index.plt = index.plt
              )
         )
}
###############################################################################
rare.subsampling <- function(X,
                             Y,
                             n.ssp,
                             alpha,
                             b,
                             criterion,
                             likelihood,
                             p.plt,
                             ddL.plt.correction,
                             P.plt,
                             index.plt) {
  N <- nrow(X)
  N1 <- sum(Y)
  N0 <- N - N1
  n.plt <- length(index.plt)
  ## length(index.plt) might be smaller than n.plt, so here we reset n.plt.
  w.ssp <- offset <- NA
  nm <- rare.calculate.nm(X, Y, ddL.plt.correction, P.plt, criterion)
  ## Currently only Poisson sampling method has been implemented.
  if(criterion %in% c('optL', 'optA')){
    H <- quantile(nm, 1 - n.ssp / (b * N)) # if consider threshold
    nm[nm > H] <- H
    NPhi <- (N0 / N) * sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- n.ssp * ((1 - alpha) * nm / NPhi + alpha / N)
    index.ssp <- poisson.index(N, Y + (1 - Y) * p.ssp)
    p.ssp <- pmin(p.ssp[index.ssp], 1)
    ## calculate offset or weights:
    if (likelihood == 'logOddsCorrection') {
      offset <- -log(p.ssp)
    } else if (likelihood == 'weighted') {
      w.ssp <- 1 / (Y[index.ssp] + (1 - Y[index.ssp]) * p.ssp)
    }
  } else if (criterion == 'LCC'){
    dm <- sum(nm[index.plt] / p.plt) / n.plt
    p.ssp <- (n.ssp + N1) * nm / dm
    index.ssp <- poisson.index(N, p.ssp)
    p.ssp <- pmin(p.ssp[index.ssp], 1)
    ## calculate offset or weights:
    if (likelihood == 'logOddsCorrection') {
      nm.1 <- abs(1 - P.plt[index.ssp])
      nm.0 <- abs(P.plt[index.ssp])
      pi.1 <- pmin((n.ssp + N1) * nm.1 / dm, 1)
      pi.0 <- pmin((n.ssp + N1) * nm.0 / dm, 1)
      offset <- log(pi.1 / pi.0)
    } else if (likelihood == 'weighted') {
      w.ssp <- 1 / p.ssp
    }
  }
  return (list(index.ssp = index.ssp,
               w.ssp = w.ssp,
               offset = offset
               )
          )
}
###############################################################################
rare.subsample.estimate <- function(x.ssp,
                                    y.ssp,
                                    n.ssp,
                                    N,
                                    w.ssp,
                                    offset,
                                    beta.plt,
                                    likelihood) {
  if (likelihood == "weighted"){
    results.ssp <- rare.coef.estimate(x.ssp, y.ssp, weights = w.ssp)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta
    # P.ssp <- pbeta(x.ssp, beta.ssp) # as same as results.ssp$pbeta
    var.ssp <- results.ssp$cov
    ddL.ssp <- rare.ddL(x.ssp, P.ssp, w = w.ssp * n.ssp / N)
    dL.sq.ssp <- rare.dL.sq(x.ssp, y.ssp, P.ssp, w = w.ssp^2 * n.ssp^2 / N^2)
  } else if (likelihood == 'logOddsCorrection'){
    results.ssp <- rare.coef.estimate(X = x.ssp,
                                      Y = y.ssp,
                                      start = beta.plt,
                                      offset = offset)
    beta.ssp <- results.ssp$beta
    P.ssp <- results.ssp$pbeta
    # P.ssp <- pbeta(x.ssp, beta.ssp, offset)
    var.ssp <- results.ssp$cov
    ddL.ssp <- rare.ddL(x.ssp, P.ssp, w = 1 / n.ssp)
    dL.sq.ssp <- rare.dL.sq(x.ssp, y.ssp, P.ssp, w = 1 / n.ssp ^ 2)
  }
  return (list(beta.ssp = beta.ssp,
               ddL.ssp = ddL.ssp,
               dL.sq.ssp = dL.sq.ssp,
               var.ssp = var.ssp
               )
          )
}
###############################################################################
rare.combining <- function(ddL.plt,
                           ddL.ssp,
                           dL.sq.plt,
                           dL.sq.ssp,
                           n.plt,
                           n.ssp,
                           beta.plt,
                           beta.ssp) {
  ddL.plt <- n.plt * ddL.plt
  ddL.ssp <- n.ssp * ddL.ssp
  dL.sq.plt <- n.plt ^ 2 * dL.sq.plt
  dL.sq.ssp <- n.ssp ^ 2 * dL.sq.ssp
  MNsolve <- solve(ddL.plt + ddL.ssp)
  beta.cmb <- c(MNsolve %*% (ddL.plt %*% beta.plt + ddL.ssp %*% beta.ssp))
  var.cmb <- MNsolve %*% (dL.sq.plt + dL.sq.ssp) %*% MNsolve
  return(list(beta.cmb = beta.cmb,
              var.cmb = var.cmb
              )
  )
}
###############################################################################
#' relogit Main results summary
#'
#' @param object A list object output by the main function, which contains the
#' @param ... Additional arguments passed to the summary function.
#'  results of the estimation of the parameters, the estimation of the
#'  covariance matrix, subsample index, etc.
#'
#' @return A series of data.frame will be printed.
#' @export
#'
#' @examples
#' #TBD
summary.ssp.relogit <- function(object, ...) {
  coef <- object$beta
  se <- sqrt(diag(object$var))
  N <- object$N
  n.ssp.expect <- object$subsample.size.expect
  n.ssp.actual <- length(object$index)
  n.ssp.unique <- length(unique(object$index))
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
  coef_table <- data.frame(
    Estimate = round(coef, digits = 4),
    `Std. Error` = round(se, digits = 4),
    `z value` = round(coef / se, digits = 4),
    `Pr(>|z|)` = format.p.values(2 * (1 - pnorm(abs(coef / se))),
                                 threshold = 0.0001),
    check.names = FALSE
  )
  rownames(coef_table) <- names(coef)
  print(coef_table)
  # Add more summary information as needed
}
