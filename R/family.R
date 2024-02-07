binomial.expand <- function (){
  psi <- function(eta, offset = NULL) {
    ifelse(is.null(offset),
           p <- log(1 + exp(eta)),
           p <- log(1 + exp(eta + offset)))
    return(as.vector(p))
  }
  d.psi <- function(eta, offset = NULL) {
    ifelse(is.null(offset),
           p <- 1 - 1 / (1 + exp(eta)),
           p <- 1 - 1 / (1 + exp(eta + offset)))
    return(as.vector(p))
  }
  dd.psi <- function(eta, offset = NULL) {
    ifelse(is.null(offset),
           p <- exp(eta) / (1 + exp(eta)) ^ 2,
           p <- exp(eta + offset) / (1 + exp(eta + offset)) ^ 2)
    return(as.vector(p))
  }
  u <- function(eta) as.vector(eta)
  d.u <- function(eta) 1
  dd.u <- function(eta) 0
  structure(list(family.name = "binomial",
                 canonical = TRUE,
                 psi = psi,
                 d.psi = d.psi,
                 dd.psi = dd.psi,
                 u = u,
                 d.u = d.u,
                 dd.u = dd.u),
            class = "family")
}
###############################################################################
poisson.expand <- function (){
  psi <- function(eta, offset = NULL) exp(eta)
  d.psi <- function(eta, offset = NULL) exp(eta)
  dd.psi <- function(eta, offset = NULL) exp(eta)
  u <- function(eta) eta
  d.u <- function(eta) 1
  dd.u <- function(eta) 0
  structure(list(family.name = "poisson",
                 canonical = TRUE,
                 psi = psi,
                 d.psi = d.psi,
                 dd.psi = dd.psi,
                 u = u,
                 d.u = d.u,
                 dd.u = dd.u),
            class = "family")
}
###############################################################################
gamma.expand <- function (){
  psi <- function(eta, offset = NULL) -log(eta)
  d.psi <- function(eta, offset = NULL) -1/eta
  dd.psi <- function(eta, offset = NULL) 1/eta^2
  u <- function(eta) -eta
  d.u <- function(eta) -1
  dd.u <- function(eta) 0
  structure(list(family.name = "gamma",
                 canonical = TRUE,
                 psi = psi,
                 d.psi = d.psi,
                 dd.psi = dd.psi,
                 u = u,
                 d.u = d.u,
                 dd.u = dd.u),
            class = "family")
}
