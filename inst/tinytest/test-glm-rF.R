set.seed(1)

N <- 5e3
d_rare <- 3
d_cont <- 2
p_rare <- c(0.01, 0.02, 0.05)
beta0 <- c(0.5, rep(0.5, d_rare), rep(0.5, d_cont))
corr <- 0.5
sigmax <- matrix(corr, d_cont, d_cont) + diag(1 - corr, d_cont)
X_cont <- mvtnorm::rmvnorm(N, rep(0, d_cont), sigmax)
Z <- do.call(cbind, lapply(seq_along(p_rare), function(i) {
  rbinom(N, 1, p_rare[i])
}))
X <- cbind(Z, X_cont)
d <- ncol(X)
rareFeature.index <- 1:d_rare
P <- 1 / (1 + exp(-(beta0[1] + X %*% beta0[-1])))
Y <- as.integer(rbinom(N, 1, P))

colnames(X) <- paste0("X", 1:d)
data <- data.frame(Y = Y, X)
formula <- Y ~ .
n.plt <- 300
n.ssp <- 600
family <- "quasibinomial"


expect_silent(
  rf1 <- ssp.glm.rF(
    formula = formula,
    data = data,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family,
    criterion = "BL-Uni",
    sampling.method = "poisson",
    likelihood = "weighted",
    balance.plt = TRUE,
    balance.Y = FALSE,
    rareFeature.index = rareFeature.index
  ),
  info = "ssp.glm.rF should run without error for BL-Uni"
)

expect_true(inherits(rf1, "list"),
            info = "Output should be a list.")

expect_true(inherits(rf1, "ssp.glm.rF"),
            info = "Output should have class 'ssp.glm.rF'.")


criteria_to_test <- c("Lopt", "Aopt", "R-Lopt", "BL-Lopt", "Uni")

for (crit in criteria_to_test) {
  expect_silent(
    ssp.glm.rF(
      formula = formula,
      data = data,
      n.plt = n.plt,
      n.ssp = n.ssp,
      family = family,
      criterion = crit,
      sampling.method = "poisson",
      likelihood = "weighted",
      balance.plt = TRUE,
      balance.Y = FALSE,
      rareFeature.index = rareFeature.index
    ),
    info = paste("ssp.glm.rF should run for criterion =", crit)
  )
}


expect_silent(
  rf_subset <- ssp.glm.rF(
    formula = formula,
    data = data,
    subset = 1:(N/2),
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family,
    criterion = "BL-Uni",
    sampling.method = "poisson",
    likelihood = "weighted",
    rareFeature.index = rareFeature.index
  ),
  info = "ssp.glm.rF should run with subset argument"
)


## ------------------------------------------------------------
## Test passing arguments through ...
## ------------------------------------------------------------

expect_silent(
  rf_ctrl <- ssp.glm.rF(
    formula = formula,
    data = data,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family,
    criterion = "Lopt",
    sampling.method = "poisson",
    likelihood = "weighted",
    maxit = 30,
    rareFeature.index = rareFeature.index
  ),
  info = "ssp.glm.rF should accept ... and pass to svyglm()"
)


## ------------------------------------------------------------
## Test contrasts with a factor variable
## ------------------------------------------------------------

data2 <- data
data2$F1 <- sample(c("A", "B", "C"), N, replace = TRUE)
formula2 <- Y ~ . + F1

expect_silent(
  rf_contrast <- ssp.glm.rF(
    formula = formula2,
    data = data2,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family,
    criterion = "BL-Uni",
    sampling.method = "poisson",
    likelihood = "weighted",
    contrasts = list(F1 = "contr.sum"),
    rareFeature.index = rareFeature.index
  ),
  info = "ssp.glm.rF should run with contrasts argument"
)


##############################################################
# Cleanup
##############################################################
rm(list = ls())
gc()
