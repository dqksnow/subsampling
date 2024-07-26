set.seed(1)
N <- 2 * 1e4
beta0 <- c(-5, -rep(0.7, 6))
d <- length(beta0) - 1
X <- matrix(0, N, d)
corr <- 0.5
sigmax <- corr ^ abs(outer(1:d, 1:d, "-"))
sigmax <- sigmax / 4
X <- MASS::mvrnorm(n = N, mu = rep(0, d), Sigma = sigmax)
Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ .
n.plt <- 200
n.ssp <- 1000


expect_silent(subsampling.results <- ssp.relogit(formula = formula,
                                             data = data,
                                             n.plt = n.plt,
                                             n.ssp = n.ssp), 
              info = "It should run without errors on valid input.")
# expect_true(inherits(subsampling.results, "list"),
#             info = "Output should be a list.")
# expect_true(inherits(subsampling.results, "ssp.relogit"), 
#             info = "Output should be of class 'ssp.relogit'")

# Cleanup
rm(list = ls())
gc()