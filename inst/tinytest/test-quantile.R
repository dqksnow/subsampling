add_quantile <- function(a, b) {
  return(a + b)
}
expect_equal(add_quantile(1, 1), 2,
             info = "1+1=2")

set.seed(1)
N <- 1e4
B <- 5
tau <- 0.75
beta.true <- rep(1, 7)
d <- length(beta.true) - 1
corr  <- 0.5
sigmax  <- matrix(0, d, d)
for (i in 1:d) for (j in 1:d) sigmax[i, j] <- corr^(abs(i-j))
X <- MASS::mvrnorm(N, rep(0, d), sigmax)
err <- rnorm(N, 0, 1) - qnorm(tau)
Y <- beta.true[1] + X %*% beta.true[-1] + err * rowMeans(abs(X))
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ X
n.plt <- 100
n.ssp <- 100

expect_silent(optL.results <- 
                subsampling.quantile(formula,
                                     data,
                                     tau = tau,
                                     n.plt = n.plt,
                                     n.ssp = n.ssp,
                                     B = B,
                                     boot = TRUE,
                                     criterion = 'OptL',
                                     sampling.method = 'WithReplacement',
                                     likelihood = 'Weighted'), 
              info = "It should run without errors on valid input")

expect_coef <- c(0.9239319, 0.9841376, 1.0964170, 0.9649661, 1.0107086,
                 1.0012242, 0.9452326)

expect_true(inherits(optL.results, "list"), info = "Output should be a list")
expect_true(inherits(optL.results, "subsampling.quantile"), 
            info = "Output should be of class 'subsampling.quantile'")


expect_true(all(c("beta", "index") %in% names(optL.results)), 
            info = "Output list should contain 'beta' and 'index'")

expect_equivalent(optL.results$beta, expect_coef, 
                  tolerance = 1e-1, 
                  info = "Coefficients should match expected values")

# Cleanup
rm(list = ls())
gc()