set.seed(1)
N <- 2e4 # toy example
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
                ssp.quantreg(formula = formula,
                             data = data,
                             tau = tau,
                             n.plt = n.plt,
                             n.ssp = n.ssp,
                             B = B,
                             boot = TRUE,
                             criterion = 'optL',
                             sampling.method = 'withReplacement',
                             likelihood = 'weighted'), 
              info = "It should run without errors on valid input.")

expect_true(inherits(optL.results, "list"), info = "Output should be a list.")
expect_true(inherits(optL.results, "ssp.quantreg"), 
            info = "Output should be of class 'ssp.quantreg.'")

expect_equivalent(length(optL.results$index), 
                  B, 
                  info = "Subsamples should be divided into B lists.")

expect_warning(ssp.quantreg(formula = formula,
                            data = data,
                            tau = tau,
                            n.plt = n.plt,
                            n.ssp = 1000,
                            B = B,
                            boot = TRUE,
                            criterion = 'optL',
                            sampling.method = 'withReplacement',
                            likelihood = 'weighted'))

expect_silent(optL.results <- 
                ssp.quantreg(formula = formula,
                             data = data,
                             tau = tau,
                             n.plt = n.plt,
                             n.ssp = n.ssp,
                             B = B,
                             boot = FALSE,
                             criterion = 'optL',
                             sampling.method = 'withReplacement',
                             likelihood = 'weighted'), 
              info = "It should run without errors on valid input.")

expect_equivalent(length(optL.results$index), 
                  n.ssp*B, 
                  info = "When boot=F, Subsamples should not be divided into 
                  groups.")

expect_silent(optL.results <- 
                ssp.quantreg(formula = formula,
                             data = data,
                             subset = c(1:(N/2)), 
                             tau = tau,
                             n.plt = n.plt,
                             n.ssp = n.ssp,
                             B = B,
                             boot = TRUE,
                             criterion = 'optL',
                             sampling.method = 'withReplacement',
                             likelihood = 'weighted'),
              info = "It should run without errors when use subset argument.")

# Cleanup
rm(list = ls())
gc()