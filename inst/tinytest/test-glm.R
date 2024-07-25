set.seed(1)
N <- 1e4
beta0 <- rep(-0.5, 7)
d <- length(beta0) - 1
X <- matrix(0, N, d)
generate_rexp <- function(x) x <- rexp(N, rate = 2)
X <- apply(X, 2, generate_rexp)
Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
print(paste('N: ', N))
print(paste('sum(Y): ', sum(Y)))
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ .
n.plt <- 500
n.ssp <- 1000
subsampling.results <- ssp.glm(formula, data, n.plt, n.ssp,
                               family = 'binomial', criterion = "optL", sampling.method = 'poisson',
                               likelihood = "logOddsCorrection")
summary(subsampling.results)



expect_silent(subsampling.results <- ssp.glm(formula,
                                             data,
                                             n.plt,
                                             n.ssp,
                                             family = 'binomial', 
                                             criterion = "optL",
                                             sampling.method = 'poisson',
                                             likelihood = "logOddsCorrection"), 
              info = "It should run without errors on valid input.")
expect_true(inherits(subsampling.results, "list"),
            info = "Output should be a list.")
expect_true(inherits(subsampling.results, "ssp.glm"), 
            info = "Output should be of class 'ssp.quantreg.'")

# Cleanup
rm(list = ls())
gc()