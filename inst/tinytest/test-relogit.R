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
expect_true(inherits(subsampling.results, "list"),
            info = "Output should be a list.")
expect_true(inherits(subsampling.results, "ssp.relogit"),
            info = "Output should be of class 'ssp.relogit'")

expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'LCC',
                            likelihood = 'logOddsCorrection'), 
              info = "It should run without errors on valid input.")

expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'optA',
                            likelihood = 'weighted'), 
              info = "It should run without errors on valid input.")

expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'optA',
                            likelihood = 'weighted',
                            epsilon = 1e-6),
              info = "It should run without errors when pass
              arguments to svyglm() through '...' .")

expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'optA',
                            likelihood = 'weighted',
                            control = list(alpha=0.2)),
              info = "It should run without errors when use control argument.")

expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            subset = c(1:(N/2)), 
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'optA',
                            likelihood = 'logOddsCorrection'),
             info = "It should run without errors when use subset argument.")

data$F1 <- sample(c("A", "B", "C"), N, replace=TRUE)
colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""), "F1")
expect_silent(subsampling.results <- 
                ssp.relogit(formula = formula,
                            data = data,
                            n.plt = n.plt,
                            n.ssp = n.ssp,
                            criterion = 'optA',
                            likelihood = 'logOddsCorrection',
                            contrasts = list(F1="contr.treatment")), 
              info = "It should run without errors when use contrasts argument.")


# Cleanup
rm(list = ls())
gc()