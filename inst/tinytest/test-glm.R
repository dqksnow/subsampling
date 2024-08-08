set.seed(1)
N <- 1e4
beta0 <- rep(-0.5, 7)
d <- length(beta0) - 1
X <- matrix(0, N, d)
generate_rexp <- function(x) x <- rexp(N, rate = 2)
X <- apply(X, 2, generate_rexp)
Y <- rbinom(N, 1, 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1])))
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ .
n.plt <- 500
n.ssp <- 1000
family <- 'quasibinomial'
expect_silent(ssp.results <- ssp.glm(formula = formula,
                                             data = data,
                                             n.plt = n.plt,
                                             n.ssp = n.ssp,
                                             family = family, 
                                             criterion = "optL",
                                             sampling.method = 'poisson',
                                             likelihood = "logOddsCorrection"), 
              info = "It should run without errors on valid input.")
expect_true(inherits(ssp.results, "list"),
            info = "Output should be a list.")
expect_true(inherits(ssp.results, "ssp.glm"), 
            info = "Output should be of class 'ssp.glm'")

expect_silent(ssp.results <- ssp.glm(formula = formula,
                                     data = data,
                                     n.plt = n.plt,
                                     n.ssp = n.ssp,
                                     family = family, 
                                     criterion = "optA",
                                     sampling.method = 'poisson',
                                     likelihood = "logOddsCorrection"), 
              info = "It should run without errors on valid input.")
expect_silent(ssp.results <- ssp.glm(formula = formula,
                                     data = data,
                                     n.plt = n.plt,
                                     n.ssp = n.ssp,
                                     family = family, 
                                     criterion = "LCC",
                                     sampling.method = 'poisson',
                                     likelihood = "weighted"), 
              info = "It should run without errors on valid input.")
expect_silent(ssp.results <- ssp.glm(formula = formula,
                                     data = data,
                                     n.plt = n.plt,
                                     n.ssp = n.ssp,
                                     family = family, 
                                     criterion = "uniform",
                                     sampling.method = 'poisson',
                                     likelihood = "logOddsCorrection"), 
              info = "It should run without errors on valid input.")
expect_silent(ssp.results <- 
                ssp.glm(formula = formula,
                        data = data,
                        subset = c(1:(N/2)), 
                        n.plt = n.plt,
                        n.ssp = n.ssp,
                        family = family, 
                        criterion = "optL",
                        sampling.method = 'poisson',
                        likelihood = "weighted"),
              info = "It should run without errors when use subset")

expect_silent(ssp.results <- 
                ssp.glm(formula = formula,
                        data = data,
                        subset = c(1:(N/2)), 
                        n.plt = n.plt,
                        n.ssp = n.ssp,
                        family = family, 
                        criterion = "optA",
                        sampling.method = 'poisson',
                        likelihood = "logOddsCorrection",
                        maxit = 30),
              info = "It should run without errors when pass
              arguments to svyglm() through '...' .")

expect_silent(ssp.results <- 
                ssp.glm(formula = formula,
                        data = data,
                        subset = c(1:(N/2)), 
                        n.plt = n.plt,
                        n.ssp = n.ssp,
                        family = family, 
                        criterion = "LCC",
                        sampling.method = 'poisson',
                        likelihood = "logOddsCorrection",
                        control = list(alpha=0.1)),
              info = "It should run without errors when use control argument.")


data$F1 <- sample(c("A", "B", "C"), N, replace=TRUE)
colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""), "F1")
expect_silent(ssp.results <- 
                ssp.glm(formula = formula,
                        data = data,
                        n.plt = n.plt,
                        n.ssp = n.ssp,
                        family = family, 
                        criterion = "optL",
                        sampling.method = 'poisson',
                        likelihood = "logOddsCorrection",
                        contrasts = list(F1 = 'contr.sum')),
              
              info = "It should run without errors when use contrasts.")











# Cleanup
rm(list = ls())
gc()