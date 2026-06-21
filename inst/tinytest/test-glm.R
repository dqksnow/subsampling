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

set.seed(101)
uniform.results <- ssp.glm(formula = formula,
                           data = data,
                           n.plt = n.plt,
                           n.ssp = n.ssp,
                           family = family,
                           criterion = "uniform",
                           sampling.method = "withReplacement",
                           likelihood = "weighted")

expect_equal(uniform.results$subsample.size.expect,
             n.plt + n.ssp,
             info = "For criterion = 'uniform', the expected subsample size should be n.plt + n.ssp.")

expect_true(!is.null(uniform.results$index),
            info = "Returned object should store the drawn subsample indices in 'index'.")

expect_equal(length(uniform.results$index),
             n.plt + n.ssp,
             info = "With replacement uniform sampling should draw exactly n.plt + n.ssp observations.")

expect_true(all(is.finite(uniform.results$coef)),
            info = "Coefficient estimates should be finite for a valid uniform subsample fit.")

expect_true(all(is.finite(diag(uniform.results$cov))),
            info = "Estimated variances should be finite for a valid uniform subsample fit.")

expect_error(ssp.glm(formula = formula,
                     data = data,
                     n.plt = n.plt,
                     n.ssp = n.ssp,
                     family = family,
                     criterion = "optL",
                     sampling.method = "withReplacement",
                     likelihood = "logOddsCorrection"),
             info = paste(
               "'logOddsCorrection' should error with sampling.method =",
               "'withReplacement'."
             ))

expect_error(ssp.glm(formula = formula,
                     data = data,
                     n.plt = n.plt,
                     n.ssp = n.ssp,
                     family = "poisson",
                     criterion = "optL",
                     sampling.method = "poisson",
                     likelihood = "logOddsCorrection"),
             info = paste(
               "'logOddsCorrection' should error for non-binomial families."
             ))










# Cleanup
rm(list = ls())
gc()
