set.seed(1)
# softmax regression
d <- 3 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 1e4
beta.true <- 0.2 * matrix(-1, d, K)
beta.true.sum <- cbind(rep(1, d), beta.true)
set.seed(1)
mu <- rep(0, d)
sigma <- matrix(0.5, nrow = d, ncol = d)
diag(sigma) <- rep(1, d)
X <- MASS::mvrnorm(N, mu, sigma)
prob <- exp( X %*% beta.true.sum)
prob <- prob / rowSums(prob)
Y <- apply(prob, 1, function(row) sample(0:K, size = 1, prob = row))
n.plt <- 300
n.ssp <- 500
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ . -1
expect_silent(results <- ssp.softmax(formula = formula,
                                          data = data, 
                                          n.plt = n.plt,
                                          n.ssp = n.ssp, 
                                          criterion = 'MSPE', 
                                          sampling.method = 'withReplacement',
                                          likelihood = 'weighted',
                                          constraint = 'baseline'), 
              info = "It should run without errors on valid input.")

expect_silent(results <- ssp.softmax(formula = formula,
                                          data = data, 
                                          n.plt = n.plt,
                                          n.ssp = n.ssp, 
                                          criterion = 'optL', 
                                          sampling.method = 'withReplacement',
                                          likelihood = 'weighted',
                                          constraint = 'baseline'), 
              info = "It should run without errors on valid input.")

expect_true(inherits(results, "list"), info = "Output should be a list.")
expect_true(inherits(results, "ssp.softmax"), 
            info = "Output should be of class 'ssp.softmax'")

expect_silent(results <- ssp.softmax(formula = formula,
                                          data = data,
                                          subset = c(1:(N/2)),
                                          n.plt = n.plt, 
                                          n.ssp = n.ssp, 
                                          criterion = 'MSPE', 
                                          sampling.method = 'withReplacement',
                                          likelihood = 'weighted',
                                          constraint = 'baseline'), 
              info = "It should run without errors when use subset argument.")

data$F1 <- sample(c("A", "B", "C"), N, replace=TRUE)
colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""), "F1")

expect_silent(results <- 
                ssp.softmax(formula = formula,
                            data = data,
                            n.plt = n.plt, 
                            n.ssp = n.ssp, 
                            criterion = 'MSPE', 
                            sampling.method = 'withReplacement',
                            likelihood = 'weighted',
                            constraint = 'baseline',
                            contrasts = list(F1="contr.treatment")),
              info = "It should run without errors when use contrast argument.")
## Cleanup
rm(list = ls())
gc()