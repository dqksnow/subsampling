set.seed(1)
# softmax regression
d <- 3 # dim of covariates
K <- 5 # K + 1 classes
G <- rbind(rep(-1/(K+1), K), diag(K) - 1/(K+1)) %x% diag(d)
N <- 2e4
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
n.plt <- 500
n.ssp <- 1000
data <- as.data.frame(cbind(Y, X))
formula <- Y ~ X
expect_silent(WithRep.MSPE <- ssp.softmax(formula, data, n.plt, n.ssp, 
                                    criterion = 'MSPE', 
                                    sampling.method = 'withReplacement',
                                    likelihood = 'weighted',
                                    constraint = 'baseline'), 
              info = "It should run without errors on valid input.")
## Cleanup
rm(list = ls())
gc()