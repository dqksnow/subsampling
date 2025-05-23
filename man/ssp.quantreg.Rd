% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_main_function.R
\name{ssp.quantreg}
\alias{ssp.quantreg}
\title{Optimal Subsampling Methods for Quantile Regression Model}
\usage{
ssp.quantreg(
  formula,
  data,
  subset = NULL,
  tau = 0.5,
  n.plt,
  n.ssp,
  B = 5,
  boot = TRUE,
  criterion = "optL",
  sampling.method = "withReplacement",
  likelihood = c("weighted"),
  control = list(...),
  contrasts = NULL,
  ...
)
}
\arguments{
\item{formula}{A model formula object of class "formula" that describes the model to be fitted.}

\item{data}{A data frame containing the variables in the model. Denote \eqn{N} as the number of observations in \code{data}.}

\item{subset}{An optional vector specifying a subset of observations from \code{data} to use for the analysis. This subset will be viewed as the full data.}

\item{tau}{The interested quantile.}

\item{n.plt}{The pilot subsample size (first-step subsample size).
This subsample is used to compute the pilot estimator and estimate the optimal subsampling probabilities.}

\item{n.ssp}{The expected size of the optimal subsample (second-step subsample). For \code{sampling.method = 'withReplacement'}, The exact subsample size is \code{n.ssp}. For \code{sampling.method = 'poisson'}, \code{n.ssp} is the expected subsample size.}

\item{B}{The number of subsamples for the iterative sampling algorithm. Each subsample contains \code{n.ssp} observations. This allows us to estimate the covariance matrix.}

\item{boot}{If TRUE then perform iterative sampling algorithm and estimate the covariance matrix. If FALSE then only one subsample with size \code{B*n.ssp} is returned.}

\item{criterion}{It determines how subsampling probabilities are computed.
Choices include \code{optL}(default) and \code{uniform}.
\itemize{
\item \code{optL} Minimizes the trace of a transformation of the asymptotic covariance matrix of the subsample estimator.
\item \code{uniform} Assigns equal subsampling probability
\eqn{\frac{1}{N}} to each observation, serving as a baseline subsampling strategy.
}}

\item{sampling.method}{The sampling method for drawing the optimal subsample.
Choices include \code{withReplacement} and \code{poisson}(default). \code{withReplacement} draws exactly \code{n.ssp}
subsamples from size \eqn{N} full dataset with replacement, using the specified
subsampling probabilities. \code{poisson} draws observations independently by
comparing each subsampling probability with a realization of uniform random
variable  \eqn{U(0,1)}.}

\item{likelihood}{The type of the maximum likelihood function used to
calculate the optimal subsampling estimator. Currently \code{weighted} is implemented which applies a weighted likelihood function where each observation is weighted by the inverse of its subsampling probability.}

\item{control}{The argument \code{control} contains two tuning parameters \code{alpha} and \code{b}.
\itemize{
\item \code{alpha} \eqn{\in [0,1]} is the mixture weight of the user-assigned subsampling
probability and uniform subsampling probability. The actual subsample
probability is \eqn{\pi = (1-\alpha)\pi^{opt} + \alpha \pi^{uni}}. This protects the estimator from extreme small
subsampling probability. The default value is 0.
\item \code{b} is a positive number which is used to constaint the poisson subsampling probability. \code{b} close to 0 results in subsampling probabilities closer to uniform probability \eqn{\frac{1}{N}}. \code{b=2} is the default value.
See relevant references for further details.
}}

\item{contrasts}{An optional list. It specifies how categorical variables are represented in the design matrix. For example, \code{contrasts = list(v1 = 'contr.treatment', v2 = 'contr.sum')}.}

\item{...}{A list of parameters which will be passed to \code{quantreg::rq()}.}
}
\value{
\code{ssp.quantreg} returns an object of class "ssp.quantreg" containing the following components (some are optional):

\describe{
\item{model.call}{The original function call.}
\item{coef.plt}{The pilot estimator. See Details for more information.}
\item{coef}{The estimator obtained from the optimal subsample.}
\item{cov}{The covariance matrix of \code{coef}}
\item{index.plt}{Row indices of pilot subsample in the full dataset.}
\item{index.ssp}{Row indices of of optimal subsample in the full dataset.}
\item{N}{The number of observations in the full dataset.}
\item{subsample.size.expect}{The expected subsample size}
\item{terms}{The terms object for the fitted model.}
}
}
\description{
Draw subsample from full dataset and fit quantile regression model. For a quick start, refer to the \href{https://dqksnow.github.io/subsampling/articles/ssp-quantreg.html}{vignette}.
}
\details{
Most of the arguments and returned variables have the same meaning with \link{ssp.glm}. Refer to \href{https://dqksnow.github.io/subsampling/articles/ssp-logit.html}{vignette}

A pilot estimator for the unknown parameter \eqn{\beta} is required because
optL subsampling probabilities depend on \eqn{\beta}. There is no "free lunch" when determining optimal subsampling probabilities. For quantile regression, this
is achieved by drawing a size \code{n.plt} subsample with replacement from full
dataset, using uniform sampling probability.

If \code{boot}=TRUE, the returned value \code{subsample.size.expect} equals to \code{B*n.ssp}, and the covariance matrix for \code{coef} would be calculated.
If \code{boot}=FALSE, the returned value \code{subsample.size.expect} equals to \code{B*n.ssp}, but the covariance matrix won't be estimated.
}
\examples{
#quantile regression
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
Y <- beta.true[1] + X \%*\% beta.true[-1] + 
err * rowMeans(abs(X))
data <- as.data.frame(cbind(Y, X))
colnames(data) <- c("Y", paste("V", 1:ncol(X), sep=""))
formula <- Y ~ .
n.plt <- 200
n.ssp <- 100
optL.results <- ssp.quantreg(formula,data,tau = tau,n.plt = n.plt,
n.ssp = n.ssp,B = B,boot = TRUE,criterion = 'optL',
sampling.method = 'withReplacement',likelihood = 'weighted')
summary(optL.results)
uni.results <- ssp.quantreg(formula,data,tau = tau,n.plt = n.plt,
n.ssp = n.ssp,B = B,boot = TRUE,criterion = 'uniform',
sampling.method = 'withReplacement', likelihood = 'weighted')
summary(uni.results)
}
\references{
Wang, H., & Ma, Y. (2021). Optimal subsampling for quantile regression in big data. \emph{Biometrika}, \strong{108}(1), 99-112.
}
