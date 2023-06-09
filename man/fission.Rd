% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fission.R
\name{fission}
\alias{fission}
\title{Perform post-clustering inference with data fission}
\usage{
fission(X, Sigma = NULL, tau = 0.4, cl_fun, test, ...)
}
\arguments{
\item{X}{The data matrix of size n x p where post-clustering inference must be applied}

\item{Sigma}{An estimator of the covariance matrix of size p x p of X}

\item{tau}{The tuning parameter used for data fission}

\item{cl_fun}{A clustering function that returns a factor vector containing the clusters}

\item{test}{A test function that inputs a numeric matrix, a partition, and the two clusters to test, and outputs p-values}

\item{...}{Additional arguments that can be passed to the test function}
}
\value{
A list with the following elements:
\itemize{
  \item \code{Cluster}: The estimated 2-clusters partition using \code{cl_fun} on f(X)
  \item \code{p.value}: The p-values computed on g(X) using labels estimated on f(X)
}
}
\description{
Perform post-clustering inference with data fission
}
\examples{
X <- matrix(rnorm(200), ncol = 2)
fiss <- data_fission(X, tau = 0.4)
cl_fun <- function(x, K) {
  km <- kmeans(x, center = K)
  return(km$cluster)
}
par(mfrow = c(1, 3))
plot(X)
cl <- cl_fun(X, K = 2)
plot(fiss$fX, col = cl)
plot(fiss$gX, col = cl)
dev.off()
pval <- fission(X, tau = 0.4, cl_fun = cl_fun, k = 2, k1 = 1, k2 = 2, test = t_test.fission)

}
