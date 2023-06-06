#' Perform post-clustering inference with data fission
#'
#' @param X The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p of X
#' @param tau The tuning parameter used for data fission
#' @param cl_fun A clustering function that returns a factor vector containing the clusters
#' @param k The number of clusters built using cl_fun
#' @param k1 The first cluster of interest among the k
#' @param k2 The second cluster of interest among the k
#' @param test A test function that inputs a numeric matrix, a partition, and the two clusters to test, and outputs p-values
#' @param ... Additional arguments that can be passed to the test function
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{Cluster}: The estimated partition using cl_fun on f(X)
#'   \item \code{p.value}: The p-values computed on g(X) using labels estimated on f(X)
#' }
#'
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' fiss <- data_fission(X, tau = 0.4)
#' cl_fun <- function(x, K) {
#'   km <- kmeans(x, center = K)
#'   return(km$cluster)
#' }
#' par(mfrow = c(1, 3))
#' plot(X)
#' cl <- cl_fun(X, K = 2)
#' plot(fiss$fX, col = cl)
#' plot(fiss$gX, col = cl)
#' dev.off()
#' pval <- fission(X, tau = 0.4, cl_fun = cl_fun, k = 2, k1 = 1, k2 = 2, test = t_test.fission)
#'
#' @export
fission <- function(X, Sigma = NULL, tau = 0.4, cl_fun, k, k1, k2, test, ...) {
  fiss <- data_fission(X, Sigma, tau)
  cl <- cl_fun(fiss$fX, K = k)
  p.value <- test(fiss$gX, cl, k1, k2, ...)
  return(list(Cluster = cl, p.value = p.value))
}
