#' Post Clustering inference with Data Fission
#'
#' @param X The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p  of X
#' @param tau The tuning parameter used of data fission
#' @param cl_fun The clustering function that will applied on f(X)
#' @param k The number of clusters built using cl_fun
#' @param k1 The first cluster of interest among the k
#' @param k2 The second cluster of interest among the k
#' @param test The statistical test to use for inference
#' @param ... further arguments that could be parsed in test
#'
#' @return
#' A list with the following elements \itemize{
#' \item \code{Cluster}: The estimated partition using cl_fun on f(X)
#' \item \code{p.value}: The t-test p-values computed on g(X) using labels estimated on f(X)
#' }
#' @export
#'
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' fission <- data_fission(X, tau = 0.4)
#' cl_fun <- function(x, K){
#'   km <- kmeans(x, center  = K)
#'   return(km$cluster)
#' }
#' par(mfrow = c(1,3))
#' plot(X)
#' cl <- cl_fun(X, K=2)
#'plot(fission$fX, col = cl)
#'plot(fission$gX, col = cl)
#'dev.off()
#'pval <- apply_fission(X, tau = 0.4, cl_fun = cl_fun, k=2, k1=1, k2=2, test=t.test)

apply_fission <- function(X, Sigma = NULL, tau = 0.4, cl_fun, k, k1, k2, test,...){
  fission <- data_fission(X, Sigma, tau)
  cl <- cl_fun(fission$fX, K = k)
  p.value <- apply(fission$gX, 2, function(x){test(x=x[cl==k1], y=x[cl==k2])$p.value})
  return(list(Cluster = cl, p.value = p.value))
}

