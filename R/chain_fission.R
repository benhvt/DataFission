#' Perform post-clustering inference using chain data fission
#'
#' @param X The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p of X. If NULL, sample covariance estimator is used
#' @param tau The tuning parameter used for data fission
#' @param nFission The number of repeated data fissions to perform
#' @param merge_function The merging function used to derive chain fission p-values. It can be one of "KS", "Fisher", "geometric", or "harmonic"
#' @param cl_fun A clustering function that returns a factor vector containing the clusters
#' @param k The number of clusters to be built using cl_fun
#' @param k1 The first cluster of interest among the k
#' @param k2 The second cluster of interest among the k
#' @param parallel A logical flag indicating whether parallel computation should be enabled
#' @param ncores An integer indicating the number of cores to be used when parallel is TRUE
#' @param test A test function that inputs a numeric matrix, a partition and the two clusters to test and that outputs p-values
#' @param ... Additional arguments that can be passed to the test function
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{Cluster}: The nFission estimated partitions
#'   \item \code{p.value}: A numeric vector of length p containing the nFission merging p-values or a list of M numeric vectors of length p if more than one merging method is used
#' }
#'
#' @examples
#' cl_fun <- function(x, K) {
#'   km <- kmeans(x, center = K)
#'   return(km$cluster)
#' }
#'
#' X1_cl1 <- matrix(rnorm(100 * 25, mean = 0), ncol = 25)
#' X1_cl2 <- matrix(rnorm(100 * 25, mean = 1.5), ncol = 25)
#' X1 <- rbind(X1_cl1, X1_cl2)
#' X2 <- matrix(rnorm(200 * 50), ncol = 50)
#' X <- cbind(X1, X2)
#' plot(X, col = rep(1:2, each = 100))
#' chain_fission(X, nFission = 10, cl_fun = cl_fun, k = 2, k1 = 1, k2 = 2, test = t_test.fission)
#' fission(X, tau = 0.4, cl_fun = cl_fun, k = 2, k1 = 1, k2 = 2, test = t_test.fission)
#'
#' @import pbapply
#' @importFrom stats cov
#'
#' @export

chain_fission <- function(X,
                          Sigma = NULL,
                          tau = 0.4,
                          nFission = 1000,
                          merge_function = c("KS",
                                                   "Fisher",
                                                   "geometric",
                                                   "harmonic"),
                          cl_fun,
                          k,
                          k1,
                          k2,
                          parallel = FALSE,
                          ncores = NULL,
                          test,
                          ...){
  if(is.null(Sigma)){
    Sigma <- stats::cov(X)
  }

  if (parallel){
    multi_fission <- pbapply::pblapply(1:nFission, function(x){
      fission <- fission(X,
                         Sigma = Sigma,
                         tau = tau,
                         cl_fun = cl_fun,
                         k = k,
                         k1 = k1,
                         k2 = k2,
                         test = test,
                         ...)
      return(fission)
    }, cl = ncores)
  }

  if (!parallel){
    multi_fission <- lapply(1:nFission, function(x){
      fission <- fission(X,
                         Sigma = Sigma,
                         tau = tau,
                         cl_fun = cl_fun,
                         k = k,
                         k1 = k1,
                         k2 = k2,
                         test = test,
                         ...)
      return(fission)
    })
  }
  multi_fission_p.value <- lapply(multi_fission, function(x){x$p.value})
  multi_fission_clustering <- lapply(multi_fission, function(x){x$Cluster})
  all_pvalue_multifission <- do.call("rbind", multi_fission_p.value)
  if (length(merge_function) == 1){
    pval <- merge_pvalue(all_pvalue_multifission, method = merge_function)
  } else {
    pval <- lapply(merge_function, function(x){merge_pvalue(all_pvalue_multifission, method = x)})
    names(pval) <- merge_function

  }
  return(list(p.value = pval,
              Cluster = multi_fission_clustering))
}

