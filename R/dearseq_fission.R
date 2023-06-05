#' Post Clustering inference with Data Fission using dearseq
#'
#' @param X The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p  of X
#' @param tau The tuning parameter used of data fission
#' @param cl_fun The clustering function that will applied on f(X)
#' @param k The number of clusters built using cl_fun
#' @param k1 The first cluster of interest among the k
#' @param k2 The second cluster of interest among the k
#' @param preprocessed A logical flag indicating whether the expression data
#' have already been preprocessed (e.g. log2 transformed).
#' Default is FALSE, in which case y is assumed to contain raw counts and is
#' normalized into log(counts) per million.
#'
#' @import dearseq
#'
#' @return
#' A list with the following elements \itemize{
#' \item \code{Cluster}: The estimated partition using cl_fun on f(X)
#' \item \code{p.value}: A data.frame with p rows containing the row and the
#' adjusted p-values returned by \code{dear_seq} when it applies on g(X)
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
#'pval <- dearseq_fission(X, tau = 0.4, cl_fun = cl_fun, k=2, k1=1, k2=2, preprocessed = TRUE)

dearseq_fission <- function(X, Sigma = NULL, tau = 0.4, cl_fun, k, k1, k2, preprocessed = FALSE){
  if (!preprocessed){
    X <- sistmr::counts_normalization(raw_counts = X,
                                      MARGIN = "row",
                                      plot = FALSE)
  }
  fission <- data_fission(X, Sigma, tau)
  cl <- cl_fun(fission$fX, K = k)
  ToTest <- as.matrix(as.numeric(cl), ncol = 1)
  dear_res <- dear_seq(exprmat = t(fission$gX),
                      variables2test = ToTest,
                      which_test = "asymptotic",
                      preprocessed = TRUE)
  return(list(Cluster = cl, p.value = dear_res$pvals))
}
