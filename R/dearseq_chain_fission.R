#' Post clustering inference using chain data fission and dearseq
#'
#' @param X  The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p  of X
#' @param tau The tuning parameter used of data fission
#' @param nFission The number of repeated data fission to perform
#' @param preprocessed A logical flag indicating whether the expression data
#' have already been preprocessed (e.g. log2 transformed).
#' Default is FALSE, in which case y is assumed to contain raw counts and is
#' normalized into log(counts) per million.
#' @param merge_function The merging function used to derive chain fission p values
#' @param cl_fun The clustering function that will applied on f(X)
#' @param k The number of clusters built using cl_fun
#' @param k1 The first cluster of interest among the k
#' @param k2 The second cluster of interest among the k
#' @param parallel A logical flag indicating whether parallel computation should be enabled
#' @param ncores 	An integer indicating the number of cores to be used when parallel is TRUE
#'
#' @import pbapply
#' @import dearseq
#'
#' @return
#' #' A list with the following elements \itemize{
#' \item \code{Cluster}: the nFission estimated partition on f(X)
#' \item \code{p.value}: a numeric vector of length p containing the nFission merging p-values or a list
#' of M numeric vectors of length p if more than one merging method are used}
#' @export
#'
#' @examples
#' cl_fun <- function(x, K){
#' km <- kmeans(x, center  = K)
#'   return(km$cluster)
#' }
#' X1 <- c(rnorm(100, sd = 1), rnorm(100, mean = 2, sd = 1))
#' X2 <- c(rnorm(200))
#' X <- cbind(X1,X2)
#' plot(X, col = rep(1:2, each = 100))
#' dearseq_chain_fission(X, nFission = 10, cl_fun = cl_fun, k=2,  k1=1, k2=2, preprocessed = TRUE)
#' dearseq_fission(X, tau = 0.4, cl_fun = cl_fun, k=2, k1=1, k2=2, preprocessed = TRUE)

dearseq_chain_fission <- function(X,
                                Sigma = NULL,
                                tau = 0.4,
                                nFission = 1000,
                                preprocessed = FALSE,
                                merge_function = c("KS",
                                                   "Fisher",
                                                   "geometric",
                                                   "harmonic"),
                                cl_fun,
                                k,
                                k1,
                                k2,
                                parallel = FALSE,
                                ncores = NULL){
  if(is.null(Sigma)){
    Sigma <- stats::cov(X)
  }

  if (!preprocessed){
    X <- sistmr::counts_normalization(raw_counts = X,
                                      MARGIN = "row",
                                      plot = FALSE)
    }


  if (parallel){
    multi_fission <- pbapply::pblapply(1:nFission, function(x){
      fission <- dearseq_fission(X,
                               Sigma = Sigma,
                               tau = tau,
                               cl_fun = cl_fun,
                               k = k,
                               k1 = k1,
                               k2 = k2,
                               preprocessed = TRUE)
      return(fission)
    }, cl = ncores)
  }

  if (!parallel){
    multi_fission <- lapply(1:nFission, function(x){
      fission <- dearseq_fission(X,
                                 Sigma = Sigma,
                                 tau = tau,
                                 cl_fun = cl_fun,
                                 k = k,
                                 k1 = k1,
                                 k2 = k2,
                                 preprocessed = TRUE)
      return(fission)
    })
  }
  multi_fission_p.value <- lapply(multi_fission, function(x){x$p.value$rawPval})
  multi_fission_clustering <- lapply(multi_fission, function(x){x$Cluster})
  all_pvalue_multifission <- do.call("rbind", multi_fission_p.value)
  if (length(merge_function) == 1){
    pval <- merge_pvalue(all_pvalue_multifission, method = merge_function)
  } else {
    # browser()
    pval <- lapply(merge_function, function(x){merge_pvalue(all_pvalue_multifission, method = x)})
    names(pval) <- merge_function

  }
  return(list(p.value = pval,
              Cluster = multi_fission_clustering))
}

