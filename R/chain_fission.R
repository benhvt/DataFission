#' Perform post-clustering inference using chain data fission
#'
#' @param X The data matrix of size n x p where post-clustering inference must be applied
#' @param Sigma An estimator of the covariance matrix of size p x p of X. If NULL, sample covariance estimator is used
#' @param tau The tuning parameter used for data fission
#' @param nFission The number of repeated data fissions to perform
#' @param merge_function The merging function used to derive chain fission p-values. It can be one of "KS", "shapiro", "Fisher", "geometric", or "harmonic"
#' @param K The number of clusters to build (must be an argument of the cl_fun)
#' @param cl_fun A clustering function that returns a factor vector containing the clusters
#' @param cl_ref A reference clustering that could be used to map labels obtained on $f(X)$. Default is NULL indicating that mapping is done according to the labels obtained on the first fission
#' @param parallel A logical flag indicating whether parallel computation should be enabled
#' @param ncores An integer indicating the number of cores to be used when parallel is TRUE
#' @param test A test function that inputs a numeric matrix, a partition and the two clusters to test and that outputs p-values
#' @param ... Additional arguments that can be passed to the test function
#'
#' @import diceR
#' @import dplyr
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{Cluster}: The nFission estimated 2-clusters partitions
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
#' chain_fission(X, nFission = 10, cl_fun = cl_fun, K=2, test = t_test.fission)
#' fission(X, tau = 0.4, cl_fun = cl_fun, K=2, test = t_test.fission)
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
                          K,
                          cl_ref = NULL,
                          parallel = FALSE,
                          ncores = NULL,
                          test,
                          ...){
  if(is.null(Sigma)){
    Sigma <- stats::cov(X)
  }

  if (parallel){
    multi_fission <- pbapply::pblapply(1:nFission, function(x){
      fission <- data_fission(X,
                         Sigma = Sigma,
                         tau = tau)
      return(fission)
    }, cl = ncores)
  }

  if (!parallel){
    multi_fission <- lapply(1:nFission, function(x){
      fission <- data_fission(X,
                              Sigma = Sigma,
                              tau = tau)
      return(fission)
    })
  }
  multi_fission_fX <- lapply(multi_fission, function(l){l$fX})
  multi_fission_gX <- lapply(multi_fission, function(l){l$gX})

  if (parallel){
    multi_fission_clustering <- pblapply(multi_fission_fX,
                                         cl_fun,
                                         K,
                                         cl = ncores)
  }
  if (!parallel){
    multi_fission_clustering <- lapply(multi_fission_fX,
                                       cl_fun,
                                       K)
  }


  # browser()

  if (is.null(cl_ref)){
    cl_ref <- multi_fission_clustering[[1]]
  }

  multi_fission_clustering_relabel <- lapply(multi_fission_clustering, function(c){
    diceR::relabel_class(pred.cl = c, ref.cl = cl_ref)
  })

  if (parallel){
    multi_fission_p.value <- pblapply(1:nFission, function(n){
      test(multi_fission_gX[[n]], cl = multi_fission_clustering_relabel[[n]], ...)
    }, cl = ncores)
  }

  if(!parallel){
    multi_fission_p.value <- lapply(1:nFission, function(n){
      test(multi_fission_gX[[n]], cl = multi_fission_clustering_relabel[[n]], ...)
    })
  }

  all_pvalue_multifission <- do.call("rbind", multi_fission_p.value)

    if (length(merge_function)==1){
    aggregate_pvalues <- all_pvalue_multifission %>%
      group_by(.data$Cluster, .data$Variable) %>%
      summarise(AggregateP = merge_pvalue(.data$pvalues, method = merge_function)) %>%
      mutate(MergingFunction = merge_function)
  }
  else{
    aggregate_pvalues <- all_pvalue_multifission %>%
      group_by(.data$Cluster, .data$Variable) %>%
      reframe(AggregateP = sapply(merge_function, function(m){merge_pvalue(.data$pvalues, method = m)})) %>%
      mutate(MergingFunction = rep(merge_function, ncol(X)*length(table(all_pvalue_multifission$Cluster))))
  }

  return(list(p.value = aggregate_pvalues,
              Cluster = multi_fission_clustering_relabel,
              allfissions_pvals = all_pvalue_multifission))
}

