#' Performs t-test for K >= 2 clusters
#'
#' @param x a numerical vectors of size n
#' @param cl A vector of length n specifying the class labels for each observation in x
#' @param ... Additional arguments to be passed to the \code{t.test} function.
#'
#' @return A dataframe containing p-values for the t-test between each clusters pairs
#
#' @export
#'
#' @examples
#' n <- 150
#' cl <- sample(1:5, size = n, replace = T)
#' x <- rnorm(n)
#' t_test_multiclusters(x, cl)

t_test_multiclusters <- function(x, cl, ...) {

  if(!is.factor(cl)){
    cl <- as.factor(cl)
  }

  K <- length(levels(cl))
  all_pairs <- utils::combn(1:K, 2)

  all_tests <- lapply(1:ncol(all_pairs), function(k){
    cl2test <- all_pairs[,k]
    pval <- stats::t.test(x[cl == cl2test[1]], x[cl == cl2test[2]],...)$p.value
    return(data.frame(Cluster = paste0("Cluster ", cl2test[1], " vs Cluster ", cl2test[2]),
                      pvalues = pval))
  })
  return(do.call("rbind.data.frame", all_tests))
}
