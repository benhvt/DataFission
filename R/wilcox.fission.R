#' Perform Wilcoxon rank-sum tests on columns of a matrix or data frame.
#'
#' @param X A matrix of size n x p containing the data.
#' @param cl A vector of length n specifying the class labels for each observation in X.
#' @param k1 The first class label for the Wilcoxon rank-sum test.
#' @param k2 The second class label for the Wilcoxon rank-sum test.
#' @param ... Additional arguments to be passed to the \code{wilcox.test} function.
#'
#' @return A numeric vector of length p containing the p-values obtained from the Wilcoxon rank-sum tests.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' cl <- c(rep(1, 5), rep(2, 5))
#'
#' # Perform Wilcoxon rank-sum test using wilcox.fission
#' p_values <- wilcox.fission(X, cl, 1, 2)
#'
#' # Display the results
#' p_values
#'
#' @references
#' For more information on the wilcox.test function, refer to its documentation: ?wilcox.test
#'
#' @export
wilcox.fission <- function(X, cl, k1, k2, ...) {
  wilcox_res <- apply(X, 2, function(x) {
    stats::wilcox.test(x[cl == k1], x[cl == k2], ...)
  })
  return(wilcox_res$p.value)
}
