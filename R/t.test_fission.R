#' Perform independent two-sample t-tests on columns of a matrix or data frame.
#'
#' @param X A matrix of size n x p containing the data.
#' @param cl A vector of length n specifying the class labels for each observation in X.
#' @param k1 The first class label for the t-test.
#' @param k2 The second class label for the t-test.
#' @param ... Additional arguments to be passed to the \code{t.test} function.
#'
#' @return A numeric vector of length p containing the p-values obtained from the t-tests.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' cl <- c(rep(1, 5), rep(2, 5))
#'
#' # Perform t-test using t_test.fission
#' p_values <- t_test.fission(X, cl, 1, 2)
#'
#' # Display the results
#' p_values
#'
#' @references
#' For more information on the t.test function, refer to its documentation: ?t.test
#'
#' @export
t_test.fission <- function(X, cl, k1, k2, ...) {
  ttest_res <- apply(X, 2, function(x) {
    stats::t.test(x[cl == k1], x[cl == k2], ...)$p.value
  })
  return(ttest_res)
}
