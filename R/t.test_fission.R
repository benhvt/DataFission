#' Perform independent two-sample t-tests on columns of a matrix or data frame.
#'
#' @param X A matrix of size n x p containing the data.
#' @param cl A vector of length n specifying the class labels for each observation in X.
#' @param ... Additional arguments to be passed to the \code{t.test} function.
#'
#' @return A dataframe containing p-values for the t-test between each clusters pairs for each variables
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' cl <- c(rep(1, 5), rep(2, 5))
#'
#' # Perform t-test using t_test.fission
#' p_values <- t_test.fission(X, cl)
#'
#' # Display the results
#' p_values
#'
#' @references
#' For more information on the t.test function, refer to its documentation: ?t.test
#'
#' @export
t_test.fission <- function(X, cl, ...) {
  if(!is.factor(cl)){
    cl <- as.factor(cl)
  }

  K <- length(levels(cl))

  # browser()

  ttest_res <- apply(X, 2, function(x) {
    t_test_multiclusters(x = x, cl = cl, ...)
  })
  results <- do.call("rbind.data.frame", ttest_res)
  if (is.null(colnames(X))){
    results$Variable <-  rep(paste0("X", 1:ncol(X)), each = choose(K, 2))
  }
  else{
    results$Variable <- rep(colnames(X), each = choose(K, 2))
  }

  return(results)
}
