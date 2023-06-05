#' Perform Differential Expression Analysis usgin dearseq.
#'
#' @param X A matrix of size n x p containing the expression data.
#' @param cl  A numeric design matrix of size n x K containing the K variables to be tested.
#' @param ... Additional arguments to be passed to the \code{dear_seq} function.
#'
#'@import dearseq
#'
#' @return A numeric vector containing the p-values obtained from the dearseq analysis.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' cl <- c(rep(1, 5), rep(2, 5))
#'
#' # Perform DEAR-Seq analysis using dear_seq.fission
#' p_values <- dear_seq.fission(X, cl, preprocessed = T)
#'
#' # Display the results
#' p_values
#'
#' @references
#' For more information on the dear_seq function, refer to its documentation: ?dear_seq
#'
#' @export
dear_seq.fission <- function(X, cl, ...) {
  dear_res <- dear_seq(exprmat = t(X), variables2test = as.matrix(cl, ncol = 1), ...)
  return(dear_res$pvals$rawPval)
}
