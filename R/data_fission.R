#' Data fission of Gaussian data
#'
#' @param X The data matrix of size n x p for fission.
#' @param Sigma  An estimator of the covariance matrix of size p x p  of X.
#' @param tau The tuning parameter used of data fission
#'
#' @import mvtnorm
#' @return A list with the following elements \itemize{
#' \item \code{fX}: The dataset that will be used for clustering
#' \item \code{gX}: The dataset that will be used for inference}
#' @export
#'
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' fission <- data_fission(X, tau = 0.4)
#' cl <- kmeans(fission$fX, centers = 2)
#' par(mfrow = c(1,3))
#' plot(X)
#'plot(fission$fX, col = cl$cluster)
#'plot(fission$gX, col = cl$cluster)


data_fission <- function(X, Sigma = NULL, tau = 0.4){
  if (!is.matrix(X))
    stop("X should be a matrix")
  if(!is.null(Sigma)){
    if(!is.matrix(Sigma))
      stop("Sigma must be a matrix")
    if (ncol(Sigma)!= ncol(X) | nrow(Sigma)!= ncol(X))
      stop("Sigma must be of dimension p x p")
  }

  if(is.null(Sigma)){
    Sigma <- stats::cov(X)
  }
  n <- nrow(X)
  p <- ncol(X)
  Z <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = Sigma)
  fX <- X + tau*Z
  gX <- X - ((1/tau)*Z)
  return(list(fX = fX,
              gX = gX))
}
