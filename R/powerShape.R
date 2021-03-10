#' M-estimators of Shape from the Power Family.
#'
#' Power M-estimators of shape and location were recently suggested in
#' Frahm et al. (2020). They have a tuning parameter `alpha` taking values in
#' \code{[0,1]}. The extreme case `alpha` = 1 corresponds to Tyler's shape
#' matrix and `alpha` = 0 to the classical covariance matrix. These special
#' cases have their own, more efficient functions \code{\link{tylerShape}} and
#' \code{\link{classicShape}}, respectively.
#' If the true location is known, it should be supplied as `center`, otherwise
#' it is estimated simultaneously with the shape.
#'
#' These functions assume that the data were generated from an
#' elliptical distribution, for Tyler's estimate this can be relaxed to
#' generalized elliptical distributions.
#'
#' For multivariate normally distributed data, \code{\link{classicShape}} is the maximum
#' likelihood estimator of location and scale. It is a special case of the
#' power M-estimator with tail index `alpha` = 0, which returns the empirical
#' covariance matrix and the empirical mean vector.
#'
#' The function \code{\link{tylerShape}} maximizes the likelihood function after projecting
#' the observed data of each individual onto the unit hypersphere, in which case
#' we obtain an angular central Gaussian distribution. It is a special case of
#' the power M-estimator with tail index `alpha` = 1, which returns Tyler's
#' M-estimator of scatter and an affine equivariant multivariate median
#' according to Hettmansperger and Randles (2002).
# TODO: Hettmansperger and randles
#'
#' The function \code{\link{powerShape}} requires an additional parameter, the so-called
#' tail index `alpha`. For heavy tailed data, the index should be chosen closer
#' to 1, whereas for light tailed data the index should be chosen closer to 0.
#'
#' @aliases powerShape
#' @aliases tylerShape
#' @aliases classicShape
#'
#' @usage powerShape(x, alpha, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage tylerShape(x, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage classicShape(x, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#'
#' @param x A numeric data matrix or data.frame without missing data.
#' @param alpha Tail index, a numeric value in the interval \code{[0, 1]}.
#'    Determines the power function. For more information see 'Details'.
#' @param center An optional vector of the data's center. If `NULL` the center
#'   will be estimated simultaneously with the shape.
#' @param normalization A string determining how the shape matrix is standardized.
#' The possible values are
#' \itemize{
#'     \item{`'det'`}{such that the returned shape estimate has determinant 1.}
#'     \item{`'trace'`}{such that the returned shape estimate has trace \code{ncol(x)}.}
#'     \item{`'one'`}{such that the returned shape estimate's top left entry
#'     (`S[1, 1]`) is 1.}
#'   }
#' @param maxiter A positive integer, restricting the maximum number of iterations.
#' @param eps A numeric, specifying the tolerance level of when the iteration stops.
#'
#' @return A list with class 'shapeNA' containing the following components:
#'   \item{S}{The estimated shape matrix.}
#'   \item{scale}{The scale with which the shape matrix may be scaled to obtain
#'   a scatter estimate. If `alpha` = 1, then this value is `NA`, as Tyler's
#'   shape matrix has no natural scale.}
#'   \item{mu}{The location parameter, either provided by the user or estimated.}
#'   \item{alpha}{The tail index with which the Power M-estimator has been called.}
#'   \item{naBlocks}{`NULL`, since \code{\link{powerShape}} operates only on complete data.}
#'   \item{iterations}{Number of computed iterations before convergence.}
#'   \item{call}{The matched call.}
#'
#' @export
#'
#' @references Tyler, D.E. (1987). A Distribution-Free M-Estimator of Multivariate Scatter. The Annals of Statistics, 15, 234.251. \doi{10.1214/aos/1176350263}.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569. \doi{10.1016/j.jmva.2019.104569}.
#' @references Hettmansperger, T. P., & Randles, R. H. (2002). A practical affine equivariant multivariate median. Biometrika, 89(4), 851-860. \doi{10.1093/biomet/89.4.851}
#'
#' @seealso [powerShapeNA], [tylerShapeNA] and [classicShapeNA] for the
#'     corresponding functions for data with missing values.
#'
#' @examples
#'     ## Generate example data
#'     S <- toeplitz(c(1, 0.1))
#'     x <- mvtnorm::rmvt(100, S)
#'     ## Compute some M-estimators
#'     res0 <- classicShape(x, center = c(0, 0))
#'     res1 <- powerShape(x, alpha = 0.67, normalization = 'one')
#'     res2 <- tylerShape(x, normalization = 'trace')
#'     ## Get location estimates
#'     res1$mu
#'     res2$mu
#'     ## Get shape estimates
#'     res0$S
#'     res1$S
#'     res2$S
#'     ## Print summary
#'     summary(res0)
powerShape <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (any(is.na(x))) {
    stop("Missing values found. Use powerShapeNA()")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter, eps)
    } else {
      xCentered <- sweep(x, 2, center)
      zeroEntry <- rowSums(xCentered) == 0
      if (any(zeroEntry)) {
        xCentered <- xCentered[!zeroEntry, ]
        z <- sum(zeroEntry)
        if (z==1) {
          message("Found ", z, " observation coinciding with given center.")
        } else {
          message("Found ", z, " observations coinciding with given center.")
        }
      }
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter, eps)
      res$mu <- center
    }
    res$alpha <- alpha
    res$call <- fctCall
    if (alpha == 1) {
      res$scale <- NA
    }
    return(res)
    }
  )
}
