#' M-estimators of the Shape from the Power Family when Data is Missing
#'
#' Power M-estimators of shape and location were recently suggested in
#' Frahm et al. (2020). They have a tuning parameter `alpha` taking values in
#' \code{[0,1]}. The extreme case `alpha` = 1 corresponds to Tyler's shape
#' matrix and `alpha` = 0 to the classical covariance matrix. These special
#' cases have their own, more efficient functions \code{\link{tylerShapeNA}} and
#' \code{\link{classicShapeNA}}, respectively.
#' If the true location is known, it should be supplied as `center`, otherwise
#' it is estimated simultaneously with the shape.
#'
#' These functions assume that the data were generated from an
#' elliptical distribution, for Tyler's estimate this can be relaxed to
#' generalized elliptical distributions
#' The missingness mechanism should be MCAR or, under
#' stricter distributional assumptions, MAR. See the references for details.
# TODO MAR also okay?
#'
#' For multivariate normally distributed data, \code{\link{classicShapeNA}} is the maximum
#' likelihood estimator of the location and scale. It is a special case of the
#' power M-estimator with tail index `alpha` = 0, which returns the
#' empirical covariance matrix and the empirical mean vector.
#'
#' The function \code{\link{tylerShapeNA}} maximizes the likelihood function after projecting
#' the observed data of each individual onto the unit hypersphere, in which case
#' we obtain an angular central Gaussian distribution. It is a special case of
#' the power M-estimator with tail index `alpha` = 1, which returns Tyler's
#' M-estimator of scatter and an affine equivariant multivariate median
#' according to Hettmansperger and Randles (2002).
#'
#' The function \code{\link{powerShapeNA}} requires an additional parameter, the so-called
#' tail index `alpha`. For heavy tailed data, the index should be chosen closer
#' to 1, whereas for light tailed data the index should be chosen closer to 0.
#'
#' @aliases powerShapeNA
#' @aliases classicShapeNA
#' @aliases tylerShapeNA
#'
#' @usage powerShapeNA(x, alpha, center = NULL, normalization = c("det", "trace", "one"),
#'          maxiter = 1e4, eps = 1e-6)
#' @usage tylerShapeNA(x, center = NULL, normalization = c("det", "trace", "one"),
#'           maxiter = 1e4, eps = 1e-6)
#' @usage classicShapeNA(x, center = NULL, normalization = c("det", "trace", "one"),
#'          maxiter = 1e4, eps = 1e-6)
#'
#' @param x A data matrix or data.frame with missing data and `p` > 2 columns.
#' @param alpha Tail index, a numeric value in the interval \code{[0, 1]}.
#'     Determines the power function. For more information see 'Details'.
#' @param center An optional vector of the data's center, if `NULL` the center
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
#' @param eps A numeric, specifying tolerance level of when the iteration stops.
#'
#' @return A list with class 'shapeNA' containing the following components:
#' \describe{
#'   \item{S}{The estimated shape matrix.}
#'   \item{scale}{The scale with which the shape matrix may be scaled to obtain
#'       a scatter estimate. If `alpha` = 1, then this value will be `NA`, as
#'       Tyler's shape matrix has no natural scale.}
#'   \item{mu}{The location parameter, either provided by the user or estimated.}
#'   \item{alpha}{The tail index with which the Power M-estimator has been called.}
#'   \item{naBlocks}{An `naBlocks` object, with information about the missingness
#'       of the data.}
#'   \item{iterations}{Number of computed iterations before convergence.}
#'   \item{call}{The matched call.}
#' }
#'
#' @export
#'
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler's M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54, 374-393. \doi{10.1016/j.csda.2009.08.019}.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569. \doi{10.1016/j.jmva.2019.104569}.
#' @references Hettmansperger, T. P., & Randles, R. H. (2002). A practical affine equivariant multivariate median. Biometrika, 89(4), 851-860. \doi{10.1093/biomet/89.4.851}
#'
#' @seealso [powerShape], [tylerShape] and [classicShape] for the
#'     corresponding functions for data without missing values.
#'
#' @examples
#'     ## Generate a data set with missing values
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech = 'MCAR')$amp
#'
#'     ## Compute some M-estimators
#'     res0 <- classicShapeNA(y, center = c(0, 0, 0))
#'     res1 <- powerShapeNA(y, alpha = 0.67, normalization = 'one')
#'     res2 <- tylerShapeNA(y, normalization = 'trace')
#'
#'     ## Get location estimates
#'     res1$mu
#'     res2$mu
#'     ## Get shape estimates
#'     res0$S
#'     res1$S
#'     res2$S
#'
#'     ## Print summary
#'     summary(res0)
#'     ## Inspect missingness pattern
#'     plot(res0$naBlocks)
#'     barplot(res0)
powerShapeNA <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (!any(is.na(x))) {
    stop("No missing values found. Use powerShape().")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter,eps)
    } else {
      xCentered <- sweep(x, 2, center)
      zeroEntry <- rowSums(xCentered, na.rm = TRUE) == 0
      if (any(zeroEntry)) {
        xCentered <- xCentered[!zeroEntry, ]
        z <- sum(zeroEntry)
        if (z==1) {
          message("Found ", z, " observation coinciding with given center.")
        } else {
          message("Found ", z, " observations coinciding with given center.")
        }
      }
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter,eps)
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
