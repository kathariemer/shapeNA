#' Compute shape estimate for full data
#'
#' Given a data matrix \code{x} from a continuous distribution, return a shape estimate and, if not supplied as \code{center}, a location estimate.
#'
# aliases could also be in a comma separated list
#' @aliases powerShape
#' @aliases tylerShape
#' @aliases classicShape
#'
#' @usage powerShape(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6)
#' @usage tylerShape(x, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6)
#' @usage classicShape(x, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6)
#'
#' @param x numeric data matrix or data.frame without missing data. Representing sample from continuous distribution
#' @param alpha numeric, determines power function
#' @param center optional vector of center, if NULL the center will be estimated simultaneously to the shape estimate
#' @param normalization string, determines scale of returned shape estimate
#' @param maxiter integer, maximum number of itreations
#' @param tol numeric, tolerance level
#'
#' @return a shapeNA object, which contains a shape and center estimate
#' @export
#'
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler’s M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54(2), 374-393.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569.
#'
#' @seealso powerShapeNA
#' @seealso tylerShapeNA
#' @seealso classicShapeNA
#'
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out=5)))
#'     res <- powerShape(x, alpha=0.67, normalization='one')
#'
powerShape <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6) {
  if (any(is.na(x))) {
    stop("Missing values found. Use powerShapeNA()")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter, tol)
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
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter, tol)
      res$mu <- center
    }
    res$alpha <- alpha
    res$call <- fctCall
    return(res)
    }
  )
}
