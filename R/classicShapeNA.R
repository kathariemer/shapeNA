#' @export
classicShapeNA <- function(x, center = NULL,  normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (!any(is.na(x))) {
    stop("No missing values found. Use classicShape().")
  }
  try(
    res <- powerShapeNA(x, 0, center, normalization, maxiter, eps)
  )
  res$call <- match.call()
  return(res)
}
