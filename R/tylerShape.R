#' @export
tylerShape <- function(x, center = NULL,  normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6) {
  if (any(is.na(x))) {
    stop("Missing values found. Use tylerShapeNA().")
  }
  try(
    res <- powerShape(x, 1, center, normalization, maxiter, tol)
  )
  res$call <- match.call()
  return(res)
}
