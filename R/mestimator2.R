mestimator_mean_cov <- function(x,  tol=1e-6, ...) {
  UseMethod("mestimator_mean_cov")
}

mestimator_mean_cov.default <- function(x, ...) {
  stop("No implementation for object of provided class. Please supply a data matrix of observations.")
}

mestimator_mean_cov.data.frame <- function(x, ...) {
  mestimator_mean_cov(as.matrix(x), ...)
}

mestimator_mean_cov.matrix <- function(x, powerfct, normalization, maxiter=1e4, tol=1e-6, ...) {
  if (tol <= 0 || maxiter <= 0) {
    stop("Nonpositive arguments maxiter or tol.")
  }
  if (any(is.na(x))) {
    return(mestimator_mean_cov.naBlocks(naBlocks(x), powerfct, normalization, maxiter, tol, ...))
  }

  x_nonCentered <- x
  centerAndKeepNonzeroObs <- function(mu) {
    x <- sweep(x_nonCentered, 2, mu)
    zeroObservation <- rowSums(x) == 0
    if (any(zeroObservation)) {
      x <- x[!zeroObservation,, drop=FALSE]
      message("Remove observation at center")
    }
    return(x)
  }

  i <- 0
  dist <- 2*tol
  n <- nrow(x_nonCentered)
  p <- ncol(x_nonCentered)
  mu <- numeric(p)
  S0 <- diag(p)
  while (i < maxiter && dist > tol) {
    x <- centerAndKeepNonzeroObs(mu)
    n <- nrow(x)
    # squared mahalanobis distances xi
    xi <- rowSums(x * t(solve(S0, t(x))))
    try(w <- powerfct(xi, p = p))
    S <- (crossprod(x, w$w*x))/n
    mu2 <-  colMeans(sweep(x, 1, w$v, "*"))

    dist <- max(norm(S-S0, "M"), max(abs(mu2)))
    S0 <- S
    mu <- mu+mu2
    i <- i+1
  }

  if (i >= maxiter) {
    warning(paste("No convergence in", i, "steps."))
  }

  shape <- normalization(S0)
  res <- list(S=shape$S, scale=shape$scale, mu=mu, alpha=NULL, iterations=i, naBlocks=NULL)
  class(res) <- "shapeNA"
  return(res)
}

mestimator_mean_cov.naBlocks <- function(x, powerfct, normalization, maxiter, tol, ...) {
  if (tol <= 0 || maxiter <= 0) {
    stop("Nonpositive arguments maxiter or tol.")
  }

  y_nonCentered <- x$data
  centerAndKeepNonzeroObs <- function(mu) {
    x <- sweep(y_nonCentered, 2, mu)
    zeroObservation <- rowSums(x, na.rm = TRUE) == 0
    if (any(zeroObservation)) {
      x <- x[!zeroObservation,, drop=FALSE]
      message("Remove observation at center")
    }
    return(x)
  }

  i <- 0
  dist <- 2*tol
  n <- nrow(y_nonCentered)
  p <- ncol(y_nonCentered)

  # use closure: access powerfct and dimension p without supplying as argument
  covAndMeanOfSubset <- function(x, S, varCount, n) {
      Sinv <- solve(S)
      if (!isSymmetric(Sinv)) {
        Sinv <- (Sinv + t(Sinv))/2
      }
      rootS <- matroot(Sinv)
      xi <- rowSums((x%*%Sinv)* x) # squared mahalanobis
      try(
        w <- powerfct(xi, p, varCount)
      )
      a <- (crossprod(x, w$w*x)) # term for Sigma
      b <- sweep(x%*%rootS, 1, w$v, FUN = "*")
      return(list(S=Sinv %*% a %*% Sinv - n * Sinv, mu=colSums(b)))
  }


  S0 <- diag(p)
  blockIdx <- x$N
  blockPattern <- x$P
  mu <- numeric(p)
  while (i < maxiter && dist > tol) {
    y <- centerAndKeepNonzeroObs(mu)
    n <- nrow(y)
    # split data into blocks of same missingness
    rows <- 1:blockIdx[1]
    tryCatch(
    a <- covAndMeanOfSubset(
      y[rows, , drop = FALSE], S0, n=length(rows), varCount=p)
      , error = function(e) {
        message(paste("Iteration", i))
        message(e)
        stop("Matrix not positive definite")
      }
    )
    a_Sigma <- a$S
    a_mu <- a$mu
    for (j in 2:length(blockIdx)) {
      rows <- (blockIdx[j-1]+1):blockIdx[j]
      cols <- asBinaryVector(blockPattern[j], p)
      tryCatch(
      a <- covAndMeanOfSubset(
        y[rows, cols, drop = FALSE], S0[cols, cols], n=length(rows), varCount=sum(cols))
      , error = function(e) {
        message(paste("Iteration", i, "block", j))
        message(e)
        stop("Matrix not positive definite")
      }
      )
      a_Sigma[cols, cols] <- a_Sigma[cols, cols] + a$S
      a_mu[cols] <- a_mu[cols] + a$mu
    }
    # adjust scatter estimate
    a_Sigma <- S0 %*% (a_Sigma/n) %*% S0
    a_mu <- (a_mu %*% matroot(S0))/n
    dist <- max(norm(a_Sigma, "M"), max(abs(a_mu)))
    S0 <- S0 + a_Sigma
    mu <- mu + a_mu
    i <- i+1
  }

  if (i >= maxiter) {
    warning(paste("No convergence in", i, "steps."))
  }

  ogOrder <- order(x$permutation)
  S0 <- S0[ogOrder, ogOrder]
  mu <- mu[ogOrder]

  shape <- normalization(S0)
  res <- list(S=shape$S, scale=shape$scale, mu=mu, alpha=NULL, iterations=i, naBlocks=x)
  class(res) <- "shapeNA"
  return(res)
}

