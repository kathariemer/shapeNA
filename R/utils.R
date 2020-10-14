# format data / preprocessing for Mestimator
# s.t. first row and first column have most responses
# and last row and last column has least number of responses
# with monotone increase
# @param data nxp matrix
# @param cleanup logical flag, if true, observations with less than 2 responses
#    are discarded
# @param plot logical flag
# @return list with
#    data, which has been reordered and cleaned up
#    permutation, which indicates new ordering of the variables
#    N - signifying the rows in which the missingness pattern changes
#    D - signifying the number of responses per observation in the corresponding block
#    P - signifying missingness pattern
naBlocks <- function(data, cleanup=TRUE, plot=FALSE) {
  ## remove observations without responses or only 1 response ##
  x <- data
  if (cleanup) {
    hasValue <- !is.na(x) # R matrix
    #rmIdx <- apply(hasValue, 1, function(x) sum(x) <=1)
    rmIdx <- rowSums(hasValue)<=1
    if (any(rmIdx)) {
      x <- data[!rmIdx, ]
      message(paste("Removed", sum(rmIdx), "observations with less than 2 responses each."))
    }
  }
  ## compute dimensions of cleaned up data ##
  n <- nrow(x)
  p <- ncol(x)
  if (n < 2) {
    stop(paste("Only", n, "observations."))
  } else if (n <= p) {
    stop('Dimension exceeds sample size.')
  }
  ## reorder data to form contiguous blocks with the same pattern ##
  hasValue <- !is.na(x) # R matrix
  respPerVar <- colSums(hasValue)
  xIdx <- order(respPerVar, decreasing = TRUE)
  x <- x[, xIdx]
  hasValue <- hasValue[, xIdx]

  ## get indices where patterns change and the number of observations/the corresponding pattern ##
  blockInfo <- getMissingnessBlocks(hasValue)
  rle <- blockInfo$rle
  N <- cumsum(rle$lengths)
  if (length(N) == 1) {
    warning("No missing data.\n")
    return(x)
  }
  D <- rowSums(hasValue[blockInfo$order[N],])
  P <- rle$values  # save pattern info

  x <- x[blockInfo$order, ]

  if (N[1] < D[1]) {
    # first block has fewer responses than rows
    cat("N: ", N, "\nD: ", D, "\n")
     stop('Missingness is too strong! Estimator cannot converge.')
  }
  res <- list(data=as.matrix(x), permutation=xIdx, rowPermutation=blockInfo$order,
              N=N, D=D, P=P, kn=respPerVar/n, kp=respPerVar/p)
  class(res) <- "naBlocks"
  if (plot) {
    plot(res)
  }
  return(res)
}


#' plot missingness pattern of data
#'
#' @export
plot.naBlocks <- function(x, orderProp=TRUE, ...) {
  idx <- x$N
  bprop <- (idx - c(0, idx[-length(idx)]))/nrow(x$data)
  if (orderProp){
    bIdx <- order(bprop, decreasing = TRUE)
    bprop <- bprop[bIdx]
    plotAsBoolMat(x$data[idx[bIdx],], round(bprop,3))
  } else {
    bIdx <- 1:length(idx)
    plotAsBoolMat(x$data[idx,], round(bprop,3))
  }
  invisible(list(blockProps=bprop, blockIdx=bIdx))
}

# given a number n, return its binary representation as a vector
# @param n decimal number
# @param d length of vector
# @param logical whether to return a numeric or logical vector
asBinaryVector <- function(n, d, logical = TRUE) {
  b <- numeric(d)
  for (i in d:1) {
    b[i] <- n %% 2
    n <- n %/% 2
  }
  if (n > 0) {
    warning(paste("Number was larger than 2^",d,"-1. Allow for ~",floor(log2(n))+1," more digits.\n", sep=''))
  }
  if (logical) {
    b <- as.logical(b)
  }
  return(b)
}

# Given a integer, which represents a binary pattern, return a vector, which defines the corresponding indices
#
# @param x integer representing a binary pattern
# @param d maximal exponent, i.e. 2^d > x
#
# @return subset of 1:d, which selects the entries corresponding to coefficient 1 of binary representation of x
binaryToIndex <- function(x, d) {
  b <- numeric(d)
  for (i in d:1) {
    b[i] <- x %% 2
    x <- x %/% 2
  }
  if (x > 0) {
    warning(paste("Number was larger than 2^",d,"-1. Allow for ~",floor(log2(x))+1," more digits.\n", sep=''))
  }
  b <- (1:d)[b==1]
  return(b)
}

# plot boolean matrix as a grid of colored rectangles
# blue fields are TRUE, white fields are FALSE, yellow fields are NA
# @param R boolean matrix
plotAsBoolMat <- function(R, blockNames, ...) {
  # R <- obj$pattern
  color2 <- colors["white"]
  if (any(is.na(R))) {
    color2 <- colors["yellow"]
    B <- is.na(R)
  } else {
    B <- R == 0
  }
  r <- nrow(R)
  c <- ncol(R)
  plot(1, type='n', xlim=c(0,c+1), ylim = c(r+1,0), axes=FALSE, xlab=NA, ylab=NA)
  xleft <- rep(1, r)
  ybottom <- 0:(r-1)
  for (i in 1:c) {
    rowColor <- ifelse(B[,i], color2, colors["blue"])
    rect(xleft, ybottom, xleft+1, ybottom+1, col = rowColor, ...)
    xleft <- xleft+1
  }
  text(0.5+(1:c),r+1, colnames(R))
  if (!missing(blockNames)){
    text(0, (1:r)-0.5, blockNames)
  }
}

# transform numeric input to have values between 0 and 1
# the largest value gets mapped to 0 (black), the smallest to 1 (white)
grayscale <- function(A) {
  m <- min(A)
  M <- max(A)
  return((A-m)/(m-M)+1)
}

#' print method for elements of class shapeNA
#'
#' Only print M-estimates and alpha level
#'
#' @export
print.shapeNA <- function(obj) {
  if (is.null(obj$mu)) {
    print(list(alpha=obj$alpha, S=obj$S))
  } else {
    print(list(alpha=obj$alpha, S=obj$S, mu=obj$mu))
  }
}

#' Crude visualization of shape estimate
#'
#' If estimate comes from missing data, additionally the columns are
#' marked with a colored bar, indicating their missingness proportion
#'
#' @export
plot.shapeNA <- function(obj, legend=TRUE, message=TRUE) {
  S <- grayscale(obj$S)
  p <- nrow(S)
  plot(1, type = "n", xlim = c(0, p*1.25), ylim = c(0, p+1),
       axes=FALSE, xlab=NA, ylab=NA)
  for (i in 1:p) {
    for (j in 1:p) {
      rect(i-1, j-1, i, j, col=gray(S[i,j]))
    }
  }
  if (!is.null(obj$naBlocks)) {
    k <- obj$naBlocks$kn
  topocolors <- rainbow(12)[1:10]
  for (i in 1:p) {
    rect(i-k[i], p+0.5, i, p+1, col=topocolors[ceiling(10*k[i])])
  }
  if (legend) {
    legend(p+0.3, p+0.8, (10:1)/10, col=rev(topocolors),
           horiz = FALSE, pch=20, ncol=2)
  }
  if (!is.null(obj$naBlocks) && message) {
    message(paste("% of total obs:", paste(round(obj$naBlocks$kn,2), collapse = "\t")))
  }
  }
}

tr <- function(X) {
  return(sum(diag(X)))
}

# return a scaled scale estimator S, which satisfies trace(S) = p
scaleTr <- function(Sigma) {
  #s <- nrow(Sigma)/sum(diag(Sigma))
  s <- tr(Sigma)/nrow(Sigma)
  if (s == 0) {
    stop("Cannot scale matrix, trace equals 0")
  }
  return(list(S = Sigma/s, scale=s))
}

scaleDet <- function(Sigma) {
  d <- det(Sigma)^(1/ncol(Sigma))
  if (d == 0) {
    stop("Cannot scale matrix, determinant equals 0")
  }
  return(list(S = Sigma/d, scale=d))
}

scaleOne <- function(Sigma) {
  s <- Sigma[1,1]
  if (s == 0) {
    stop("Cannot scale matrix, first entry equals 0")
  }
  return(list(S = Sigma/Sigma[1,1], scale = s))
}

colors <- c(red = "#d55e00",
            purple = "#cc79a7",
            blue = "#0072b2",
            yellow = "#FFCD70",
            green = "#009e73",
            black = "#000000",
            white = "#ffffff")

# use eigen decomposition to compute square root of matrix
# @param Sigma symmetric, pos-definite matrix
matroot <- function(Sigma) {
  A <- eigen(Sigma)
  Gamma <- A$vectors
  D <- diag(sqrt(A$values))
  return(tcrossprod(Gamma%*%D, Gamma))
}

# runlength encoding for missingness matrix, where
# the rows get transformed to an integer corresponding
# to it's missingness pattern (big endian)
# Example: c(1,1,0) <-> 4\*1 + 2\*1 + 1\*0
# @param R reordered missingness matrix
# @return permutation for reordering of data and runlength encoded integer pattern,
#   i.e a missingness of (0,1,0,1) translates to 8*0 + 4*1 + 2*0 * 1*1 = 5
getMissingnessBlocks <- function(R) {
  p2 <- 2^((ncol(R)-1):0)
  pattern <- rowSums(sweep(R, 2, p2, "*"))
  lev <- summary(as.factor(pattern))
  patternF <- factor(pattern, levels = names(sort(lev, decreasing = TRUE)))
  #rowOrder <- order(rowSums(R), pattern, decreasing = c(TRUE, FALSE))
  rowOrder <- order(patternF)
  return(list(order=rowOrder, rle = rle(pattern[rowOrder])))
}

#' Print missingness pattern
#'
#' @export
print.naBlocks <- function(obj) {
  p <- length(obj$permutation)
  pattern <- sapply(obj$P, asBinaryVector, d=p, logical=FALSE)
  res <- matrix(c(t(pattern), D=p-obj$D), nrow=length(obj$N))
  rCount <- c(obj$N, 0) - c(0, obj$N)
  row.names(res) <- rCount[1:length(obj$N)]
  colNames <- c(paste("V", obj$permutation, sep=""), " #")
  colnames(res) <- colNames
  print(res)
}

#' Barplot showcasing missingness proportion of the original data
#'
#' @export
barplot.naBlocks <- function(obj, sortNA=FALSE) {
  mprop <- colSums(is.na(obj$data))/nrow(obj$data)
  if (!sortNA) {
    mprop <- mprop[order(obj$permutation)]
  }
  barplot(mprop)
  invisible(mprop)
}

#' summary method for class shapeNA
#'
#' @param obj an object of class shapeNA, usually from a call to powerShape or similar functions
#' @param ... further arguments
#'
#' @return object of class shapeNA
#' @export
#'
#' @examples
#'     obj <- tylerShape(mvtnorm::rmvt(100, diag(3)))
#'     summary(obj)
summary.shapeNA <- function(obj, ...) {
  class(obj) <- 'summary.shapeNA'
  return(obj)
}

#' print method for class summary.shapeNA
#'
#' @param obj object returned from summary.shapeNA
#' @param ... further arguments
#'
#' @return invisibly return NULL
#' @export
#'
print.summary.shapeNA <- function(obj, ...) {
  cat('Call: ')
  print(obj$call)
  cat('---\nalpha:', obj$alpha)
  cat('\nCenter:\n')
  print(obj$mu, digits=3)
  cat('\nShape:\n')
  print(obj$S, digits=3)
  if (!is.null(obj$naBlocks)) {
    cat('Missingness pattern:\n')
    print(obj$naBlocks)
  }
  cat('---\nComputed', obj$iterations, 'iterations.\n')
  invisible(NULL)
}

# get power function to pass to powerShape/powerShapeNA
powerFunction <- function(alpha) {
  if (alpha > 1 || alpha < 0) {
    stop("Ensure 0 <= alpha <= 1.")
  }
  else if (alpha == 1) {
    powerfct <- power_weight_tyler
  } else if (alpha==0) {
    powerfct <- power_weight_gauss
  } else {
    powerfct <-
      function(xi, p, m)
        power_weight(xi, p, m, alpha = alpha)
  }
  return(powerfct)
}

# get normalization function to pass to powerShape/powerShapeNA
normalizationFunction <- function(normalization) {
  ## determine normalization of scatter estimate ##
  if (length(normalization) >= 1) {
    normalization <- normalization[1]
  }
  normalization <- pmatch(tolower(normalization), c("det", "trace", "one"))
  if (normalization == 1) {
    scatterNormFct <- scaleDet
  } else if (normalization == 2) {
    scatterNormFct <- scaleTr
  } else if (normalization == 3) {
    scatterNormFct <- scaleOne
  } else {
    stop("Unknown method of normalization")
  }
  return(scatterNormFct)
}

asCov <- function(obj) {
  if (obj$alpha == 1) {
    stop("Covariance matrix for Tyler's M-estimate undefined")
  }
  return(obj$scale * obj$S)
}
