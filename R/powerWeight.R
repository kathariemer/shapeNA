## Weight functions ##

# Power weight function for complete and incomplete data
# Generalization of Gaussian and tyler's weight function
#
# @param xi vector of mahalanobis distances
# @param p dimension of data
# @param m dimension of available data, if missing the complete case is assumed
#
# @return vector of weights for scatter (w) and location estimate (v)
#
# @examples
power_weight <- function(xi, p, m, alpha) {
  if (missing(m) || p == m) {
    ## complete data case ##
    return(list(w=(p/xi)^alpha, v=xi^(-alpha/2)))
  } else {
    ## missing data case ##
    s1 <- (m/2)+1
    s2 <- (p-m)/2
    return(list(w=beta(s1,s2)/beta(s1-alpha,s2)*(p/xi)^alpha, v=xi^(-alpha/2)))
  }
}

# Gaussian weight function for complete and incomplete data
# special case of power weight function with alpha=0
#
# @param xi vector of mahalanobis distances
# @param p dimension of data
# @param m dimension of available data, if missing the complete case is assumed
#
# @return vector of weights for scatter (w) and location estimate (v)
#
# @examples
power_weight_gauss <- function(xi, p, m) {
  nu <- rep(1, length(xi))
  return(list(w=nu, v=nu))
}

# Tyler's weight function for complete and incomplete data
# special case of power weight function with alpha=1
#
# @param xi vector of squared mahalanobis distances
# @param p dimension of data
# @param m dimension of available data, if missing the complete case is assumed
#
# @return vector of weights for scatter (w) and location estimate (v)
#
# @examples
power_weight_tyler <- function(xi, p, m) {
  if (any(xi <= 0)) {
    stop(paste("Nonpositive weights", paste(round(xi,3), collapse = " ")))
  }
  nu <- 1/xi
  if (missing(m) || p == m) {
    return(list(w=p*nu, v=sqrt(nu)))
  } else {
    s1 <- m/2
    s2 <- (p-m)/2
    return(list(w=beta(s1+1, s2)/beta(s1, s2) * (p*nu), v=sqrt(nu)))
  }
}

