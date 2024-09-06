# library(Rcpp)
# sourceCpp("AGammaProcess_util.cpp")

e1 <- function(x, d=1e-9) {
  return(2/d * (1-pchisq(2*x, d)))
}
e1inv <- function(y, d=1e-9) {
  return(qchisq(1-(d/2)*y, d) / 2)
}
inverse_gamma_measure <- function(x, alpha, beta) {
  1 / beta * e1inv(x / alpha)
}
Sigma_fun_eye <- function(x, d=2) { # x in [0,1]
  n <- length(x)
  if (n > 1) {
    array(diag(d), dim=c(d,d,n))
  } else {
    diag(d) * (x>=0 & x <= 1) 
  }
  }
Sigma_fun_times_omega_inv_fun <- function(x, d, Sigma_fun, omega_fun) {
  n <- length(x)
  if (n > 1) {
    S <- d * Sigma_fun(x)
    for (i in 1:n) {
      S[,,i] <- S[,,i] / omega_fun(x[i])
    }
    S
  } else {
    d * Sigma_fun(x) / omega_fun(x)
  }
}
Sigma_fun_varma <- function(x, d=2) {
  stopifnot(d==2)
  # Sigma <- matrix(c(1,0.5,-0.5,1),ncol=2,nrow=2)
  Sigma <- matrix(c(1,-0.5,-0.5,1),ncol=2,nrow=2)
  n <- length(x)
  if (n > 1) {
    psd_varma(pi*x, ar=diag(c(0.5,-0.7)), sigma=Sigma)$psd
  } else {
    psd_varma(pi*x, ar=diag(c(0.5,-0.7)), sigma=Sigma)$psd[,,1]
  }
}

# returns a function taking one argument -- scaled beta density
create_omega_fun_from_beta_density <- function(C_alpha, g0.alpha, g0.beta) {
  ret <- function(x) { C_alpha * dbeta(x,g0.alpha,g0.beta) }
  ret
}
