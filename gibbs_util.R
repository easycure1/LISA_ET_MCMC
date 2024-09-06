library(Rcpp)
# Note: Need to source 'gibbs_util.cpp'
library(MASS)
library(compiler)

print_warn <- function(msg) {
  print(msg)
  warning(msg)
}

fast_mean <- function(x) {
  sum(x) / length(x)
}

# Fourier frequencies on unit interval
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

#####
# Store list for beta mixtures - WORK IN PROGRESS
# TO DO: A more elegant solution to the storage problem
# Can come across huge memory/storage problems with large kmax and n
coarsened_bernstein <- function(omega, l) {
  res <- matrix(NA, nrow=l, ncol=length(omega))
  for (i in 1:l) {
    res[i,] <- coarsened_bernstein_i(omega, l, i)
  }
  res
}
coarsened_bernstein_i <- function(omega, l, i) {
  k <- l^2
  b_tmp <- 0 * omega
  for (j in ((i-1)*l+1):(i*l)) {
    b_tmp <- b_tmp + dbeta(omega, j, k+1-j)
  }
  b_tmp <- b_tmp / l
  b_tmp
}
dbList <- function(n, kmax, normalized=F, bernstein_l=0, bernstein_r=1, coarsened=F) {
  db.list <- vector("list", kmax)
  omega <- omegaFreq(n); NN <- length(omega)
  omega_for_dblist <- seq(bernstein_l, bernstein_r, length.out=NN)
  if (coarsened) {
    stopifnot(!normalized)
    cat("Using coarsened Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    for (kk in 1:kmax) {
      print(kk)
      db.list[[kk]] <- coarsened_bernstein(omega_for_dblist, kk)
    }  
  } else {
    cat("Using standard Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    for (kk in 1:kmax) {
      db.list[[kk]] <- matrix(dbeta(omega_for_dblist,                                                 
                                    rep(1:kk, each = NN),
                                    rep(kk:1, each = NN)),
                              ncol = NN,
                              byrow = TRUE)
      if (normalized) {
        db.list[[kk]] <- db.list[[kk]] / kk
      }
    }  
  }
  return(db.list)
}
#####

#' Compute F_n X_n with the real-valued Fourier matrix F_n
fast_ft <- compiler::cmpfun(function(x, real=T) {
  # Function computes FZ (i.e. fast Fourier transformed data)
  # Outputs coefficients in correct order and rescaled
  n <- length(x)
  sqrt2 <- sqrt(2)
  sqrtn <- sqrt(n)
  # Cyclically shift so last observation becomes first
  x <- c(x[n], x[-n])  # Important since fft() uses 0:(n-1) but we use 1:n
  # FFT
  fourier <- fft(x)
  if (real) {
    # Extract non-redundant real and imaginary coefficients in correct order and rescale
    FZ <- rep(NA, n)
    FZ[1] <- Re(fourier[1]) # First coefficient
    if (n %% 2) {
      N <- (n-1)/2
      FZ[2*(1:N)] <- sqrt2 * Re(fourier[2:(N+1)]) # Real coefficients
      FZ[2*(1:N)+1] <- sqrt2 * Im(fourier[2:(N+1)]) # Imaginary coefficients
    } else {
      FZ[n] <- Re(fourier[n / 2 + 1]) # Last coefficient
      FZ[2 * 1:(n / 2 - 1)] <- sqrt2 * Re(fourier[2:(n / 2)]) # Real coefficients
      FZ[2 * 1:(n / 2 - 1) + 1] <- sqrt2 * Im(fourier[2:(n / 2)]) # Imaginary coefficients
    }
  } else {
    N <- ifelse(n %% 2, (n+1)/2, n/2+1)
    FZ <- fourier[1:N]
  }
  return(FZ / sqrtn)
})

#' Compute F_n^t X_n with the real-valued Fourier matrix F_n
fast_ift <- compiler::cmpfun(function(x, real=T, TOL=1e-15) {
  # Function computes inverse Fourier transform
  # Can be used for finding FCFZ
  if (real) {
    n <- length(x)
    sqrtn <- sqrt(n)
    sqrtn2 <- sqrt(n / 2)
    # Construct complex vector
    CFZ <- rep(NA, n)
    CFZ[1] <- x[1] * sqrtn
    if (n %% 2) {
      N <- (n-1)/2
      CFZ[2:(N+1)] <- (x[2 * (1:N)] + 1i * x[2 * (1:N)+1] ) * sqrtn2
      CFZ[(N+2):n] <- rev(Conj(CFZ[2:(N+1)])) # Include complex complex conjugates
    } else {
      CFZ[n / 2 + 1] <- x[n] * sqrtn
      CFZ[2:(n / 2)] <- (x[2 * (1:(n / 2 - 1))] + x[2 * (1:(n / 2 - 1)) + 1] * 1i) * sqrtn2
      CFZ[(n / 2 + 2):n] <- rev(Conj(CFZ[2:(n / 2)])) # Include complex complex conjugates
    }
  } else {
    N <- length(x)
    n_is_even <- (abs(Im(x[N])) < TOL)
    n <- ifelse(n_is_even, 2*(N-1), 2*N-1)
    CFZ <- c(x, rev(Conj(x[-c(1,N)]))) * sqrt(n)
  }
  # Inverse FFT (normalised)
  Z <- fft(CFZ, inverse = TRUE) / n
  # Cyclically shift
  Z <- c(Z[-1], Z[1])
  if (real) {
    Z <- Re(Z)
  }
  return(Z)
})

#define density of t-distribution in terms of excess kurtosis
dtex.kurt <- function(x, ex.kurt) {
  nu <- 6 / ex.kurt + 4
  dt(x, nu)
}

uniformmax <- function(sample) {
  max(abs(sample - median(sample)) / mad(sample), na.rm=T)
}

logfuller<-function(x, xi = 0.001){
  log(x + xi) - xi / (x + xi)  
}

logrosenthal <- function(x, inverse=F) {
  if (!inverse) return(sign(x) * log(1 + abs(x)))
  else return(sign(x) * (exp(abs(x)) - 1))
}

se_kurt <- function(n) {
  
  # Standard error for kurtosis
  return(2 * sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3))) * 
           sqrt((n ^ 2 - 1) / ((n - 3) * (n + 5))))
  
}

# negative log likelihood of iid standard normal observations [unit variance]
nll_norm <- function(epsilon_t, ...) {
  m <- length(epsilon_t)
  cll <-  1 / 2 * (sum(epsilon_t ^ 2) + m * log(2*pi))
  return(cll)
}

# unnormalized negative log likelihood of iid standard normal observations [unit variance]
nll_norm_unnormalized <- function(epsilon_t, ...) {
  cll <-  1 / 2 * sum(epsilon_t ^ 2)
  return(cll)
}

# negative log likelihood of iid t observations with given excess kurtosis [unit variance]
nll_t <- function(epsilon_t, ex.kurt, ...) {
  nu <- (6 / ex.kurt) + 4
  sigma2 <- nu / (nu - 2)
  cll <- -sum(dt(epsilon_t * sqrt(sigma2),df=nu,log=T)) - log(sigma2)/2
  return(cll)
}

generalizedGaussian.alpha <- function(beta) { # for unit variance
  lalpha <- 1 / 2 * (lgamma(1/beta) - lgamma(3/beta))
  alpha <- exp(lalpha)
  return(alpha)
}
generalizedGaussian.kurtosis <- function(beta) {
  log.kurt <- lgamma(5 / beta) + lgamma(1 / beta) - 2*lgamma(3 / beta)
  ex.kurt <- exp(log.kurt) - 3
  return(ex.kurt)
}
l_generalizedGaussian <- function(x, beta) {
  alpha <- generalizedGaussian.alpha(beta)
  llik <- log(beta) - log(2*alpha) - lgamma(1/beta) - (abs(x) / alpha)^beta
  return(llik)
}
nll_generalizedGaussian <- function(epsilon_t, beta) {
  cll <- -sum(l_generalizedGaussian(epsilon_t, beta))
  return(cll)
}

plotPsdEstimate <- function(mcmc, lambda, psd.true) {
  plot(lambda, psd.true, type = "l", ylim = c(0, 2*max(psd.true)))
  lines(lambda, mcmc$fpsd.s, col = 2, lty = 2)
  lines(lambda, mcmc$fpsd.s05, col = 2, lty = 2)  # PW CI
  lines(lambda, mcmc$fpsd.s95, col = 2, lty = 2)  # PW CI
  lines(lambda, mcmc$log.confupper, col = 4, lty = 2)  # Uniform CI
  lines(lambda, mcmc$log.conflower, col = 4, lty = 2)  # Uniform CI
  title(mcmc$algorithm)
}

reduceMemoryStorageMCMC <- function(mcmc) { # Discard memory intensive traces
    ret <- mcmc
    ret$fpsd.sample <- NULL
    
    ret$pdgrm <- NULL
    ret$V <- NULL
    ret$W <- NULL
    ret$df <- NULL
    
    ret$cy <- NULL
    ret$mixtureTrace <- NULL
    ret$phi_mu <- NULL
    ret$phi_sigma2 <- NULL
    
    return(ret)
  }

logDet_stickBreaking <- function(v) {
  L <- length(v)
  sum(log(1-v)*((L-1):0))
}

my_rdirichlet <- function(alpha) {
  ga <- rgamma(length(alpha), alpha, rep(1,length(alpha)))
  ga / sum(ga)
}

my_ddirichlet_unnormalized <- function(x, alpha, log) {
  ld <- sum(x*(1-alpha))
  if (log) {
    res <- ld
  } else {
    res <- exp(ld)
  }
  res
}
