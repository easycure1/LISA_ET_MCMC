## TODO's
# - Implement real-valued periodogram, psd and Whittle likelihood (see Section 0.1.2. in draft)
# - remove MTS package dependency (VARMAcov_muted)
##

# + model selection with scree plot
# + multivariate periodogram from data
# + multivariate ARMA psd and transfer function from coefficients
# + Likelihoods
# + bayesian VAR(p)
# + coherency, phase spectrum etc
# + encapsulate the nonparametric prior (see 2d_matrixDirichlet)
# + correction matrix in frequency domain

# Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # lambda function using C++11
# library(Rcpp)
# source("/home/alexander/Code/beyondWhittle/gibbs_util.R")
# sourceCpp("gibbs_multivariate_util.cpp")

# Fourier frequencies on unit interval
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

Adj <- function(m) {
  return(t(Conj(m)))
}

center <- function(x) {
  return(x - mean(x))
}

eye <- function(d) {
  return(diag(rep(1,d)))
}

sqrt_hpd <- function(A) {
  # assert: A is hpsd
  ee <- eigen(A)
  U <- ee$vectors
  Lambda_half <- diag(sqrt(ee$values))
  U %*% Lambda_half %*% Adj(U)
}

##
## Convert vector parametrization (beta) to matrix-parametrization (phi),
## as e.g. used in MTS::VAR()$ar
##
phiFromBeta_normalInverseWishart <- function(beta, K, p) {
  return(t(matrix(data=beta, nrow=K*p, ncol=K)))
}

##
## Prior covariance matrix of beta (autoregression coefficent vector) 
## according to Minnesota prior
##
V_prior_minnesota <- function(Z, p, a1, a2, a3=NULL) {
  K <- ncol(Z)
  if (!is.null(a3)) {
    stop("No support for exogeneous variables yet.")
  }
  sigma_ii <- apply(Z, 2, function(ts) {
    return(as.numeric(ar(ts,order.max=p,aic=F,method="ols")$var.pred))
  })
  v_diag <- rep(NA, K*K*p)
  i <- 1
  for (m in 1:K) {
    for (q in 1:p) {
      for (mm in 1:K) {
        if (mm==m) {
          v_diag[i] <- a1 / (q^2) # TODO p^2 here?
        } else {
          v_diag[i] <- a2 / (q^2) * sigma_ii[m] / sigma_ii[mm] # TODO p^2 here?
        }
        i <- i+1
      }
    }
  }
  V_prior <- diag(v_diag)
  return(V_prior)
}
##
## Regressor matrix of VAR model
##
VAR_regressor_matrix <- function(data, var.order) {
  K <- ncol(data)
  n <- nrow(data)
  p <- var.order
  ZZ <- NULL 
  for (t in (p+1):n) { 
    z_mt <- c(t(data[(t-1):(t-p),]))
    for (m in 1:K) {
      ZZ <- rbind(ZZ, c(
        rep(0, (m-1)*K*p),
        z_mt,
        rep(0, (K-m)*K*p)))
    }
  }
  ZZ
}


##
## 1D DFT (real or complex valued formulation)
##
fast_ft_1d <- function(z, real=F) {
  stop("Deprecated --> Use fast_ift in gibbs_util.R instead")
  stopifnot(!(length(z)%%2)) # Uneven lengths not supported yet. TODO: Also adjust llike_whittle to uneven lengths
  n <- length(z)
  N <- n / 2 + 1
  # Cyclically shift
  z <- c(z[n], z[-n])
  fourier <- fft(z) / sqrt(n)
  if (real) {
    FZ <- rep(NA, n)
    FZ[1] <- Re(fourier[1])
    FZ[n] <- Re(fourier[N])
    FZ[2 * 1:(n / 2 - 1)] <- sqrt(2) * Re(fourier[2:(n / 2)])
    FZ[2 * 1:(n / 2 - 1) + 1] <- sqrt(2) * Im(fourier[2:(n / 2)])
  } else {
    FZ <- fourier[1:N]
  }
  return(FZ)
}
fast_ift_1d <- function(FZ, real=F) {
  stop("Deprecated --> Use fast_ift in gibbs_util.R instead")
  if (real) {
    n <- length(FZ)
    N <- n / 2 + 1
    # Construct complex vector
    CFZ <- rep(NA, n)
    CFZ[1] <- FZ[1] * sqrt(n)
    CFZ[N] <- FZ[n] * sqrt(n)
    CFZ[2:(N-1)] <- (FZ[2 * (1:(n / 2 - 1))] + FZ[2 * (1:(n / 2 - 1)) + 1] * 1i) * sqrt(n / 2)
    CFZ[(n / 2 + 2):n] <- rev(Conj(CFZ[2:(n / 2)]))
  } else {
    N <- length(FZ)
    n <- 2 * (N-1)
    CFZ <- c(FZ, rev(Conj(FZ[-c(1,N)]))) * sqrt(n)
  }
  Z <- fft(CFZ, inverse = TRUE) / n
  # Cyclically shift
  Z <- c(Z[-1], Z[1])
  if (real) {
    Z <- Re(Z)
  }
  return(Z)
}

##
## Multivariate discrete (fast) Fourier Transform
##
mdft <- function(Z, real=F) {
  FZ <- apply(Z, 2, fast_ft, real=real)
  return(FZ)
}
midft <- function(FZ, real=F) {
  Z <- apply(FZ, 2, fast_ift, real=real)
  return(Z)
}
vech <- function(A) { 
  # vectorize \mathbb R^{n \times d} to \mathbb R^{nd}
  # (see \tilde x_n in Proof of Theorem 42 in draft_20171127)
  c(t(A))
}
ivech <- function(v, d) {
  # Invese transformation of vech (dimension d needs to be given)
  stopifnot(length(v) %% d == 0)
  n <- length(v) / d
  matrix(v, nrow=n, ncol=d, byrow=T) 
}

##
## Compute Periodgram matrix (given the DFT)
##
mpdgrm <- function(FZ) {
  N <- nrow(FZ)
  d <- ncol(FZ)
  res <- array(data=NA, dim=c(d,d,N))
  for (j in 1:N) {
    res[,,j] <- FZ[j,] %*% Adj(FZ[j,])
  }
  res <- res / 2 / pi
  return(res)
}

##
## VARMA spectral density
## 
## lambda: frequency to evaluate
## Ar: coeffient matrix (d times p*d)
## Ma: coeffient matrix (d times q*d)
## Sigma: innovation variance (d times d)
## (See Brockwell/Davis, $11.5)
psd_varma <- function(lambda, 
                      ar=matrix(nrow=nrow(sigma),ncol=0), 
                      ma=matrix(nrow=nrow(sigma),ncol=0), 
                      sigma) { # TODO: CPP
  d <- nrow(sigma)
  N <- length(lambda)
  stopifnot(nrow(ar)==d && !(ncol(ar)%%d))
  stopifnot(nrow(ma)==d && !(ncol(ma)%%d))
  transfer_ar <- transfer_polynomial(lambda, -ar) # note the minus
  transfer_ma <- transfer_polynomial(lambda, ma)
  psd <- varma_transfer2psd(transfer_ar, transfer_ma, sigma)
  return(list(psd=psd,
              transfer_ar=transfer_ar,
              transfer_ma=transfer_ma))
}

##
## Whittle likelihood (unnormalized)
## See Theorem 42 in draft_20171127
##
#llike_whittle <- function(FZ, fpsd, excludeBoundary=T) { # Complex multivariate Whittle likelihood
#  N <- nrow(FZ)
#  boundaryFrequecies <- c(1,N)
#  ll <- sldcmvnorm_t(FZ[-boundaryFrequecies,], 2*pi*fpsd[,,-boundaryFrequecies]) # ignore boundary frequencies
#  if (!excludeBoundary) {
#    ll <- ll + sldmvnorm_t(Re(FZ[boundaryFrequecies,]),
#                           2*pi*Re(fpsd[,,boundaryFrequecies])) # add real-valued boundary part
#  }
#  return(ll)
#}

##
## VAR(p) partial likelihood (unnormalized)
## Note: Fine for fixed p, but not suited for model comparison
##
llike_var_partial <- function(zt, ar, sigma) {
  epsilon_t <- epsilon_var(zt, ar)
  #ll <- sum(mvtnorm::dmvnorm(epsilon_t,sigma=sigma,log=T))
  ll <- sldmvnorm(epsilon_t, sigma)
  return(ll)
}

##
## VAR(p) full log likelihood (including marginals, normalized)
##
llike_var_full <- function(zt, ar, sigma) {
  n <- nrow(zt)
  d <- ncol(zt)
  p <- ncol(ar) / nrow(ar)
  epsilon_t <- epsilon_var(zt, ar)
  cll <- sldmvnorm(epsilon_t, sigma)
  zt_p <- c(t(zt[1:p,]))
  gamma_p <- VARMAcov_muted(Phi=ar, Sigma=sigma, lag=p-1)$autocov[,(1:(d*p))]
  Gamma_p <- acvToeplitz(gamma_p)
  # while(det(Gamma_p) < 1e-6) {
  #   Gamma_p <- Gamma_p + 1e-6*diag(d*p)
  #   warning("Gamma_p singular")
  # }
  Gamma_p_inv <- solve(Gamma_p)
  mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
  mll <- -log(det(Gamma_p)) / 2 + mll_unnormalized
  return(cll + mll)
}
VARMAcov_muted <- function (Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120) 
{ # code taken from the MTS package, but muted (removed output)
  warning("Please reduce dependency of MTS package entirely")
  m1 = MTS::PSIwgt(Phi = Phi, Theta = Theta, lag = trun, plot = FALSE)
  Psi = m1$psi.weight
  nc = dim(Psi)[2]
  k = dim(Psi)[1]
  if (is.null(Sigma)) {
    wk = Psi
  }
  else {
    wk = NULL
    for (i in 0:trun) {
      ist = i * k
      wk = cbind(wk, Psi[, (ist + 1):(ist + k)] %*% Sigma)
    }
  }
  Gam0 = wk %*% t(Psi)
  SE = diag(1/sqrt(diag(Gam0)))
  covmtx = Gam0
  cormtx = SE %*% Gam0 %*% SE
  for (i in 1:lag) {
    ist = i * k
    Gami = wk[, (ist + 1):nc] %*% t(Psi[, 1:(nc - ist)])
    covmtx = cbind(covmtx, Gami)
    cormtx = cbind(cormtx, SE %*% Gami %*% SE)
  }
  # for (i in 0:lag) {
  #   ist = i * k
  #   cat("Auto-Covariance matrix of lag: ", i, "\n")
  #   print(round(covmtx[, (ist + 1):(ist + k)], 5))
  # }
  # for (i in 0:lag) {
  #   ist = i * k
  #   cat("cross correlation matrix of lag: ", i, "\n")
  #   print(round(cormtx[, (ist + 1):(ist + k)], 4))
  # }
  VARMAcov <- list(autocov = covmtx, ccm = cormtx)
}

##
## Corrected VAR likelihood (frequency domain): f- and q-parametrization (see maths below)
##
llike_var_corrected <- function(FZ, ar, f_param_half, sigma, f, excludeBoundary=T) { # exclude lambda_{0,N} for mean-centered TS
  d <- ncol(FZ)
  f_half <- chol_cube(f, excludeBoundary) 
  f_half_inv <- inv_cube(f_half, excludeBoundary)
  CFZ <- get_CFZ(FZ, f_half_inv, f_param_half, excludeBoundary)
  FCFZ <- Re(midft(CFZ))
  ll <- 2*sum(logdet_cube(f_half_inv, excludeBoundary)) + # times 2, because of functional determinant in real-valued formulation
    2*sum(logdet_cube(f_param_half, excludeBoundary)) + # times 2, because of functional determinant in real-valued formulation
    llike_var_partial(FCFZ, ar, sigma=sigma)
  return(ll)
}
##
## --- Some maths: ---
## CFZ = L_param L_f^{-1} FZ     // term in corrected likelihood
## L_param Q L_param^* = f       // prior specification -- see Jentsch/Kreiss
## L_f = L_param L_Q             // closedness-under-multiplication of lower triangular + uniqueness of Chol
## L_param L_f^{-1} = L_param (L_param L_Q)^{-1}
## --------------------
##
llike_var_corrected_q <- function(FZ, ar, f_param_half, f_param_half_trans, sigma, q, excludeBoundary=T) { # see notes 20170223 for q(=tildeQ)
  d <- ncol(FZ)
  ll_t <- 0
  for(i in 1:10){
    CFZ <- get_CFZ_q(FZ[,,i], q, f_param_half, excludeBoundary)
    FCFZ <- Re(midft(CFZ))
    ll <- 2 * sum(-logdet_cube(q,excludeBoundary)/2) + # times 2, because of functional determinant in real-valued formulation
      llike_var_partial(FCFZ, ar, sigma=sigma)
    ll_t <- ll_t + ll
  }
  return(ll_t)
}


##
## Complex inverse Wishart distribution
##
rciWishart <- function(nu, Psi_half) {
  Sigma_inv_half <- solve(Psi_half)
  return(solve(rcWishart(nu, Sigma_inv_half)))
}

##
## Generate hermitian (component-wise) normal matrices, with mean 0
## Note: Not to be confused with a Gaussian ensemble!!
## Used as MH-proposal for positive definite weights
##
rcomplexMatrixNorm <- function(d, sd) {
  innov <- (rnorm(d*d, mean=0, sd=sd) + 1i*rnorm(d*d, mean=0, sd=sd)) / sqrt(2)
  M <- matrix(data=innov, nrow=d, ncol=d)
  return(M + Adj(M))
}

##
## Visualization
##
plotMPsd <- function(f, g=NULL, 
                     lty=rep(1, 1+length(g)), # line type
                     col=rep(1, 1+length(g)), # color
                     type=rep("l", 1+length(g)), # line or points (e.g. for pdgrm)?
                     pch=rep("0", 1+length(g)), # point type for points
                     log=F, ylim.compound=T, 
                     mar=rep(2,4), 
                     mfrow=c(dim(f)[1],dim(f)[1]),
                     lambda_scaling=pi,
                     ylab_prefix="f_",
                     xlab="lambda",
                     excludeBoundary=T,
                     ...) {
  # f: main psd to plot
  # g: list of additional psds (CI's, grount truth, etc)
  # log: plot diagonals on log scale
  # ylim.compound: ylim according to c(f,g) instead of only g
  stopifnot(length(lty) == 1+length(g))
  stopifnot(length(col) == 1+length(g))
  stopifnot(length(type) == 1+length(g))
  N <- dim(f)[3]
  if (excludeBoundary) {
    lambda <- (1:(N-2))/(N-1)*lambda_scaling
    lambda_ind <- 2:(N-1)
  } else {
    lambda <- (0:(N-1))/(N-1)*lambda_scaling
    lambda_ind <- 1:N
  }

  if (any(is.complex(f))) {
    f_plot <- realValuedPsd(f)
  } else {
    f_plot <- f
  }
  if (any(sapply(g, is.complex))) {
    stopifnot(all(sapply(g, is.complex)))
    g_plot <- lapply(g, realValuedPsd)
  } else {
    g_plot <- g
  }
  d <- dim(f)[1]
  par(mfrow=mfrow, mar=mar)
  for (i in 1:d) {
    for (j in 1:d) {
      ylab <- paste0(ylab_prefix, i, j)
      if (i==j && log) {
        if (log) {
          ylab <- paste0("log(", ylab, ")")
        }
        if (ylim.compound) {
          ylim_min <- min(min(log(f_plot[i,j,lambda_ind])), as.numeric(sapply(g_plot, function(z) { if (is.null(z)) return(Inf) else return(min(log(z[i,j,lambda_ind])))})))
          ylim_max <- max(max(log(f_plot[i,j,lambda_ind])), as.numeric(sapply(g_plot, function(z) { if (is.null(z)) return(-Inf) else return(max(log(z[i,j,lambda_ind])))})))
          ylim=c(ylim_min, ylim_max)
        } else {
          ylim <- NULL
        }
        plot(x=lambda, 
             y=log(f_plot[i,j,lambda_ind]), 
             col=col[1], 
             lty=lty[1], 
             ylim=ylim, 
             type=type[1], 
             pch=pch[1], 
             xlab=xlab, 
             ylab=ylab, 
             ...)
        for (ll in seq_len(length(g))) {
          lines(x=lambda, y=log((g_plot[[ll]])[i,j,lambda_ind]), col=col[1+ll], lty=lty[1+ll], type=type[1+ll], pch=pch[1+ll], ...)
        }
      } else {
        if (i<j) {
          ylab <- paste0("Re(", ylab, ")")
        } else {
          if (i==j) {
            ylab <- ylab
          } else {
            ylab <- paste0("Im(", ylab, ")")
          }
        }
        if (ylim.compound) {
          ylim_min <- min(min((f_plot[i,j,lambda_ind])), as.numeric(sapply(g_plot, function(z) { if (is.null(z)) return(Inf) else return(min((z[i,j,lambda_ind])))})))
          ylim_max <- max(max((f_plot[i,j,lambda_ind])), as.numeric(sapply(g_plot, function(z) { if (is.null(z)) return(-Inf) else return(max((z[i,j,lambda_ind])))})))
          ylim=c(ylim_min, ylim_max)
        } else {
          ylim <- NULL
        }
        plot(x=lambda, 
             y=f_plot[i,j,lambda_ind], 
             col=col[1], 
             lty=lty[1], 
             ylim=ylim, 
             type=type[1], 
             pch=pch[1], 
             xlab=xlab,
             ylab=ylab, 
             ...)
        for (ll in seq_len(length(g))) {
          lines(x=lambda, y=(g_plot[[ll]])[i,j,lambda_ind], col=col[1+ll], lty=lty[1+ll], type=type[1+ll], pch=pch[1+ll], ...)
        }
      }
    }
  }
  par(mfcol=c(1,1))
}

##
## Generate real time series with a given PSD (see theorem 2 in Dai/Guo04)
##
tsFromSpectrum <- function(f) {
  d <- dim(f)[1]
  N <- dim(f)[3]
  n <- 2*N - 2
  FZ <- matrix(data=NA, nrow=N, ncol=2)
  FZ[1,] <- FZ[N,] <- rep(0,d)
  f_half <- chol_cube(f, T)
  for (j in 2:(N-1)) {
    innov <- (rnorm(d) + 1i*rnorm(d)) / sqrt(2)
    FZ[j,] <- sqrt(2*pi) * f_half[,,j] %*% innov
  }
  Z <- Re(midft(FZ))
  return(Z)
}

is_hpd <- function(A, tol=1e-15) {
  (A==Adj(A))[1] && (!hasEigenValueSmallerZero(A, tol))
}
is_spd <- function(A, tol=1e-5) {
  (A==t(A))[1] && (!hasEigenValueSmallerZero(A, tol))
}

##
## Draw from VARMA model
## Use "..." to parse arguments to rand.gen
## Canonical choice for rand.gen: "MASS::mvrnorm"
##
varma_sim <- function(n, d, 
                      ar=matrix(nrow=d,ncol=0), 
                      ma=matrix(nrow=d,ncol=0), 
                      rand.gen,
                      burnin=2e4,
                      ...)
{ 
  stopifnot(nrow(ar) == d && !(ncol(ar) %% d))
  stopifnot(nrow(ma) == d && !(ncol(ma) %% d))
  p <- ncol(ar) / d
  q <- ncol(ma) / d
  stopifnot(burnin >= max(p,q))
  X_sim <- epsilon_sim <- matrix(NA, nrow=n+burnin, ncol=d)
  if (max(p,q) > 0) {
    X_sim[1:max(p,q),] <- rand.gen(n=max(p,q), ...)
    epsilon_sim[1:max(p,q),] <- rand.gen(n=max(p,q), ...)
  }
  for (t in (max(p,q)+1):nrow(X_sim)) {
    epsilon_sim[t,] <- rand.gen(n=1, ...)
    X_sim[t,] <- epsilon_sim[t,]
    for (j in seq_len(p)) {
      X_sim[t,] <- X_sim[t,] + t(ar[,((j-1)*d+1):(j*d)] %*% matrix(X_sim[t-j,], nrow=d, ncol=1))
    }
    for (j in seq_len(q)) {
      X_sim[t,] <- X_sim[t,] + t(ma[,((j-1)*d+1):(j*d)] %*% matrix(epsilon_sim[t-j,], nrow=d, ncol=1))
    }
  }
  if (burnin > 0) {
    # # Un-comment if innovations are needed, too
    #return(list(X=X_sim[-(1:burnin),], epsilon=epsilon_sim[-(1:burnin),]))
    return(X_sim[-(1:burnin),])
  } else {
    # # Un-comment if innovations are needed
    #return(list(X=X_sim, epsilon=epsilon_sim))
    return(X_sim)
  }
}

##
## Some postprocessing stuff
##
uci_matrix <- function(fpsd.sample, alpha, uniform_among_components=F) {
  d <- dim(fpsd.sample)[1]
  N <- dim(fpsd.sample)[3]
  fpsd.uci05 <- fpsd.uci95 <- array(NA, dim=c(d, d, N))
  if (uniform_among_components) {
    # Use the same C_\alpha^* value for all components
    for (i in 1:d) {
      fpsd.sample[i,i,,] <- log(fpsd.sample[i,i,,])
    }
    fpsd.s <- apply(fpsd.sample, c(1,2,3), median)
    fpsd.mad <- apply(fpsd.sample, c(1,2,3), mad)
    fpsd.help <- uniformmax_multi(fpsd.sample)
    Cvalue <- quantile(fpsd.help, 1-alpha)
    fpsd.uci05 <- fpsd.s - Cvalue * fpsd.mad
    fpsd.uci95 <- fpsd.s + Cvalue * fpsd.mad
    for (i in 1:d) {
      fpsd.uci05[i,i,] <- exp(fpsd.uci05[i,i,])
      fpsd.uci95[i,i,] <- exp(fpsd.uci95[i,i,])
    }
  } else {
    # Use individual C_\alpha^* among each component
    for (i in 1:d) {
      for (j in 1:d) {
        uci_tmp <- uci_help(fpsd.sample[i,j,,], alpha, log=(i==j))
        fpsd.uci05[i,j,] <- uci_tmp$conflower
        fpsd.uci95[i,j,] <- uci_tmp$confupper
      }
    }
  }
  return(list(fpsd.uci05=fpsd.uci05, 
              fpsd.uci95=fpsd.uci95))
}
uci_help <- function(fpsd.sample, alpha, log=F) {
  if (log) {
    fpsd.sample <- log(fpsd.sample) #logfuller(fpsd.sample)
  }
  fpsd.s <- apply(fpsd.sample, 1, median)
  fpsd.mad <- apply(fpsd.sample, 1, mad)
  fpsd.help <- apply(fpsd.sample, 1, uniformmax)
  Cvalue <- quantile(fpsd.help, 1-alpha)
  conflower <- fpsd.s - Cvalue * fpsd.mad
  confupper <- fpsd.s + Cvalue * fpsd.mad
  if (log) {
    conflower <- exp(conflower)
    confupper <- exp(confupper)
  }
  return(list(conflower=conflower,
              confupper=confupper))
}
uniformmax <- function(sample) {
  max(abs(sample - median(sample)) / mad(sample), na.rm=T)
}
uniformmax_multi <- function(mSample) {
  d <- dim(mSample)[1]
  N <- dim(mSample)[3]
  N_sample <- dim(mSample)[4]
  C_help <- array(NA, dim=c(d,d,N,N_sample))
  for (j in 1:N) {
    for (r in 1:d) {
      for (s in 1:d) {
        C_help[r,s,j,] <- uniformmax_help(mSample[r,s,j,])
      }
    }
  }
  apply(C_help, 4, max, na.rm=T)
}
uniformmax_help <- function(sample) {
  abs(sample - median(sample)) / mad(sample)
}

##
## Some further util stuff that might be useful
##
acfMatrixFromPsd <- function(fpsd, h, excludeBoundary=T) { # TODO: better approximation for small n?
  stopifnot(h >= 0)
  N <- dim(fpsd)[3]
  n <- 2 * (N-1)
  d <- dim(fpsd)[1]
  res <- matrix(data=0, ncol=d, nrow=d)
  lambda <- pi*omegaFreq(n)
  stopifnot(length(lambda) == N)
  for (j in (1+excludeBoundary):(N-excludeBoundary)) {
    res <- res + fpsd[,,j] * exp(1i * h * lambda[j])
    if (j > 1 && j < N) {
      res <- res + Conj(fpsd[,,j]) * exp(-1i * h * lambda[j])
    }
  }
  res <- res / N * pi
}
squaredCoherency <- function(f) {
  d <- dim(f)[1]
  N <- dim(f)[3]
  k <- array(NA, dim=c(d,d,N))
  for (i in 1:d) {
    for (j in 1:d) {
      k[i,j,] <- abs(f[i,j,])^2 / f[i,i,] / f[j,j,]
    }
  }
  return(k)
}
coSpectrum <- function(f) {
  return(Re(f))
}
quadratureSpectrum <- function(f) {
  return(Im(f))
}
amplitudeSpectrum <- function(f) {
  return(abs(f))
}
phaseSpectrum <- function(f) {
  return(Arg(f))
}
