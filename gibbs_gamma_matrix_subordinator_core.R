# library(Rcpp)
# source("AGammaProcess_util.R")
# sourceCpp("gibbs_gamma_matrix_subordinator_core.cpp")

# Note: no adjoint but just transpose here
# TODO: Cpp
trans_cube_no_adj <- function(f) {
  d <- dim(f)[1]
  N <- dim(f)[3]
  res <- array(NA, dim=c(d,d,N))
  for (i in 1:N) {
    res[,,i] <- t(f[,,i])
  }
  res
}

vnp_forecast <- function(noise, f, k) {
  N <- dim(f)[3]
  f[,,1] <- Re(f[,,1]) # note: ignore imaginary part of f(0)
  if (!(n%%2)) {
    f[,,N] <- Re(f[,,N]) # note: ignore imaginary part of f(0)
  }
  d <- ncol(noise)
  stopifnot(d==dim(f)[1] && d==dim(f)[2])
  n <- nrow(noise)+k
  ff <- inv_cube(f, F)

  # Construct ACV vector under f^{-1} psd
  if (!(n%%2))
    ff_rev <- rev_cube(ff[,,-c(1,N)])
  else
    ff_rev <- rev_cube(ff[,,-1])
  ff_r <- trans_cube_no_adj(ff_rev)
  ff_for_fft <- c_cube(ff, ff_r)
  stopifnot(dim(ff_for_fft)[3]==n)
  # TODO Solve this more efficient, without loops
  gamma <- array(NA, dim=c(d,d,n))
  for (i in 1:d) {
    for (j in 1:d) {
      gamma[i,j,] <- Re(fft(ff_for_fft[i,j,], inverse=T)) / 2 / pi / n
    }
  }
  fc <- matrix(NA, nrow=k, ncol=d)
  gamma_0_inv <- solve(gamma[,,1])
  for (i in 1:k) {
    # TODO Solve this more efficient, without loops
    mu_forecast_tmp <- rep(0,d)
    for (j in 1:nrow(noise)) {
      mu_forecast_tmp <- mu_forecast_tmp +
        noise[1,] %*% gamma[,,nrow(noise)+2-j]
    }
    mu_forecast <- -gamma_0_inv %*% t(mu_forecast_tmp)
    Sigma_forecast <- gamma_0_inv
    fc[i,] <- MASS::mvrnorm(1, mu=mu_forecast, Sigma=Sigma_forecast)
    noise <- rbind(noise, fc[i,])
  }
  stopifnot(nrow(noise)==n)
  fc
}

##
## Complementary incomplete gamma function (approx.)
##
e_1 <- function(x, d=1e-12) {
  return(2/d * (1-pchisq(2*x, d)))
}
e_ab <- function(x, a, b, d=1e-12) {
  return(a*e_1(b*x, d))
}
e_1_inv <- function(y, d=1e-12) {
  return(qchisq(1-(d/2)*y, d) / 2)
}
e_ab_inv <- function(y, a, b, d=1e-12) {
  return(e_1_inv(y/a, d) / b)
}

##
## MULTIVARIATE WHITTLE LIKELIHOOD (WHITTLE AND CORRECTED PARAMETRIC)
##
llike_matrixGamma <- function(omega,
                              FZ,
                              r,
                              U,
                              Z,
                              k,
                              db.list,
                              corrected,
                              phi,
                              sigma_ar,
                              prior.q,
                              prior.cholesky,
                              excludeBoundary,
                              verbose) {
  f <- get_f_matrix(U, r, Z, k, db.list, prior.cholesky)
  # Stay numerically stable
  if (numericalUnstable(f, excludeBoundary=F)) {
    if (verbose) print_warn("Discarding f in likelihood, because of numerical instablity")
    return(-Inf)
  }
  if (corrected) {
    # Corrected parametric likelihood
    if (prior.q) {
      ll <- llike_var_corrected_q(FZ=FZ,
                                  ar=phi$ar,
                                  f_param_half=phi$f_param_half,
                                  f_param_half_trans=phi$f_param_half_trans,
                                  sigma=sigma_ar,
                                  q=f,
                                  excludeBoundary=excludeBoundary)
    } else {
      ll <- llike_var_corrected(FZ=FZ,
                                ar=phi$ar,
                                f_param_half=phi$f_param_half,
                                sigma=sigma_ar,
                                f=f,
                                excludeBoundary=excludeBoundary)
    }
  } else {
    # Whittle's likelihood
    ll <- llike_whittle_sum(FZ=FZ,
                            f=f)
  }
  return(ll)
}
get_f_matrix <- function(U, r, Z, k, db.list, prior.cholesky) {
  W <- cubeTimesVector(U, r) # get ordinate matrices W from polar decomposition (r,U)
  if (prior.cholesky) {
    # smooth Cholesky components of Gamma Measure with d*d Bernstein polynomials
    f <- fFromCholeskySmoothing(W=W, Z=Z, k=k, db.list=db.list)
  } else {
    # smooth the whole Matrix Gamma Measure with a single Bernstein polynomial
    w <- get_w_rcpp(W, Z, k)   # accumulate weight matrices w for the k mixtures from (Z,W)-tuples [abscissa Z, ordinate W]
    w[,,1] <- Re(w[,,1]) # ensure that f(0) is spd (not complex)
    w[,,k] <- Re(w[,,k]) # ensure that f(pi) is spd (not complex)
    f <- get_mix_rcpp(w, db.list[[k]]) # get Bernstein-mixture with weight matrices w
  }
  return(f)
}

fFromCholeskySmoothing <- function(W, Z, k, db.list) {
  d <- dim(W)[1]
  stopifnot(length(k)==d*d)
  L <- array(0, dim=c(d,d,length(db.list[[1]][1,]))) # Cholesky composition of f: to be constructed with k's and Gamma process
  jj <- 1
  for (i in 1:d) {
    w_ii_half <- chol_cube_semidefinite(get_w_rcpp(W, Z, k[jj]))[i,i,]
    L[i,i,] <- get_mix_rcpp_1d(w_ii_half, db.list[[k[jj]]])
    jj <- jj+1
    for (j in seq_len(i-1)) {
      w_re_half <- Re(chol_cube_semidefinite(get_w_rcpp(W, Z, k[jj]))[i,j,])
      w_im_half <- Im(chol_cube_semidefinite(get_w_rcpp(W, Z, k[jj+1]))[i,j,])
      L[i,j,] <- get_mix_rcpp_1d(w_re_half, db.list[[k[jj]]]) +
        1i * get_mix_rcpp_1d(w_im_half, db.list[[k[jj+1]]], excludeBoundary=T)
      jj <- jj+2
    }
  }
  f <- mult_cube(L, trans_cube(L))
  return(f)
  }


##
## GAMMA PROCESS PRIOR: BASED ON THE TRACE NORM TR((X^*X)^{1/2})=TR(X) FOR HPD X
##
lprior_matrixGamma <- function(r,
                               U,
                               Z,
                               k,
                               C_alpha,
                               omega_fun,
                               k.theta,
                               eta,
                               Sigma_fun,
                               phi,
                               verbose) {

  # # remove boundary stuff
  # L <- length(r)
  # r <- r[-c(1,L)]
  # U <- U[,,-c(1,L)]
  # Z <- Z[-c(1,L)]
  # #
  if (min(Z) <= 0 || max(Z) >= 1) {
    print(Z)
    stop()
  }
  omega_vals <- do.call(omega_fun, list(Z)) #eval(omega_fun(Z))
  Sigma_vals <- do.call(Sigma_fun, list(Z)) #eval(Sigma_fun(Z))
  #tryCatch({
  beta_vals <- beta_fun_AGamma_process_cube(U, Sigma_vals)
  # }, error=function(e) {
  #   #print(U)
  #   #print(Sigma_vals)
  #   stop("Error in beta_fun")
  # })
  # TODO Employ log determinant of this transformation
  W <- e_ab(sort(r, decreasing=T), C_alpha, beta_vals)
  if (is.unsorted(W)) { # that might happen because of numerical integral approximation
    if (verbose) print_warn("W unsorted -- sorting it")
    W <- sort(W)
  }
  V <- c(W[1], diff(W)) # V ~iid~ exp(1)
  #tryCatch({
  lp <- sum(dexp(V, 1, log=T)) +
    lalphaStar_AGamma_process(U, eta, omega_vals, Sigma_vals) -
    sum(k.theta * k * log(k)) +
    lprior_parametricPart(phi) + # Only needed for corrected likelihood (==0 otherwise)
    logdet_radialJacobian(C_alpha, beta_vals, r) # Prior is for v, sampling is for r -- need logDeterminant of transformation Jacobian
  # }, error=function(e){
  #   #print(Z)
  #   #print(Sigma_vals)
  #   stop("Error in alpha_fun")
  # })

  return(lp)
}
##
## Unnormalized prior for the parametric (VAR) part of the semiparametric procedure
## (Normal-Inverse-Wishart for beta)
## Note: The covariance matrix is NOT modeled here, but within the Gamma process!
##
## phi$beta: beta vector in Normal-Inverse-Wishart representation ('rolled out' version of phi$ar)
## phi$mu_beta: prior mean for beta vector
## phi$V_beta: prior covariance matrix for beta vector
##
lprior_parametricPart <- function(phi) {
  if (is.null(phi$mu_beta)) {
    # dummy fallback, for non-toggle
    return(0)
  } else {
    # asume V_beta_half and mu_beta are fixed
    return(-.5* t(phi$beta-phi$mu_beta) %*% phi$V_beta_inv %*% (phi$beta-phi$mu_beta))
  }
}
# log determinant of Jacobian of mapping (r_1,...,r_L) -> (v_1,...,v_L)
# unnormalized -- see (62) in draft_20171127
logdet_radialJacobian <- function(C_alpha, beta_vals, r_vals) {
  sum(-beta_vals * r_vals -log(r_vals))
}


##
## Log Posterior, unnormalized
##
lpost_matrixGamma <- function(omega,
                              FZ,
                              r,
                              U,
                              Z,
                              k,
                              C_alpha,
                              omega_fun,
                              k.theta,
                              pdgrm,
                              db.list,
                              eta,
                              Sigma_fun, # prior parameter for AGamma process
                              corrected,
                              phi,
                              sigma_ar, # corresponding to AR fit
                              prior.q,
                              prior.cholesky,
                              excludeBoundary,
                              verbose) {

    ll <- llike_matrixGamma(omega=omega,
                                 FZ=FZ,
                                 r=r,
                                 U=U,
                                 Z=Z,
                                 k=k,
                                 db.list=db.list,
                                 corrected=corrected,
                                 phi=phi,
                                 sigma_ar=sigma_ar,
                                 prior.q=prior.q,
                                 prior.cholesky=prior.cholesky,
                                 excludeBoundary=excludeBoundary,
                                 verbose=verbose)

  lp <- lprior_matrixGamma(r=r,
                           U=U,
                           Z=Z,
                           k=k,
                           C_alpha=C_alpha,
                           omega_fun=omega_fun,
                           k.theta=k.theta,
                           eta=eta,
                           Sigma_fun=Sigma_fun,
                           phi=phi,
                           verbose=verbose)
  ##
  ## BEGIN DEBUG
  ##
  if (is.na(ll)) {
    likeDump <- list(r=r,
                     U=U,
                     Z=Z,
                     k=k,
                     corrected=corrected,
                     phi=phi,
                     sigma_ar=sigma_ar,
                     prior.q=prior.q,
                     prior.cholesky=prior.cholesky)
    save(list=c("likeDump"), file="likeDump")
    print_warn("Likelihood NA computed, dumped output")
  }
  if (is.na(lp)) {
    priorDump <- list(r=r,
                      U=U,
                      Z=Z,
                      k=k,
                      C_alpha=C_alpha,
                      k.theta=k.theta,
                      eta=eta,
                      phi=phi)
    save(list=c("priorDump"), file="priorDump")
    print_warn("Prior NA computed, dumped output")
  }
  ##
  ## END DEBUG
  ##
  return(ll + lp)
}

reduceMemoryStorage_matrixGamma <- function(mcmc) {
  ## Delete memory-intensive traces
  ret <- (list(data=mcmc$data,
               fpsd.s=mcmc$fpsd.s,
               fpsd.mean=mcmc$fpsd.mean,
               fpsd.s05=mcmc$fpsd.s05,
               fpsd.s95=mcmc$fpsd.s95,
               # fpsd.uci05=mcmc$fpsd.uci05,
               # fpsd.uci95=mcmc$fpsd.uci95,
               fpsd.uuci05=mcmc$fpsd.uuci05,
               fpsd.uuci95=mcmc$fpsd.uuci95,
               lpostTrace=mcmc$lpostTrace,
               theta=mcmc$theta,
               k=mcmc$k))
  return(ret)
}

reduceMemoryStorage_vnp <- function(mcmc, data) {
  ## Delete memory-intensive traces
  ret <- (list(data=data,
               fpsd.median=mcmc$psd.median,
               fpsd.mean=mcmc$psd.mean,
               fpsd.s05=mcmc$psd.p05,
               fpsd.s95=mcmc$psd.p95,
               # fpsd.uci05=mcmc$fpsd.uci05,
               # fpsd.uci95=mcmc$fpsd.uci95,
               fpsd.uuci05=mcmc$psd.u05,
               fpsd.uuci95=mcmc$psd.u95,
               lpostTrace=mcmc$lpost,
               r=mcmc$r,
               x=mcmc$x,
               U=mcmc$U,
               k=mcmc$k))
  return(ret)
}

lalphaStar_AGamma_process <- function(U, # Cube in \mathbb S_d^+, U(x)'s
                                      eta, # constant > d-1
                                      omega_vals, # vector of positive reals, omega(x)'s
                                      Sigma_vals # cube in \mathcal S_d^+, Sigma(x)'s
) {
  # Log density of A-Gamma(eta,omega,Sigma) process
  # measure alpha_star on [0,1] \times \mathbb S_d^+
  # See Lemma 12 and (35) in draft_20171020
  #tryCatch({
  res <- sum(
    log(omega_vals) +
      lalphaStar_AGamma(U, eta, Sigma_vals)
  )
  return(res)
  # }, error=function(e) {
  #   # print(U)
  #   # print(eta)
  #   # print(omega_vals)
  #   # print(Sigma_vals)
  #   # print("----------")
  #   stop(e)
  # })
}
