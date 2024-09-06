##
## Whittle/Corrected Likelihood /w AGamma(eta_0,omega,Sigma) prior
## eta_0 fixed, omega strict weight function, Sigma matrix function
##

#
normalizedBernsteins <- F # normalize Bernstein polynomials in maximum norm (may or may not imrpove mixing)
#

# source("../gibbs_util.R")
# source("gibbs_multivariate_util.R")
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # lambda function using C++11
# library(Rcpp)
# source("mittelbach_util.R")
# source("gibbs_gamma_matrix_subordinator_core.R")

##
## New MCMC algorithm
##
gibbs_m_nuisance <- function(data,
                             mcmc_params,
                             seg_n=1,
                             truncation=FALSE,
                             trunc_N=NULL,
                             corrected,
                             prior_params,
                             model_params) {
  
  warning("Using conjugate sampling for theta. TODO: Get that consistent with MH!!")
  
  # --- STUFF TODO ----------------------------
  # -------------------------------------------
  # + implement pure parametric VAR nuisance model
  # + implement corrected VAR nuisance model (+ double steps?)
  # + include model_params, similar to gibbs_NP(C)_nuisance.R
  # + fix lprior function for inhomogeneous hpd AGamma [C_alpha, omega_fun, Sigma_fun]
  # + add print_interval
  # + adjust gibbs_gamma_matrix_subordinator_test.R to new interface
  # + add numerical_thresh
  # + add omega_fun, Sigma_fun
  # + remove M, g0.alpha, g0.beta 
  #   M replaced by C_alpha, {g0.alpha, g0.beta} by omega_fun
  # + remove tau
  # + add bernstein_{l,r}
  # + employ Jacobians for rho
  # -------------------------------------------
  # - implement alpha.toggle? (damping function stuff)
  # -------------------------------------------
  
  
  
  # Note: "Assuming symmetric theta proposals in MH steps"
  
  time_b <- Sys.time()
  
  
  stopifnot(ncol(data) > 1) # use 1d version instead
  d <- ncol(data)
  
  #if ((n%%2) != 0) {
  #  TXT <- "TODO: Adjust llike_whittle in gibbs_multivariate_util.R for odd lengths. 
  #  There may also be more adjustments needed in the code framework." 
  #  stop(TXT)
  #}
  
  # Segments
  stopifnot(!is.null(seg_n)); stopifnot(seg_n>0) 
  seg_n <- seg_n
  n <- nrow(data) / seg_n
  n_init <- n ## Store the original full length of each segment
  data_new <- array(rep(NA, n*d*seg_n), dim = c(n, d, seg_n))
  for (ii in 1:seg_n) {
    data_new[,,ii] <- data[(n*(ii-1)+1):(n*ii),]
  }
  
  # Trunction of the Fourier transforms
  if (truncation) {
    stopifnot(!is.null(trunc_N))
    if (trunc_N <= 0 || trunc_N > (n+2)/2) {
      stop("The number of trunction should be an integer that is positive and smaller than the length of the Fourier transform of each segment")
    }
    truncation <- truncation
    trunc_N <- trunc_N
    n <- trunc_N*2 + 2 ##? need to check it again
  }
  
  # MCMC parameters
  stopifnot(!is.null(mcmc_params$Ntotal)); stopifnot(mcmc_params$Ntotal>0)
  Ntotal <- mcmc_params$Ntotal
  stopifnot(!is.null(mcmc_params$burnin)); stopifnot(mcmc_params$burnin>=0 && mcmc_params$burnin<Ntotal)
  burnin <- mcmc_params$burnin
  stopifnot(!is.null(mcmc_params$thin)); stopifnot(mcmc_params$thin>=1)
  thin <- mcmc_params$thin
  stopifnot(!is.null(mcmc_params$print_interval)); stopifnot(mcmc_params$print_interval>0)
  print_interval <- mcmc_params$print_interval
  stopifnot(!is.null(mcmc_params$numerical_thresh)); stopifnot(mcmc_params$numerical_thresh>0)
  NUMERICAL_THRESH <- mcmc_params$numerical_thresh # >= 1e-12 recommended
  stopifnot(!is.null(mcmc_params$verbose))
  verbose <- mcmc_params$verbose
  
  # Adaptive MCMC parameters
  stopifnot(!is.null(mcmc_params$Nadaptive)); stopifnot(mcmc_params$Nadaptive>0)
  Nadaptive <- mcmc_params$Nadaptive
  stopifnot(!is.null(mcmc_params$adaption.batchSize)); stopifnot(mcmc_params$adaption.batchSize>0)
  adaption.batchSize <- mcmc_params$adaption.batchSize
  stopifnot(!is.null(mcmc_params$adaption.targetAcceptanceRate)); stopifnot(mcmc_params$adaption.targetAcceptanceRate>0 && mcmc_params$adaption.targetAcceptanceRate<1)
  adaption.targetAcceptanceRate <- mcmc_params$adaption.targetAcceptanceRate
  
  # AGAMMA PRIOR PAREMETERS
  stopifnot(!is.null(prior_params$prior.cholesky))
  prior.cholesky <- prior_params$prior.cholesky
  if (prior.cholesky) {
    stop('Development of prior.cholesky has been stalled, 
         because prior knowledge is difficult to interpret -- 
         see the prior visualization in 
         inverse_matrix_levy_measure_algorithm.R
         with prior.cholesky=T')
  }
  stopifnot(!is.null(prior_params$eta)); stopifnot(prior_params$eta > d-1) # eta fixed
  eta <- prior_params$eta
  stopifnot(!is.null(prior_params$omega_fun)); stopifnot(class(prior_params$omega_fun)=="function")
  omega_test_val <- prior_params$omega_fun(0.5); stopifnot(omega_test_val > 0); rm(omega_test_val)
  omega_fun <- prior_params$omega_fun
  stopifnot(!is.null(prior_params$Sigma_fun)); stopifnot(class(prior_params$Sigma_fun)=="function")
  Sigma_test_val <- prior_params$Sigma_fun(0.5); stopifnot(ncol(Sigma_test_val)==d && nrow(Sigma_test_val)==d); stopifnot(is_hpd(Sigma_test_val)); rm(Sigma_test_val)
  Sigma_fun <- prior_params$Sigma_fun
  stopifnot(!is.null(prior_params$k.theta)); stopifnot(all(prior_params$k.theta > 0)); stopifnot((prior.cholesky && length(prior_params$k.theta)==d*d) || ((!prior.cholesky) && length(prior_params$k.theta)==1)) # may be a vector (iff prior.cholesky)
  k.theta <- prior_params$k.theta
  stopifnot(!is.null(prior_params$kmax)); stopifnot(prior_params$kmax > 0)
  kmax <- prior_params$kmax
  stopifnot(!is.null(prior_params$bernstein_l) && !is.null(prior_params$bernstein_r)); stopifnot(prior_params$bernstein_l >= 0 && prior_params$bernstein_r <= 1)
  bernstein_l <- prior_params$bernstein_l
  bernstein_r <- prior_params$bernstein_r
  stopifnot(!is.null(prior_params$coarsened))
  coarsened <- prior_params$coarsened
  stopifnot(!is.null(prior_params$L)); stopifnot(prior_params$L > 0)
  L <- prior_params$L
  if (corrected) {
    stopifnot(!is.null(prior_params$toggle))
    toggle <- prior_params$toggle
    stopifnot(!is.null(prior_params$prior.q))
    prior.q <- prior_params$prior.q
    stopifnot(!is.null(prior_params$var.order)); stopifnot(prior_params$var.order>0)
    var.order <- prior_params$var.order
    # not implmeneted yet
    # stopifnot(!is.null(prior_params$alpha.toggle)) 
    # alpha.toggle <- prior_params$alpha.toggle
    if (toggle) {
      phi.fit <- NULL
      sigma.fit <- NULL
      stopifnot(!is.null(prior_params$mu_beta)) # prior parameters for parametric part (Normal Inverse Wishart)
      mu_beta <- prior_params$mu_beta
      stopifnot(!is.null(prior_params$V_beta))  # prior parameters for parametric part (Normal Inverse Wishart)
      V_beta <- prior_params$V_beta
    } else {
      stopifnot(!is.null(prior_params$phi.fit) && !is.null(prior_params$sigma.fit))
      phi.fit <- prior_params$phi.fit
      sigma.fit <- prior_params$sigma.fit
      mu_beta <- NULL
      V_beta <- NULL
    }
  } else {
    toggle <- alpha.toggle <- prior.q <- FALSE
    phi.fit <- sigma.fit <- NULL
  }
  # Integral of Levy measure part on [0,1] x \mathbb S_d^+
  C_alpha <- integrate(omega_fun, lower=0, upper=1)$value
  
  # Model paramaters
  stopifnot(!is.null(model_params$theta_dim)); stopifnot(model_params$theta_dim >= 0)
  theta_dim <- model_params$theta_dim
  stopifnot(!is.null(model_params$get_noise)); stopifnot(class(model_params$get_noise)=="function")
  get_noise <- model_params$get_noise
  stopifnot(!is.null(model_params$get_data)); stopifnot(class(model_params$get_data)=="function")
  get_data <- model_params$get_data
  stopifnot(!is.null(model_params$initialize_theta)); stopifnot(class(model_params$initialize_theta)=="function")
  initialize_theta <- model_params$initialize_theta
  stopifnot(!is.null(model_params$lprior_theta)); stopifnot(class(model_params$lprior_theta)=="function")
  lprior_theta <- model_params$lprior_theta
  stopifnot(!is.null(model_params$propose_next_theta)); stopifnot(class(model_params$propose_next_theta)=="function")
  propose_next_theta <- model_params$propose_next_theta
  # stopifnot(!is.null(model_params$excludeBoundary))
  # excludeBoundary <- model_params$excludeBoundary # (!is.null(model_params$excludeBoundary) && model_params$excludeBoundary)
  if (n %% 2) {
    boundaryFrequecies <- c(1,n)
  } else {
    boundaryFrequecies <- 1
  }
  
  omega <- omegaFreq(n)
  N <- length(omega)
  lambda <- pi * omega
  pdgrm_scaling <- c(pi, rep(2*pi, n-1))
  if (!(n%%2)) pdgrm_scaling[n] <- pi
  
  # Pre-computed beta mixtures up to kmax as list
  db.list <- dbList(n, kmax, normalized=normalizedBernsteins,
                    bernstein_l, bernstein_r, coarsened) 
  
  # initielize storage: parameter of interest (if nuisance model)
  theta <- matrix(NA, nrow=theta_dim, ncol=Ntotal)
  
  # initialize storage: nonparametric part
  lpostTrace <- llikeTrace <- lpriorTrace <- rep(NA, Ntotal)
  U <- array(NA, dim=c(d,d,L,Ntotal)) # parametrization: W = rU
  U__phi <- array(NA, dim=c(d*d-1, L, Ntotal)) # Mittelbach et al
  U__IND_ALL <- 1:(d*d-1)
  U__IND_PIHALF <- (1:(d-1))^2
  U__IND_PI <- setdiff(U__IND_ALL, U__IND_PIHALF)
  U__SCALING <- U__IND_ALL
  U__SCALING[U__IND_PI] <- pi
  U__SCALING[U__IND_PIHALF] <- pi / 2
  r <- matrix(NA, nrow=L, ncol=Ntotal)
  Z <- matrix(NA, nrow=L, ncol=Ntotal)
  k__DIM <- d*d - ((d*d-1)*(1-prior.cholesky)) # d*d k's for prior.cholesky, 1 k for !prior.choelsky
  k <- array(NA, dim=c(k__DIM, Ntotal))
  
  # initialize storage: parametric VAR part
  if (corrected && toggle) {
    param__beta <- array(NA, dim=c(d*d*var.order, Ntotal))
    param__phi <- array(NA, dim=c(d,d*var.order,Ntotal)) # redundant, for convenience
  }
  
  # starting values
  nu_U <- rep(10, L)
  k[,1] <- rep(kmax / 2, k__DIM) # TODO This needs tuning
  #Z[,1] <- rbeta(L, 1, 1) # TODO Draw from omega_star ?
  #r[,1] <- rexp(L, 1)
  Z[,1] <- seq(1/(2*max(k[,1])), 1-1/(2*max(k[,1])), length.out=L)
  r[,1] <- rep(1/L, L)
  U_start <- cholesky_runif(L, d)
  U__phi[,,1] <- U_start$phi
  U[,,,1] <- U_start$U
  # for (j in 1:L) {
  #   # TODO: Can we improve this to draw uniformly U, as in Mittelbach ?
  #   U__phi[U__IND_PIHALF,j,1] <- runif(length(U__IND_PIHALF), 0, 1) * pi / 2 
  #   U__phi[U__IND_PI,j,1] <- runif(length(U__IND_PI), 0, 1) * pi 
  #   U[,,j,1] <- cholesky_UFromPhi(U__phi[,j,1]) 
  # }
  f_start <- get_f_matrix(U[,,,1], r[,1], Z[,1], k[,1], db.list, prior.cholesky)
  # Make sure to start numerically stable
  if (numericalUnstable(f_start, excludeBoundary=F)) {
    cat("Re-drawing starting value due to numerical stability..")
    unstable_start <- T
  } else {
    unstable_start <- F
  }
  while(numericalUnstable(f_start, excludeBoundary=F)) {
    cat(".")
    Z[,1] <- rbeta(L, 1, 1) # TODO Draw from omega_star ?
    r[,1] <- rexp(L, 1)
    f_start <- get_f_matrix(U[,,,1], r[,1], Z[,1], k[,1], db.list, prior.cholesky)
  }
  if (unstable_start) cat("DONE!\n")
  if (corrected) {
    if (toggle) {
      a1 <- ar(data, order.max=var.order, aic=F)
      beta_start <- NULL
      for (jj in 1:d) {
        beta_start <- c(beta_start, c(t(a1$ar[,jj,])))
      }
      param__beta[,1] <- beta_start
      param__phi[,,1] <- phiFromBeta_normalInverseWishart(param__beta[,1], d, var.order)
      sigma.fit <- a1$var.pred
      V_beta_inv <- solve(V_beta)
      f_param_half <- chol_cube(psd_varma(lambda, param__phi[,,1], sigma=sigma.fit)$psd, excludeBoundary=F)
      ##
      ## NOTE: a bit abuse of notation for phi: 
      ## - Parse f_param_half to likelihood to save redundant computations
      ## - Parse beta parametrization to the Normal-Inverse-Wishart prior
      ##
      phi.fit <- list(ar=param__phi[,,1],
                      f_param_half=f_param_half,
                      f_param_half_trans=trans_cube(f_param_half),
                      beta=param__beta[,1],  ## 
                      mu_beta=mu_beta,       ## include stuff for prior computation, iff toggle
                      V_beta_inv=V_beta_inv) ##
    } else {
      ##
      ## NOTE: a bit abuse of notation for phi: 
      ## Parse f_param_half to likelihood to save redundant computations
      ##
      f_param_half <- chol_cube(psd_varma(lambda, ar=phi.fit$ar, sigma=sigma.fit)$psd, excludeBoundary=F)
      phi.fit <- list(ar=phi.fit$ar,
                      f_param_half=f_param_half,
                      f_param_half_trans=trans_cube(f_param_half))
    }
  }
  theta[,1] <- initialize_theta(data)
  
  # proposal variances: nonparametric part
  eps_r <- eps_Z <- eps_U <- seq(1, L) / (seq(1, L) + 2 * sqrt(n))
  lsd_r <- log(eps_r) / 2
  lsd_Z <- log(eps_Z) / 2
  lsd_U <- log(eps_U) / 2
  
  # proposal variances: parametric part
  if (corrected && toggle) {
    eps_param <- rep(1/n/sqrt(n), var.order)
    lsd_param <- log(eps_param) / 2
  }
  
  ##
  ## MH-within Gibbs sampler
  ##
  for (i in 1:(Ntotal-1)) {
    
    if (!(i%%print_interval)) {
      cat("iteration ", i, "/", Ntotal, "\n", sep="")
    }
    
    #storage segmented data
    if (seg_n == 1) {
      noise <- array(rep(NA, n_init*d*1), dim = c(n_init, d, 1))
      FZ <- array(rep(NA, N*d*1), dim = c(N, d, 1))
      noise[,,1] <- get_noise(data_new, theta[,i])
      FZ[,,1] <- mdft(noise)
    } else {
      noise <- array(rep(NA, n_init*d*seg_n), dim = c(n_init, d, seg_n))
      FZ <- array(rep(NA, N*d*seg_n), dim = c(N, d, seg_n))
      for (ii in 1:seg_n) {
        noise[,,ii] <- get_noise(data_new[,,ii], theta[,i])
        FZ[,,ii] <- mdft(noise[,,ii])[1:N,]
      }
    }
    
    #noise <- get_noise(data, theta[,i])  # noise = data - signal
    #FZ <- mdft(noise)  # Frequency domain
    
    if (i==1) {
      ##
      ## f.store: previous lpost value to save some computation time in MH steps
      ## Needs to be updated for every proposal acceptance
      ##
      f.store <- lpost_matrixGamma(omega=omega,
                                   FZ=FZ,
                                   r=r[,i],
                                   U=U[,,,i],
                                   Z=Z[,i],
                                   k[,i],
                                   C_alpha=C_alpha,
                                   omega_fun=omega_fun,
                                   k.theta=k.theta,
                                   db.list=db.list,
                                   eta=eta,
                                   Sigma_fun=Sigma_fun,
                                   corrected=corrected,
                                   phi=phi.fit,
                                   sigma_ar=sigma.fit,
                                   prior.q=prior.q,
                                   prior.cholesky=prior.cholesky,
                                   excludeBoundary=T, # note
                                   verbose=verbose)
      
      lpostTrace[1] <- f.store + lprior_theta(theta[,1])
      lpriorTrace[1] <- lprior_matrixGamma(r=r[,1],
                                           U=U[,,,1],
                                           Z=Z[,1],
                                           k=k[,1],
                                           C_alpha=C_alpha,
                                           omega_fun=omega_fun,
                                           k.theta=k.theta,
                                           eta=eta,
                                           Sigma_fun=Sigma_fun,
                                           phi=phi.fit,
                                           verbose=verbose) +
        lprior_theta(theta[,1])
      llikeTrace[1] <- lpostTrace[1] - lpriorTrace[1]
    }
    
    ##
    ## Step 0: Adjust propsal variances
    ##
    #print("Step 0")
    if ((i < Nadaptive) && (i > 1) && (i %% adaption.batchSize == 1)) {
      batch <- (i-adaption.batchSize):(i-1)
      adaption.delta <- min(0.05, 1/(i^(1/2))) # c.f. Rosenthal
      ### r
      batch.r <- r[, batch]
      batch.r.acceptanceRate <- apply(batch.r, 1, acceptanceRate) 
      lsd_r <- lsd_r + ((batch.r.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      eps_r <- exp(2*lsd_r)
      # ### Z (Note: not adapted, to improve mixing)
      # print_warn("!!If Z proposals are adapted, take care to make them stay in [0,1]!!")
      # batch.Z <- Z[, batch]
      # batch.Z.acceptanceRate <- apply(batch.Z, 1, acceptanceRate)
      # lsd_Z <- lsd_Z + ((batch.Z.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      # eps_Z <- exp(2*lsd_Z)
      ### U (first component suffices)
      batch.U <- U[1,1,,batch]
      batch.U.acceptanceRate <- apply(batch.U, 1, acceptanceRate)
      lsd_U <- lsd_U + ((batch.U.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
      eps_U <- exp(2*lsd_U)
      ### parametric part (first component of phi matrices suffices)
      if (corrected && toggle) {
        batch.param <- param__phi[1,(1:var.order)*(d-1)+1,batch,drop=F]
        batch.param.acceptanceRate <- apply(batch.param, 1, acceptanceRate)
        lsd_param <- lsd_param + ((batch.param.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
        eps_param <- exp(2*lsd_param)
      }
    }
    
    ##
    ## Step 1: Sample k from MH
    ##
    #print("Step 1")
    k.old <- k[, i]
    for (j in 1:k__DIM) {
      k.star <- k.old
      k.star[j] <- round(rt(1, 1)) + k.old[j]  # Cauchy distribution discretized
      while (k.star[j] < 3 || k.star[j] > kmax) {  # A bit hacky
        k.star[j] <- round(rt(1, 1)) + k.old[j]
      }
      f.k.star <- lpost_matrixGamma(omega=omega,
                                    FZ=FZ,
                                    r=r[,i],
                                    U=U[,,,i],
                                    Z=Z[,i],
                                    k.star,
                                    C_alpha=C_alpha,
                                    omega_fun=omega_fun,
                                    k.theta=k.theta,
                                    db.list=db.list,
                                    eta=eta,
                                    Sigma_fun=Sigma_fun,
                                    corrected=corrected,
                                    phi=phi.fit,
                                    sigma_ar=sigma.fit,
                                    prior.q=prior.q,
                                    prior.cholesky=prior.cholesky,
                                    excludeBoundary=T, # note
                                    verbose=verbose)
      f.k <- f.store
      # Accept/reject
      alpha1 <- min(0, f.k.star - f.k)  # log acceptance ratio
      if (log(runif(1, 0, 1)) < alpha1) {
        k[j, i+1] <- k.star[j]  # Accept k.star
        k.old <- k.star
        f.store <- f.k.star
      } else {
        k[j, i+1] <- k[j, i]  # Reject and use previous
      }
    }
    
    ##
    ## Step 2: Sample Z's
    ##
    #print("Step 2")
    Z.old <- Z[, i] # very previous Z: update in case of acceptance!
    for (l in 1:L) {
      Z.star <- Z.old
      Z.star[l] <- Z.star[l] + runif(1, -eps_Z[l], eps_Z[l])
      Z.star[l][Z.star[l] < 0] <- Z.star[l] + 1 # shift in [0,1]
      Z.star[l][Z.star[l] > 1] <- Z.star[l] - 1 # shift in [0,1]
      f.Z.star <- lpost_matrixGamma(omega=omega,
                                    FZ=FZ,
                                    r=r[,i],
                                    U=U[,,,i],
                                    Z=Z.star,
                                    k=k[,i+1],
                                    C_alpha=C_alpha,
                                    omega_fun=omega_fun,
                                    k.theta=k.theta,
                                    db.list=db.list,
                                    eta=eta,
                                    Sigma_fun=Sigma_fun,
                                    corrected=corrected,
                                    phi=phi.fit,
                                    sigma_ar=sigma.fit,
                                    prior.q=prior.q,
                                    prior.cholesky=prior.cholesky,
                                    excludeBoundary=T, # note
                                    verbose=verbose)
      # Accept / reject
      f.Z <- f.store
      alpha2 <- min(0, f.Z.star - f.Z) # Note: symmetric proposals
      if (log(runif(1,0,1)) < alpha2) {
        Z[l, i+1] <- Z.star[l] # accept
        f.store <- f.Z.star
        Z.old <- Z.star
      } else {
        Z[l,i+1] <- Z[l, i] # reject
      }
    }
    
    ##
    ## Step 3: Sample r's on log scale
    ##
    #print("Step 3")
    r.old <- r[,i]
    for (l in 1:L) {
      r.star <- r.old
      r.star[l] <- rlnorm(1, meanlog=log(r.old[l]), sdlog=sqrt(eps_r[l]))
      if (r.star[l] < NUMERICAL_THRESH) { # keep numerically stable (incomplete gamma function approximation)
        if (verbose) print_warn(paste0("Discaring r_", l, " prosal with value ", r.star[l], " due to numerics"))
        r[l, i+1] <- r[l, i]
      } else {
        f.r.star <- lpost_matrixGamma(omega=omega,
                                      FZ=FZ,
                                      r=r.star,
                                      U=U[,,,i],
                                      Z=Z[,i+1],
                                      k=k[,i+1],
                                      C_alpha=C_alpha,
                                      omega_fun=omega_fun,
                                      k.theta=k.theta,
                                      db.list=db.list,
                                      eta=eta,
                                      Sigma_fun=Sigma_fun,
                                      corrected=corrected,
                                      phi=phi.fit,
                                      sigma_ar=sigma.fit,
                                      prior.q=prior.q,
                                      prior.cholesky=prior.cholesky,
                                      excludeBoundary=T, # note
                                      verbose=verbose)
        f.r <- f.store
        # Accept / reject
        alpha3 <- min(0, f.r.star - 
                        f.r + 
                        dlnorm(r.old[l], meanlog=log(r.star[l]), sd=sqrt(eps_r[l]), log=T) -
                        dlnorm(r.star[l], meanlog=log(r.old[l]), sd=sqrt(eps_r[l]), log=T))
        #)
        if (log(runif(1, 0, 1)) < alpha3) {
          r[l,i+1] <- r.star[l] # accept
          r.old <- r.star
          f.store <- f.r.star
        } else {
          r[l,i+1] <- r[l,i] # reject
        }
      }
    }
    
    ##
    ## Step 4: Sample U's
    ##
    #print("Step 4")
    U.old <- U[,,,i]
    for (l in 1:L) {
      
      # MH proposal with boundary inverse-reflections (similar to Z's -- see Choudhuri 2004)
      rejectedU <- F
      U__phi.star <- U__phi[,l,i] + runif(length(U__SCALING), -eps_U[l], eps_U[l]) * U__SCALING
      U__phi.star[U__phi.star < 0] <- U__phi.star[U__phi.star < 0] + U__SCALING[U__phi.star < 0] # put in interval
      U__phi.star[U__phi.star > U__SCALING] <- U__phi.star[U__phi.star > U__SCALING] - U__SCALING[U__phi.star > U__SCALING] # put in interval
      U.star <- U.old
      U.star[,,l] <- cholesky_UFromPhi(U__phi.star)
      if (hasEigenValueSmallerZero(U.star[,,l], TOL=NUMERICAL_THRESH)) { # stay positive definite
        if (verbose) print_warn(paste0("Discaring U_", l, " prosal",  #"with value ", U.star[,,l], 
                                       " because of eigenvalues"))
        U[,,l,i+1] <- U[,,l,i] # reject and use previous 
        U__phi[,l,i+1] <- U__phi[,l,i]
        rejectedU <- T
      } 
      if (matCond(U.star[,,l]) < NUMERICAL_THRESH) { # stay numerically stable
        if (verbose) print_warn(paste0("Discaring U_", l, " prosal ", # "with value ", U.star[,,l], 
                                       " because of numerics"))
        U[,,l,i+1] <- U[,,l,i] # reject and use previous 
        U__phi[,l,i+1] <- U__phi[,l,i]
        rejectedU <- T
      }
      if (!rejectedU) {
        f.U.star <- lpost_matrixGamma(omega=omega,
                                      FZ=FZ,
                                      r=r[,i+1],
                                      U=U.star,
                                      Z=Z[,i+1],
                                      k=k[,i+1],
                                      C_alpha=C_alpha,
                                      omega_fun=omega_fun,
                                      k.theta=k.theta,
                                      db.list=db.list,
                                      eta=eta,
                                      Sigma_fun=Sigma_fun,
                                      corrected=corrected,
                                      phi=phi.fit,
                                      sigma_ar=sigma.fit,
                                      prior.q=prior.q,
                                      prior.cholesky=prior.cholesky,
                                      excludeBoundary=T, # note
                                      verbose=verbose)
        # accept/reject
        f.U <- f.store
        # Note: symmetric Uniform proposals -- but need jacobian of transformation [phi_U |-> U]
        alpha4 <- min(0, f.U.star + 
                        cholesky_jacobianLogDeterminant(U__phi.star) - 
                        f.U -
                        cholesky_jacobianLogDeterminant(U__phi[,l,i])) 
        if (log(runif(1,0,1)) < alpha4) {
          U[,,l,i+1] <- U.star[,,l] # accept
          U__phi[,l,i+1] <- U__phi.star
          U.old <- U.star
          f.store <- f.U.star
        } else {
          U[,,l,i+1] <- U[,,l,i] # reject
          U__phi[,l,i+1] <- U__phi[,l,i]
        }
      }
    }
    
    ##
    ## Step 5: Sample parametric part from full conjugate conditional
    ##
    #print("Step 5")
    if (corrected && toggle) {
      # TODO This proposal needs improval
      param__beta.old <- param__beta[,i]
      for (jj in 1:var.order) {
        indices_jj <- ((jj-1)*d*d+1):(jj*d*d)
        param__beta.star <- param__beta.old
        param__beta.star[indices_jj] <- param__beta.star[indices_jj] + 
          MASS::mvrnorm(1, mu=rep(0, d*d), Sigma=diag(eps_param[jj], d*d, d*d))
        param__phi.star <- phiFromBeta_normalInverseWishart(param__beta.star, d, var.order)
        f_param.star <- psd_varma(lambda, param__phi.star, sigma=sigma.fit)$psd
        
        # plotMPsd(f_param.star, main="proposed")
        rejectedPhi <- F
        if (any(apply(f_param.star, 3, hasEigenValueSmallerZero, TOL=NUMERICAL_THRESH))) { # stay positive definite
          if (verbose) print_warn("Discarding f_param proposal, because of positive definiteness")
          param__beta[,i+1] <- param__beta[,i]
          param__phi[,,i+1] <- param__phi[,,i]
          rejectedPhi <- T
        } 
        if (numericalUnstable(f_param.star, excludeBoundary=F, TOL=NUMERICAL_THRESH)) { # stay numerically stable
          if (verbose) print_warn("Discarding f_param proposal, because of numerical instablity")
          param__beta[,i+1] <- param__beta[,i]
          param__phi[,,i+1] <- param__phi[,,i]
          rejectedPhi <- T
        }
        if (!rejectedPhi) {
          f_param_half.star <- chol_cube(f_param.star, excludeBoundary=F)
          phi.fit.star <- list(ar=param__phi.star,
                               f_param_half=f_param_half.star,
                               f_param_half_trans=trans_cube(f_param_half.star),
                               beta=param__beta.star,  ## 
                               mu_beta=mu_beta,        ## include stuff for prior computation, too
                               V_beta_inv=V_beta_inv)  ##
          f.phi.star <- lpost_matrixGamma(omega=omega,
                                          FZ=FZ,
                                          r=r[,i+1],
                                          U=U[,,,i+1],
                                          Z=Z[,i+1],
                                          k=k[,i+1],
                                          C_alpha=C_alpha,
                                          omega_fun=omega_fun,
                                          k.theta=k.theta,
                                          db.list=db.list,
                                          eta=eta,
                                          Sigma_fun=Sigma_fun,
                                          corrected=corrected,
                                          phi=phi.fit.star, #
                                          sigma_ar=sigma.fit,
                                          prior.q=prior.q,
                                          prior.cholesky=prior.cholesky,
                                          excludeBoundary=T, # note
                                          verbose=verbose)
          f.phi <- f.store
          alpha5 <- min(0, f.phi.star - f.phi) # Note: Normal proposal symmetric for beta
          if (log(runif(1,0,1)) < alpha5) {
            # accept
            param__beta[,i+1] <- param__beta.star
            param__beta.old <- param__beta.star
            param__phi[,,i+1] <- param__phi.star
            phi.fit <- phi.fit.star
            f.store <- f.phi.star
          } else {
            # reject
            param__beta[,i+1] <- param__beta[,i]
            param__phi[,,i+1] <- param__phi[,,i]
          }
        }
      }
    }
    
    
    ##
    ## Step 6: Sample parameter of interest
    ##
    
    # MH Step for theta
    if (theta_dim > 0) {
      if (corrected) {
        q_for_theta <- get_f_matrix(U[,,,i+1], r[,i+1], Z[,i+1], k[,i+1], db.list, prior.cholesky)
        f_for_theta <- mult_cube(mult_cube(phi.fit$f_param_half, q_for_theta), phi.fit$f_param_half_trans)
        previous_theta <- theta[,i] # might change for corrected in case of double steps
      } else {
        f_for_theta <- get_f_matrix(U[,,,i+1], r[,i+1], Z[,i+1], k[,i+1], db.list, prior.cholesky)
        previous_theta <- theta[,i] # might change for corrected in case of double steps
      }
      theta_prop <- propose_next_theta(data=data, f=f_for_theta, previous_theta=previous_theta, NULL)
      theta_star <- theta_prop$theta_star
      noise_star <- get_noise(data, theta_star)  # noise = data - signal
      FZ_star <- mdft(noise_star)  # Frequency domain
      
      ##
      f.theta.star <- lpost_matrixGamma(omega=omega,
                                        FZ=FZ_star, # note
                                        r=r[,i+1],
                                        U=U[,,,i+1],
                                        Z=Z[,i+1],
                                        k=k[,i+1],
                                        C_alpha=C_alpha,
                                        omega_fun=omega_fun,
                                        k.theta=k.theta,
                                        db.list=db.list,
                                        eta=eta,
                                        Sigma_fun=Sigma_fun,
                                        corrected=corrected,
                                        phi=phi.fit,
                                        sigma_ar=sigma.fit,
                                        prior.q=prior.q,
                                        prior.cholesky=prior.cholesky,
                                        excludeBoundary=F, # note
                                        verbose=verbose)
      f.theta <- lpost_matrixGamma(omega=omega,
                                   FZ=FZ, # note
                                   r=r[,i+1],
                                   U=U[,,,i+1],
                                   Z=Z[,i+1],
                                   k=k[,i+1],
                                   C_alpha=C_alpha,
                                   omega_fun=omega_fun,
                                   k.theta=k.theta,
                                   db.list=db.list,
                                   eta=eta,
                                   Sigma_fun=Sigma_fun,
                                   corrected=corrected,
                                   phi=phi.fit,
                                   sigma_ar=sigma.fit,
                                   prior.q=prior.q,
                                   prior.cholesky=prior.cholesky,
                                   excludeBoundary=F, # note
                                   verbose=verbose)
      ##
      
      alphaTheta <- min(0, f.theta.star + lprior_theta(theta_star) - 
                          f.theta - lprior_theta(theta[,i]) +
                          theta_prop$lprop_previous_theta - 
                          theta_prop$lprop_theta_star)
      if (log(runif(1,0,1)) < alphaTheta) {
        theta[,i+1] <- theta_star
        f.store <- lpost_matrixGamma(omega=omega,
                                     FZ=FZ_star, # note
                                     r=r[,i+1],
                                     U=U[,,,i+1],
                                     Z=Z[,i+1],
                                     k=k[,i+1],
                                     C_alpha=C_alpha,
                                     omega_fun=omega_fun,
                                     k.theta=k.theta,
                                     db.list=db.list,
                                     eta=eta,
                                     Sigma_fun=Sigma_fun,
                                     corrected=corrected,
                                     phi=phi.fit,
                                     sigma_ar=sigma.fit,
                                     prior.q=prior.q,
                                     prior.cholesky=prior.cholesky,
                                     excludeBoundary=T, # note
                                     verbose=verbose)
      } else {
        theta[,i+1] <- theta[,i]
      }
    }
    
    
    ##
    ## Final step: Update some traces
    ##
    #print("Step final")
    lpostTrace[i+1] <- f.store + lprior_theta(theta[,i+1])
    lpriorTrace[i+1] <- lprior_matrixGamma(r=r[,i+1],
                                           U=U[,,,i+1],
                                           Z=Z[,i+1],
                                           k=k[,i+1],
                                           C_alpha=C_alpha,
                                           omega_fun=omega_fun,
                                           k.theta=k.theta,
                                           eta=eta,
                                           Sigma_fun=Sigma_fun,
                                           phi=phi.fit,
                                           verbose=verbose) +
      lprior_theta(theta[,i+1])
    llikeTrace[i+1] <- lpostTrace[i+1] - lpriorTrace[i+1]
    
  } # END MCMC LOOP
  
  ##
  ## Post processing
  ##
  #print("Postprocessing results ...")
  keep <- seq(burnin+1, Ntotal, by=thin)
  
  if (prior.cholesky) {
    k <- k[,keep]
  } else {
    k <- array(data=k[,keep], dim=c(1, length(keep))) 
  }
  r <- r[,keep]
  Z <- Z[,keep]
  U <- U[,,,keep]
  U__phi <- U__phi[,,keep]
  if (corrected && toggle) {
    param__beta <- param__beta[,keep]
    param__phi <- param__phi[,,keep]
  } else {
    param__beta <- NULL
    param__phi <- NULL
  }
  theta <- theta[,keep,drop=F]
  
  W <- array(dim=c(d,d,L,length(keep))) # for convenience
  fpsd.sample <- array(NA, dim=c(d, d, N, length(keep)))
  
  
  # store the coherence
  comb_t_n <- ncol(combn(d, 2))
  comb_t <- combn(d, 2)
  coherence.sample <- array(NA, dim = c(N, comb_t_n, length(keep)))
  
  # Corrected (not toggled yet)
  if (corrected && prior.q && (!toggle)) {
    f_param_half <- chol_cube(psd_varma(lambda, phi.fit$ar, sigma=sigma.fit)$psd, excludeBoundary=F)
    f_param_half_trans <- trans_cube(f_param_half)
  }
  for (isample in 1:length(keep)) {
    if (corrected && prior.q) {
      if (toggle) {
        f_param_half <- chol_cube(psd_varma(lambda, param__phi[,,isample], sigma=sigma.fit)$psd, excludeBoundary=F)
        f_param_half_trans <- trans_cube(f_param_half)
      }
      q_sample <- get_f_matrix(U[,,,isample], r[,isample], Z[,isample], k[,isample], db.list, prior.cholesky)
      f_sample <- mult_cube(mult_cube(f_param_half, q_sample), f_param_half_trans) # "prior on q=f/f_param"
    } else {
      f_sample <- get_f_matrix(U[,,,isample], r[,isample], Z[,isample], k[,isample], db.list, prior.cholesky)
    }
    fpsd.sample[,,,isample] <- realValuedPsd(f_sample)
    for (ipair in 1:comb_t_n){
      index1 <- comb_t[1, ipair]
      index2 <- comb_t[2, ipair]
      f1 <- fpsd.sample[index1, index1,,isample]
      f2 <- fpsd.sample[index2, index2,,isample]
      f_Re <- fpsd.sample[index1, index2,,isample]
      f_Im <- fpsd.sample[index2, index1,,isample]
      coherence.sample[,ipair, isample] <- sqrt(f_Re^2+f_Im^2)/sqrt(f1*f2)
    }
  }
  fpsd.s <- apply(fpsd.sample, c(1,2,3), median)
  fpsd.mean <- apply(fpsd.sample, c(1,2,3), mean)
  # pointwise 90%
  fpsd.s05 <- apply(fpsd.sample, c(1,2,3), quantile, 0.05)
  fpsd.s95 <- apply(fpsd.sample, c(1,2,3), quantile, 0.95)
  # # pointwise 95%
  # fpsd.s025 <- apply(fpsd.sample, c(1,2,3), quantile, 0.025)
  # fpsd.s975 <- apply(fpsd.sample, c(1,2,3), quantile, 0.975)
  # # pointwise 99%
  # fpsd.s005 <- apply(fpsd.sample, c(1,2,3), quantile, 0.005)
  # fpsd.s995 <- apply(fpsd.sample, c(1,2,3), quantile, 0.995)
  
  alpha_uci <- 0.1 # same as in 1D
  uci_tmp <- uci_matrix(fpsd.sample, alpha=alpha_uci)
  fpsd.uci05 <- uci_tmp$fpsd.uci05
  fpsd.uci95 <- uci_tmp$fpsd.uci95
  uuci_tmp <- uci_matrix(fpsd.sample, alpha=alpha_uci, uniform_among_components=T)
  fpsd.uuci05 <- uuci_tmp$fpsd.uci05
  fpsd.uuci95 <- uuci_tmp$fpsd.uci95
  
  
  #squared coherence (prior squared) (median and 90% pointwise region)
  coherence.s <- apply(coherence.sample, c(1, 2), median)
  coherence.s05 <- apply(coherence.sample, c(1, 2), quantile, 0.05)
  coherence.s95 <- apply(coherence.sample, c(1, 2), quantile, 0.95)
  
  # # Construct forecasts
  # N_FORECAST <- 5
  # rm("db.list") # clear some memory
  # db.list_forecast <- dbList(n+N_FORECAST, kmax, normalized=normalizedBernsteins, bernstein_l, bernstein_r) 
  # N_MCMC_IT <- length(keep)
  # data_forecast <- array(data=NA, dim=c(N_FORECAST, d, N_MCMC_IT))
  # for (isample in 1:N_MCMC_IT) {
  #   if (!(i%%print_interval)) {
  #     cat("forecasting ", isample, "/", N_MCMC_IT, "\n", sep="")
  #   }
  #   noise <- get_noise(data, theta[,isample])
  #   if (corrected && prior.q) {
  #     stop("TODO: Include forecast to corrected branch")
  #   }
  #   else {
  #     f_forecast <- get_f_matrix(U[,,,isample], 
  #                                r[,isample], 
  #                                Z[,isample], 
  #                                k[,isample], 
  #                                db.list_forecast, 
  #                                prior.cholesky)
  #   }
  #   noise_forecast <- vnp_forecast(noise, f_forecast, N_FORECAST)
  #   noise_all <- rbind(noise, noise_forecast)
  #   data_forecast[,,isample] <- get_data(noise_all, theta[,isample])[-(1:nrow(noise)),]
  # }
  
  time_e <- Sys.time()
  time_d <- time_e - time_b
  
  ##
  ## Return stuff
  ##
  return(list(#data=data,
    #k=k,
    #r=r,
    #Z=Z,
    #U=U,
    #U__phi=U__phi,
    #W=W,
    fpsd.s=fpsd.s,
    fpsd.mean=fpsd.mean,
    fpsd.s05=fpsd.s05,
    fpsd.s95=fpsd.s95,
    # fpsd.s025=fpsd.s025,
    # fpsd.s975=fpsd.s975,
    # fpsd.s005=fpsd.s005,
    # fpsd.s995=fpsd.s995,
    # fpsd.uci05=fpsd.uci05,
    # fpsd.uci95=fpsd.uci95,
    fpsd.uuci05=fpsd.uuci05,
    fpsd.uuci95=fpsd.uuci95,
    coherence.s = coherence.s,
    coherence.s05 = coherence.s05,
    coherence.s95 = coherence.s95,
    llikeTrace=llikeTrace,
    lpostTrace=lpostTrace, # log posterior: don't discard burnin to investigate convergence
    lpriorTrace=lpriorTrace,
    #param__phi=param__phi,
    #theta=theta,
    time_d))#,
  #data_forecast=data_forecast))
}

