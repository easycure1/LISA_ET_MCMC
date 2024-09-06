##
## The transformation \mathcal B \colon \mathbb C^{d\times d} \to \mathbb R^{2d\times 2d}
## from Hannan1970, page 224.
##
B_matrix_trans <- function(A) {
  rbind(cbind(Re(A), -Im(A)), cbind(Im(A), Re(A)))
}

##
## Time series model X_t=e_t, E[e_t]=0, without additional parameters
##
psd_dummy_model <- function() {
  theta_dim <- 0
  warning("excludeBoundary parameter is deprecated")
  warning("mean centering not taken into account for forcasting!")
  # excludeBoundary <- T 
  get_noise <- function(data, theta, ...) {
    # mean centered version
    t(t(data)-apply(data,2,mean))
  }
  get_data <- function(noise, theta, ...) {
    noise
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    # dummy
    numeric(0)
  }
  initialize_theta <- function(data, ...) {
    # dummy
    numeric(0)
  }
  lprior_theta <- function(theta, ...) {
    # dummy
    0
  }
  model_params <- list(theta_dim=theta_dim,
                       get_noise=get_noise,
                       get_data=get_data,
                       propose_next_theta=propose_next_theta,
                       lprior_theta=lprior_theta,
                       initialize_theta=initialize_theta)
  return(model_params)
}

##
## Multivariate linear model, with design matrix X \in \mathbb R^{nd\times r},
## where r is the dimension of theta.
## 
## Model: \mathbb R^{nd} Y = X\theta + e, (Y data, X design matrix, e nuisance TS)
##
## Note that Y is obtained from the vech tansformation (see gibbs_multivariate_util.R)
##
## X_fun: function to compute the design matrix, depending on the sample size n and d
## ...  : further arguments to  be parsed to call of X_fun
##
normal_linear_model <- function(X_fun, n, d, lprior_theta=function(theta){1},
                                     prop.scaling=1) {
  X <- X_fun(n, d)
  theta_dim <- ncol(X) # = r
  stopifnot(nrow(X) %% d == 0)  # nrow(X) = n*d
  stopifnot(n == nrow(X) / d) # n <- nrow(X) / d
  warning("excludeBoundary parameter is deprecated")
  print("Warning: Using random walk scheme!")
  # excludeBoundary <- F 
  # construct the FT'ed design matrix X_tilde from X
  X_tilde <- matrix(NA, nrow=nrow(X), ncol=ncol(X))
  for (j in 1:ncol(X)) {
    X_tilde[,j] <- vech(mdft(ivech(X[,j], d), real=T))
  }
  get_noise <- function(data, theta, ...) {
    stopifnot(nrow(data)==n)
    stopifnot(ncol(data)==d)
    ivech(vech(data) - X %*% theta, d)
  }
  get_data <- function(noise, theta, ...) {
    nn <- nrow(noise)
    stopifnot(ncol(noise)==d)
    X_nn <- X_fun(nn, d)
    ivech(vech(noise) + X_nn %*% theta, d)
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    #
    # Propose from conjugate under Whittle's likelihood with flat prior on theta
    #
    # Idea: \mathbb R^{nd} Y = X\theta + e (Y data, X design matrix, e nuisance TS)
    #   --FZ--> Y_tilde=X_tilde\theta+e_tilde, e_tilde ~ Whittle(f)
    #   --CFZ-> Y_tilde_c=X_tilde_c\theta+e_tilde_c, e_tilde ~ N(0, I_{nd})
    # and the marginal posterior for theta under ~iid~N(0,I_{nd}) with flat prior is known
    # (see Section 9.5.2 in Christensen et al)
    #
    stopifnot(nrow(data)==n)
    stopifnot(ncol(data)==d)
    Y_tilde <- vech(mdft(data, real=T))
    N <- floor(n/2)-1
    # Apply Whittle correction to get white noise residuals
    Y_tilde_c <- rep(NA, n*d)
    X_tilde_c <- matrix(NA, nrow=nrow(X), ncol=ncol(X))
    # . lambda=0 --> real-valued
    D_0 <- 2 * pi * f[,,1] # TODO Is this the right scaling?
    C_0 <- solve(Re(chol_cpp(D_0)))
    ind_set_0 <- 1:d
    Y_tilde_c[ind_set_0] <- C_0 %*% Y_tilde[ind_set_0]
    X_tilde_c[ind_set_0,] <- C_0 %*% X_tilde[ind_set_0,]
    # . lambda \in (0,\pi) --> hpd-valued
    for (j in 1:N) {
      ind_set_j <- (d+1+(j-1)*2*d):(d+j*2*d)
      stopifnot(all(is.na(Y_tilde_c[ind_set_j])))
      stopifnot(all(is.na(X_tilde_c[ind_set_j,])))
      D_j <- B_matrix_trans(2*pi*f[,,j+1]) # TODO Is this the right scaling?
      C_j <- solve(Re(chol_cpp(D_j)))
      Y_tilde_c[ind_set_j] <- C_j %*% Y_tilde[ind_set_j]
      X_tilde_c[ind_set_j,] <- C_j %*% X_tilde[ind_set_j,]
    }
    # . lambda=pi --> real-valued
    if (!(n%%2)) {
      ind_set_pi <- ((n-1)*d+1):(n*d)
      stopifnot(all(is.na(Y_tilde_c[ind_set_pi])))
      stopifnot(all(is.na(X_tilde_c[ind_set_pi,])))
      stopifnot(dim(f)[[3]]==N+2)
      D_pi <- 2*pi*f[,,N+2] # TODO Is this the right scaling?
      C_pi <- solve(Re(chol_cpp(D_pi)))
      Y_tilde_c[ind_set_pi] <- C_pi %*% Y_tilde[ind_set_pi]
      X_tilde_c[ind_set_pi,] <- C_pi %*% X_tilde[ind_set_pi,]
    } else {
      stopifnot(dim(f)[[3]]==N+1)
    }
    stopifnot(all(!is.na(Y_tilde_c)))
    stopifnot(all(!is.na(X_tilde_c)))
    # Propose from joint marginal posterior (under Whittle + flat prior)
    Sigma.inv <- t(X_tilde_c) %*% X_tilde_c
    # Sigma <- ginv(Sigma.inv)
    Sigma <- solve(Sigma.inv)
    mu <- previous_theta # Sigma %*% t(X_tilde_c) %*% Y_tilde_c
    theta_star <- mvrnorm(1, mu, prop.scaling*Sigma)
    # Compute log proposal density at proposal and previous (for MH acceptance probability)
    # Note: unnormalized multivariate normal
    lprop_theta_star <- 0 #as.numeric(
      #-1/2 * t(theta_star-mu) %*% (Sigma.inv/prop.scaling) %*% (theta_star-mu))
    lprop_previous_theta <- 0 #as.numeric(
      #-1/2 * t(previous_theta-mu) %*% (Sigma.inv/prop.scaling) %*% (previous_theta-mu))
    list(theta_star=theta_star,
         lprop_theta_star=lprop_theta_star,
         lprop_previous_theta=lprop_previous_theta)
  }
  initialize_theta <- function(data, ...) {
    stopifnot(nrow(data)==n)
    stopifnot(ncol(data)==d)
    Y <- vech(data)
    ginv(t(X) %*% X) %*% t(X) %*% Y # ML
  }
  model_params <- list(theta_dim=theta_dim,
                       get_noise=get_noise,
                       get_data=get_data,
                       propose_next_theta=propose_next_theta,
                       lprior_theta=lprior_theta,
                       initialize_theta=initialize_theta)
  return(model_params)
}

##
## Design matrix for Multivariate linear trend model *with common slope*: 
##                     X_t=bt + mu + e_t, e_t~Nuisance TS
## with mu=(mu_1,...,mu_d) and b (scalar!).
designMatrix_linearTrend_commonSlope <- function(n,d) {
  #
  # theta[1,...,d] = mu
  # theta[d+1] = b
  #
  X <- matrix(NA, nrow=n*d, ncol=d+1)
  X[,1:d] <- designMatrix_meanModel(n,d)
  X_dp1_tmp <- matrix(rep(1:n,d), ncol=d, nrow=n)
  X[,d+1] <- vech(X_dp1_tmp)
  X
}

##
## Design matrix for Multivariate mean model:
##                     X_t= mu + e_t, e_t~Nuisance TS
## with mu=(mu_1,...,mu_d).
designMatrix_meanModel <- function(n,d) {
  #
  # theta[1,...,d] = mu
  #
  X <- matrix(NA, nrow=n*d, ncol=d)
  for (j in 1:d) {
    e_j <- rep(0,d)
    e_j[j] <- 1
    X[,j] <- rep(e_j, n)
  }
  X
}

##
## Design matrix for linear trend temperature model
## Two time series, splitted seasonal (d=8)
##
## Assuming 8 individual intercepts and 4 slopes,
## that vary among season but NOT among subjects
##
designMatrix_seasonalTemperatureModel <- function(n, d) {
  stopifnot(d==8)
  # 
  # d = 8, r=12
  # theta[1,...,8] = mu
  # theta[9,...,12] = b
  #
  r <- 12
  d <- 8
  X_tmp <- suppressWarnings(matrix(c(1,1,0,0,0,0,0,0,0,0), nrow=8, ncol=4)) # a bit hacky
  X <- matrix(NA, nrow=0, ncol=r)
  for (t in 1:n) {
    X <- rbind(X, cbind(diag(d), t*X_tmp))
  }
  X
}

##
## Design matrix for Multivariate linear trend model *with individual slopes*: 
##                     X_t=bt + mu + e_t, e_t~Nuisance TS
## with mu=(mu_1,...,mu_d) and b=(b_1,...,b_d).
designMatrix_linearTrend_individualSlope <- function(n,d) {
  #
  # theta[1,...,d] = mu
  # theta[d+1,...,2*d] = b
  #
  r <- 2*d
  X <- matrix(NA, nrow=0, ncol=r)
  for (t in 1:n) {
    X <- rbind(X, cbind(diag(d), t*diag(d)))
  }
  X
}

# ##
# ## Simple mean model: X_t=mu+e_t, mu~N_d(mu_0,Sigma_0), e_t~Nuisance TS
# ##
# normal_mean_model <- function(d,
#                               mu.mu0=rep(0,d), # = mu_0, a priori mean of mu
#                               mu.Sigma0=diag(1e6,d), # = sigma_0, a priori covariance of mu
#                               mu.prop.scaling=1 # proposal scaling parameter
# ) {
#   stop("Deprecated -- Please use normal_linear_model with designMatrix_meanModel")
#   theta_dim <- d
#   warning("excludeBoundary parameter is deprecated")
#   # excludeBoundary <- F
#   mu.Sigma0_inv <- solve(mu.Sigma0)
#   get_noise <- function(data, theta, ...) {
#     t(t(data)-theta)
#   }
#   propose_next_theta <- function(data, f, previous_theta, ...) {
#     n <- nrow(data)
#     d <- ncol(data)
#     # sigma_prop <- mu.prop.scaling / n * (5*Re(f[,,1]) + diag(1e-4/n,d)) # some shrinkage in proposal covariance
#     sigma_prop <- mu.prop.scaling / n * 5 * Re(f[,,1])
#     MASS::mvrnorm(n=1, mu=apply(data, 2, mean), Sigma=sigma_prop)
#   }
#   initialize_theta <- function(data, ...) {
#     apply(data, 2, mean)
#   }
#   lprior_theta <- function(theta, ...) { # careful: unnormalized log prior density
#     theta_c <- theta - mu.mu0
#     -as.numeric(t(theta_c) %*% mu.Sigma0_inv %*% theta_c)
#   }
#   model_params <- list(theta_dim=theta_dim,
#                        get_noise=get_noise,
#                        propose_next_theta=propose_next_theta,
#                        lprior_theta=lprior_theta,
#                        initialize_theta=initialize_theta)
#   return(model_params)
# }

# ##
# ## Multivariate linear trend model *with common slope*:
# ##                     X_t=bt + mu + e_t, mu~N_d(mu_0,Sigma_0),
# ##                     b~N(mu_1,sigma_1^2), e_t~Nuisance TS
# ## with mu=(mu_1,...,mu_d) and b (scalar!) end e_t a priori independent
# ##
# normal_linearTrend_model <- function(d, mu.mu0=0, mu.Sigma0=diag(1e6,d),
#                                 b.mu1=0, b.sd=1,
#                                 prop.scaling=1) {
#   stop("Deprecated: Please use normal_linear_model with designMatrix_linearTrend_commonSlope")
#   #
#   # theta[1,...,d] = mu
#   # theta[d+1] = b
#   #
#   theta_dim <- d+1
#   warning("excludeBoundary parameter is deprecated")
#   # excludeBoundary <- F
#   mu.Sigma0_inv <- solve(mu.Sigma0)
#   lm_res <- function(x, b) {
#     n <- length(x)
#     Yt <- 1:n
#     x - b*Yt
#   }
#   get_noise <- function(data, theta, ...) {
#     data_m <- t(t(data)-theta[1:d]) # subtract means
#     apply(data_m, 2, lm_res, theta[d+1]) # subtract slope
#   }
#   propose_next_theta <- function(data, f, previous_theta, ...) {
#     #
#     # Propose from conjugate under Whittle's likelihood with flat prior on theta
#     #
#     # Idea: \mathbb R^{nd} Y = X\theta + e (Y data, X design matrix, e nuisance TS)
#     #   --FZ--> Y_tilde=X_tilde\theta+e_tilde, e_tilde ~ Whittle(f)
#     #   --CFZ-> Y_tilde_c=X_tilde_c\theta+e_tilde_c, e_tilde ~ N(0, I_{nd})
#     # and the marginal posterior for theta under ~iid~N(0,I_{nd}) with flat prior is known
#     # (see Section 9.5.2 in Christensen et al)
#     #
#     Y_tilde <- vech(mdft(data, real=T))
#     n <- nrow(data)
#     N <- floor(n/2)-1
#     d <- ncol(data)
#     # construct design matrix X
#     X <- matrix(NA, nrow=n*d, ncol=d+1)
#     for (j in 1:d) {
#       e_j <- rep(0,d)
#       e_j[j] <- 1
#       X[,j] <- rep(e_j, n)
#     }
#     X_dp1_tmp <- matrix(rep(1:n,d), ncol=d, nrow=n)
#     X[,d+1] <- vech(X_dp1_tmp)
#     # construct the FT'ed design matrix X_tilde from X
#     X_tilde <- matrix(NA, nrow=nrow(X), ncol=ncol(X))
#     for (j in 1:ncol(X)) {
#       X_tilde[,j] <- vech(mdft(ivech(X[,j], d), real=T))
#     }
#     # for (j in 1:d) {
#     #   X_j_tmp <- matrix(0, nrow=n, ncol=d)
#     #   X_j_tmp[,j] <- rep(1,n)
#     #   X_tilde[,j] <- vech(mdft(X_j_tmp, real=T))
#     # }
#     # X_tilde[,d+1] <- vech(mdft(X_dp1_tmp, real=T))
#     # Apply Whittle correction to get white noise residuals
#     Y_tilde_c <- rep(NA, n*d)
#     X_tilde_c <- matrix(NA, nrow=n*d, ncol=d+1)
#     # . lambda=0 --> real-valued
#     D_0 <- pi * f[,,1] # TODO Is this the right scaling?
#     C_0 <- solve(Re(chol_cpp(D_0)))
#     ind_set_0 <- 1:d
#     Y_tilde_c[ind_set_0] <- C_0 %*% Y_tilde[ind_set_0]
#     X_tilde_c[ind_set_0,] <- C_0 %*% X_tilde[ind_set_0,]
#     # . lambda \in (0,\pi) --> hpd-valued
#     for (j in 1:N) {
#       ind_set_j <- (d+1+(j-1)*2*d):(d+j*2*d)
#       stopifnot(all(is.na(Y_tilde_c[ind_set_j])))
#       stopifnot(all(is.na(X_tilde_c[ind_set_j,])))
#       D_j <- B_matrix_trans(2*pi*f[,,j+1]) # TODO Is this the right scaling?
#       C_j <- solve(Re(chol_cpp(D_j)))
#       Y_tilde_c[ind_set_j] <- C_j %*% Y_tilde[ind_set_j]
#       X_tilde_c[ind_set_j,] <- C_j %*% X_tilde[ind_set_j,]
#     }
#     # . lambda=pi --> real-valued
#     if (!(n%%2)) {
#       ind_set_pi <- ((n-1)*d+1):(n*d)
#       stopifnot(all(is.na(Y_tilde_c[ind_set_pi])))
#       stopifnot(all(is.na(X_tilde_c[ind_set_pi,])))
#       stopifnot(dim(f)[[3]]==N+2)
#       D_pi <- pi*f[,,N+2] # TODO Is this the right scaling?
#       C_pi <- solve(Re(chol_cpp(D_pi)))
#       Y_tilde_c[ind_set_pi] <- C_pi %*% Y_tilde[ind_set_pi]
#       X_tilde_c[ind_set_pi,] <- C_pi %*% X_tilde[ind_set_pi,]
#     } else {
#       stopifnot(dim(f)[[3]]==N+1)
#     }
#     stopifnot(all(!is.na(Y_tilde_c)))
#     stopifnot(all(!is.na(X_tilde_c)))
#     # Propose from joint marginal posterior (under Whittle + flat prior)
#     Sigma.inv <- t(X_tilde_c) %*% X_tilde_c
#     Sigma <- ginv(Sigma.inv)
#     mu <- Sigma %*% t(X_tilde_c) %*% Y_tilde_c
#     mvrnorm(1, mu, prop.scaling*Sigma)
#   }
#   initialize_theta <- function(data, ...) {
#     x <- data
#     n <- nrow(x)
#     Yt <- 1:n
#     lm_coefs <- apply(data, 2, function(x){lm(x~Yt)$coef})
#     mu_start <- lm_coefs[1,]
#     b_start <- mean(lm_coefs[2,])
#     c(mu_start, b_start)
#   }
#   lprior_theta <- function(theta, ...) { # careful: unnormalized log prior density
#     # mu_c <- theta[1:d] - mu.mu0
#     # b_c <- theta[d] - b.mu1
#     # -as.numeric(t(mu_c) %*% mu.Sigma0_inv %*% mu_c)/2 -
#     #   -(b_c/b.sd)^2
#     1
#   }
#   model_params <- list(theta_dim=theta_dim,
#                        get_noise=get_noise,
#                        propose_next_theta=propose_next_theta,
#                        lprior_theta=lprior_theta,
#                        initialize_theta=initialize_theta)
#   return(model_params)
# }