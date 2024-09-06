n <- 819200

d <- 2
param_ar <- rbind(c(.1, 0, 0, 0), c(0, -.3, 0, -.5))
param_ma <- matrix(nrow=2, ncol=0)
param_sigma <- matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2)

data <- beyondWhittle::sim_varma(model=list(ar=param_ar, ma=param_ma, sigma=param_sigma),
                                 n=n, d=2)



source('_include_multivariate_cluster.R')

Ntotal <- 3000 # TODO 80000
burnin <- 1000  # TODO 30000
Nadaptive <- burnin
thin <- 1 # TODO 5
print_interval <- 500 # TODO 10000


my_Sigma_fun <- function(x) {
  Sigma_fun_eye(x,d=2) * 1e4
}



mcmc_params <- list(Ntotal=Ntotal,
                    burnin=burnin,
                    thin=thin,
                    print_interval=print_interval,
                    numerical_thresh=1e-12,
                    Nadaptive=Nadaptive,
                    adaption.batchSize=50,
                    adaption.targetAcceptanceRate=0.44,
                    verbose=F)
prior_params <- list(prior.cholesky=F,
                     #var.order=p_order,
                     eta=2,
                     omega_fun=create_omega_fun_from_beta_density(1,1,1),
                     Sigma_fun=my_Sigma_fun,
                     k.theta=0.01,
                     kmax=300,
                     bernstein_l=0.1, # note
                     bernstein_r=0.9, # note
                     coarsened=F, # note
                     L=32)
model_params <- psd_dummy_model()

mcmc_vnp_avg <- gibbs_m_nuisance(data=data,
                                 mcmc_params=mcmc_params,
                                 seg_n=200,
                                 corrected=F,
                                 prior_params=prior_params,
                                 model_params=model_params)

save(mcmc_vnp_avg, file = 'mcmc_vnp_avg.RData')



