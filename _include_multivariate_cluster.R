##
## Dummy file to include all files that may be needed
## for multivariate spectral inference
##

Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # lambda function using C++11
library(Rcpp)

sourceCpp("gibbs_util.cpp")
source("gibbs_util.R")

sourceCpp("AGammaProcess_util.cpp")
source("AGammaProcess_util.R")

source("gibbs_gamma_matrix_subordinator.R")
sourceCpp("gibbs_gamma_matrix_subordinator_core.cpp")
source("gibbs_gamma_matrix_subordinator_core.R")

sourceCpp("gibbs_multivariate_util.cpp")
source("gibbs_multivariate_util.R")

source("gibbs_VAR_nuisance.R")

sourceCpp("mittelbach_util.cpp")
source("mittelbach_util.R")

source("mNuisanceModels.R")
