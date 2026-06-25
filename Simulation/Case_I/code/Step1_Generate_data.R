##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


B <- 3
K <- 3
G <- 1000

## Y: dim B * (n_b (b=1,...,B) * G) matrix list
## mu: dim B * (G * K) matrix list
## n_vec: dim B vector
## pi: dim B * K matrix
## Z: dim B * n_b (b=1,...,B) vector list
## alpha: dim G vector
## xi: dim G * K matrix with first column being zeros
## gamma: dim B * G matrix with first row being zeros
## sigma: dim B * G matrix
## tau: scalar

n_vec <- c(100, 110, 120)

pi_real <- matrix(c(0.2,0.2,0.6,
                    0.1,0.8,0.1,
                    0.6,0.1,0.3), nrow = B, ncol = K, byrow = T)


Z_real <- list()
for(b in 1:B){
  Z_real[[b]] <- c(rep(1, pi_real[b,1]*n_vec[b]), 
                   rep(2, pi_real[b,2]*n_vec[b]), 
                   rep(3, pi_real[b,3]*n_vec[b]))
}

alpha_real <- rep(2, times=G)
xi_real <- matrix(c(rep(0,G), rep(2,G), rep(3,G)), nrow = G, ncol = K)
gamma_real <- matrix(c(rep(0,G), rep(1,G), rep(2,G)), nrow = B, ncol = G, byrow = T)
sigma_real <- matrix(c(rep(sqrt(0.1),G), rep(sqrt(0.2),G), rep(sqrt(0.15),G)), nrow = B, ncol = G, byrow = T)



## Set seed for data generation in simulation case I
set.seed(1)


## Data generation
mu <- list()
Y <- list()
for(b in 1:B){
  mu[[b]] <- matrix(rep(alpha_real, K), nrow =G, ncol = K) + xi_real + matrix(rep(gamma_real[b, ], K), nrow = G, ncol = K)
  Y[[b]] <- matrix(NA, nrow = n_vec[b], ncol = G)
  for(i in 1:n_vec[b]){
    Y[[b]][i, ] = sapply(1:G, function(g){rnorm(1, mean = mu[[b]][g, Z_real[[b]][i]], sd = sigma_real[b, g])})
  }
}


## True subgroup mean
truemean_by_subgroup <- list()
for(k in 1:K) truemean_by_subgroup[[k]] = matrix(rep(alpha_real+xi_real[, k], times = length(which(unlist(Z_real)==k))), ncol=G, byrow = TRUE)

## Save true subgroup mean
save(truemean_by_subgroup, file = "../input_data/truemean_by_subgroup.RData")


## Save raw simulated data and subgroup label
save(Y, file = "../input_data/Y.RData")
save(Z_real, file = "../input_data/Z_real.RData")




