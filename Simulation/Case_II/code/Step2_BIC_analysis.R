##############################################################
#################     Simulation Case II    ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(matrixStats)


##################### Hyper-parameters ##########################

c <- 2
lambda <- 1
eta <- 0.1
p <- 2
q <- 2
delta <- 20
r <- 2
s <- 2


############### Update Variational Parameters ###################

## Z - varphi[[b]][i, k]: dim B * (n_b * K) matrix list
update_varphi <- function(c_til, lambda_til, omega, tau_til, rho, r_til, s_til){
  dig_c <- digamma(c_til)
  penalty <- -0.5 * (omega^2 + tau_til^2)
  
  temp <- vector("list", B)
  
  for (b in 1:B) {
    a_b <- r_til[b, ] / s_til[b, ]
    Y_centered <- sweep(Y[[b]], 2, lambda_til + rho[b, ], "-")
    X_b <- sweep(Y_centered, 2, a_b, "*")
    
    log_prob <- X_b %*% omega
    const_k <- colSums(sweep(penalty, 1, a_b, "*"))
    log_prob <- sweep(log_prob, 2, dig_c[b, ] + const_k, "+")
    
    row_max <- apply(log_prob, 1, max)
    exp_log_prob <- exp(log_prob - row_max)
    temp[[b]] <- exp_log_prob / rowSums(exp_log_prob)
  }
  
  return(temp)
}

## pi - c_til[b, k]: dim B * K matrix
update_c_til <- function(varphi){
  return( t(sapply(1:B, function(b){c + colSums(varphi[[b]])})) )
}

## alpha - lambda_til[g]: dim G vector
update_lambda_til <- function(varphi, omega, rho, r_til, s_til){
  temp1_vec <- numeric(G)
  for(b in 1:B){
    a_b <- r_til[b, ] / s_til[b, ]              # length G
    phi_b <- varphi[[b]]                        # n_b x K
    phi_omega <- phi_b %*% t(omega)
    resid_b <- sweep(Y[[b]] - phi_omega, 2, rho[b, ], FUN = "-")
    temp1_vec <- temp1_vec + colSums(resid_b) * a_b
  }
  temp1_vec <- temp1_vec + lambda / eta^2
  
  temp2_vec <- 1 / (colSums(n_vec * (r_til / s_til)) + 1 / eta^2)
  
  return(temp1_vec * temp2_vec)
}

## alpha - eta_til[g]: dim G vector
update_eta_til <- function(r_til, s_til){
  return( sqrt(1 / (colSums(n_vec * (r_til / s_til)) + 1 / eta^2)) )
}

## gamma - rho[b, g]: dim B * G matrix
update_rho <- function(varphi, lambda_til, omega, r_til, s_til){
  temp1_mat <- matrix(0, nrow = B - 1, ncol = G)
  for (b in 2:B) {
    a_b <- r_til[b, ] / s_til[b, ]     # length G
    phi_b <- varphi[[b]]               # n_b x K
    phi_omega <- phi_b %*% t(omega)
    resid_b <- sweep(Y[[b]] - phi_omega, 2, lambda_til, FUN = "-")
    temp1_mat[b - 1, ] <- colSums(resid_b) * a_b
  }
  temp2_mat <- 1 / ((r_til[-1, ] / s_til[-1, ]) * n_vec[-1] + 1 / delta^2)
  
  return(rbind(rep(0, G), temp1_mat * temp2_mat))
}

## gamma - delta_til[b, g]: dim B * G matrix
update_delta_til <- function(r_til, s_til){
  return( rbind(rep(0, G), sqrt(1 / (n_vec[-1] * r_til[-1, ] / s_til[-1, ] + 1/delta^2))) )
}

## sigma^2 - r_til[b, g]: dim B * G matrix
update_r_til <- function(){
  return( r + matrix(rep(n_vec/2, times=G), nrow = B, ncol = G) ) 
}

## sigma^2 - s_til[b, g]: dim B * G matrix
update_s_til <- function(varphi, lambda_til, omega, rho, eta_til, tau_til, delta_til){
  omega2_tau2 <- omega^2 + tau_til^2
  temp <- matrix(0, nrow = B, ncol = G)
  for (b in 1:B) {
    phi_b <- varphi[[b]]   # n_b x K
    Y_centered <- sweep(Y[[b]], 2, lambda_til + rho[b, ], FUN = "-")
    phi_omega <- phi_b %*% t(omega)
    phi_quad <- phi_b %*% t(omega2_tau2)
    contrib_b <- Y_centered^2 - 2 * Y_centered * phi_omega + phi_quad
    contrib_b <- sweep(contrib_b, 2, eta_til^2 + delta_til[b, ]^2, FUN = "+")
    temp[b, ] <- s + 0.5 * colSums(contrib_b)
  }
  return(temp)
}

## xi - omega[g, k]: dim G * K matrix
update_omega <- function(varphi, lambda_til, rho, r_til, s_til, p_til, q_til){
  temp1_mat <- matrix(0, nrow = G, ncol = K - 1)
  temp2_denom <- matrix(0, nrow = G, ncol = K - 1)
  for (b in 1:B) {
    a_b <- r_til[b, ] / s_til[b, ]   # length G
    Y_centered <- sweep(Y[[b]], 2, lambda_til + rho[b, ], FUN = "-")   # n_b x G
    phi_b <- varphi[[b]][, 2:K, drop = FALSE]   # n_b x (K-1)
    temp1_mat <- temp1_mat + sweep(t(Y_centered) %*% phi_b, 1, a_b, FUN = "*")
    temp2_denom <- temp2_denom + outer(a_b, colSums(phi_b))
  }
  temp2_mat <- 1 / (temp2_denom + p_til / q_til)
  
  return(cbind(rep(0, times = G), temp1_mat * temp2_mat))
}

## xi - tau_til[g, k]: dim G * K matrix
update_tau_til <- function(varphi, r_til, s_til, p_til, q_til){
  temp_denom <- matrix(0, nrow = G, ncol = K - 1)
  for (b in 1:B) {
    a_b <- r_til[b, ] / s_til[b, ]
    phi_b <- varphi[[b]][, 2:K, drop = FALSE]
    temp_denom <- temp_denom + outer(a_b, colSums(phi_b))
  }
  temp <- 1 / (temp_denom + p_til / q_til)
  return(cbind(rep(0, times = G), sqrt(temp)))
}

## tau^2 - p_til: scalar
update_p_til <- function(){
  return( p + G * (K - 1) / 2 ) 
}

## tau^2 - q_til: scalar
update_q_til <- function(omega, tau_til){
  return( q + sum(omega^2 + tau_til^2) / 2 )
}


################# Calculate the ELBO Function ##################### 

ELBO <- function(varphi, c_til, lambda_til, eta_til, rho, delta_til, r_til, s_til, omega, tau_til, p_til, q_til){
  
  term_Z <- 0
  for (b in 1:B) {
    term_Z <- term_Z - sum(ifelse(varphi[[b]] == 0, 0, varphi[[b]] * log(varphi[[b]])))
  }
  
  term_pi <- sum(lgamma(c_til))
  
  eta2 <- eta_til^2
  delta2 <- delta_til[-1, , drop = FALSE]^2
  tau2 <- tau_til[, -1, drop = FALSE]^2
  
  term_alpha <- 0.5 * sum(log(eta2) - ((lambda_til - lambda)^2 + eta2) / eta^2)
  term_gamma <- 0.5 * sum(log(delta2) - (rho[-1, , drop = FALSE]^2 + delta2) / delta^2)
  term_xi <- 0.5 * sum(log(tau2))
  term_sigma <- -sum(r_til * log(s_til))
  term_tau <- -p_til * log(q_til)
  
  return(term_Z + term_pi + term_alpha + term_gamma + term_xi + term_sigma + term_tau)
}


################# BIC Function ##################### 

CALCULATE_BIC <- function(K){
  
  ## Initiation for parameters
  varphi_0 <- list()
  Z_0 <- list()
  for(b in 1:B){
    Z_0[[b]] <- sample(1:K, size = n_vec[b], replace = TRUE)
    
    # one-hot
    varphi_0[[b]] <- matrix(0, nrow = n_vec[b], ncol = K)
    for(i in 1:n_vec[b]){
      varphi_0[[b]][i, Z_0[[b]][i]] <- 1
    }
  }
  
  ## omega[g,1]==0, tau_til[g,1]==0, rho[1,g]==0, delta_til[1,g]==0
  c_til_0 <- matrix(c, nrow = B, ncol = K)
  lambda_til_0 <- rep(lambda, times = G)
  eta_til_0 <- rep(eta, times = G)
  omega_0 <- matrix(0, nrow = G, ncol = K)
  tau_til_0 <- matrix(c(rep(0, times=G), rep(sqrt(q/(p-1)), times=G*(K-1))), nrow = G, ncol = K)
  rho_0 <- matrix(0, nrow = B, ncol = G)
  delta_til_0 <- matrix(c(rep(0, times=G), rep(delta, times=G*(B-1))), nrow = B, ncol = G, byrow=T) 
  r_til_0 <- matrix(r, nrow = B, ncol = G)
  s_til_0 <- matrix(s, nrow = B, ncol = G)
  p_til_0 <- p
  q_til_0 <- q
  
  
  varphi_t = varphi_0
  c_til_t = c_til_0
  lambda_til_t = lambda_til_0
  eta_til_t = eta_til_0
  omega_t = omega_0
  tau_til_t = tau_til_0
  rho_t = rho_0
  delta_til_t = delta_til_0
  r_til_t = r_til_0
  s_til_t = s_til_0 
  p_til_t = p_til_0
  q_til_t = q_til_0
  elbo_value = numeric(0)
  
  
  tol = 1e-5
  cat("######## THE CAVI STARTS for K =", K, "########", "\n")
  t = 1
  while(TRUE){
    # pi
    c_til_t = update_c_til(varphi_t)
    # sigma^2
    r_til_t = update_r_til()
    s_til_t = update_s_til(varphi_t, lambda_til_t, omega_t, rho_t, eta_til_t, tau_til_t, delta_til_t)
    
    # calculate the ELBO
    elbo_value[t] = ELBO(varphi_t, c_til_t, lambda_til_t, eta_til_t, rho_t, delta_til_t, r_til_t, s_til_t, omega_t, tau_til_t, p_til_t, q_til_t)
    if(t>=2){
      if(abs(elbo_value[t] / elbo_value[t-1] - 1) < tol){
        cat("######## THE CAVI IS COMPLETED WHEN t = ", t, "########", "\n")
        break
      }
    }
    
    # xi
    omega_t = update_omega(varphi_t, lambda_til_t, rho_t, r_til_t, s_til_t, p_til_t, q_til_t)
    tau_til_t = update_tau_til(varphi_t, r_til_t, s_til_t, p_til_t, q_til_t)
    # gamma
    rho_t = update_rho(varphi_t, lambda_til_t, omega_t, r_til_t, s_til_t)
    delta_til_t = update_delta_til(r_til_t, s_til_t)
    # tau^2
    p_til_t = update_p_til()
    q_til_t = update_q_til(omega_t, tau_til_t)
    # alpha
    lambda_til_t = update_lambda_til(varphi_t, omega_t, rho_t, r_til_t, s_til_t)
    eta_til_t = update_eta_til(r_til_t, s_til_t)
    # Z
    varphi_t = update_varphi(c_til_t, lambda_til_t, omega_t, tau_til_t, rho_t, r_til_t, s_til_t)
    
    t = t + 1
  }
  
  ## point estimation for parameters
  alpha_hat <- lambda_til_t
  xi_hat <- omega_t
  gamma_hat <- rho_t
  sigma_hat <- sqrt(s_til_t/(r_til_t-1))
  pi_hat <- c_til_t/rowSums(c_til_t)
  
  ## calculate the BIC
  BIC_value = -2 * sum(sapply(1:B, function(b){sum(sapply(1:n_vec[b], function(i){
    temp = sapply(1:K, function(k){
      log(pi_hat[b, k])+sum(-log(2*pi)/2-log(sigma_hat[b, ])-(Y[[b]][i, ]-alpha_hat-xi_hat[, k]-gamma_hat[b, ])^2/(2*sigma_hat[b, ]^2))
    })
    log(sum(exp(temp - max(temp))))+max(temp)
  }))})) + (K*G+(2*B-1)*G)*log(sum(n_vec)*G)
  
  return(BIC_value)
}



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
B <- length(Y)
G <- ncol(Y[[1]])
n_vec <- sapply(Y, nrow)
 

## Set seed for BUSVI in simulation case II
seed = 2004
set.seed(seed)


## BIC analysis from K = 2 to 10
BIC_rec = c(NA)
for(K in 2 : 10){
  set.seed(seed)
  BIC_rec[K] = CALCULATE_BIC(K)
}


# Save BIC values for K = 2 to 10
# save(BIC_rec, file = "../result_data/BIC_rec.RData")
# BIC_rec = c(6255717, 4523408, 5565690, 5055219, 4673600, 5836180, 4973165, 6128242, 6277377)


## Save BIC plot
png(file="../figures/Figure2(c).png", width = 500, height = 400)
options(scipen = 999)
par(mar = c(4, 4.5, 2, 2), family = "Arial", font.lab = 1, font.axis = 1)  
BIC_plot = plot(x = 2:10, y = BIC_rec[-1], pch = 16, cex = 0.8, cex.axis = 24 / par("ps"), ylab = "BIC", cex.lab = 24 / par("ps"), xlab = "")
lines(x = 2:10, y = BIC_rec[-1])
dev.off()


