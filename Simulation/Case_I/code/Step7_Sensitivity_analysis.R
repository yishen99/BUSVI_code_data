##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

library(invgamma)
library(aricode)
library(ggplot2)


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
update_c_til <- function(varphi, c_value){
  return( t(sapply(1:B, function(b){c_value + colSums(varphi[[b]])})) )
}

## alpha - lambda_til[g]: dim G vector
update_lambda_til <- function(varphi, omega, rho, r_til, s_til, lambda_value, eta_value) {
  temp1_vec <- numeric(G)
  for(b in 1:B){
    a_b <- r_til[b, ] / s_til[b, ]              # length G
    phi_b <- varphi[[b]]                        # n_b x K
    phi_omega <- phi_b %*% t(omega)
    resid_b <- sweep(Y[[b]] - phi_omega, 2, rho[b, ], FUN = "-")
    temp1_vec <- temp1_vec + colSums(resid_b) * a_b
  }
  temp1_vec <- temp1_vec + lambda_value / eta_value^2
  
  temp2_vec <- 1 / (colSums(n_vec * (r_til / s_til)) + 1 / eta_value^2)
  
  return(temp1_vec * temp2_vec)
}

## alpha - eta_til[g]: dim G vector
update_eta_til <- function(r_til, s_til, eta_value) {
  return( sqrt(1 / (colSums(n_vec * (r_til / s_til)) + 1 / eta_value^2)) )
}

## gamma - rho[b, g]: dim B * G matrix
update_rho <- function(varphi, lambda_til, omega, r_til, s_til, delta_value) {
  temp1_mat <- matrix(0, nrow = B - 1, ncol = G)
  for (b in 2:B) {
    a_b <- r_til[b, ] / s_til[b, ]
    phi_b <- varphi[[b]]
    phi_omega <- phi_b %*% t(omega)
    resid_b <- sweep(Y[[b]] - phi_omega, 2, lambda_til, FUN = "-")
    temp1_mat[b - 1, ] <- colSums(resid_b) * a_b
  }
  temp2_mat <- 1 / (
    (r_til[-1, ] / s_til[-1, ]) * n_vec[-1] + 1 / delta_value^2
  )
  return(rbind(rep(0, G), temp1_mat * temp2_mat))
}

## gamma - delta_til[b, g]: dim B * G matrix
update_delta_til <- function(r_til, s_til, delta_value){
  return( rbind(rep(0, G), sqrt(1 / (n_vec[-1] * r_til[-1, ] / s_til[-1, ] + 1/delta_value^2))) )
}

## sigma^2 - r_til[b, g]: dim B * G matrix
update_r_til <- function(r_value){
  return(
    r_value + matrix(rep(n_vec/2, times=G), nrow = B, ncol = G)
  )
}

## sigma^2 - s_til[b, g]: dim B * G matrix
update_s_til <- function(varphi, lambda_til, omega, rho,
                         eta_til, tau_til, delta_til, s_value){
  omega2_tau2 <- omega^2 + tau_til^2
  temp <- matrix(0, nrow = B, ncol = G)
  for (b in 1:B) {
    phi_b <- varphi[[b]]
    Y_centered <- sweep(Y[[b]], 2, lambda_til + rho[b, ], FUN = "-")
    phi_omega <- phi_b %*% t(omega)
    phi_quad <- phi_b %*% t(omega2_tau2)
    contrib_b <- Y_centered^2 - 2 * Y_centered * phi_omega + phi_quad
    contrib_b <- sweep(contrib_b, 2, eta_til^2 + delta_til[b, ]^2, FUN = "+")
    temp[b, ] <- s_value + 0.5 * colSums(contrib_b)
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
update_p_til <- function(p_value) {
  return(p_value + G * (K - 1) / 2)
}

## tau^2 - q_til: scalar
update_q_til <- function(omega, tau_til, q_value) {
  return(q_value + sum(omega^2 + tau_til^2) / 2)
}


################# Calculate the ELBO Function ##################### 

ELBO <- function(varphi, c_til, lambda_til, eta_til, rho,
                 delta_til, r_til, s_til, omega, tau_til,
                 p_til, q_til, lambda_value, eta_value, delta_value) {
  
  term_Z <- 0
  for (b in 1:B) {
    term_Z <- term_Z - sum(
      ifelse(varphi[[b]] == 0, 0, varphi[[b]] * log(varphi[[b]]))
    )
  }
  
  term_pi <- sum(lgamma(c_til))
  
  eta2 <- eta_til^2
  delta2 <- delta_til[-1, , drop = FALSE]^2
  tau2 <- tau_til[, -1, drop = FALSE]^2
  
  term_alpha <- 0.5 * sum(
    log(eta2) -
      ((lambda_til - lambda_value)^2 + eta2) / eta_value^2
  )
  
  term_gamma <- 0.5 * sum(
    log(delta2) -
      (rho[-1, , drop = FALSE]^2 + delta2) / delta_value^2
  )
  
  term_xi <- 0.5 * sum(log(tau2))
  
  term_sigma <- -sum(r_til * log(s_til))
  
  term_tau <- -p_til * log(q_til)
  
  return(term_Z + term_pi + term_alpha + term_gamma + term_xi + term_sigma + term_tau)
}



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
load("../input_data/Z_real.RData")
B <- length(Y)
G <- ncol(Y[[1]])
K <- 3 ### From BIC analysis
n_vec <- sapply(Y, nrow)





## Function of sensitivity analysis for each hyper-parameter
hyperparameter_sensitivity_analysis <- function(hyperparam_name, value) {
  ## Hyper-parameters
  c <- 2
  lambda <- 1
  eta <- 0.1
  p <- 2
  q <- 2
  delta <- 20
  r <- 2
  s <- 2
  
  
  if (hyperparam_name == "c") {
    c <- value
  }
  else if (hyperparam_name == "lambda") {
    lambda <- value
  }
  else if (hyperparam_name == "eta") {
    eta <- value
  }
  else if (hyperparam_name == "p") {
    p <- value
  }
  else if (hyperparam_name == "q") {
    q <- value
  }
  else if (hyperparam_name == "delta") {
    delta <- value
  }
  else if (hyperparam_name == "r") {
    r <- value
  }
  else if (hyperparam_name == "s") {
    s <- value
  }
  else
    return("ERROR!")
  
  
  ## Set seed for BUSVI in simulation case I
  seed = 5
  set.seed(seed)
  
  ### Initiation for parameters
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
  
  
  ## CAVI
  tol = 1e-5
  cat("######## THE CAVI STARTS ########", "\n")
  t = 1
  while(TRUE){
    # pi
    c_til_t = update_c_til(varphi_t, c)
    # sigma^2
    r_til_t = update_r_til(r)
    s_til_t = update_s_til(varphi_t, lambda_til_t, omega_t, rho_t, eta_til_t, tau_til_t, delta_til_t, s)
    
    # calculate the ELBO
    elbo_value[t] = ELBO(
      varphi_t, c_til_t, lambda_til_t, eta_til_t, rho_t,
      delta_til_t, r_til_t, s_til_t, omega_t, tau_til_t,
      p_til_t, q_til_t,
      lambda_value = lambda,
      eta_value = eta,
      delta_value = delta
    )
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
    rho_t = update_rho(varphi_t, lambda_til_t, omega_t, r_til_t, s_til_t, delta_value = delta)
    delta_til_t = update_delta_til(r_til_t, s_til_t, delta_value = delta)
    # tau^2
    p_til_t = update_p_til(p_value = p)
    q_til_t = update_q_til(omega_t, tau_til_t, q_value = q)
    # alpha
    lambda_til_t = update_lambda_til(varphi_t, omega_t, rho_t, r_til_t, s_til_t, lambda_value = lambda, eta_value = eta)
    eta_til_t = update_eta_til(r_til_t, s_til_t, eta_value = eta)
    # Z
    varphi_t = update_varphi(c_til_t, lambda_til_t, omega_t, tau_til_t, rho_t, r_til_t, s_til_t)
    
    
    if(t %% 50 == 0) cat("######## THE CAVI IS RUNNING WHEN t = ", t, "########", "\n")
    t = t + 1
  }
  
  
  ## Estimation for parameters
  alpha_hat <- lambda_til_t
  xi_hat <- omega_t
  gamma_hat <- rho_t
  sigma_hat <- sqrt(s_til_t/(r_til_t-1))
  pi_hat <- c_til_t/rowSums(c_til_t)
  
  ## BUSVI clustering labels
  Z_BUSVI <- list()
  for(b in 1:B){
    Z_BUSVI[[b]] = rep(NA, times = n_vec[b])
    for(i in 1:n_vec[b]){
      Z_BUSVI[[b]][i] = which.max(varphi_t[[b]][i, ])
    }
  }
  
  
  ## Compute ARI
  ARI_BUSVI <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_BUSVI[[1]]), 3), 
                          Batch2 = round(ARI(Z_real[[2]], Z_BUSVI[[2]]), 3), 
                          Batch3 = round(ARI(Z_real[[3]], Z_BUSVI[[3]]), 3),
                          Overall = round(ARI(unlist(Z_real), unlist(Z_BUSVI)), 3))
  
  return(ARI_BUSVI)
  
}


hyperparam_names <- c(
  "c",
  "lambda",
  "eta",
  "p",
  "q",
  "delta",
  "r",
  "s"
)


## Set the range for varying hyper-parameters
hyperparam_value_list <- list(
  c = seq(1, 5, by = 0.4),
  lambda = seq(0, 5, by = 0.5),
  eta = seq(0.05, 0.25, by = 0.02),
  p = seq(1.1, 5.1, by = 0.4),
  q = seq(1, 5, by = 0.4),
  delta = seq(5, 55, by = 5),
  r = seq(1.1, 5.1, by = 0.4),
  s = seq(1, 5, by = 0.4)
)




## Run sensitivity analysis and create 12 data frames
for (hyperparam_name in names(hyperparam_value_list)){
  
  hyperparam_value <- hyperparam_value_list[[hyperparam_name]]
  value_ARI_list <- list()
  for(ind in 1:11){
    new_row = hyperparameter_sensitivity_analysis(hyperparam_name, hyperparam_value[ind])
    value_ARI_list[[ind]] = new_row
  }
  
  df_name <- paste0("df_", hyperparam_name, "_ARI")
  assign(
    df_name,
    cbind(hyperparameter = hyperparam_name, value = hyperparam_value, do.call(rbind, value_ARI_list))
  )
}


# df_c_ARI = data.frame(hyperparameter = "c", 
#                      value = seq(0.5, 1.5, by = 0.05),
#                      Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_eta_alpha_ARI = data.frame(hyperparameter = "eta_alpha", 
#                               value = seq(-0.2, 0.2, by = 0.02),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_tau_alpha_ARI = data.frame(hyperparameter = "tau_alpha", 
#                               value = seq(0.05, 0.15, by = 0.005),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_p_tau_ARI = data.frame(hyperparameter = "p_tau", 
#                           value = seq(2.5, 12.5, by = 0.5),
#                           Batch1 = c(0.331,0.854,0.854,0.854,0.854,rep(1,16)), 
#                           Batch2 = c(0.677,0.895,0.895,0.895,0.895,rep(1,16)), 
#                           Batch3 = c(0.538,0.839,0.839,0.839,0.839,rep(1,16)), 
#                           Overall = c(0.455,0.761,0.761,0.761,0.761,rep(1,16)))

# df_q_tau_ARI = data.frame(hyperparameter = "q_tau", 
#                           value = seq(0.5, 1.5, by = 0.05),
#                           Batch1 = c(rep(1,20),0.854), 
#                           Batch2 = c(rep(1,20),0.895), 
#                           Batch3 = c(rep(1,20),0.839), 
#                           Overall = c(rep(1,20),0.761))

# df_tau_gamma_ARI = data.frame(hyperparameter = "tau_gamma", 
#                           value = seq(0.005, 0.105, by = 0.005),
#                           Batch1 = 1, 
#                           Batch2 = c(1,1,0.981,0.993,rep(1,17)), 
#                           Batch3 = c(0.839,0.825,1,1,rep(1,17)), 
#                           Overall = c(0.733,0.735,0.989,0.996,rep(1,17)))

# df_eta_psi_ARI = data.frame(hyperparameter = "eta_psi", 
#                               value = seq(0.25, 1.25, by = 0.05),
#                               Batch1 = c(0.854,0.854,0.854,0.854,1,1,0.844,rep(1,14)), 
#                               Batch2 = c(0.895,0.895,0.895,0.895,1,1,0.895,rep(1,14)), 
#                               Batch3 = c(0.839,0.839,0.839,0.839,1,1,0.839,rep(1,14)), 
#                               Overall = c(0.761,0.761,0.761,0.761,1,1,0.758,rep(1,14)))

# df_tau_psi_ARI = data.frame(hyperparameter = "tau_psi", 
#                             value = seq(0.4, 1.4, by = 0.05),
#                             Batch1 = c(0.816,0.816,0.816,rep(1,18)), 
#                             Batch2 = c(0.895,0.895,0.895,rep(1,18)), 
#                             Batch3 = c(0.839,0.839,0.839,rep(1,18)), 
#                             Overall = c(0.781,0.781,0.781,rep(1,18)))

# df_phi_alpha_ARI = data.frame(hyperparameter = "phi_alpha", 
#                       value = seq(0.01, 0.11, by = 0.005),
#                       Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_xi_ARI = data.frame(hyperparameter = "phi_xi", 
#                               value = seq(0.05, 0.15, by = 0.005),
#                               Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_gamma_ARI = data.frame(hyperparameter = "phi_gamma", 
#                            value = seq(0.01, 0.05, by = 0.002),
#                            Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)

# df_phi_psi_ARI = data.frame(hyperparameter = "phi_psi", 
#                              value = seq(0.5, 3, by = 0.125),
#                              Batch1 = 1, Batch2 = 1, Batch3 = 1, Overall = 1)


df_SA_ARI = rbind(df_c_ARI, df_lambda_ARI, df_eta_ARI, df_p_ARI, df_q_ARI, df_delta_ARI,df_r_ARI, df_s_ARI)



## Draw ARI plots for each hyper-parameter
c_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "c"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "c", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(a).png", plot = c_SY_plot, width=10, height=5)


lambda_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "lambda"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "lambda", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(b).png", plot = lambda_SY_plot, width=10, height=5)


eta_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "eta"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "eta", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(c).png", plot = eta_SY_plot, width=10, height=5)


p_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "p"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "p", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(d).png", plot = p_SY_plot, width=10, height=5)


q_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "q"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "q", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(e).png", plot = q_SY_plot, width=10, height=5)


delta_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "delta"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "delta", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(f).png", plot = delta_SY_plot, width=10, height=5)


r_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "r"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "r", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(g).png", plot = r_SY_plot, width=10, height=5)


s_SY_plot <- ggplot(subset(df_SA_ARI, hyperparameter == "s"), aes(x = value, y = Overall)) +
  geom_line(linetype = "longdash") + 
  geom_point(size = 3) +
  theme_gray() +
  labs(x = "s", y = "ARI") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text = element_text(size = 24), axis.title = element_text(size = 24))
ggsave("../figures/FigureS1(h).png", plot = s_SY_plot, width=10, height=5)





