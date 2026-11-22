##############################################################
################    breast cancer data    ####################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################


library(BUScorrect)
library(aricode)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())
t1 = Sys.time()


## Read data Y
load("../input_data/integrativeGEdata_microarray.RData")
Y <- Data_list
Y <- lapply(Y, t)
B <- length(Y)
G <- ncol(Y[[1]])
K <- 4 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## set seed the same of BUSVI 
seed = 11
set.seed(seed)



## BUS clustering labels
Y_trs = list(t(Y[[1]]), t(Y[[2]]), t(Y[[3]]))
BUSfits = BUSgibbs(Y_trs, n.subtypes = 4, showIteration = FALSE)
Z_BUS = BUSfits $ Subtypes


###################### Batch Correction #########################

alpha_BUS = BUSfits $ alpha # G
mu_BUS = BUSfits $ mu # G * K
gamma_BUS = BUSfits $ gamma # G * B
sigma_BUS = BUSfits $ sigma_sq # G * B

Y_BUS <- list()
Y_BUS[[1]] = Y[[1]]
for(b in 2:B){
  Y_BUS[[b]] <- matrix(NA, nrow = n_vec[b], ncol = G)
  for(i in 1:n_vec[b]){
    Y_BUS[[b]][i, ] = alpha_BUS + mu_BUS[, Z_BUS[[b]][i]] + (Y[[b]][i, ]-alpha_BUS-mu_BUS[, Z_BUS[[b]][i]]-gamma_BUS[, b]) / (sigma_BUS[, b]/sigma_BUS[, 1])
  }
}


###################### Timing #########################

t2 = Sys.time()
Running_Time = round(difftime( t2, t1, units = "secs"), 3)
cat(paste0("  The BUS procedure takes: ", Running_Time, " seconds", "\n"))



## Save Corrected data and clustering label
save(Y_BUS, file = "../result_data/Y_BUS.RData")
save(Z_BUS, file = "../result_data/Z_BUS.RData")



