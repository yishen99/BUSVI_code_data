##############################################################
#################     Simulation Case II    ##################
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
load("../input_data/Y.RData")
load("../input_data/Z_real.RData")
B <- length(Y)
G <- ncol(Y[[1]])
K <- 3 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## set the same seed as BUSVI in simulation case II
seed = 2004
set.seed(seed)


## BUS clustering labels
Y_trs = list(t(Y[[1]]), t(Y[[2]]), t(Y[[3]]))
BUSfits = BUSgibbs(Y_trs, n.subtypes = 3, showIteration = FALSE)
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


###################### ARI Computing #########################

ARI_BUS <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_BUS[[1]]), 3), 
                      Batch2 = round(ARI(Z_real[[2]], Z_BUS[[2]]), 3), 
                      Batch3 = round(ARI(Z_real[[3]], Z_BUS[[3]]), 3),
                      Overall = round(ARI(unlist(Z_real), unlist(Z_BUS)), 3))

## Save ARI
write.csv(ARI_BUS, row.names = "BUS ARI", file = "../result_data/ARI_BUS.csv")




##############  Save Corrected data and clustering labels #################

save(Y_BUS, file = "../result_data/Y_BUS.RData")
save(Z_BUS, file = "../result_data/Z_BUS.RData")



