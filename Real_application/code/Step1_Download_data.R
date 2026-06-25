##############################################################
################    breast cancer data    ####################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

########################## Note ##############################
## Please first install required R packages.
## R: System_preparation.R

## Please download the dataset to the "input_data" folder. 
## integrativeGEdata_microarray.RData
## https://tandf.figshare.com/ndownloader/articles/6826961/versions/1



## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


load("../input_data/integrativeGEdata_microarray.RData")






