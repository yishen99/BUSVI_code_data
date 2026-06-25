##############################################################
################## System preparation ########################
##############################################################
rm(list=ls())

### Install basic packages, if necessary
if (!require("aricode", quietly = TRUE)) 
  install.packages("aricode", version="1.0.3")
if (!require("circlize", quietly = TRUE))
  install.packages("circlize", version="0.4.16")
if (!require("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2", version="3.5.1")
if (!require("invgamma", quietly = TRUE)) 
  install.packages("invgamma", version="1.1")
if (!require("matrixStats", quietly = TRUE)) 
  install.packages("matrixStats", version="1.5.0")
if (!require("remotes", quietly = TRUE))
  install.packages("remotes", version="2.5.0")


### install the Bioconductor package manager, if necessary
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")


### BUS
if (!require("BUScorrect", quietly = TRUE))
  BiocManager::install("BUScorrect")


### Attention!
# For Windows users, please note that the version of Rtools 
# needs to be compatible with the version of R!

