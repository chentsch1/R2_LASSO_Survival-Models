##############################################################################
############# Main: CCCC of R2-Measures  #####################################
##############################################################################

rm(list=ls())
wddummy <- getwd()
setwd("U:/BIOMETRIE_Heidelberg/Masterarbeit/Gütemaße_Überlebenszeit/
                                            R-Code Gütemaße/analysis")

library(fields)
library(ggplot2)
source("ccc.r")

setwd("U:/BIOMETRIE_Heidelberg/Masterarbeit/Gütemaße_Überlebenszeit/
                                        R-Code Gütemaße/sim results")

##-----------------------------------------------------##
##  Function: make_list()                              ##
##-----------------------------------------------------##

load_list <- function(files){
  R2.list <- list()
	for (i in 1:length(files)){
    file <- files[i]
    load(file)
    R2.list[[i]] <- R2
  }
  return(R2.list)
}

##-----------------------------------------------------##
## load and merge data for analysis and set factors    ##
##-----------------------------------------------------##


files <- list.files()
R2.all <- load_list(files)

# corlist <- list()
# for (i in 1:length(R2.all)){
#    corlist[[i]] <- ccc(R2.all[[i]][,3:8])}

corarray <- array(rep(0,6*6*486), dim=c(6,6,486))
for (i in 1:length(R2.all)){
    corarray[,,i] <- ccc(R2.all[[i]][,3:8])
}
dimnames(corarray)[1] <- list(colnames(R2.all[[1]][,3:8]))
dimnames(corarray)[2] <- list(colnames(R2.all[[1]][,3:8]))
dimnames(corarray)[3] <- list(files)


meanccc <- apply(X=corarray, MARGIN=c(1,2), FUN=mean)
colnames(meanccc) <- colnames(R2.all[[1]][,3:8])
rownames(meanccc) <- colnames(R2.all[[1]][,3:8])
 
#########################
setwd(wddummy)
######################### 
                                                                                                                                                           