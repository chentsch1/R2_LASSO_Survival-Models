##############################################################################
############# Sim-Function: LASSO-Fitting and Results R^2-Measures  ##########
##############################################################################

# rm(list=ls())
library(survival)
library(penalized)
source("datasplit.r")
source("R2.r")
source("datagenfun.r")

##-----------------------------------------------------##
## Function for complete simulaton call                ##
##-----------------------------------------------------##

sim <- function(n, sigma, betac, lambda=1e-10, PHviol=FALSE, inform=FALSE, 
                        pnoninf=10, covnoninf=0, nb=1, lambdaI=NULL, iter=1){
   
  R2 <- matrix(rep(0,7*iter), nrow=iter, ncol=7)
  coeflist <- list()
  datasets <- list()
  survobject <- list()
  linpred <- list() 
  cens <- rep(0,iter)
  evtime <- list()
  ind <- rep(0, iter)
  i <- 0

  for (i in 1:iter){
    seed <- (10*i)
    set.seed(seed)
    #evt <- matrix(rep(0,12), ncol=3, nrow=4)

  ##--------------------------------------------------##
  ##  generate data with 2*n Observations             ##
  ##--------------------------------------------------##
    
   datagen <- dataGenFun(n=n, sigma=sigma, betac=betac, lambda=lambda, 
              PHviol=PHviol, inform=inform, pnoninf=pnoninf, lambdaI=lambdaI)
    	
  ##--------------------------------------------------##
  ##  splitting dataset (2/3 to 1/3)                  ##
  ##--------------------------------------------------##

    nstatus0 <- (length(datagen$status)-sum(datagen$status))
    cens[i] <- nstatus0/nrow(datagen)
    names(cens[i]) <- paste("censoring iteraton ",i)
    datasets[[i]] <- datasplit(datagen, 0.33, seed)
    
  ##----------------------------------------------------##
  ## LASSO fiting and computing prediction on test-data ##
  ##----------------------------------------------------##

    optfit <- optL1(Surv(datasets[[i]]@train$time, datasets[[i]]@train$status), 
              penalized=datasets[[i]]@train[,3:(ncol(datasets[[i]]@train)-2)], 
              data=datasets[[i]]@train, standardize=T) 
 	 survobject[[i]] <- Surv(datasets[[i]]@test$time, datasets[[i]]@test$status)
 	 system.time(linpred[[i]] <- predict(optfit$fullfit,                   
                penalized=datasets[[i]]@test[,3:(ncol(datasets[[i]]@test)-2)], 
				                 data=datasets[[i]]@test, lperrg=TRUE))

  ##----------------------------------------------------##
  ## coputing R^2 measures and fitted beta-coefficients ##
  ##----------------------------------------------------##
    
    coeflist[[i]] <- coefficients(optfit$fullfit)
  
 }   # end loop 
        
    R2 <- r2(survobject, linpred) 
    
  ##----------------------------------------------------##
  ## save R2 and beta-coefficients                      ##
  ##----------------------------------------------------##
    
    meancens <- mean(cens)
    save(R2, coeflist, cens, datasets, evtime,  file = paste("n", n, "sigma", 
         sigma, "noninf", pnoninf, "cens", meancens, "viol", PHviol, "inform", 
          inform, ".rda", sep = "_") )

}


###### input #################################################################
betac <- c(-0.25,0.25,0.5)  # constant true beta vector of covariates
iter  <- 100    # Anzahl der Gesamtdurchläufe  -> 1
##############################################################################



#### example call ###

sim(n=600, sigma=0.5, betac=betac, lambda=0.365, PHviol=FALSE, inform=FALSE, 
                 pnoninf=5000, covnoninf=0.5, nb=10, lambdaI=NULL, iter=iter)
