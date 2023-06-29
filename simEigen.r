################################################################################
############# Main: LASSO-Fitting and Results R^2-Measures  ####################
################################################################################

rm(list=ls())
setwd("/net/192.168.245.150/export/cluster-home/zu_meinen_Dateien")
library(survival)
library(penalized)
source("datasplit.r")
source("R2.r")
source("datagenfun.r")
source("Lambda_set.R")


##-----------------------------------------------------##
## Function for complete simulaton call                ##
##-----------------------------------------------------##

#sim(n=75, sigma=0.5, lambda=1e-10, PHviol=TRUE, inform=TRUE, pnoninf=10, lambdaI=1.126, iter=5)


sim <- function(sim.par)
{
	n <- sim.par[1]
	sigma <- sim.par[2]
	lambda <- sim.par[3]
	PHviol <- as.logical(sim.par[4])
	inform <- as.logical(sim.par[5])
	pnoninf <- sim.par[6]
	iter <- sim.par[7]
	inform_var <- as.logical(sim.par[8])
	covnoninf <- 0.5
	nb <- 10
	betac <- c(-0.25,0.25,0.5)
	
	
	setwd("/net/192.168.245.150/export/cluster-home/zu_meinen_Dateien")
	R2 <- matrix(rep(0,8*iter), nrow=iter, ncol=8)
							
	coeflist <- vector("list",length=iter)
	survobject <- vector("list",length=iter)
	linpred <- vector("list",length=iter)
	cens <- rep(0,iter)
	i <- 0
	
	for (i in 1:iter){
		seed <- (10*i)
		set.seed(seed)
		#evt <- matrix(rep(0,12), ncol=3, nrow=4)
		
		##--------------------------------------------------##
		##  generate data with 2*n Observations             ##
		##--------------------------------------------------##
		datagen <- dataGenFun(n=n, sigma=sigma, betac=betac, lambda=lambda,
						  PHviol=PHviol, inform=inform, pnoninf=pnoninf, lambdaI=lambda)


		##--------------------------------------------------##
		##  splitting dataset (2/3 to 1/3)                  ##
		##--------------------------------------------------##

		
		nstatus0 <- (length(datagen$status)-sum(datagen$status))
		cens[i] <- nstatus0/nrow(datagen)
		names(cens[i]) <- paste("censoring iteraton ",i)
		datasets <- datasplit(datagen, 0.33, seed)
		
		##----------------------------------------------------##
		## LASSO fiting and computing prediction on test-data ##
		##----------------------------------------------------##
		if(inform_var){
			evt3 <- optfit <- optL1(Surv(datasets@train$time, datasets@train$status),
								penalized=datasets@train[,3:(ncol(datasets@train)-2)],
								data=datasets@train, standardize=TRUE, fold=10, trace=FALSE)
			survobject[[i]] <- Surv(datasets@test$time, datasets@test$status)
			linpred[[i]] <- predict(optfit$fullfit, penalized=datasets@test[,3:(ncol(datasets@test)-2)],
								data=datasets@test, lperrg=TRUE)
		}else{
			evt3 <- optfit <- optL1(Surv(datasets@train$time, datasets@train$status),
									penalized=datasets@train[,4:(ncol(datasets@train)-2)],
									data=datasets@train, standardize=TRUE, fold=10, trace=FALSE)
			survobject[[i]] <- Surv(datasets@test$time, datasets@test$status)
			linpred[[i]] <- predict(optfit$fullfit, penalized=datasets@test[,4:(ncol(datasets@test)-2)],
									data=datasets@test, lperrg=TRUE)
		}
	
		##----------------------------------------------------##
		## coputing R^2 measures and fitted beta-coefficients ##
		##----------------------------------------------------##
		
		coeflist[[i]] <- coefficients(optfit$fullfit)
	}   
	# end loop
	
	R2 <- r2(survobject, linpred)
	##----------------------------------------------------##
	## save R2 and beta-coefficients                      ##
	##----------------------------------------------------##
	
    meancens <- mean(cens)
    save(R2, coeflist, cens,  file = paste("n", n, "sigma", sigma, "noninf", lambda,
		pnoninf, "cens", round(meancens,4), "viol", PHviol, "inform", inform,
		"inform_var", inform_var, ".rda", sep = "_") )
}



###### input ###################################################################
#betac <- c(-0.25,0.25,0.5)  # constant true beta vector of covariates
#iter  <- 100    # Anzahl der Gesamtdurchläufe  -> 100
################################################################################

##-----------------------------------------------------------##
## Call for simulaton Function with extrem parameters        ##
##-----------------------------------------------------------##

#memory.limit(3000)
#
#

#sim(n=75, sigma=0.5, lambda=1e-10, PHviol=TRUE, inform=TRUE, pnoninf=10, lambdaI=1.126, iter=5)




###############################################################################


Sim.Param <- expand.grid(n=c(300,450,600), sigma=c(0.3,0.5,0.8), 
				lambda=c(0.3,0.7,1.25), PHviol=c(TRUE),
				inform=c(FALSE,TRUE), pnoninf=c(100,1000,5000),
				iter=100, inform_var= c(FALSE,TRUE))



Sim.Param <- data.frame(apply(Sim.Param,1,Lambda.set))


library(snow)
cl <- makePVMcluster(40)

clusterEvalQ(cl, eval(parse(text="library(survival);library(penalized); 
                      library(splines);library(MASS);library(flexmix)")))
clusterExport(cl,c("dataGenFun","datasplit","plf","plfIndiv","r2", 
                  "r2ko","r2l","r2n","r2oxs","r2r","rextrval","xo"))
ERG <- clusterApplyLB(cl,Sim.Param,sim)







