
################################################################################
########### Function: Compute R^2 Measures #####################################
################################################################################

library(survival)
library(flexmix)
#source("koq3.r")

##-----------------------------------------------------------##
## help function: partial likelihood function for Cox-Modell ##
##-----------------------------------------------------------##

# y = Surv Model Object with 2 columns: survival time, censoring status
# f = prediction function of Cox model

plf <- function(y, f) {
        time <- y[, 1]
        event <- y[, 2]
        time2 <- sort(unique(time))
        event2 <- tapply(event, time, FUN=sum)
        f2 <- tapply(f*event, time, FUN=sum)
        n <- length(time)
        if (length(f) == 1)
            f <- rep(f, n)
        ef <- exp(f)
        risk <- rep(0, n)
        for (i in 1:n) risk[i] <- sum((time >= time[i]) * ef)
        risk2 <- tapply(risk, time, FUN=unique)
        sum(f2 - event2 * log(risk2))
    }


##----------------------------------------------------------##
## help function: caculation of R^2{L} by Cox Snell & Magee ##
##----------------------------------------------------------##

r2l <- function(object, linpred, linpred0=0){
		n <- length(object[,1])
		if((length(linpred0)==1)&&(linpred0=1)){linpred0 <- rep(0,n)}
		return(1-exp((-2/n)*(plf(object,linpred)-plf(object,linpred0))))}



##----------------------------------------------------------##
## help function: caculation of R^2{N} by Nagelkerke        ##
##----------------------------------------------------------##

r2n <- function(object, linpred, linpred0=0){
		n <- length(object[,1])
		if((length(linpred0)==1)&(linpred0=1)){linpred0 <- rep(0,n)}
		return( (1-exp((-2/n)*(plf(object,linpred)-plf(object,linpred0))))/
                                  (1-exp((2/n)*plf(object,linpred0))) )}



##-------------------------------------------------------------##
## help function: caculation of R^2{OXS} by O'Quigley Xu Stare ##
##-------------------------------------------------------------##

r2oxs <- function(object, linpred, linpred0=0){
		n <- length(object[,1])
		events <- sum(object[,2]==TRUE)
		if((length(linpred0)==1)&(linpred0=1)){linpred0 <- rep(0,n)}
		return( 1-exp((-2/events)*( plf(object,linpred)-plf(object,linpred0))) )}



##-----------------------------------------------------------##
## help function: caculation of R^2{R} by Roysten            ##
##-----------------------------------------------------------##

r2r <- function(object, linpred, linpred0=0){
		n <- length(object[,1])
		events <- sum(object[,2]==TRUE)
		if((length(linpred0)==1)&(linpred0=1)){linpred0 <- rep(0,n)}
		r2oxs <- 1-exp((-2/events)*( plf(object,linpred)-plf(object,linpred0)))
		return(r2oxs/(r2oxs+(1-r2oxs)*(pi^2)/6))}



##-----------------------------------------------------------##
## help function: caculation of R^2{XO} Xu & Q'Qigley        ##
##-----------------------------------------------------------##

# Hilfsfunktion: individual contributions to GammaHat,
# Formula (2) in O'Quigley, Xu & Stare
plfIndiv <- function(y, t, f) {
    stime <- y[,1]
    denom <- sum( (stime>=t) * exp(f))
    numer <- (stime>=t) * exp(f)
    numer / denom
    }

xo <- function(object, linpred, linpred0=0){

time <- object[,1]
n <- length(time)
if(length(linpred0==1)&&(linpred0==1)){linpred0 <- rep(0,n)}

# second factor for fixed i
secondFactor <- rep(0, n)
for (i in 1:n){
  pijbeta <- plfIndiv(y=object, t=time[i], f=linpred)
  pij0 <- plfIndiv(y=object, t=time[i], f=linpred0)
  indvec <- pij0>0
  secondFactor[i] <- sum(pijbeta[indvec] * log(pijbeta[indvec]/pij0[indvec]))
  }

# KM jumps
STE <- survfit(object~1)
# print("STE:")
# print(STE)
kmj <- c(1,STE$surv[1:(n-1)]) - STE$surv

GammaHat <- 2 * sum( kmj * secondFactor )
return( 1-exp(-GammaHat) )
}

##-----------------------------------------------------------##
##  help function: caculation of R^2{R} by Kent & O'Quigley  ##
##-----------------------------------------------------------##

r2ko <- function(object, linpred1, linpred0){ 

#--------------------------------------------------------
    ELL <- function(theta, theta1){
    
    #	Expected Log-Likelihood function for the Weibull regression model
    #
    #	Note: 	negative Log-likelihood value is returned
    #		to facilitate finding of extreme value of ELL in
    #		find.mu.alfa

      n <- length(theta)-2
      linpred1 <- (theta1[1:n])
      linpred <- (theta[1:n])
      mu1 <- (theta1[(n+1)])
      mu <- (theta[(n+1)])
      alfa1 <- (theta1[(n+2)])
      alfa <- (theta[(n+2)])
		  a <- alfa/alfa1
    
		#print(linpred1)
		#print(length(linpred1))

		b.linpred <- as.matrix(linpred - a * linpred1)
		b <- (mu - a * mu1) + b.linpred
		ga1 <- gamma(a + 1)
		
    #	negative value of ELL is returned !
    return( -(log(alfa) - 0.57721566 * a + mean(b - exp(b) * ga1)) )}

#-- koq main function ------------------------------------------

	n <- length(linpred1)
	which <- rep(F, n)
	
	theta1 <- c(linpred1, 0, 1)
	theta <- c(linpred0, 0, 1)
	
  #	Set lower and upper bounds for mu (-Inf,Inf)
  #	and alfa (0,Inf)

		lower <- c(which, T, 0) * theta
		lower[c(which, T, 0) == T] <-  - Inf
		upper <- c(which, T, T) * theta
		upper[c(which, T, T) == T] <- Inf	#

  #	find theta=c(linpred,mu,alfa) which maximize
  #	Expected Log-Likelihood given by ell
	
	 b0 <- nlminb(start = theta, objective = ELL, lower = lower, upper = upper, 
              theta1 = theta1)
	  theta0 <- b0$par
  
    GAMMA <- 2 * (ELL(theta = theta0, theta1 = theta1) - ELL(theta = theta1, theta1 = theta1))
	 return(1 - exp( - GAMMA))
}
	 

##----------------------------------------------------------##
## main function: coputation of R^2 Measures and output     ##
##----------------------------------------------------------##


r2 <- function(object, linpred, linpred0=0){
	
	if(is.list(object) && is.list(linpred)){
		if(!is.list(linpred0)){ 
			linpred00 <- list()
			h <- 0
			for (h in 1:length(linpred)){
			linpred00[[h]] <- linpred0}
		}
		else{
			linpred00 <- linpred0
		}
		R2 <- matrix(rep(0,8), ncol=8, nrow=length(linpred))
		rnames <- rep("",length(object))
		part.loglik0 <- rep(0, length(object))
		part.loglik <- rep(0,length(object))
		i <- 0
		for (i in 1:length(object)){
				if(length(linpred00[[i]])<length(linpred[[i]])){linpred00[[i]] <- rep(0,length(linpred[[i]]))}
			R2[i,1] <- plf(object[[i]], linpred00[[i]])
			R2[i,2] <- plf(object[[i]], linpred[[i]])
			R2[i,3] <- r2l(object[[i]], linpred[[i]], linpred00[[i]])
			R2[i,4] <- r2n(object[[i]], linpred[[i]], linpred00[[i]])
			R2[i,5] <- r2oxs(object[[i]], linpred[[i]], linpred00[[i]])
			R2[i,6] <- r2r(object[[i]], linpred[[i]], linpred00[[i]])
			R2[i,7] <- xo(object[[i]], linpred[[i]], linpred00[[i]])
			R2[i,8] <- r2ko(object[[i]], linpred[[i]], linpred00[[i]])
			rnames[i] <- paste("R-Square iteration",i)
		}
		colnames(R2) <- c("part.loglik0","part.loglik","R2L","R2N","R2OXS","R2R","R2XO","R2KO")
		rownames(R2) <- rnames
		meanR2 <- c(mean(R2[,1]),mean(R2[,2]),mean(R2[,3]),mean(R2[,4]),mean(R2[,5]), mean(R2[,6]),mean(R2[,7]),mean(R2[,8]))
		R2 <- rbind(R2, meanR2)
	}
	else{
		if((length(linpred0)==1) && (linpred0==1)){linpred0 <- rep(1,length(linpred))}
		R2 <- matrix(rep(0,8), ncol=8)
		R2[,1] <- plf(object, rep(1,length(linpred)))
		R2[,2] <- plf(object, linpred)
		R2[,3] <- r2l(object, linpred, linpred0)
		R2[,4] <- r2n(object, linpred, linpred0)
		R2[,5] <- r2oxs(object, linpred, linpred0)
		R2[,6] <- r2r(object, linpred,linpred0)
		R2[,7] <- xo(object, linpred, linpred0)
		R2[,8] <- r2ko(object, linpred, linpred0)
	
		colnames(R2) <- c("part.loglik0","part.loglik","R2L","R2N","R2OXS","R2R","R2XO","R2KO")
		rownames(R2) <- "R-Square"}

	return(R2)
}

