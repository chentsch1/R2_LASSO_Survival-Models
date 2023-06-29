##############################################################################
################## Kent O'Quigley R2 #########################################
##############################################################################

koq <- function(linpred1, linpred0, np, verbose = F)
{
#
#----------------------------------------------------------------
	ELL <- function(theta, theta1)
	{
#	Expected Log-Likelihood function for the Weibull regression model
#	(see reference 1 in help file)
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
		if(verbose) {
    print("a in ELL:")
		print(a)
		print("alfa:")
		print(alfa)
		print("mu:")
		print(mu)
		print("länge linpred1 in ELL:")
		#print(linpred1)
		#print(length(linpred1))
		}
		b.linpred <- as.matrix(linpred - a * linpred1)
		b <- (mu - a * mu1) + b.linpred
		ga1 <- gamma(a + 1)
		
#	negative value of ELL is returned !
  return( -(log(alfa) - 0.57721566 * a + mean(b - exp(b) * ga1)) )
  }

#
#-- koq function ---------------------------------------------------------
#
	n <- length(linpred1)
	which <- rep(F, n)
	
	if(verbose) {
		cat("linpred  = ", linpred, "\n")
		cat("linpred1 = ", linpred1, "\n")
	}
	
	theta1 <- c(linpred1, 0, 1)
	theta <- c(linpred0, 0, 1)
	
#	find theta=c(beta,mu,alfa) which maximize
#	Expected Log-Likelihood given by ell
#
#	Set lower and upper bounds for mu (-Inf,Inf)
#	and alfa (0,Inf)

		lower <- c(which, T, 0) * theta
		lower[c(which, T, 0) == T] <-  - Inf
		upper <- c(which, T, T) * theta
		upper[c(which, T, T) == T] <- Inf	#
	
	 b0 <- nlminb(start = theta, objective = ELL, lower = lower, upper = upper, 
              theta1 = theta1)
	  theta0 <- b0$par
  if(verbose) {
  print("b0:")
  print(b0)
	print("theta0:")
	print(theta0)
  }	
  GAMMA <- 2 * (ELL(theta = theta0, theta1 = theta1) - ELL(theta = theta1, 
                theta1 = theta1))
	KOQ <- (1 - exp( - GAMMA))
	if(verbose) {
  print("KOQ:")
	print(KOQ)
	}
  return(KOQ)
}
