# packages

library(survival)
library(Design)
library(MASS)

rextrval <- function(x) log( -log(1-x) )


### data generating function ###

# Weibull outcome, three informative predictor variables
# sigma = scale parameter of Weibull distribution
# betac = parameter vector
# n = sample size
# covar = Kovarianz der informativen Kovariablen
# lambda = parameter for exponential (censoring) distribution
# PHviol = TRUE/FALSE, Verletzung der PH-Annahme ja/nein?
# inform = TRUE/FALSE informative Censoring ja/nein?
#    falls ja: lambdaI = Anpassungsparameter
# pnoninf = Anzahl noise-Variablen
# covnoninf = Kovarianz der Noise-Variablen
# nb = Blockgröße der korrelierten Noise-Variablen

dataGenFun <- function(n, sigma, betac, covar=0.5, lambda, PHviol=FALSE, inform=FALSE,
                       pnoninf, covnoninf = 0, nb = 1, lambdaI = NULL){
                       
    #set.seed(seed)
    # Erzeugen der Daten
    u <- runif(n)
    # Noise-Variable fuer Weibull-Modell
    w <- rextrval(u)

    # informative Kovariablen
    pinf <- length(betac)
    Sigma <- matrix(covar, nrow = pinf, ncol = pinf)
    diag(Sigma) <- 1
    X <- mvrnorm(n, mu = rep(0, pinf), Sigma = Sigma)
    
    # Namen der informativen Kovariablen
    nameinf <- rep("",pinf)
    for(jj in 1:pinf){
      nameinf[jj] <- paste("x",jj,sep = "")}
    colnames(X) <- nameinf
    
    # falls PH-Annahme verletzt:
    # verwende verschiedene (zufällige) sigma-Werte fuer jede Beobachtung
    if(PHviol) sigma <- sigma * exp(X[,1])

    # Blockweise Korrelation von Xnoise
    Sigma <- matrix(covnoninf, nrow = nb, ncol = nb)
    diag(Sigma) <- 1
    Xnoise <- mvrnorm(n, mu = rep(0, nb), Sigma = Sigma)
    for (kk in 1:(pnoninf/nb-1))
		Xnoise <- cbind(Xnoise, mvrnorm(n, mu = rep(0, nb), Sigma = Sigma))

    # wahre survival-Zeiten
    survtime <- exp(X%*%betac + sigma*w)

    # Zensierungs-Zeiten
    # falls random censoring:
    # unabhängige exponentialverteilte Zensierungszeiten
    # falls informatives censoring:
    # exponentialverteilte Zensierungszeiten,
    # die mit den survival-Zeiten korreliert sind.
    # lambdaI steuert die Korrelation
    if (inform==TRUE)
        censtime <- rexp(n, rate=lambdaI*mean(survtime)/survtime) 
    else
        censtime <- rexp(n, rate=lambda)
    
    # Zensierungs-Indikator
    status <- survtime<censtime
    # beobachtete Survival-Zeiten
    time <- pmin(survtime, censtime)

    DF <- data.frame(time, status, X, Xnoise, survtime, censtime)
	DF$time <- pmin(DF$time, 30)
    return(DF)
}






