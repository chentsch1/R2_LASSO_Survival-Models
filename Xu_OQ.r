library(survival)

##############################################################################
################## measure by Xu & O'Quigley (R^2_{XO}) ######################
##############################################################################

# Hilfsfunktion: individual contributions to GammaHat,
# Formula (2) in O'Quigley, Xu & Stare
plfIndiv <- function(y, t, f) {
    stime <- y[,1]
    denom <- sum( (stime>=t) * exp(f))
    numer <- (stime>=t) * exp(f)
    numer / denom
    }

# object = Surv-Objekt (survival_time, event_status) der Testdaten
# linpred = linearer Praediktor
# linpred0 = linearer Praediktor aus Null-Modell

xo <- function(object, linpred, linpred0=1){

time <- object[,1]
n <- length(time)
if(length(linpred0==1)&(linpred0==1)){linpred0 <- rep(1,n)}

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
kmj <- c(1,STE$surv[1:(n-1)]) - STE$surv

GammaHat <- 2 * sum( kmj * secondFactor )
return( 1-exp(-GammaHat) )

}


############################ Beginn Beispiel #################################

## Beispiel, ohne penalized package
ovarianTR <- ovarian[1:16,]
ovarianTE <- ovarian[17:26,]
model0 <- coxph(Surv(futime, fustat)~1, data=ovarianTR)
model1 <- coxph(Surv(futime, fustat)~age, data=ovarianTR)
f0 <- predict(model0)[1:nrow(ovarianTE)]
f1 <- predict(model1, newdata=ovarianTE)

xo(Surv(ovarianTE$futime, ovarianTE$fustat), f1, f0)

############################# Ende Beispiel ##################################