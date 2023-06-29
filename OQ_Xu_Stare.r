library(survival)

##############################################################################
################## measure by O'Qugley, Xu & Stare (R^2_{OXS}) ###############
##############################################################################

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

##-----------------------------------------------------------##
## main function: caculation of R^2{OXS} by O'Quigley Xu Stare ##
##-----------------------------------------------------------##


r2oxs <- function(object, linpred, linpred0=1){
		n <- length(object[,1])
		events <- sum(object[,2]==TRUE)
		if((length(linpred0)==1)&(linpred0=1)){linpred0 <- rep(1,n)}
		return(1-(exp((-2/events)*( plf(object,linpred)-plf(object,linpred0)))))}


