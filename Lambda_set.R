
#Nicht informatives Censoring (-> falls inform=FALSE):
# bei sigma = 0.3
#     bei lambda = 0.365: 30% censoring
#     bei lambda = 0.770: 50% censoring
#     bei lambda = 1.525: 70% censoring
# bei sigma = 0.5
#     bei lambda = 0.385: 30% censoring
#     bei lambda = 0.845: 50% censoring
#     bei lambda = 1.74: 70% censoring
# bei sigma = 0.8
#     bei lambda = 0.4: 30% censoring
#     bei lambda = 0.935: 50% censoring
#     bei lambda = 2.15: 70% censoring

#bei obigen Settings kannst Du lambdaI mit irgendwas belegen, spielt keine Rolle.

##########

#Informatives Censoring (-> falls inform=TRUE):

# bei sigma = 0.3
#     bei lambdaI = 0.33: 30% censoring
#     bei lambdaI = 0.64: 50% censoring
#     bei lambdaI = 1.115: 70% censoring
# bei sigma = 0.5
# bei lambdaI = 0.332: 30% censoring
# bei lambdaI = 0.65: 50% censoring
# bei lambdaI = 1.126: 70% censoring
# bei sigma = 0.8
# bei lambdaI = 0.317: 30% censoring
# bei lambdaI = 0.620: 50% censoring
# bei lambdaI = 1.075: 70% censoring


Lambda.set  <-  function(x){
	if(!x[5]){
		if(x[2] == 0.3 && x[3] == 0.3){
			x[3] <- 0.365
		}
		if(x[2] == 0.3 && x[3] == 0.7){
			x[3] <- 0.77
		}
		if(x[2] == 0.3 && x[3] == 1.25){
			x[3] <- 1.525
		}
		##############################
		#
		##############################
		if(x[2] == 0.5 && x[3] == 0.3){
			x[3] <- 0.385
		}
		if(x[2] == 0.5 && x[3] == 0.7){
			x[3] <- 0.845
		}
		if(x[2] == 0.5 && x[3] == 1.25){
			x[3] <- 1.74
		}
		##############################
		#
		##############################
		if(x[2] == 0.8 && x[3] == 0.3){
			x[3] <- 0.4
		}
		if(x[2] == 0.8 && x[3] == 0.7){
			x[3] <- 0.935
		}
		if(x[2] == 0.8 && x[3] == 1.25){
			x[3] <- 2.15
		}
	}
	else{
		if(x[2] == 0.3 && x[3] == 0.3){
			x[3] <- 0.33
		}
		if(x[2] == 0.3 && x[3] == 0.7){
			x[3] <- 0.64
		}
		if(x[2] == 0.3 && x[3] == 1.25){
			x[3] <- 1.115
		}
		##############################
		#
		##############################
		if(x[2] == 0.5 && x[3] == 0.3){
			x[3] <- 0.332
		}
		if(x[2] == 0.5 && x[3] == 0.7){
			x[3] <- 0.65
		}
		if(x[2] == 0.5 && x[3] == 1.25){
			x[3] <- 1.126
		}
		##############################
		#
		##############################
		if(x[2] == 0.8 && x[3] == 0.3){
			x[3] <- 0.317
		}
		if(x[2] == 0.8 && x[3] == 0.7){
			x[3] <- 0.620
		}
		if(x[2] == 0.8 && x[3] == 1.25){
			x[3] <- 1.075
		}
	}
	x
}

