##############################################################################
############# Main: LASSO-Fitting and Results R^2-Measures  ##################
##############################################################################

rm(list=ls())
wddummy <- getwd()
setwd("U:/BIOMETRIE_Heidelberg/Masterarbeit/Gütemaße_Überlebenszeit/R-Code 
      Gütemaße/sim results")

library(fields)
library(ggplot2)


##-----------------------------------------------------##
##  Function: set_par()                                ##
##-----------------------------------------------------##

set_par <- function(x){
	N <- if(grepl("n_300",x)){ 300 }else{ if(grepl("n_450",x)){ 450 }
            else{ if(grepl("n_600",x)) {600}}}
	sigma <- if(grepl("sigma_0.3",x)){ 0.3 }else{ if(grepl("sigma_0.5",x)){0.5}
            else{ if(grepl("sigma_0.8",x)) {0.8}}}
  pnoninf <- if(grepl("5000",x)){ 5000 }else{ if(grepl("1000",x)){ 1000 }
            else{ if(grepl("100",x)) {100}}}
  inform <- if(grepl("inform_TRUE",x)){ TRUE }
            else{ if(grepl("inform_FALSE",x)){ FALSE }}
  viol <- if(grepl("viol_TRUE",x)){ TRUE }
            else{ if(grepl("viol_FALSE",x)){ FALSE }}
  misspred <- if(grepl("misspred_TRUE",x)){ TRUE } 
              else{ if(grepl("misspred_FALSE",x)){ FALSE }}
  
  Erg <- data.frame(N,sigma,pnoninf,inform,viol,misspred)
  names(Erg) <- c("N","sigma","pnoninf","inform","viol","misspred")
	
	load(x)
	R2 <- R2[-100,-(1:2)]
	rownames(R2)<-NULL
	cens <- round(cens,3)
  Erg <- data.frame(Erg,cens,R2)
  Erg
}

##-----------------------------------------------------##
## load and merge data for analysis and set factors    ##
##-----------------------------------------------------##


files <- list.files()  #[c(45,47,49)]

ERG <- lapply(files,set_par)
ERG_1 <- NULL
for(i in 1:length(ERG)){ 
	ERG_1 <- rbind(ERG_1,ERG[[i]])}

ERG_2 <- melt(ERG_1, id=c("N", "sigma", "pnoninf", "inform", "viol", "cens", 
              "misspred"), na.rm=TRUE) 
names(ERG_2)

ERG_2$cens[ERG_2$cens < 0.4] <- 0.3
ERG_2$cens[(0.4 < ERG_2$cens) & (ERG_2$cens < 0.6)] <- 0.5 
ERG_2$cens[ERG_2$cens > 0.6] <- 0.7
ergrow <- row(ERG_2)[,1]
ERG_2 <- cbind(ergrow, ERG_2)
rows <- ERG_2$ergrow[ERG_2$cens == 0.6]
ifelse(((ERG_2$cens[rows-1] < 0.6) & (ERG_2$cens[rows+1] < 0.6)), 
              (ERG_2$cens[rows] <- 0.5), (ERG_2$cens[rows] <- 0.7))
ifelse(((ERG_2$cens[rows-1] < 0.4) & (ERG_2$cens[rows+1] < 0.4)), 
              (ERG_2$cens[rows] <- 0.3), (ERG_2$cens[rows] <- 0.5))

ERG_2$variable <- factor(ERG_2$variable, levels = c("R2L", "R2N", "R2OXS", 
                                                    "R2R", "R2XO", "R2KO"))
ERG_2$N <- factor(ERG_2$N, levels = c("600","450","300"), 
                  labels=c("N=600", "N=450", "N=300"))
ERG_2$sigma <- factor(ERG_2$sigma, levels = c("0.3", "0.5", "0.8"), 
                  labels = c("sigma=0.3", "sigma=0.5", "sigma=0.8"))
ERG_2$pnoninf <- factor(ERG_2$pnoninf, levels = c("100", "1000", "5000"), 
                  labels = c("pnoninf=100", "pnoninf=1000", "pnoninf=5000"))
ERG_2$inform <- factor(ERG_2$inform, levels = c("FALSE", "TRUE"), 
                  labels = c("noninformative", "informative"))
ERG_2$viol <- factor(ERG_2$viol, levels = c("FALSE", "TRUE"), 
                  labels = c("pHviol=FALSE", "pHviol=TRUE"))
ERG_2$cens <- factor(ERG_2$cens, levels = c(0.3, 0.5, 0.7),
                  labels = c("30%", "50%", "70%"))
ERG_2$misspred <- factor(ERG_2$misspred, levels = c("TRUE", "FALSE"), 
                  labels = c("no missing predictor", "missing predictor"))


ERG_sub1 <- subset(x=ERG_2, subset=(ERG_2$inform=="noninformative" & 
                    ERG_2$viol=="pHviol=FALSE" & ERG_2$pnoninf=="pnoninf=100" 
                                        & ERG_2$sigma=="sigma=0.5"))
ERG_sub2 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$viol=="pHviol=FALSE" & ERG_2$pnoninf=="pnoninf=100" 
                  & ERG_2$sigma=="sigma=0.5"))
ERG_sub3 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$inform=="noninformative" 
                  & ERG_2$pnoninf=="pnoninf=100" 
                  & ERG_2$sigma=="sigma=0.5"))

ERG_sub4 <- subset(x=ERG_2, subset=(ERG_2$inform=="noninformative" 
                  & ERG_2$viol=="pHviol=FALSE" & ERG_2$N=="N=600" 
                  & ERG_2$sigma=="sigma=0.5"))
ERG_sub5 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$viol=="pHviol=FALSE" & ERG_2$N=="N=600" 
                  & ERG_2$sigma=="sigma=0.5"))
ERG_sub6 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$inform=="noninformative" & ERG_2$N=="N=600" 
                  & ERG_2$sigma=="sigma=0.5"))

ERG_sub7 <- subset(x=ERG_2, subset=(ERG_2$inform=="noninformative" 
                  & ERG_2$viol=="pHviol=FALSE" & ERG_2$pnoninf=="pnoninf=100" 
                  & ERG_2$N=="N=600"))
ERG_sub8 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$viol=="pHviol=FALSE" & ERG_2$pnoninf=="pnoninf=100" 
                  & ERG_2$N=="N=600"))
ERG_sub9 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$inform=="noninformative" 
                  & ERG_2$pnoninf=="pnoninf=100" & ERG_2$N=="N=600"))

ERG_sub10 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$inform=="noninformative" & ERG_2$viol=="pHviol=FALSE" 
                  & ERG_2$pnoninf=="pnoninf=100" & ERG_2$N=="N=600"))
ERG_sub11 <- subset(x=ERG_2, subset=(ERG_2$misspred=="no missing predictor" 
                  & ERG_2$inform=="noninformative" & ERG_2$viol=="pHviol=FALSE" 
                  & ERG_2$pnoninf=="pnoninf=100" & ERG_2$N=="N=600" 
                  & ERG_2$sigma=="sigma=0.5"))

                              
##-----------------------------------------------------##
## Draw Boxplots                                       ##
##-----------------------------------------------------##

######### Boxplot 1: R2 by cens, mispred & N ###########

attach(ERG_sub1)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366,
                            0.430763,0.430763,0.430763), 
                            misspred = c("no missing predictor",
                              "no missing predictor","no missing predictor",
                              "missing predictor","missing predictor",
                              "missing predictor"), 
                            N = c("N=600","N=450","N=300","N=600",
                                                          "N=450","N=300"))  

# Draw boxplot
p <- ggplot(ERG_sub1,aes(factor(cens),value), facets=(N ~misspred))
p + geom_boxplot(aes(fill=variable)) + facet_grid(N ~misspred) + 
                              ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), 
                              hline.data) + 
opts(title = expression("R2-Measures by Censoring, Missing 
                                          Predictor Covariates & N"), 
            labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
                            noninformative cens, p noninf=100, sigma=0.5", 
                            y = "R2", colour = "Measure")) 

detach(ERG_sub1)


######### Boxplot 2: R2 by cens, inform & N ###########

attach(ERG_sub2)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

devAskNewPage(ask = TRUE)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366,
                                  0.5085366,0.5085366,0.5085366), 
                            inform = c("noninformative", "noninformative",
                                        "noninformative", "informative", 
                                        "informative", "informative"), 
                            N = c("N=600","N=450","N=300",
                                    "N=600","N=450","N=300"))  

# Draw boxplot
p <- ggplot(ERG_sub2,aes(factor(cens),value), facets=(N ~inform))
p + geom_boxplot(aes(fill=variable)) + facet_grid(N ~inform) + ylim(c(-0.2,1)) 
                              + geom_hline(aes(yintercept = z), hline.data) + 
opts(title = expression("R2-Measures by Censoring (Degree / Informative) & N"), 
            labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
            no missing predictor, p noninf=100, sigma=0.5", 
                                y = "R2", colour = "Measure"))

detach(ERG_sub2)


######### Boxplot 3: R2 by cens, viol & N ###########

attach(ERG_sub3)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

devAskNewPage(ask = TRUE)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366), 
                            viol = c("pHviol=FALSE","pHviol=FALSE",
                                                      "pHviol=FALSE"), 
                            N = c("N=600","N=450","N=300"))  

# Draw boxplot
p <- ggplot(ERG_sub3,aes(factor(cens),value), facets=(N ~viol))
p + geom_boxplot(aes(fill=variable)) + facet_grid(N ~viol) + ylim(c(-0.2,1)) 
                            + geom_hline(aes(yintercept = z), hline.data) + 
opts(title = expression("R2-Measures by Censoring, pH Violation & N"), 
    labels = c(x = "censoring \n \n all R2-Measures at no missing predictor, 
          noninformative cens, p noninf=100, sigma=0.5", 
                            y = "R2", colour = "Measure")) 

detach(ERG_sub3)


####### Boxplot 4: R2 by cens, mispred & pnoninf ###########

attach(ERG_sub4)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

devAskNewPage(ask = TRUE)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366,
                                      0.430763,0.430763,0.430763), 
                            misspred = c("no missing predictor",
                                    "no missing predictor",
                                    "no missing predictor",
                                    "missing predictor",
                                    "missing predictor",
                                    "missing predictor"), 
                            pnoninf = c("pnoninf=100", "pnoninf=1000", 
                            "pnoninf=5000", "pnoninf=100", 
                            "pnoninf=1000", "pnoninf=5000"))  

# Draw boxplot
p <- ggplot(ERG_sub4,aes(factor(cens),value), facets=(pnoninf ~misspred))
p + geom_boxplot(aes(fill=variable)) + facet_grid(pnoninf ~misspred) 
            + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) + 
opts(title = expression("R2-Measures by Censoring, 
      Missing Predictor Covariates & p of Noninformative Covariates"), 
            labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
            noninformative cens, N=600, sigma=0.5", 
                    y = "R2", colour = "Measure")) 

detach(ERG_sub4)


######### Boxplot 5: R2 by cens, inform & pnoninf ###########

attach(ERG_sub5)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

devAskNewPage(ask = TRUE)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366,
                                    0.5085366,0.5085366,0.5085366), 
                            inform = c("noninformative", "noninformative",
                            "noninformative", "informative", "informative", 
                            "informative"), 
                            pnoninf = c("pnoninf=100", "pnoninf=1000", 
                            "pnoninf=5000", "pnoninf=100", "pnoninf=1000", 
                            "pnoninf=5000"))     

# Draw boxplot
p <- ggplot(ERG_sub5,aes(factor(cens),value), facets=(inform ~pnoninf))
p + geom_boxplot(aes(fill=variable)) + facet_grid(pnoninf ~inform) 
          + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) + 
opts(title = expression("R2-Measures by Censoring (Degree / Informative) 
                              & p of Noninformative Covariates"), 
            labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
            no missing predictor, N=600, sigma=0.5", 
                          y = "R2", colour = "Measure")) 

detach(ERG_sub5)


######### Boxplot 6: R2 by cens, viol & pnoninf ###########

attach(ERG_sub6)

# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)

devAskNewPage(ask = TRUE)

# Set Line Values 
hline.data <- data.frame(z = c(0.5085366,0.5085366,0.5085366), 
                            viol = c("pHviol=FALSE","pHviol=FALSE",
                                                    "pHviol=FALSE"), 
                            pnoninf = c("pnoninf=100", "pnoninf=1000", 
                                                        "pnoninf=5000"))

# Draw boxplot
p <- ggplot(ERG_sub6,aes(factor(cens),value), facets=(pnoninf ~viol))
p + geom_boxplot(aes(fill=variable)) + facet_grid(pnoninf ~viol) 
            + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) + 
opts(title = expression("R2-Measures by Censoring, pH violation 
                                    & p of Noninformative Covariates"), 
      labels = c(x = "censoring \n \n all R2-Measures at no missing predictor, 
      noninformative cens, N=600, sigma=0.5", y = "R2", colour = "Measure"))

detach(ERG_sub6)

####### Boxplot 7: R2 by cens, mispred & sigma ###########                                                                                                
                                                                                                                                                            
attach(ERG_sub7)                                                                                                                                       
# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)                                                                                    
devAskNewPage(ask = TRUE)                                                                                                                                   
                                                                                                                                                            
# Set Line Values                                                                                                                                           
x <- 0.11
hline.data <- data.frame(z = c(0.729463,0.5085366,0.3016534,
                                0.6121161,0.430763,0.2585898),                                                                                     
                            misspred = c("no missing predictor",
                            "no missing predictor","no missing predictor",
                            "missing predictor","missing predictor",
                            "missing predictor"),
                            sigma = c("sigma=0.3", "sigma=0.5", "sigma=0.8",
                            "sigma=0.3", "sigma=0.5", "sigma=0.8"))                      
                                                                                                                                                            
# Draw boxplot                                                                                                                                              
p <- ggplot(ERG_sub7,aes(factor(cens),value), facets=(sigma ~misspred))                                                                                   
p + geom_boxplot(aes(fill=variable)) + facet_grid(sigma ~misspred) 
          + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) +                      
      opts(title = expression("R2-Measures by Censoring, 
                          Missing Predictor Covariates & Sigma"),                                      
            labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
            noninformative cens, N=600, p noninf=100",
                            y = "R2", colour = "Measure"))       
                                                                                                                                                            
detach(ERG_sub7)                                                                                                                                            
                                                                                                                                                            
                                                                                                                                                            
######### Boxplot 8: R2 by cens, inform & sigma ###########                                                                                               
                                                                                                                                                            
attach(ERG_sub8)                                                                                                                                            
                                                                                                                                                            
# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)                                                                                                   
                                                                                                                                                            
devAskNewPage(ask = TRUE)                                                                                                                                   
                                                                                                                                                            
# Set Line Values                                                                                                                                           
hline.data <- data.frame(z = c(0.729463,0.5085366,0.3016534,
                                  0.729463,0.5085366,0.3016534),                                                                      
                            inform = c("noninformative", "noninformative",
                            "noninformative", "informative", 
                            "informative", "informative"),                   
                            sigma = c("sigma=0.3", "sigma=0.5", "sigma=0.8",
                            "sigma=0.3", "sigma=0.5", "sigma=0.8"))                      
                                                                                                                                                            
# Draw boxplot                                                                                                                                              
p <- ggplot(ERG_sub8,aes(factor(cens),value), facets=(sigma ~inform))                                                                                     
p + geom_boxplot(aes(fill=variable)) + facet_grid(sigma ~inform) 
            + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) +                        
opts(title = expression("R2-Measures by Censoring 
                              (Degree / Informative) & Sigma"),                                                                 
labels = c(x = "censoring \n \n all R2-Measures at pHviol=FALSE, 
no missing predictor, N=600, p noninf=100", y = "R2", colour = "Measure"))     
                                                                                                                                                            
detach(ERG_sub8)                                                                                                                                           
                                                                                                                                                                                                                                                       
                                                                                                                                                            
######### Boxplot 9: R2 by cens, viol & sigma ###########                                                                                                 
                                                                                                                                                            
attach(ERG_sub9)                                                                                                                                            
                                                                                                                                                            
# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)                                                                                                   
                                                                                                                                                            
devAskNewPage(ask = TRUE)                                                                                                                                   
                                                                                                                                                            
# Set Line Values                                                                                                                                           
hline.data <- data.frame(z = c(0.729463,0.5085366,0.3016534),                                                                      
                            viol = c("pHviol=FALSE","pHviol=FALSE",
                            "pHviol=FALSE"),                               
                            sigma = c("sigma=0.3", "sigma=0.5", "sigma=0.8"))                      
                                                                                                                                                            
# Draw boxplot                                                                                                                                              
p <- ggplot(ERG_sub9,aes(factor(cens),value), facets=(sigma ~viol))                                                                                       
p + geom_boxplot(aes(fill=variable)) + facet_grid(sigma ~viol) 
          + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) +                          
opts(title = expression("R2-Measures by Censoring, pH violation & Sigma"),                                                                                
labels = c(x = "censoring \n \n all R2-Measures at no missing predictor, 
      noninformative cens, N=600, p noninf=100", y = "R2", colour = "Measure")) 
                                                                                                                                                            
detach(ERG_sub9)                                                                                                                                            


######### Boxplot 10: R2 reference plot N=600, pnoninf=100 ###########                                                                                                 
                                                                                                                                                            
ERG_sub10$one <- rep(1, nrow(ERG_sub10))

attach(ERG_sub10)                                                                                                                                            
                                                                                                                                                            
# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)                                                                                                   
                                                                                                                                                            
devAskNewPage(ask = TRUE)                                                                                                                                   
                                                                                                                                                            
# Set Line Values                                                                                                                                           
hline.data <- data.frame(z = c(0.729463,0.5085366,0.3016534),                                                                      
                    viol = c("pHviol=FALSE","pHviol=FALSE","pHviol=FALSE"),                               
                          sigma = c("sigma=0.3", "sigma=0.5", "sigma=0.8"))                      
                                                                                                                                                            
# Draw boxplot                                                                                                                                              

p <- ggplot(ERG_sub10,aes(factor(cens),value), facets=(sigma ~one))                                                                                       
p + geom_boxplot(aes(fill=variable)) + facet_grid(sigma ~one) 
              + ylim(c(-0.2,1)) + geom_hline(aes(yintercept = z), hline.data) +                          
opts(title = expression("R2-Measures by Censoring & Sigma"),                                                                                
labels = c(x = "censoring \n \n all R2-Measures at no missing predictor, 
      noninformative cens, N=600, p noninf=100", y = "R2", colour = "Measure")) 
                                                                                                                                                            
detach(ERG_sub10)                                                                                                                                            

######### Boxplot 11: R2 reference plot at sigma=0.5, N=600, pnoninf=100 #####                                                                                                 
                                                                                                                                                            

attach(ERG_sub11)                                                                                                                                            
                                                                                                                                                            
# table(ERG_sub2$cens, ERG_sub2$sigma, ERG_sub2$misspred)                                                                                                   
                                                                                                                                                            
devAskNewPage(ask = TRUE)                                                                                                                                   
                                                                                                                                                            
# Set Line Values                                                                                                                                           
hline.data <- data.frame(z = c(0.5085366),                                                                      
                            viol = c("pHviol=FALSE"),                               
                            sigma = c("sigma=0.5"))                      
                                                                                                                                                            
# Draw boxplot                                                                                                                                              

p <- ggplot(ERG_sub11,aes(factor(cens),value))                                                                                       
p + geom_boxplot(aes(fill=variable)) + ylim(c(-0.2,1)) 
                      + geom_hline(aes(yintercept = z), hline.data) +                          
opts(title = expression("R2-Measures by Censoring"),                                                                                
      labels = c(x = "censoring \n \n all R2-Measures at no missing predictor, 
      noninformative cens, sigma=0.5, N=600, p noninf=100", 
                                y = "R2", colour = "Measure")) 
                                                                                                                                                            
detach(ERG_sub11)                                                                                                                                            


#########################
setwd(wddummy)
#########################                                                                                                                                                            