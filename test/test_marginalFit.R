
#testing effective implementation of looped lm:
library(normtools)
#Generating example dataset:
set.seed(1)
p = 1000 #number of sites
n = 100 #number of observations
BIGY = matrix(rnorm(p*n),nrow=n) #Samples given per row (sites per column)
x = rnorm(n) + 0.9*BIGY[,sample(1:p,1)] #create a response variable (with true correlation to 1 site)

#helpfunction to check if two numbers are equal
numequal = function(x,y) all.equal(as.numeric(x),as.numeric(y))

#TESTING P-values of slope param
#Test model of type Y_i~(1+x) for i=1,...,p
ret1 = margResponseModel(BIGY,cbind(1,x))
pvals1 = rep(NA,p) #manual calculation of pvalues
for(i in 1:p) pvals1[i] = coef(summary(lm(BIGY[,i]~x)))[2,4]
numequal(ret1$pval[2,],pvals1)
#plot(-log10(pvals1))

#Test model of type x~(1+Y_i) for i=1,...,p
ret2 = margCovariateModel(BIGY,x)
pvals2 = rep(NA,p) #manual calculation of pvalues
for(i in 1:p) pvals2[i] = coef(summary(lm(x~BIGY[,i])))[2,4]
numequal(ret2$pval,pvals2)
#plot(-log10(pvals2))


#Other tests (sigmasq,stderr,tval,residuals)
ret = margResponseModel(BIGY,cbind(1,x),T,T,T,T,T,T) #returning all variables
i = 1 #site to check
lmfit = summary(lm(BIGY[,i]~x))

#coefficients
numequal(ret$coef[,i], coef(lmfit)[,1])

#Sigmasq
numequal(ret$sigmasq[i], lmfit$sigma^2)

#Std.err^2 (variance of coefhat)
numequal(ret$varCoefhat[,i], coef(lmfit)[,2]^2)

#t-values
numequal(ret$tval[,i], coef(lmfit)[,3])

#Residuals
numequal(ret$resid[,i], lmfit$residuals)

#CALC FULL EXAMPLE
if(0) {
  set.seed(1)
  p = 400000 #number of sites
  n = 1000 #number of observations
  BIGY = matrix(rnorm(p*n),nrow=n) #Samples given per row (sites per column)
  x = rnorm(n) + 0.8*BIGY[,sample(1:p,1)] #create a response variable (with true correlation to 1 site)
  
  system.time({ #FAST WAY TO OBTAIN P-values across all sites
    pvals = margResponseModel(BIGY,cbind(1,x),coefReturn = FALSE)$pval[2,] #obtain p-value for slope param
  })[3] #6s
  pvals[pvals==0] = 1e-20 #avoid -inf
# summary(-log10(pvals))
}

