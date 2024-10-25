

#VERY EFFICIENT INFERENCE OF A MULTIRESPONSE VARIABLE (assume big matrix is (samples x covariates))
#Fitting model BIGY~X for each col of BIGY
margResponseModel = function(BIGY,X,coefReturn=TRUE, pvalReturn=TRUE, sigmasqReturn=FALSE, tvalReturn=FALSE, varCoefhatReturn=FALSE, residReturn=FALSE) {
  #note does not handle NAs
  p = ncol(X) #number of covariates
  if(is.null(p)) {
    X = cbind(X)
    p = ncol(X)
  }
  n = nrow(X) #number of samples
  
  #CHeck:
  if(nrow(BIGY)!=n) stop("Number of samples in BIGY was the not same as in X!")
  
  Xt <- t(X) 
  XtX <- Xt%*%X
  condNumber = kappa(XtX) #calculate conditional number of matrix to invert
  invXtX = solve(XtX)# + lambda0) #get inverted matrix (must be possible)
  XtY <- Xt%*%BIGY #dim(XtY)=17x400k, is crossprod() faster?
  bhat = invXtX%*%XtY #get bhat for each outcomes

  RES = BIGY - X%*%bhat #get residuals 
  SE = colSums(RES^2) #get squared errors
  sigmasq = SE/(n-p) #get estimated variance (per site), adjust for number of param (unbiased)
  
  if(residReturn) {
    colnames(RES) = colnames(BIGY)
    rownames(RES) = rownames(BIGY)
  } else {
    RES = NULL #don't return residuals
  }
  
  #Obtain p-values etc:
  pvalT <- tval <- hatVarBetaHat <-  NULL
  if(pvalReturn || tvalReturn || varCoefhatReturn) {
    hatVarBetaHat = diag(invXtX)%*%t(sigmasq) # dim(hatVarBetaHat) pxI, this is only variance!
    tval <- bhat/sqrt(hatVarBetaHat) #get z/t/f-score
    pvalT <- (1-pt(abs(tval),n-p))*2 #get p-value (two-sided test...)
    #Possible Extension: get p-values from multiple testing:
    #q = sqrt(p*qf(1-signiflvl,df1=p,df2=n-p)) #get quantile of f-distribution (simultanous testing)
    #q = sqrt(qt(1-signiflvl/2,df=n-p)) #get quantile based on t-distribution (marginal testing)
  }
  
  #Prepare returning variables (insert site name)
  siteNames = colnames(BIGY)

  if(coefReturn) {
    colnames(bhat) = siteNames
    rownames(bhat) = colnames(X)
  } else {
    bhat = NULL
  }
  if(sigmasqReturn) {
    names(sigmasq) = siteNames
  } else {
    sigmasq = NULL
  }
  if(pvalReturn) {
    names(pvalT) = siteNames
  } else {
    pvalT = NULL
  }
  if(tvalReturn) {
    names(tval) = siteNames
  } else {
    tval = NULL
  }
  if(varCoefhatReturn) {
    colnames(hatVarBetaHat) = siteNames
    rownames(hatVarBetaHat) = colnames(X)
  } else {
    hatVarBetaHat = NULL
  }
  
  #calculating the (max) log-lik of the multiresponse (normal)-model using OLS estimates
  #loglikval = -0.5*sum( n*log(2*pi*sigmasq) + (n-p) ) #sum over all sites (this is not MaximumLik)
  loglikval = -n/2*sum( log(2*pi*SE/n) + 1 ) #sum over all sites (this is MaximumLik)
  return(list( loglik=loglikval, coef=bhat, sigmasq=sigmasq, pval = pvalT, tval=tval, varCoefhat=hatVarBetaHat, resid=RES, condNumber=condNumber ) )
}

#Fitting model y~BIGX for each col of BIGX (assuming designmatrix '1 + x') 
margCovariateModel = function(BIGX,y,coefReturn=TRUE, pvalReturn=TRUE, sigmasqReturn=FALSE, tvalReturn=FALSE, varCoefhatReturn=FALSE) {
  #BIGX is a pxn matrix
  #note does not handle NAs
  p = 2 #number of covariates (intercept always assumed)
  n = length(y) #number of samples
  if(nrow(BIGX)!=n) stop("Number of samples in BIGY was the not same as in X!")
  
  XtY <- crossprod(BIGX,y)
  XtX <- colSums(BIGX*BIGX)
  Sx <- colSums(BIGX)
  Sy <- sum(y)
  YtY <- sum(y^2)
  ahat <- c(XtX*Sy - Sx*XtY)/(n*XtX - Sx^2) #intercept estimates
  bhat <- c(n*XtY - Sx*Sy)/(n*XtX - Sx^2) #slope estimates
  
  sigmasq  <- c(YtY - 2*ahat*Sy - 2*bhat*XtY + n*ahat^2 + 2*ahat*bhat*Sx + bhat^2*XtX)/(n-p) #estimated variance
  hatVarBhat <- sigmasq*n/(n*XtX - Sx^2 ) #variance of slope coefficient 
  
  tval <- bhat/sqrt(hatVarBhat) #this is t-score
  pvalT <- (1-pt(abs(tval),n-p))*2 #get p-value. assume t-distr

  #Prepare returning variables (insert site name)
  siteNames = colnames(BIGX)
  
  if(coefReturn) {
    bhat = rbind(ahat,bhat)
    colnames(bhat) = siteNames
    rownames(bhat) = c("intercept","slope")
  } else {
    bhat = NULL
  }
  if(sigmasqReturn) {
    names(sigmasq) = siteNames
  } else {
    sigmasq = NULL
  }
  if(pvalReturn) {
    names(pvalT) = siteNames
  } else {
    pvalT = NULL
  }
  if(tvalReturn) {
    names(tval) = siteNames
  } else {
    tval = NULL
  }
  if(varCoefhatReturn) {
    colnames(hatVarBetaHat) = siteNames
    rownames(hatVarBetaHat) = colnames(X)
  } else {
    hatVarBetaHat = NULL
  }
  return(list( coef=bhat, sigmasq=sigmasq, pval = pvalT, tval=tval, varCoefhat=hatVarBetaHat ) )
}



