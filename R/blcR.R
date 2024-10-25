

blc1 <- function(Y, W, maxiter=25, tol=1E-6,steptol=1e-6){
  K <- 3 #number of (latent) groups
  n <- length(Y) #number of samples
  
  #if(is.null(weights)) weights <- rep(1,n)
  mu <- a <- b <- rep(Inf,K)
  crit <- Inf
  for(i in 1:maxiter){
    eta <- colSums(W)/n 
    mu0 <- mu
    suppressWarnings({
      for(k in 1:K){
        #FUNCTION CALL
        ab <- betaEst1(y=Y, w=W[,k],steptol=steptol) 
        a[k] <- ab[1]
        b[k] <- ab[2] 
        mu[k] <- ab[1]/sum(ab)
      }
    })
    
    for(k in 1:K) W[,k] <- dbeta(Y, a[k], b[k], log=TRUE)  #Update weights
    #W <- apply(WW, c(1,2), sum, na.rm=TRUE) #(n x 3)
    
    #Obtain largest value over all classes:
    Wmax <- W[,1] #apply(W,1,max) #n long vector (max for each class)
    indLarger = W[,2]>Wmax
    Wmax[indLarger] = W[indLarger,2]
    indLarger = W[,3]>Wmax
    Wmax[indLarger] = W[indLarger,3]
    
    for(k in 1:K) W[,k] <- W[,k] - Wmax #avoiding underflow
    W <- t(eta * t(exp(W)))
    like <- rowSums(W) #apply(W,1,sum)
    W <- W/like 
    llike <- log(like) + Wmax  #UPDATING LIKELIHOOD FUNCTION
    
    crit <- max(abs(mu-mu0))
    if(crit<tol) break
  }
  return(list(a=a, b=b, eta=eta, mu=mu, w=W, llike=sum(llike)))
}

#y=Y; w=w[,k],
betaEst1 <- function(y,w,steptol=1e-6){
  #precalculations:
  N <- sum(w)
  p <- sum(w*y)/N
  v <- sum(w*y*y)/N - p*p
  
  #calculate start value:
  logab0 <- log(c(p, 1-p)) + log(pmax(1E-6,p*(1-p)/v - 1))  #obtain start values of params
  
  #Stores sufficient vars
  sumWlogY1 = sum(w*log(y))
  sumWlogY2 = sum(w*log(1-y))
  
  #beta_lpdf= (a-1)log(x) + (b-1)log(1-x) - lbeta(a,b), where lbeta(a,b)=lgam(a)+lgam(b)-lgam(a+b)
  #Minimizing the neg-loglik of beta-regression
  neg_logLik  <- function(logab){
    ab <- exp(logab)
    #-sum(w*dbeta(y,ab[1],ab[2],log=TRUE)) 
    -(ab[1]-1)*sumWlogY1 - (ab[2]-1)*sumWlogY2 + N*lbeta(ab[1],ab[2])
  }
  #FUNCTION CALL
  #opt <- try(optim(logab, betaObjf1, ydata=y, wdata=w, method="BFGS"),silent=TRUE)
  opt <- try(nlm(neg_logLik,logab0,steptol = steptol),silent=TRUE)
  
  if(inherits(opt,"try-error")) return(c(1,1))
  exp(opt$est)
  #exp(opt$param)
}