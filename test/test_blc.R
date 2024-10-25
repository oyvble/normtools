
#testing internal function RMPP::blc 

#Generate example
  n=10000
  sh1 = c(1.53,1.99,14.7)
  sh2 = c(41.6,2.54,1.08)
  Y <- NULL
  for(k in 1:length(sh1)) Y <- c(Y,rbeta(n,sh1[k],sh2[k]))# +  rbeta(n,sh1[2],sh2[2]) + rbeta(n,sh1[3],sh2[3])
  Y <- Y[sample(1:length(Y),n)]
  thresh = c(0.2,0.75) #start value thresholds
  W <- matrix(0,nrow=n,ncol=3);
  W[which(Y <= thresh[1]),1] <- 1;
  W[intersect(which(Y > thresh[1]),which(Y <= thresh[2])),2] <- 1;
  W[which(Y > thresh[2]),3] <- 1;
  #Y=as.matrix(Y);w=W
  fit1 = RPMM::blc(as.matrix(Y),W,maxiter=100,tol=1e-3)
  fit2 = blc1(Y,W,maxiter=100,tol=1e-3,steptol=1e-4)
  c(fit1$llike,fit2$llike)
  cbind(fit1$a,fit2$a)
  cbind(fit1$b,fit2$b)

  microbenchmark( #ABOUT SAME
    RPMM::blc(as.matrix(Y),W,maxiter=100,tol=1e-3),
    normtools::blc1(Y,W,maxiter=100,tol=1e-3,steptol=1e-6), 
    #blc1(Y,W,maxiter=100,tol=1e-3,steptol=1e-3), 
    times= 10
  )



if(0) {
  #COMPARING SPEED:
  library(microbenchmark)
  sum(dbeta(y,ab[1],ab[2],log=TRUE))
  
  ab = c(15.20430,1.189681)
  lbeta2 = function(a,b) lgamma(a)+lgamma(b)-lgamma(a+b)
  microbenchmark( #ABOUT SAME
    lbeta(ab[1],ab[2]),
    lbeta2(ab[1],ab[2]),
    times= 100000
  )

}

