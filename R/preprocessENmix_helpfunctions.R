#THIS FUNCTION IS ALREADY VECTORIZED
new_cm  <- function(lambda,mu2,sigma2,p1,p2,mu,sigma,s){
  #s is a long vectorr
  a <- mu+lambda*sigma^2
  C <- (sigma2^2*(s-mu)+sigma^2*mu2)/(sigma2^2+sigma^2)
  D <- sigma2*sigma/sqrt(sigma2^2+sigma^2)
  
  temp_1 <- p1*lambda*exp(lambda^2*sigma^2/2-lambda*(s-mu))*(pnorm(s,a,sigma)-pnorm(0,a,sigma))
  temp_2 <- p2*dnorm(s,mu2+mu,sqrt(sigma2^2+sigma^2))/(1-pnorm(0,mu2,sigma2))
  
  num_1 <- temp_1
  num_2 <- temp_2*(pnorm(s,C,D)-pnorm(0,C,D))
  
  denom_1 <- temp_1*(s-a+sigma*(dnorm((s-a)/sigma)-dnorm(a/sigma))/(pnorm(s,a,sigma)-pnorm(0,a,sigma)))
  denom_2 <- temp_2*(C*(pnorm(s,C,D)-pnorm(0,C,D))+D*(exp(-C^2/(2*D^2))-exp(-(s-C)^2/(2*D^2)))/sqrt(2*pi))
  
  (denom_1+denom_2)/(num_1+num_2)
}


firstpeak <- function(x,y,sn,dat)
{
  n <- length(y)
  ##sn to avoid using small peak position
  nn <- sn*2+1
  v <- matrix(NA,ncol=nn,nrow=n-nn+1)
  for(i in 1:nn){v[,i]=y[i:(n-nn+i)]}
  ix <- sn+which(apply(v<v[,(sn+1)],1,sum) == (nn-1))
  if(length(ix)>0){mu=x[ix[1]]}
  if(length(ix) == 0 | sum(dat<mu)/length(dat)>=0.15)
  {
    dat <- dat[dat<quantile(dat,p=0.10)]
    temp <- density(dat)
    flag <- temp$x>=min(dat) & temp$x<=max(dat)
    temp$x <- temp$x[flag];temp$y=temp$y[flag]
    mu <- temp$x[which.max(temp$y)]
  }
  mu
}



#Wrapper function for calling core function implemented in C
enmix_adjC <- function(meth_i=NULL,bg_i=NULL,bgParaEst) {
  #meth_i is a long integer vector (typically 10k's)
  #bg_i is a length=2 integer vector 
  if(sum(is.na(meth_i))>0) {stop("ENmix background correction does not allow missing value")}
  meth_i[meth_i <= 0]=1 #FAST
  mu <- bg_i$mu[1]
  sigma <- bg_i$sigma[1]
  
  if(bgParaEst == "est" | bgParaEst == "neg" | bgParaEst == "oob") {
    x <- (meth_i[meth_i>=mu]-mu) #FAST
    #      temp <- ENmix:::EM_estimate(x) #SLOW (1x)
    #      temp <- EM_estimateR(x) #FASTER (3x)
    temp <- EM_estimateC(x) #FASTEST (6.5x)
    lambda <- temp[[1]]
    mu2 <- temp[[2]]
    sigma2<-ifelse(temp[[3]] <= sigma, 0.1, sqrt(temp[[3]]^2-sigma^2))
    p1 <- (sum(meth_i<mu)+temp[[4]]*length(x))/length(meth_i)
    p2 <- 1-p1
    meth_adj <- new_cm(lambda,mu2,sigma2,p1,p2,mu,sigma,meth_i) #THIS IS FAST!
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values, only a few
  }else if (bgParaEst == "subtract_neg" | bgParaEst == "subtract_estBG" | bgParaEst == "subtract_q5neg" | bgParaEst == "subtract_oob"){
    meth_adj=meth_i-mu
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values
  }
  meth_adj
}

#Wrapper function for calling core function implemented in R
enmix_adjR <- function(meth_i=NULL,bg_i=NULL,bgParaEst) {
  #meth_i is a long integer vector (typically 10k's)
  #bg_i is a length=2 integer vector 
  if(sum(is.na(meth_i))>0) {stop("ENmix background correction does not allow missing value")}
  meth_i[meth_i <= 0]=1 #FAST
  mu <- bg_i$mu[1]
  sigma <- bg_i$sigma[1]
  
  if(bgParaEst == "est" | bgParaEst == "neg" | bgParaEst == "oob") {
    x <- (meth_i[meth_i>=mu]-mu) #FAST
    #      temp <- ENmix:::EM_estimate(x) #SLOW (1x)
    temp <- EM_estimateR(x) #FASTER (3x)
    #temp <- EM_estimateC(x) #FASTEST (6.5x)
    lambda <- temp[[1]]
    mu2 <- temp[[2]]
    sigma2<-ifelse(temp[[3]] <= sigma, 0.1, sqrt(temp[[3]]^2-sigma^2))
    p1 <- (sum(meth_i<mu)+temp[[4]]*length(x))/length(meth_i)
    p2 <- 1-p1
    meth_adj <- new_cm(lambda,mu2,sigma2,p1,p2,mu,sigma,meth_i) #THIS IS FAST!
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values, only a few
  }else if (bgParaEst == "subtract_neg" | bgParaEst == "subtract_estBG" | bgParaEst == "subtract_q5neg" | bgParaEst == "subtract_oob"){
    meth_adj=meth_i-mu
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values
  }
  meth_adj
}


#copy from MASS package
huber<-function (y, k = 1.5, tol = 1e-06)
{
  y <- y[!is.na(y)]
  n <- length(y)
  mu <- median(y)
  s <- mad(y)
  if (s == 0)
    stop("cannot estimate scale: MAD is zero for this sample")
  repeat {
    yy <- pmin(pmax(mu - k * s, y), mu + k * s)
    mu1 <- sum(yy)/n
    if (abs(mu - mu1) < tol * s)
      break
    mu <- mu1
  }
  list(mu = mu, s = s)
}

huber_mus <- function(x){ests <- try(huber(x)); if(class(ests)[1]=="try-error"){
  cat("Warning:Check negtive control data, or do quality control before ENmix\n");
  c(mu=median(x,na.rm=TRUE),s=sd(x,na.rm=TRUE))
}else{c(mu=ests$mu,s=ests$s)}}
huber_mu <- function(x){ests <- try(huber(x));if(class(ests)[1]=="try-error"){
  cat("Warning: Check NORM control data, or do quality control before ENmix\n");
  median(x,na.rm=TRUE)
}else{ests$mu}}

estBG  <- function(meth_i) {
  meth_i[meth_i<=0]=1e-06
  temp <- density(meth_i)
  temp <- density(meth_i[meth_i<temp$x[which.max(temp$y)]])
  flag <- temp$x>=min(meth_i) & temp$x<=max(meth_i)
  temp$x <- temp$x[flag];temp$y=temp$y[flag]
  mu <- temp$x[which.max(temp$y)]
  ##first mode
  if((sum(meth_i<mu)/length(meth_i))>=0.15){
    mu=firstpeak(temp$x,temp$y,sn=5,meth_i)
  }
  perc <- sum(meth_i<mu)/length(meth_i)
  sigma <- sqrt(sum((meth_i[meth_i<mu]-mu)^2)/sum(meth_i<mu))
  c(mu,sigma,perc)
}
