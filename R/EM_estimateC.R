#THIS IS A BOTTLENECK OF preprocessENmix

#library(fasterNormTools); print(EM_estimateC( x=rnorm(100,2000) ))
#start=c(max(density(x)$y),mean(range(x)),diff(range(x))/6,0.5);epsilon=c(0.0001,0.001,0.001,0.001)
EM_estimateC  <- function(x,start=c(max(density(x)$y),mean(range(x)),diff(range(x))/6,0.5),epsilon=c(0.0001,0.001,0.001,0.001)){
  # start=c(max(density(x)$y),mean(range(x)),diff(range(x))/6,0.5);epsilon=c(0.0001,0.001,0.001,0.001)
  lambda <- start[1]
  mu <- start[2]
  sigmaSq <- start[3]^2
  p1 <- start[4]
  n <- length(x)
  obj = .C("EM_estLoop", x, lambda, mu, sigmaSq, p1 ,n, epsilon, package="normtools") #call external C++ function 
  retVec = unlist(obj[2:5]) #obtain params
  #list(lambda,mu,sqrt(sigma.sq),p1)
  list(retVec[1],retVec[2],sqrt(retVec[3]),retVec[4])
}