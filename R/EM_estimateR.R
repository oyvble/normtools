
#THIS IS A BOTTLENECK OF preprocessENmix
EM_estimateR  <- function(x,start=c(max(density(x)$y),mean(range(x)),diff(range(x))/6,0.5),epsilon=c(0.0001,0.001,0.001,0.001)){
  lambda <- start[1]
  mu <- start[2]
  sigma.sq <- start[3]^2
  p1 <- start[4]
  n <- length(x)
  
  #Precalculations:
  x2 = x^2 #precalculate this
  log2pi = log(2*pi) #constant
  sumX = sum(x)
  #sumX2 = sum(x2) #for sigmasq_new
  
  diff <- TRUE
  while (diff) {  #is typically looped 100's of times
    #loggedVals0 = dnorm(x,mu,sqrt(sigma.sq),log=TRUE)-dexp(x,lambda,log=TRUE)
    const = 0.5*(mu^2/sigma.sq + log(sigma.sq) + log2pi) + log(lambda) #constant to update
    loggedVals = -0.5*x2/sigma.sq + x*(mu/sigma.sq + lambda) - const
    #    all.equal(loggedVals,loggedVals0)! Should be true
    z <- p1/(p1+(1-p1)*exp(loggedVals))
    sumZ = sum(z)
    sumXZ = sum(x*z)
    sum1Z = (n - sumZ) #sum(1-z)
    
    p1.new <- sumZ/n
    lambda.new <- sumZ/sumXZ
    mu.new <- (sumX - sumXZ)/sum1Z  #sum((1-z)*x)/sum(1-z)
    sigma.sq.new <- sum(((x-mu.new)^2)*(1-z))/sum1Z
    
    #mu.new2 = mu.new^2
    #sumX2Z = sum(x2*z)
    #sigma.sq.new <- (sumX2 - 2*mu.new*sumX + n*mu.new2 - sumX2Z + 2*mu.new*sumXZ - mu.new2*sumZ)/sum1Z
    
    diff <- !(all(c(abs(lambda.new-lambda),abs(mu.new-mu),abs(sigma.sq.new-sigma.sq),abs(p1.new-p1))<epsilon))
    lambda <- lambda.new; mu <- mu.new; sigma.sq <- sigma.sq.new; p1 <- p1.new
  }
  list(lambda,mu,sqrt(sigma.sq),p1)
}