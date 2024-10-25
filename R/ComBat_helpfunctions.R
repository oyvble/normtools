
#Vectorized modification of sva:::int.eprior to improve speed and memory handling
#REFERENCE: Johnsen et al 2007: Adjusting batch effects in microarray expression data using empirical Bayes methods

  int.eprior2 <- function(sdat,g.hat,d.hat,nB=2,quiet=FALSE){
   #nB - number of site segments
   bigr <- nrow(sdat) #number of sites
   if(nB > bigr) nB <- bigr
   bb <- floor(bigr/nB) #length of segments
   partid <- list()
   for(b in 1:nB) { #get indices for sites in each segments
    if(b==1) {
     partid[[b]] <- 1:bb
    } else if( b<nB) {
     partid[[b]] <- tail(partid[[b-1]],1) + 1:bb
    } else {
     partid[[b]] <- (tail(partid[[b-1]],1)+1):bigr
    }
   }
   val <- numeric()
   for(b in 1:nB) { #for each segments
    if(!quiet) print(paste0("Segment ",b," out of ",nB,"..."))
    sdat2 <- t(sdat[partid[[b]],]) #transpose from [Site x Sample] = [Sample x Site]
    g.hat2 <- g.hat[partid[[b]]]
    d.hat2 <- d.hat[partid[[b]]]
    r <- ncol(sdat2) #number of markers
    n <- nrow(sdat2) #number of samples 
    SSQz <- colSums(sdat2^2)
    Sz <- colSums(sdat2) #sum
    rmind <- (0:(r-1))*r + 1:r
    G <- rep(g.hat2,r)
    G <- t(matrix(G[-rmind],nrow=r-1,ncol=r))
    R <- replicate(r-1,Sz) #replicate matrix
    Sgz <- G*R 
    rm(R);gc()
  
    SSQg <- n*(G^2)
    sum2 <- SSQz + SSQg - 2*Sgz #at time when uses most RAM
    rm(SSQz,SSQg,Sgz);gc()
   
    D <- rep(d.hat2,r)
    D <- t(matrix(D[-rmind],nrow=r-1,ncol=r))
    LH <- exp( (-n/2)*log(2*pi*D) - sum2/(2*D) )
    rm(sum2);gc()
    CS <- rowSums(LH) 
    val2 <- rbind(rowSums(G*LH)/CS,rowSums(D*LH)/CS)
    rm(LH,G,D);gc()
    val <- cbind(val,val2)
   } 
   return(val)
  }
  
  
  #Slight modification of sva::it.sol
  it.sol2  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
    n <- rowSums(!is.na(sdat))
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while(change>conv){
      g.new <- postmean(g.hat, g.bar, n, d.old, t2)
      sum2 <- rowSums((sdat - g.new)^2, na.rm=TRUE) #rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
      d.new <- postvar(sum2, n, a, b)
      change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old)
      g.old <- g.new
      d.old <- d.new
      count <- count+1
    }
    ## cat("This batch took", count, "iterations until convergence\n")
    rbind(g.new, d.new)
    #adjust <- 
    #rownames(adjust) <- c("g.star","d.star")
    #adjust
  }
  
  
  #INTERNAL HELPFUNCTION FOR param-version
  postmean <- function(g.hat,g.bar,n,d.star,t2){
    (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
  }
  
  postvar <- function(sum2,n,a,b){
    (.5*sum2 + b) / (n/2 + a - 1)
  }
  
  
  # Inverse gamma distribution density function. (Note: does not do any bounds checking on arguments)
  dinvgamma <- function (x, shape, rate = 1/scale, scale = 1) {
    # PDF taken from https://en.wikipedia.org/wiki/Inverse-gamma_distribution
    # Note: alpha = shape, beta = rate
    stopifnot(shape > 0)
    stopifnot(rate > 0)
    ifelse(x <= 0, 0, ((rate ^ shape) / gamma(shape)) * x ^ (-shape - 1) * exp(-rate/x))
  }
  