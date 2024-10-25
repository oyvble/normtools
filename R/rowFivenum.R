
#helpfunction for fast fivenum-calculation
rowFivenum = function(X, transposed=FALSE, parallel=TRUE) {
  bool = require(Rfast)
  if(bool) {
    if(transposed) {
      X = Rfast::colSort(X, descend = FALSE, parallel = parallel)
      n = nrow(X) #number of sites
    } else {
      X = Rfast::rowSort(X, descend = FALSE, parallel = parallel)
      n = ncol(X) #number of samples
    }
    
    n4 <- floor((n + 3)/2)/2
    d <- c(n4, (n + 1)/2, n + 1 - n4)
    d1 = floor(d)
    d2 = ceiling(d)
    
    if(transposed) {
      Y = rbind(X[1,], 0.5*(X[d1,]+X[d2,]), X[n,])
    } else {
      Y = cbind(X[,1], 0.5*(X[,d1]+X[,d2]), X[,n])
    }
  } else {
    
    if(transposed) {
      Y =  apply( X, 2, stats::fivenum) #calc per col
    } else {
      Y = t(apply( X, 1, stats::fivenum)) #calc per row
    }    
  }
  return(Y)
}
