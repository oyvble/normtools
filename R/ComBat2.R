#This is a rewritten script of sva::ComBat
#The optimization of the code was done by Oyvind Bleka.
#Details: The improvment is both speed and memory (a maximum number can be given)
#Written to make both Combat param and non-param accessible to run with large dataset
#This function version inclde ref.batch: https://github.com/jtleek/sva-devel/blob/master/R/ComBat.R

#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample). Should be class type bigmemory::big.matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param maxGB Max GB to use (batched up operations)
#' @param log2transform whether log2 transformation of dat is needed
#' @param transposed whether the dat matrix is transposed (samples x probes)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#' @export

ComBat2 <- function (dat, batch, mod = NULL,  par.prior = FALSE, mean.only = FALSE, prior.plots = FALSE, maxGB=2, log2transform=FALSE, transposed=FALSE,ref.batch = NULL) {

  if(transposed) { #Data must be given as (probes x samples). Transpose if not.
    dat = t(dat)
  }
  if(log2transform) { #Recommended to transform data to real domain
    dat = log2(dat/(1-dat))
  }
  
  if (mean.only == TRUE) {
      cat("Using the 'mean only' version of ComBat\n")
  }
  if (length(dim(batch)) > 1) {
      stop("This version of ComBat only allows one batch variable")
  }
  nsites = nrow(dat) #number of sites (p)
  
  ## make batch a factor and make a set of indicators for batch
  batch <- as.factor(batch)
  if(any(table(batch)==1)){mean.only=TRUE} #only one batch observed (use Mean.only)
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  
  batchmod <- model.matrix(~-1 + batch) #no intercept, value=1 for corrsponding batch. All batches relative to each other

  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
        stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch=",ref.batch, " as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }

  cat("Found", nlevels(batch), "batches\n")
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
      batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
      mean.only = TRUE
      cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)

  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check]) #final designmatrix
  
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
      if (ncol(design) == (n.batch + 1)) {
          stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
      }
      if (ncol(design) > (n.batch + 1)) {
          if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, -c(1:n.batch)]))) {
              stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
          }
          else {
              stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
          }
      }
  }
  
  #cat("Standardizing Data across genes\n MEMORY EFFICIENT IMPLEMENTATION CARRIED OUT (no copy of dat)")
  
  #s.data is replaced with dat, where this is standardized wrt design matrix
  
  #Procedure in R:
  #1) Obtain B.hat (v x p) large matrix
  B.hat <- tcrossprod(solve(t(design) %*% design),design)
  B.hat <-  tcrossprod(B.hat,dat) #Must store '(nbatches+nvars) x p' large matrix 

  #2) Obtain grand mean (p long vector)
  if(!is.null(ref.batch)){   ## change grand.mean for ref batch
  	grand.mean <- B.hat[ref, ]
  } else {
	  #grand.mean2 <- c(t(n.batches/n.array) %*% B.hat[1:n.batch, ])
	  grand.mean <-   colSums(B.hat[1:n.batch, ]*n.batches)/n.array
  }

  #2) Obtain var.pooled (very RAM expensive)  
   if(!is.null(ref.batch)) { ### change var.pooled for ref batch
      ref.dat <- dat[, batches[[ref]]]

      tmp = design[batches[[ref]],] %*% B.hat
      tmp <- (ref.dat - t(tmp))^2 #p x n matrix
      sd.pooled <- sqrt(rowSums(tmp)/n.batches[ref])  #p long vector
      #var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
   } else {
      tmp = design %*% B.hat
      tmp <- (dat - t(tmp))^2 #p x n matrix
      sd.pooled <- sqrt(rowSums(tmp)/n.array)  #p long vector
      #var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) #=sd.pooled^2     
   }
  #max(abs(sd.pooled-sqrt(var.pooled)))
  
  #3) Obtain standardized mean (CENTERING): p x n #matrix
  if (!is.null(mod)) { #if extra variables included
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean = t(tmp %*% B.hat) #p x n matrix (update B.hat)
    stand.mean = stand.mean + grand.mean #obtain grand center
    # grand.mean %*% t(rep(1, n.array)) #p x n matrix (each sample has same mean): Each column is same 
  } else {
    stand.mean <- grand.mean #simply copy
  }
  rm(B.hat);gc()
  
  #4) CENTERING DATA  p x n #matrix
  s.data <- dat - stand.mean #stand.mean must be used further
  rm(dat);gc()
  
  #5) SCALING DATA: p long vector (be sure that s.data is p x n)
  s.data <- s.data/sd.pooled# scale

  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[, 1:n.batch]
  gamma.hat <- solve(t(batch.design) %*% batch.design) %*% t(batch.design) #v x n
  gamma.hat <- tcrossprod(gamma.hat,s.data) #gamma is v x p large vector 

  #BLOCK MATRIX USED AS INPUT:
  delta.hat <- matrix(nrow=length(batches),ncol=nsites ) #B x p
  for (i in 1:length(batches)) { #performed for each batch (doable in R)
      if (mean.only == TRUE) {
         delta.hat[i,] <- rep(1, nsites)
      }
      else {
         X <- s.data[, batches[[i]]]
         delta.hat[i,] <-  rowSums((X - rowMeans(X))^2)/(ncol(X)-1) #tmp <- apply(s.data[, batches[[i]]],1, var, na.rm = T)
      }
  }
    
 ## Find EB batch adjustments
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nsites) #init vars (n.batch x nsites)
  if (par.prior) { #PARAMETERIC
   # Following four find empirical hyper-prior values
   #precalculation:
   m <- rowSums(delta.hat)/nsites #mean per batch
   s2 <- rowSums((delta.hat - m)^2)/(nsites-1) #variance per batch
   a.prior = (2*s2 + m^2) / s2
   b.prior = (m*s2 + m^3) / s2
   
   gamma.bar <- rowSums(gamma.hat)/nsites #apply(gamma.hat, 1, mean) #average gamma per block
   t1 <- rowSums(gamma.hat)/nsites #mean per batch
   t2 <- rowSums((gamma.hat - t1)^2)/(nsites-1) #variance per batch
   
   
   ## Plot empirical and parametric priors
   if (prior.plots) {
     old_pars <- par(no.readonly = TRUE)
     on.exit(par(old_pars))
     par(mfrow=c(2,2))
     
     ## Top left
     tmp <- density(gamma.hat[1,])
     plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
     xx <- seq(min(tmp$x), max(tmp$x), length=100)
     lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
     
     ## Top Right
     qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
     qqline(gamma.hat[1,], col=2)
     
     ## Bottom Left
     tmp <- density(delta.hat[1,])
     xx <- seq(min(tmp$x), max(tmp$x), length=100)
     tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
     plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),main=expression(paste("Density Plot of First Batch ", hat(delta))))
     lines(tmp1, col=2)
     
     ## Bottom Right
     invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
     qqplot(invgam, delta.hat[1,],main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),ylab="Sample Quantiles", xlab="Theoretical Quantiles")
     lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
   }
   
   cat("Finding parametric adjustments\n")
   for(i in 1:n.batch) { #traverse each batch:
     print(paste0("Evaluting batch ",levels(batch)[i],"(",i," out of ",n.batch,")"))
     if (mean.only) {
       gamma.star[i,] <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
       delta.star[i,] <- rep(1, nsites)
     } else {
#       sdat=s.data[, batches[[i]] ];g.hat=gamma.hat[i, ];d.hat=delta.hat[i, ];g.bar=gamma.bar[i];t2=t2[i];a= a.prior[i];b=b.prior[i];conv=.0001
       temp <- normtools::it.sol2(s.data[, batches[[i]] ], gamma.hat[i, ],delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],b.prior[i])
       gamma.star[i,] <- temp[1, ]
       delta.star[i,] <- temp[2, ]
     }
   }
    
 } else { #NON-PARAMETERIC
  cat("Finding nonparametric adjustments\n")
  #required RAM size:
  getSize <- function(s) 8*s*(s-1)/1e6 #size in megabytes, 8 bit size per unit (element in matrices etc)
  reqGB <- 4*getSize(nsites)/1000 #required GB
  nB <-  ceiling(reqGB/maxGB) #number of batches to use 
   
  for (i in 1:n.batch) {
    print(paste0("Evaluting batch ",levels(batch)[i],"(",i," out of ",n.batch,")"))
    if (mean.only) {
     delta.hat[i, ] = 1
    }
    #Call internal function:
    outAdj <- normtools::int.eprior2(sdat=as.matrix(s.data[, batches[[i]]]),g.hat=gamma.hat[i, ], d.hat=delta.hat[i, ],nB=nB,quiet=TRUE)
  
    #handle if some are zero:
    isnaG <- is.na(outAdj[1,])
    isnaD <- is.na(outAdj[2,])
    if(sum(isnaG)>0 | sum(isnaD)>0) warning("Zero-weights found for adjusting gamma.star and delta.star. These are remained as original values.")
    outAdj[1,isnaG] <- gamma.hat[i,isnaG] #copy existing
    outAdj[2,isnaD] <- delta.hat[i,isnaD] #copy existing
    gamma.star[i,] <- outAdj[1, ]
    delta.star[i,] <- outAdj[2, ]
  }
 } #end choice of method
 rm(delta.hat,gamma.hat);gc() #remove data
 
 if(!is.null(ref.batch)){ #Special handling for reference batch:
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
 }
 
 cat("Adjusting the Data\n")
 for (i in 1:n.batch) { #loop through each batch
  # s.data[, i] <- (s.data[, i] - t(batch.design[i,] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
   sampleInd = batches[[i]] #index of samples in batch
   s.data[, sampleInd] <- (s.data[, sampleInd] - gamma.star[i,])/sqrt(delta.star[i, ]) #scale global
 }
 cat("Adjusting Data back to original scale\n")
 s.data <- s.data * sd.pooled #scale back data
 rm(gamma.star,delta.star,sd.pooled);gc()
 s.data <- s.data + stand.mean #centering data back
 rm(stand.mean);gc()
 
 
 ## Do not change ref batch at all in reference version
 if(!is.null(ref.batch)){
 	s.data[, batches[[ref]]] <- ref.dat #Insert back into data (unchanged)
 } 
  
 if(log2transform) {
   s.data = 1/(1+2^(-s.data)) #transform back
 }
 if(transposed) {
   s.data = t(s.data)
 }
 
 return(s.data)
}
