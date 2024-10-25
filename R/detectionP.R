#Use P-value calcs to obtain failed samples/sites:
getFailedPval = function(detPmatrix,pvalThresh = 0.05, siteFailThresh=0.05,sampleFailThresh=0.05) {
  detPmatrix[is.na(detPmatrix)] = 1 #set maximum pval if no value
  isfail = detPmatrix > pvalThresh
  dim = dim(detPmatrix) #nsites x nsamples
  failed_samples  = colnames(detPmatrix)[ colSums(detPmatrix)/dim[1] > sampleFailThresh]
  failed_sites  = rownames(detPmatrix)[ rowSums(detPmatrix)/dim[2] > siteFailThresh]
  return(list(samples=failed_samples,sites=failed_sites))
}

#SESAME: remaining Functions comes from https://github.com/zwdzwd/sesame/blob/master/R/detection.R
#Detection P-value based on ECDF of negative control OR out-of-band signal (input is flexible)
sesame_detectionPecdf = function(M,U,iR,iG,i2,bkgR,bkgG) {
  #M methylated matrix
  #U methylated matrix
  #iR index of Red channel (Type 1 probes)
  #iG index of Green channel (Type 1 probes)
  #i2 index of Type 2 probes
  #bkgR Background signal to be used for building background noise model (Red channel)
  #bkgG Background signal to be used for building background noise model (Green channel) 
  
  ## p-value is the minimium detection p-value of the 2 alleles
  detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))
  
  getMax = function(x,y) {
    z = x
    isLarger = y>x 
    isLarger[is.na(isLarger)] = FALSE
    z[isLarger] = y[isLarger]
    return(z)
  }
  for(i in 1:ncol(detP)) { #loop through each sample (easy to parallelize)
    
    #OOB CASE
    #Note that ewastools:oob to vectorizes methylated and unmethyled
    if(is.null(dim(bkgR))) { #bkgR is a list
      funcG <- ecdf(c(bkgG$M[,i],bkgG$U[,i])) #unlist both methylated/unmethylated signal (per-sample)
      funcR <- ecdf(c(bkgR$M[,i],bkgR$U[,i])) #unlist both methylated/unmethylated signal (per-sample)
      
    } else {
    #NEGATIVE CONTROL CASE
      funcG <- ecdf( bkgG[,i] ) #unlist both methylated/unmethylated signal (per-sample)
      funcR <- ecdf( bkgR[,i] ) #unlist both methylated/unmethylated signal (per-sample)
    }
    
    detP[iR,i] <- 1 - getMax( funcR(M[iR,i]),funcR(U[iR,i]) )
    detP[iG,i] <- 1 - getMax( funcG(M[iG,i]),funcG(U[iG,i]) )
    detP[i2,i] <- 1 - getMax( funcG(M[i2,i]),funcR(U[i2,i]) )
    
    detP[,i][is.na(detP[,i])] = 1 #put sites with NA to 1 (as done in SeSAMe)
  }
  return(detP)  
}

#Extracted Code from https://github.com/hhhh5/ewastools/blob/master/R/detectionP.R
#Recommended method to obtain p-value statistics based on oob measures
ewastools_detectionP <- function(M,U,iR,iG,i2){
  #raw ewastools
  
  #helpfunctions  
  summits = function (beta) {
    d <- density(beta,bw=0.01,na.rm=TRUE)
    l = which(d$x <0.4)
    u = which(d$x >0.6)
    l = l[which.max(d$y[l])]
    u = u[which.max(d$y[u])]
    d$x[c(l,u)]
  }
  
  detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))
  
  
  for(j in 1:ncol(M)){ #for each sample
    
    beta = M[,j]/(U[,j]+M[,j])
    
    # red color channel
    # locate the peaks for (un)methylated sites
    sR = summits(beta[iR])
    
    # pick 1000 CpG sites closest to the peaks
    # sR[2] is the   methylated peak, provides unmethylated background signal
    bkgU = head(order(abs(beta[iR]-sR[2])),n=1000)
    # sR[1] is the unmethylated peak, provides   methylated background signal
    bkgM = head(order(abs(beta[iR]-sR[1])),n=1000)
    
    bkgU = iR[bkgU]
    bkgM = iR[bkgM]
    
    # median and MAD for these 1000 sites
    muUR = median(U[bkgU,j],na.rm=TRUE)
    muMR = median(M[bkgM,j],na.rm=TRUE)
    
    sdUR = mad(U[bkgU,j],na.rm=TRUE)
    sdMR = mad(M[bkgM,j],na.rm=TRUE)
    
    # green color channel
    sG = summits(beta[iG])
    
    bkgU = head(order(abs(beta[iG]-sG[2])),n=1000)
    bkgM = head(order(abs(beta[iG]-sG[1])),n=1000)
    
    bkgU = iG[bkgU]
    bkgM = iG[bkgM]
    
    muUG = median(U[bkgU,j],na.rm=TRUE)
    muMG = median(M[bkgM,j],na.rm=TRUE)
    
    sdUG = mad(U[bkgU,j],na.rm=TRUE)
    sdMG = mad(M[bkgM,j],na.rm=TRUE)
    
    detP[iR,j] = pnorm(U[iR,j]+M[iR,j],mean=muUR+muMR,sd=sqrt(sdUR^2+sdMR^2),lower.tail=FALSE)
    detP[iG,j] = pnorm(U[iG,j]+M[iG,j],mean=muUG+muMG,sd=sqrt(sdUG^2+sdMG^2),lower.tail=FALSE)
    detP[i2,j] = pnorm(U[i2,j]+M[i2,j],mean=muUR+muMG,sd=sqrt(sdUR^2+sdMG^2),lower.tail=FALSE)
  }
  return(detP)
} #end function


#Copy of minfi-approach (uses negatives)
#compare with https://github.com/hansenlab/minfi/blob/master/R/detectionP.R (includes wrong formula on SD)
ewastools_detectionP.neg <- function(M,U,iR,iG,i2,bkgR,bkgG,approach="ewastools"){
  muG = apply(bkgG,2,median,na.rm=TRUE)
  sdG = apply(bkgG,2,mad   ,na.rm=TRUE)
  muR = apply(bkgR,2,median,na.rm=TRUE) 
  sdR = apply(bkgR,2,mad   ,na.rm=TRUE) 
  detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))
  for(j in 1:ncol(M)) {
    
    if(approach=="ewastools") {
      detP[iR,j] = pnorm(U[iR,j]+M[iR,j],mean=2*muR[j],sd=sqrt(2)*sdR[j],lower.tail=FALSE) #RED CHANNEL (Type1) 
      detP[iG,j] = pnorm(U[iG,j]+M[iG,j],mean=2*muG[j],sd=sqrt(2)*sdG[j],lower.tail=FALSE) #GREEN CHANNEL (Type1)  
      detP[i2,j] = pnorm(U[i2,j]+M[i2,j],mean=muR[j]+muG[j],sd=sqrt(sdR[j]^2+sdG[j]^2),lower.tail=FALSE) #BOTH CHANNEL  (Type2) 
    }
    if(approach=="minfi") {
      detP[iR,j] = pnorm(U[iR,j]+M[iR,j],mean=2*muR[j],sd=2*sdR[j],lower.tail=FALSE) #RED CHANNEL  
      detP[iG,j] = pnorm(U[iG,j]+M[iG,j],mean=2*muG[j],sd=2*sdG[j],lower.tail=FALSE) #GREEN CHANNEL  
      detP[i2,j] = pnorm(U[i2,j]+M[i2,j],mean=muR[j]+muG[j],sd=sdR[j]+sdG[j],lower.tail=FALSE) #BOTH CHANNEL  
    }
  } 
  return(detP)
}
