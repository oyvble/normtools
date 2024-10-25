
#Calculation of probability density function Z=U/V where [U,V] ~ N(mu,SIGMA)
#EXACT FORMULA provided by Oyvind Bleka (Hand-derived)
pRatioNormal = function(z, MU, SIGMA) {
  #if(length(MU)!=2 || !all(dim(SIGMA)==2)) stop("Wrong dimension of input. Dimension must be 2.")
  
  s1 = SIGMA[1,1]
  s2 = SIGMA[2,2]
  s12 = SIGMA[1,2]
  m1 = MU[1]    
  m2 = MU[2]    
  
  det = abs(s1*s2-s12^2) #obtain determinant
  a = s2*z^2 + s1 - 2*s12*z
  b = s2*m1*z + s1*m2 - s12*(m1+m2*z)
  c = s2*m1^2 + s1*m2^2 - 2*s12*m1*m2
  d = 2*det
  
  alpha = a/d
  gamma = b/a
  lambda = c/d - b^2/(a*d)
  
  #2*(alpha*(v-gamma)^2 + lambda)  #this is SSQ
  val = exp(-lambda)*gamma*sqrt(pi/alpha)/(2*pi*sqrt(det))
  return(val)
}

#Numerical calculation of quantiles
qRatioNormal = function(p, MU, SIGMA, accuracy=4, smallfold=100, npre = 100 ) {
  #if(length(MU)!=2 || !all(dim(SIGMA)==2)) stop("Wrong dimension of input. Dimension must be 2.")
  #smalltol tolerance of being small number
  #npre number of points to search in pre-search
  #accuracy indicate grid size used in calculations (must be integers)
  bw = 1/10^accuracy
  
  #Step 1: Perform rough search:
  zv = seq(0,1,l=npre) #rough search
  pv = pRatioNormal(zv, MU,SIGMA)
  
  smalltol = max(pv)/smallfold #considering fold change difference relative to max to obtain tolerance
  indOK = which(pv>smalltol) #indicate region with positive pdf
  zmin = zv[max(1,min(indOK)-1) ] #lower index
  zmax = zv[min(npre,max(indOK)+1)]#upper index

  zmin = round(zmin,accuracy) #round to closest (bw-aligned)
  zmax = round(zmax,accuracy) #round to closest (bw-aligned)

  #Step 2: Perform fine search  
  zv = seq(zmin,zmax,by=bw) #obtain grid length
  pv = pRatioNormal(zv, MU,SIGMA)
  #plot(zv,pv)
  cdf = cumsum(pv)*bw #obtain cumulative-df
  
  quantile = sapply( p, function(x) zv[max(which(cdf<x))])
 #plot(zv,cdf); abline(v=quantile); abline(h=p)
  return(quantile) #return quantiles of provided probabilities
}

#site=100;i=1;o1=0;o2=0;n=10000
#show distribution of beta for a particularr site/sample based on simulations ()
show_beta_uncertainty_sim = function(meth,site,i, o1=0, o2=0, n=10000, p = c(0.025,0.975))  {
  #input is object from ewastools::read_idats
  siteInd = which(meth$manifest$probe_id==site)
  
  #Un-Methylated  
  U_mean = meth$U[siteInd,i] #Mean
  U_sd = meth$T[siteInd,i] #SD
  U_n = meth$V[siteInd,i] # beadCounts
  U_sem =  U_sd/sqrt(U_n)
  
  #Methylated  
  M_mean = meth$M[siteInd,i] #Mean
  M_sd = meth$S[siteInd,i] #SD
  M_n = meth$N[siteInd,i] #BeadCounts
  M_sem =  M_sd/sqrt(M_n)
  
  #TOTAL (M+U)
  T_mean = M_mean + U_mean + o2
  T_sem = sqrt(M_sem^2 + U_sem^2) #standard dev of total
  
  sitename = meth$manifest$probe_id[siteInd]
  samplename = meth$meta[i,1]
  
  #obtain uncertinty about beta-values (based on uncertainty of methylated,unmethylated)
  M = rnorm(n, M_mean, M_sem)
  U = rnorm(n, U_mean, U_sem) 
  B = (M+o1)/(M+U+o2)
  hist(B,breaks=30,probability=T,main=sitename);mtext(samplename)
  
  
  #Specify COV and MEAN
  COV = diag(c(M_sem,T_sem))^2
  COV[1,2] <- COV[2,1] <- M_sem^2
  #determinant(COV,log=T)
  MU = c(M_mean, T_mean) #mean for each
  z = seq(min(B),max(B),l=1000)
  lines(z,pRatioNormal(z,MU,COV))
  
  q = qRatioNormal(p,MU,COV,acc=4) #obtain quantiles
  abline(v=q)
}


#calculation of (mean)+percentiles per site/sample
calc_beta_uncertainty = function(meth,sites=NULL, p=c(0.025,0.975), o1=0, o2=0, widthOnly=FALSE, acc=4, nCores = parallel::detectCores())  {
  #sites=NULL
  #p quantiles to calculate (q1,q2 etc) in addition to the mean (if not widthOnly)
  #widthOnly Whether only with of uncertainty (q2-q1). Should be stored (NB: require length(p)==2)
  samplenames = meth$meta[,1][[1]] #obtain sample name
  if(is.null(sites)) sites= meth$manifest$probe_id #all sites are considered
  siteInds = match(sites,meth$manifest$probe_id) #obtain index of site positions (to loop through)
  
  c1 <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(c1)
  
  #meth.o <- foreach(i=1:length(sites),.combine=rbind,.export=c("qRatioNormal","pRatioNormal")) %dopar% {
  #for(i in 1:length(sites)) { #for each site
  #  TAB <- foreach::foreach(siteInd=siteInds ,.combine=rbind,.export=c("qRatioNormal","pRatioNormal")) %dopar% {
  RESLIST <- foreach::foreach(j=1:length(samplenames),.export=c("qRatioNormal","pRatioNormal")) %dopar% {
    
    if(widthOnly) {
      ARRAY = rep(NA,length(siteInds))
    } else {
      ARRAY = matrix(NA,nrow=length(p)+1,ncol=length(siteInds))
    }
    # sitename = sites[1]
    for(i in 1:length(siteInds)) { #for each site
      #      j = 1
      siteInd = siteInds[i]
      #Un-Methylated  
      U_mean = meth$U[siteInd,j] #Mean
      U_sd = meth$T[siteInd,j] #SD
      U_n = meth$V[siteInd,j] # beadCounts
      U_sem =  U_sd/sqrt(U_n)
      
      #Methylated  
      M_mean = meth$M[siteInd,j] #Mean
      M_sd = meth$S[siteInd,j] #SD
      M_n = meth$N[siteInd,j] #BeadCounts
      M_sem =  M_sd/sqrt(M_n)
      
      #obtain MEAN and sem of TOTAL (M+U)
      MU = c(M_mean + o1, M_mean + U_mean + o2) #mean for each
      T_sem = sqrt(M_sem^2 + U_sem^2) #standard dev of total
      
      isOK = (!is.na(U_sd) && U_sd>0) && (!is.na(M_sd) && M_sd>0) #check if OK for calculation
      if(isOK) { #meaning that the number of beads=1 for either M or U
        COV = diag(c(M_sem,T_sem))^2
        COV[1,2] <- COV[2,1] <- M_sem^2
        q = qRatioNormal(p,MU,COV,acc=acc)
      } else {
        q = rep(NA,length(p)) #not possible to obtain uncertinaty
      }
      
      if(widthOnly) {
        ARRAY[i] = c(q[2]-q[1])
      } else {
        ARRAY[,i] = c(MU[1]/MU[2],q)
      }
    } 
    return(ARRAY)
  }
  parallel::stopCluster(c1)
  if(widthOnly) {
    ARRAY = array(unlist(RESLIST),dim=c(length(sites),length(samplenames)),dimnames = list(sites,samplenames)) #mean, 2.5%,97.5% per sample/site
  } else {
    ARRAY = array(unlist(RESLIST),dim=c(1+length(p),length(sites),length(samplenames)),dimnames = list(c("mean",p),sites,samplenames)) #mean, 2.5%,97.5% per sample/site
  }  
  return(ARRAY)
}

#calculation of width per site/sample
#Assume returned object calc_beta_uncertainty 
show_beta_CI = function(obj,toPDF=FALSE)  {
  #p=c(0.025,0.975); o1=0; o2=0;acc=4
  # sites = meth$manifest$probe_id[1:10*1000]
  #input is object from ewastools::read_idats
  if(length(dim(obj))!=3) stop("Wrong format on input object")
  sitenames = dimnames(obj)[[2]] #obtain sample name
  samplenames = dimnames(obj)[[3]] #obtain sample name
  
  if(toPDF)    pdf(paste0("betaUncertaintyCI_",proj,".pdf"),width=20,height=12)
  for(i in 1:length(sitenames)) {
#      i = 1 #site selection
    site = sitenames[i]
    sub = obj[,i,]
    yr = range(sub)
    
    plot(sub[1,],pch=19,main=site,ylab="",ylim=yr)
    for(j in 1:ncol(sub)) {
      lines(rep(j,2), sub[-1,j]) #show distribution
    }
    #show_beta_uncertainty1(meth,site,i=30)
  }  
  if(toPDF)    dev.off()
}


show_beta_width = function(obj,toPDF=FALSE)  {
  if(length(dim(obj))!=2) stop("Wrong format on input object")
  sitenames = dimnames(obj)[[1]] #obtain sample name
  samplenames = dimnames(obj)[[2]] #obtain sample name
}




