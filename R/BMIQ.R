#Taken from ENmix::bmiq.mc (note: this function is removed from updated versions of ENmix)
#pType=NULL;platform="450k";  nCores = detectCores(); verbose=TRUE; nfit=10000;implementation="fast";#implementation="watermelon"
bmiq.mc = function(bdat,pType=NULL,platform="450k", nfit=50000, nreq=2, nCores = detectCores(), verbose=FALSE,implementation="fast",seed=NULL) {
  #bdat: The (sites x samples) beta-matrix (values within [0,1]). 
  # pType: corresponding vector specifying probe design type (1=type1,2=type2). 
  #implementation={"watermelon"=wateRmelon (uses v1.1),"fast"=v1.4 Teschendorff-reimplement + reimplemnted blc")
  #https://aeteschendorff-lab.github.io/images/BMIQ/BMIQ_1.4.R
  #platform   
  if(is.null(pType)) {
    if(platform=="450k") require(IlluminaHumanMethylation450kmanifest)
    if(platform=="EPIC") require(IlluminaHumanMethylationEPICmanifest)
    manifestType = setNames(Manifest$Type,Manifest$Name) #obtain probe type
    #all(rownames(bdat)%in%names(manifestType)) #check that all probes has a type
    pType = rep(1,nrow(bdat)) #set probeType 1 by default
    pType[rownames(bdat)%in%names(manifestType[manifestType=="II"])] = 2 #indicate typeII probes
  }
#  table(pType)  
  nSamples = ncol(bdat) #get sampple names
  N = ceiling(nSamples/(nCores)) #segments to run
  parts = rep(1:N, each = ceiling(nSamples/N))[1:nSamples] 
  c1 <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(c1)
  if(verbose) print(paste0(N," segments to traverse.."))
  for (i in 1:N) { #for each batch
    if(verbose) print(i)
    id = which(parts == i)
    b1 = bdat[, id,drop=FALSE ]#obtaining subset of bdat to consider

#    s=4
# b2 <- BMIQfast(beta.v=b1[,s], design.v=pType, nfit=nfit,seed=seed)
    
    if(implementation=="fast") {
#beta.v=b1[,s] ; design.v=pType
      b2 <- foreach(s = 1:ncol(b1), .combine = cbind,.export = c("BMIQfast","blc1","betaEst1")  )%dopar%normtools:::BMIQfast(b1[,s], pType, nfit=nfit,nreq=nreq,seed=seed)$nbeta
    } else if(implementation=="watermelon") {
      b2 <- foreach(s = 1:ncol(b1), .combine = cbind)%dopar%wateRmelon::BMIQ(b1[,s], pType, plots = FALSE,nfit=nfit)$nbeta
    } else {
        print("Implementation not found!")
    }
    bdat[,id] = b2 #insert
  }
  parallel::stopCluster(c1)
  return(bdat)
}



###Re-implementation of "BMIQv_1.4: 23rd Sep 2014 (Teschendorff)" only exchanging
#beta.v= bdat[,s];design.v=pType; nL=3;doH=TRUE;th1.v=c(0.2,0.75);th2.v=NULL;niter=5;tol=0.001
BMIQfast <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=25,tol=0.001,nreq=2,seed=NULL,maxIter=20) {
  ### beta.v: vector consisting of beta-values for a given sample. NAs are not allowed, so these must be removed or imputed prior to running BMIQ. Beta-values that are exactly 0 or 1 will be replaced by the minimum positive above zero or maximum value below 1, respectively.
  ### design.v: corresponding vector specifying probe design type (1=type1,2=type2). This must be of the same length as beta.v and in the same order.
  #nreq Number of required prediction for each specific class (1,2,3) after Beta-mixture prediction (blc1)
  #maxIter max number of repeats until giving up
  #require(RPMM);  #inlcude blc
  if(!is.null(seed)) set.seed(seed)

  type1.idx <- which(design.v==1);
  type2.idx <- which(design.v==2);
  
  beta1.v <- beta.v[type1.idx];
  beta2.v <- beta.v[type2.idx];
  
  ### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
  if(min(beta1.v)==0){
    beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
  }
  if(min(beta2.v)==0){
    beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
  }
  if(max(beta1.v)==1){
    beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
  }
  if(max(beta2.v)==1){
    beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
  }
  
  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
  w0.m[beta1.v <= th1.v[1],1] <- 1;
  w0.m[beta1.v > th1.v[1] & beta1.v <= th1.v[2],2] <- 1;
  w0.m[beta1.v > th1.v[2],3] <- 1;
  
  
  #Alternative (faster implementation) of  apply(X,1,which.max);
  getWhichMax = function(X) {
    ncols = ncol(X)
    indMax = rep(1,nrow(X))
    for(j in 2:ncols) {
      isGreater = X[,j]>X[,j-1]
      indMax[isGreater] = j #update index
    }
    return(indMax)
  }
  
  ### fit type1: NOTE MUST RE-RUN UNTIL SUFFICIENTLY MANY OF EACH TYPE!
  isOK = FALSE
  cc = 0 #counter
  while(!isOK) {
    rand.idx <- sample(1:length(beta1.v),min(nfit,length(beta1.v)),replace=FALSE)
    em1.o <- normtools:::blc1(beta1.v[rand.idx],w0.m[rand.idx,],maxiter=niter,tol=tol); #fit beta-mix
    
    #obtain which max for each marker
    subsetclass1.v <- getWhichMax(em1.o$w) #apply(em1.o$w,1,which.max);
    predCount <- table(subsetclass1.v) #overview of predicted classess
    ##table(apply(em1.o$w,1,which.max)) Must be identical
    
    #check if sufficiently number of samples to do density
    if(sum(predCount>=nreq)==3) {
      isOK = TRUE
    } else {
      warning("Failed to fit Beta-mixture distribution to all 3 classes for Type 1 probes! Retrying...")
      cc <- cc + 1
      if(cc==maxIter) stop("Failed to apply BMIQ")
    }
  }

  #FITTING DENSITIES
  subset1 = beta1.v[rand.idx[subsetclass1.v==1]]
  subset2 = beta1.v[rand.idx[subsetclass1.v==2]]
  subset3 = beta1.v[rand.idx[subsetclass1.v==3]]
  
  subsetth1.v <- c(mean(c(max( subset1 ),min( subset2 ))),mean(c(max( subset2 ),min( subset3 ))));
  class1.v <- rep(2,length(beta1.v));
  class1.v[ beta1.v < subsetth1.v[1] ] <- 1;
  class1.v[ beta1.v > subsetth1.v[2] ] <- 3;
  nth1.v <- subsetth1.v;
  
  ### Estimate Modes 
  d1U.o <- density(beta1.v[class1.v==1])
  d1M.o <- density(beta1.v[class1.v==3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[ beta2.v<0.4 ]);
  d2M.o <- density(beta2.v[ beta2.v>0.6 ]);
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  
  
  ### now deal with type2 fit
  th2.v <- vector();
  th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
  th2.v[2] <- nth1.v[2] + (mod2M-mod1M);
  
  ### estimate initial weight matrix 
  w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
  w0.m[ beta2.v <= th2.v[1] ,1] <- 1;
  w0.m[ beta2.v > th2.v[1] & beta2.v <= th2.v[2],2] <- 1;
  w0.m[ beta2.v > th2.v[2] ,3] <- 1;
  
  
  ### Fitting EM beta mixture to type2 probes: DON*T NEED TO RE-RUN UNTIL SUFFICIENTLY MANY OF EACH TYPE!
  #hist(beta2.v)
  isOK = FALSE
  cc = 0 #counter
  while(!isOK) {
   rand.idx <- sample(1:length(beta2.v),min(nfit,length(beta2.v)),replace=FALSE)
   em2.o <- normtools:::blc1(beta2.v[rand.idx],w0.m[rand.idx,],maxiter=niter,tol=tol);

   ### for type II probes assign to state (unmethylated, hemi or full methylation)
   subsetclass2.v <- getWhichMax(em2.o$w) #apply(em2.o$w,1,which.max);
   predCount <- table(subsetclass2.v) #overview of predicted classess
   #print(predCount)
   ## table(apply(em2.o$w,1,which.max)) #Must be identical
  
   #check if sufficiently number of samples to do density
   if(sum(predCount>=nreq)==3) {
     isOK = TRUE
   } else {
     warning("Failed to fit Beta-mixture distribution to all 3 classes for Type 2 probes! Retrying...")
     cc <- cc + 1
     if(cc==maxIter) stop("Failed to apply BMIQ")
   }
  }
  #}
  
  #FITTING DENSITIES
  subset1 = beta2.v[rand.idx[subsetclass2.v==1]]
  subset2 = beta2.v[rand.idx[subsetclass2.v==2]]
  subset3 = beta2.v[rand.idx[subsetclass2.v==3]]
  
  subsetth2.v <- c(mean(c(max( subset1 ),min( subset2 ))),mean(c(max( subset2 ),min( subset3 ))));
  class2.v <- rep(2,length(beta2.v));
  class2.v[beta2.v < subsetth2.v[1]] <- 1;
  class2.v[beta2.v > subsetth2.v[2]] <- 3;
  
  classAV1.v <- vector();classAV2.v <- vector();
  for(l in 1:nL) {
    classAV1.v[l] <-  em1.o$mu[l];
    classAV2.v[l] <-  em2.o$mu[l];
  }
  
  ### start normalising type2 probes
  nbeta2.v <- beta2.v;
  ### select U probes
  lt <- 1;
  selU.idx <- which(class2.v==lt);
  selUR.idx <- selU.idx[ beta2.v[selU.idx] > classAV2.v[lt] ];
  selUL.idx <- selU.idx[ beta2.v[selU.idx] < classAV2.v[lt] ];
  
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt],em2.o$b[lt],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt],em1.o$b[lt],lower.tail=FALSE);
  nbeta2.v[selUR.idx] <- q.v;
  p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt],em2.o$b[lt],lower.tail=TRUE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt],em1.o$b[lt],lower.tail=TRUE);
  nbeta2.v[selUL.idx] <- q.v;
  
  ### select M probes
  lt <- 3;
  selM.idx <- which(class2.v==lt);
  selMR.idx <- selM.idx[ beta2.v[selM.idx] > classAV2.v[lt] ];
  selML.idx <- selM.idx[ beta2.v[selM.idx] < classAV2.v[lt] ];
  
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt],em2.o$b[lt],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt],em1.o$b[lt],lower.tail=FALSE);
  nbeta2.v[selMR.idx] <- q.v;
  
  
  if(doH){ ### if TRUE also correct type2 hemimethylated probes
    ### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
    lt <- 2;
    selH.idx <- c(which(class2.v==lt),selML.idx);
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH;
    #### need to do some patching
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    
    ## new maximum of H probes should be
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
    ## new minimum of H probes should be
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
    ndeltaH <- nmaxH - nminH;
    
    ### perform conformal transformation (shift+dilation)
    ## new_beta_H(i) = a + hf*(beta_H(i)-minH);
    hf <- ndeltaH/deltaH ;
    ### fix lower point first
    nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);
  }
  pnbeta.v <- beta.v;
  pnbeta.v[type1.idx] <- beta1.v;
  pnbeta.v[type2.idx] <- nbeta2.v;
  
  return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));
}

