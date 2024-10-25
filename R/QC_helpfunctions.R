#Helpfunction which returns the number of beads
checkBeadCounts = function(meth, nreq=3, sites=NULL) {
  
  if(is.null(sites)) sites= meth$manifest$probe_id #all sites are considered
  siteInds = match(sites,meth$manifest$probe_id) #obtain index of site positions (to loop through)
  mat = (meth$V[siteInds,] < nreq) | (meth$N[siteInds,] < nreq)
  #colnames(mat) = samples
  rownames(mat) = sites
  return(mat)
}


#helpfunction to find outlying beta-distrs based on shape
#thresh = 0.5;n=101
getBetaDistrDev = function(beta, probeType=NULL,n=101) {
  xgrid = seq(0,1,l=n)

  if(is.null(probeType)) {
   dmat = NULL #to store density values for each x
   for(i in 1:ncol(beta)) { #for each sample
     d = density(beta[,i],n = n,from=0,to=1)
     dmat = rbind(dmat,d$y) #add pdf
     #all(d$x==xgrid) Must be true
   }
   dmat_MEDIAN = robustbase::colMedians(dmat) #calculate Median density pdf
   #plot(xgrid,dmat_MEDIAN)  
   dmat_DIST = colSums(abs(t(dmat)-dmat_MEDIAN)^2)/n #obtain error distance over all x (squared sum)
   names(dmat_DIST) = colnames(beta)
  } else {

   dmat1 <- dmat2 <-  NULL #to store density values for each x
   for(i in 1:ncol(beta)) { #for each sample
	beta1 = beta[names(probeType)[probeType=="I"],i]
	beta2 = beta[names(probeType)[probeType=="II"],i]
     d1 = density(beta1,n = n,from=0,to=1)
     d2 = density(beta2,n = n,from=0,to=1)
     dmat1 = rbind(dmat1,d1$y) #add pdf
     dmat2 = rbind(dmat2,d2$y) #add pdf
     #all(d$x==xgrid) Must be true
   }
   dmat_MEDIAN1 = robustbase::colMedians(dmat1) #calculate Median density pdf
   dmat_MEDIAN2 = robustbase::colMedians(dmat2) #calculate Median density pdf
#  matplot(cbind(xgrid,xgrid),cbind(dmat_MEDIAN1,dmat_MEDIAN2),ty="l")
   dmat_DIST1 = colSums(abs(t(dmat1)-dmat_MEDIAN1)^2)/n #obtain error distance over all x (squared sum)
   dmat_DIST2 = colSums(abs(t(dmat2)-dmat_MEDIAN2)^2)/n #obtain error distance over all x (squared sum)
   dmat_DIST = cbind(dmat_DIST1,dmat_DIST2)
   colnames(dmat_DIST) = c("TypeI","TypeII") 
   rownames(dmat_DIST) = colnames(beta)
  }  
  return(dmat_DIST)
}

#helpfunction to plot control metrics (standard from Illumina)
plotControlMetrics = function(meth) {
  #meth Object from ewastools::read_idats
  ctrl = control_metrics(meth)
  
  getRng = function(x) {
    isOK = !(is.infinite(x) | is.na(x))
    return(range(x[isOK]))
  }
  
  showMetric = function(var1,var2=NULL, varNames=NULL,txt="", g=function(x) x) {
    val1 = ctrl[[var1]]
    lim1 = attr(val1,"threshold")
    rng1 = getRng(val1)
    val1[is.infinite(val1)] = rng1[2] #impute with largest observed
    
    if(!is.null(var2)) {
      val2 = ctrl[[var2]]
      lim2 = attr(val2,"threshold")
      rng2 = getRng(val2)
      val2[is.infinite(val2)] = rng2[2] #impute with largest observed
      
      plot(g(val1),g(val2),main=txt,xlab=varNames[1],ylab=varNames[2])
      abline(h=g(lim2),lty=2)
      
      outInd1 = which(val1<lim1)
      outInd2 = which(val2<lim2)
      outInd = union(outInd1,outInd2) #merge
      
      if(length(outInd)>0) {
        isNA1 = is.na(val1[outInd])
        isNA2 = is.na(val2[outInd])
        val1[outInd[isNA1]]  <- mean(val1,na.rm=T) #impute with average where isNA
        val2[outInd[isNA2]]  <- mean(val2,na.rm=T) #impute with average where isNA
        
        points(g(val1[outInd]),g(val2[outInd]),col=2,pch=19)
        text(g(val1[outInd]),g(val2[outInd]),outInd,cex=0.8,pos = 4)
        
        isNA = isNA1 |isNA2 #also Plot with other points where it was NA
        if(any(isNA)) points(g(val1[outInd[isNA]]),g(val2[outInd[isNA]]),col=2,pch=8)
      }
      
    } else {
      d = density(g(val1))
      plot(d,lty=2,xlab=varNames[1],main=txt)
      points(g(val1),rep(0,length(val1)),cex=0.8)
      outInd = which(val1<lim1)
      if(length(outInd)>0) {
        points(g(val1[outInd]),rep(0,length(outInd)),col=2,pch=19)
        text(g(val1[outInd]),rep(0,length(outInd)),outInd,cex=0.8,pos = 3)
      }
    }
    mtext(paste0("threshold=",signif(g(lim1),3)))
    abline(v=g(lim1),lty=2)
  }
  
  glog = function(x) log2(x)
  pdf(paste0("QC_controls.pdf"))
  showMetric(1,NULL,c("value (log)"),"Restoration") #REstoration
  showMetric(2,3,c("Green (log)","Red (log)"),"Staining control",glog) #Staining
  showMetric(4,5,c("Green (log)","Red (log)"),"Extension control",glog) #Extension
  showMetric(6,7,c("High/Medium","Medium/Low"),"Hybridization control",glog) #Hybridization:
  showMetric(8,9,c("TR 1 (log)","TR 2 (log)"),"Target Removal",glog) #Target removal:
  showMetric(10,11,c("Green (log)","Red (log)"),"Bisulfite Conversion I",glog) #Bisulfite Conversion I
  showMetric(12,NULL,c("value (log)"),"Bisulfite Conversion II",glog) #Bisulfite Conversion II
  showMetric(13,14,c("Green (log)","Red (log)"),"Specificity I",glog) #Specificity I
  showMetric(15,var2=NULL,c("value (log)"),"Specificity II",glog) #Specificity II
  showMetric(16,17,c("Green (log)","Red (log)"),"Non-polymorphic",glog) ##Non-polymorphic
  dev.off()
}

#CONTROL PROBE CHECK

#Input is output from 'prcomp' 
getPCAoutlier = function(pca,kk=2,prob=0.95,plott=FALSE) {
  if(kk<2) stop("Cannot be below 2 components!")
  pred1 <- predict(pca)[,1:kk]
  propV <- pca$sdev^2/sum(pca$sdev^2)
  res <- (pred1-colMeans(pred1))
  outl <- sqrt( diag((res)%*%solve(cov(pred1))%*%t(res)) ) #mahalonibis distance)
  outl <- which(pchisq(outl,kk)>prob) #get outlier
  if(plott)  {
    giveAX <- paste0("PC ",1:2," (",round(propV[1:2]*100,2),"%)")
    plot(pred1[,1],pred1[,2],xlab=giveAX[1],ylab=giveAX[2])
    if(length(outl)) {
      points(pred1[outl,1],pred1[outl,2],col=2,pch=19)
      text(pred1[outl,1],pred1[outl,2],outl,cex=0.8,pos = 2)
    }
  }
  return(outl)
}

#helpfunction to get outliers of control probes (PCA) based on ewastools object
getPCAoutlier_control = function(meth,kk=2,prob=0.95,plott=FALSE) {
  ctrl = control_metrics(meth)
  ctrlMAT = matrix(log(unlist(ctrl)),nrow=length(ctrl[[1]]))
  
  isNA = is.na(ctrlMAT)
  typeHasNA = colSums(isNA)>0
  ctrlMAT = ctrlMAT[,!typeHasNA] #remove NAs
  
  pca = prcomp(ctrlMAT, center = TRUE,scale. = TRUE)
  return( getPCAoutlier(pca,kk,prob,plott) )
}


#Slight reimplementation of minfi::plotQC
showQC = function(obj,badSampleCutoff = 10.5,title="",plott=TRUE) {
  #obj must be of class list (ewastools object) or MethylSet (includes getUnmeth,getMeth)
  
  if(class(obj)=="list") { #assuming ewastools object
    uMed <- log2(robustbase::colMedians(obj$U,na.rm = TRUE))
    mMed <- log2(robustbase::colMedians(obj$M,na.rm = TRUE))
  } else {
    uMed <- log2(robustbase::colMedians(minfi::getUnmeth(obj),na.rm = TRUE))
    mMed <- log2(robustbase::colMedians(minfi::getMeth(obj),na.rm = TRUE))
  }
  
  meds <- (mMed + uMed)/2
  whichBad <- which((meds < badSampleCutoff))
  
  #obtain distance of the bad samples (larger the worse)
  distBad = sort(setNames(badSampleCutoff - meds[whichBad],whichBad),decreasing = TRUE) 
  
  if(plott) {
    plot(mMed, uMed, xlim = c(8,14), ylim = c(8,14), xaxt = "n", yaxt = "n", main=title,
         xlab = "Meth median intensity (log2)",ylab = "Unmeth median intensity (log2)", 
         col = ifelse(1:length(meds) %in% whichBad, "red", "black"))
    axis(side = 1, at = c(9,11,13))
    axis(side = 2, at = c(9,11,13))
    abline(badSampleCutoff * 2 , -1, lty = 2)
    if (length(whichBad) > 0) {
      text(mMed[whichBad], uMed[whichBad] - 0.25,labels = whichBad, col = "red")
    }
    legend("topleft", legend = c("good", "bad, with sample index"), pch = 1,col = c("black", "red"), bty = "n")
    if (length(whichBad) > 0) text(mMed[whichBad], uMed[whichBad] - 0.25, labels = whichBad, col = "red")
  }
  return(distBad) #return index of those bad
}


#FOLLOWING CODE EXTRACTED FROM https://github.com/schalkwyk/wateRmelon/blob/master/R/pfilter.R
#getBeadcounts<-function(x){ #obtain number of bead counts for each sites (methylated/unmethylated)
 #x is ewastools object
  

############################
#wateRmelon IMPLEMENTATIONS#
############################

#  x = beta.raw;iqr=TRUE; iqrP=2; pc=1; mv=TRUE; mvP=0.15; plot=TRUE; nsites=NULL

wm_outlyx <- function(x, iqrP=2, pc=1, mvP=0.15, nsites=NULL, plot=TRUE, ...) { 
#MODIDIED CODE WAS EXTRACTED FROM https://rdrr.io/bioc/wateRmelon/src/R/outlyx.R
 ### Computes outliers within methylomic datasets
 # x    : beta matrix (sites as rows, samples as columns).M/(M+U)), no offset. Note: Sex-specific sites should be removed beforehand.
 # pc   : The desired principal component for outlier identification
 # iqr  : Logical, to indicate whether to determine outliers by interquartile ranges.
 # iqrP : The number of interquartile ranges one wishes to discriminate from the upper and lower quartiles. Default = 2, Tukey's rule suggests 1.5
 # mv   : Logical, to indicate whether to determine outliers using distance measures using a modified version of pcout from mvoutlier.
 # mvP  : The threshold bywhich one wishes to screendata based on the final weight output from the pcout functionValue between 0 and 1. Default is 0.15
 # plot : Logical, to indicate if a graphical device to displaysample outlyingness.
###

	#WITHIN HELPFUNCTIONS
	iqrFun <- function(x, pc, iqrP){ # identifying outliers based on Interquantiles
	  quantilesbx <- apply(x, 2, quantile) # fivenum also works
	  IQR <- quantilesbx[4, pc] - quantilesbx[2, pc]
	  thresh <- iqrP*IQR
	  hiOutlyx <- quantilesbx[4,pc] + thresh # Upper threshold
	  loOutlyx <- quantilesbx[2,pc] - thresh # Lower threshold
	  outhi <- x[, pc] > hiOutlyx # Upper Outliers
	  outlo <- x[, pc] < loOutlyx # Lower Outliers
	  return(list(c(rownames(x)[outhi], rownames(x)[outlo]), IQR, "hi" = hiOutlyx, "low" = loOutlyx))
	} # }}}

	mvFunFast <- function(x, mvP, ...){ # {{{
	  explvar=0.99;crit.M1=1/3;crit.c1=2.5; crit.M2=1/4;crit.c2=0.99;cs=0.25;outbound=0.25
	  
	  #Don't assume transposed x matrix    
	  p = nrow(x) #number of sites
	  n = ncol(x) #number of samples
	  
	  #This is slow:
	  #out2 <- mvFun(t(x), mvP=mvP)
	  
	  #Calculating mad per marker:
	  #x.mad=apply( t(x),2,mad) #very slow
	  madConst = 1.4826
	  MEDIAN = robustbase::rowMedians(x) #apply(xt,2,median)
	  x.cen <- x - MEDIAN #center with median (transpose back)
	  x.mad = madConst*robustbase::rowMedians(abs(x.cen))
	  
	  #Check
	  if(0) {
	    stats::mad(x[1,])==madConst*median(abs(x[1,] - median(x[1,])))
	  }
	  
	  #Remove sites which has mad=0
	  siteRemove = which(x.mad==0)
	  if(length(siteRemove>0)) {
	    MEDIAN <- MEDIAN[-siteRemove]
	    x.mad <- x.mad[-siteRemove]
	    x.cen <- x.cen[-siteRemove,,drop=FALSE]
	  }
	  
	  # PHASE 1:
	  # Step 1: robustly sphere the data:
	  x.sc <- x.cen/x.mad#scale(xt,MEDIAN ,x.mad)
	  #    sum(is.na(x.sc))
	  # Step 2: PC decomposition; compute p*, robustly sphere:
	  MEAN = rowSums(x.sc)/n #get mean per site
	  x.svd <- svd(t( x.sc - MEAN)) #center matrix (but don't scale)
	  a <- x.svd$d^2/(n-1)
	  p1 <- (1:p)[(cumsum(a)/sum(a)>explvar)][1] #get number of prinsipals to use (decided by explvar)
	  
	  x.pc <- t(x.sc)%*%x.svd$v[,1:p1]
	  xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))
	  
	  # Step 3: compute robust kurtosis weights, transform to distances:
	  wp <- abs(apply(xpc.sc^4,2,mean)-3)
	  
	  xpcw.sc <- xpc.sc%*%diag(wp/sum(wp))
	  xpc.norm <- sqrt(apply(xpcw.sc^2,1,sum))
	  x.dist1 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)
	  
	  # Step 4: determine weights according to translated biweight:
	  M1 <- quantile(x.dist1,crit.M1)
	  const1 <- median(x.dist1)+crit.c1*mad(x.dist1)
	  w1 <- (1-((x.dist1-M1)/(const1-M1))^2)^2
	  w1[x.dist1<M1] <- 1
	  w1[x.dist1>const1] <- 0
	  
	  #
	  # PHASE 2:
	  # Step 5: compute Euclidean norms of PCs and their distances:
	  xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
	  x.dist2 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)
	  
	  # Step 6: determine weight according to translated biweight:
	  M2 <- sqrt(qchisq(crit.M2,p1))
	  const2 <- sqrt(qchisq(crit.c2,p1))
	  w2 <- (1-((x.dist2-M2)/(const2-M2))^2)^2
	  w2[x.dist2<M2] <- 1
	  w2[x.dist2>const2] <- 0
	  #
	  # Combine PHASE1 and PHASE 2: compute final weights:
	  # Changed output slightly
	  pcoutbetx  <- (w1+cs)*(w2+cs)/((1+cs)^2)
	  
	  outmvbetx <- pcoutbetx < mvP #pcoutbetx$wfinal < mvP
	  return(list(rownames(as.matrix(pcoutbetx[outmvbetx])),pcoutbetx))
	} #end mvFun 

 #SCRIPT STARTS:	
	if(!is.null(nsites) && nsites<=nrow(x)) x = x[sample(1:nrow(x),nsites),,drop=FALSE]
	
  # Initialising objects to be used:
  df <- list() # Converted into dataframe later on.
  
  # Removing probes with NA values.
  x <- na.omit(x)

  #obtain normalized version of data (fast) using normtools::scaling
  pccompbetx.rot <- svd(scale(x,scale=FALSE) , nu = 0)$v #perform svd on sample-centered data (sites are repeats). Use only rotation matrix
  #pccompbetx.rot <- svd( t(x) - rowMeans(x)) , nu = 0)$v #perform svd on sample-centered data (sites are repeats). Use only rotation matrix
  #max(abs(prcomp(x, retx=FALSE)$rot - pccompbetx.rot)) #check if same calculation is obtaained
  rownames(pccompbetx.rot) = colnames(x) #insert sample name to be recognized in iqrFun
  
#BLOCK iqr (fast)
  #out1check <- iqrFun(prcomp(x, retx=FALSE)$rot, pc=pc, iqrP=iqrP) #old code
  out1 <- iqrFun(pccompbetx.rot, pc=pc, iqrP=iqrP) #input is same as prcomp$rotation
  v1 <- colnames(x) %in% out1[[1]]==TRUE #indicate sample outliers
  df[["iqr"]] <- v1
  low <- min(c(min(pccompbetx.rot[,pc]),out1[["low"]]))-out1[[2]] # Used if Plot=T
  upp <- max(c(max(pccompbetx.rot[,pc]),out1[["hi"]]))+out1[[2]]  # Used if Plot=T

#BLOCK mv: Calls to mvFun -> pcouted
  out2 <- mvFunFast(x, mvP=mvP)
  v2 <- colnames(x) %in% out2[[1]]==TRUE
  df[["mv"]] <- v2

if(plot) {
  plot(pccompbetx.rot[,pc],out2[[2]],xlim=c(low, upp),xlab="Transformed Betas",ylab="Final Weight", ...)
  abline(v=c(out1[["low"]],out1[["hi"]]),h=mvP,lty=2)
  rect(xleft=low,xright=out1[["low"]],ybottom=0,ytop=mvP,col="red",density=10)
  rect(xleft=out1[["hi"]],xright=upp,ybottom=0,ytop=mvP,col="red",density=10)
}

 # Possibly a better way to do the below.
  df[["outliers"]] <- df[["mv"]]==TRUE & df[["iqr"]]==TRUE
  df[["x"]] = pccompbetx.rot[,pc] #store x-value
  df[["y"]] = out2[[2]] #store y-value
  df <- data.frame(df)
  rownames(df) <- colnames(x)
  return(df)
}  #end wm_outlyx



wm_qual <- function(norm,raw, indPlott=NULL, sort=FALSE){ # normalized and original betas (M/(M+U))
  if(!all(dim(norm)==dim(raw))) stop("The matrices were not compatible.")
  dif  <- norm - raw
  rmsd <- sqrt(colMeans(dif^2, na.rm=TRUE)) #euclidian distance (root of mean squared difference over all sites)
  sdd  <- apply(dif, 2, sd, na.rm=TRUE)  #sd of difference over all sites (per sample)
  sadd  <- apply(abs(dif), 2, sd, na.rm=TRUE)  #sd of absolute difference over all sites (per sample)
  srms <- rmsd/sdd
  df = data.frame(rmsd,sdd,sadd,srms) 
  
  if(!is.null(indPlott)) {
    ylab = switch(indPlott, "rmsd","sdd","sadd","srms")
    xlab = "Index"
    xvals = df[,indPlott]
    if(sort) {
      xvals = sort(xvals)
      xlab = paste0(xlab," (sorted)")
    }
    plot(xvals,ylab=ylab,xlab=xlab)
  }
  return(df)
}

# Author: Tyler Gorrie-Stone, tgorri@essex.ac.uk (Revision Date: 07-01-2016)
#Modified by Oyvind BLeka for ultra-fast calculation
wm_pwod <- function(object, mul=4, verbose=TRUE){
  #pwod(getBeta(x), mul)
    # 'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges. -- probable low MAF SNP heterozygotes --
  # Arguments:
  #  object   : Matrix of Betas (sites as rows, samples as columns)
  #  mul     : Number of interquartile ranges to determine outlying probes. Default = 4. To exclude the very obvious.
  # Returns  : Matrix of Betas with outlying probes coerced to NA.
  quan = normtools::rowFivenum(object) #obtain quantiles
  #quan <- t(apply(object, 1, function(x) fivenum(x)))
      
  iqr <- quan[,4] - quan[,2] #interquantile distance per site
  bounds1 <- quan[,4] + (iqr * mul) #obtain upper boundaries (per site)
  bounds2 <- quan[,2] - (iqr * mul) #obtain lower boundaries (per site)
  d <- object > bounds1 | object < bounds2 # Upper Bound is [1], Lower Bound is [2]

  filter <- object
  filter[d] <- NA 
    
  # Calculating total number of outlying probes.
  # To give an idea of what has changed (if anything).
  tot <- sum(is.na(filter)) - sum(is.na(object))
  if(verbose) cat(tot,"probes detected.", "\n") 
  # Output
  return(filter)
}


#Extracted code from  https://github.com/schalkwyk/wateRmelon/blob/master/R/bscon_minfi.R
wm_bscon = function(meth) {
  #meth Object from ewastools::read_idats
  #Obtain site control indices:
  
  bsI_ind = subset(meth$controls, group=="BISULFITE CONVERSION I",select=index)[[1]] #get site index for BS1
  bsII_ind = subset(meth$controls, group=="BISULFITE CONVERSION II",select=index)[[1]] #get site index for BS2

  #selecting only the Bisulfite conversion I values from both green and red
  bsI.green = meth$ctrlG[bsI_ind,]
  bsI.red = meth$ctrlR[bsI_ind,]

  #selecting only the Bisulfite conversion II values from both green and red
  bsII.green = meth$ctrlG[bsII_ind,]
  bsII.red = meth$ctrlR[bsII_ind,]
  
  # calculate BS conv type I betas as an example of using an index vector
  if(nrow(bsI.green) > 11){ # 450K
    BSI.betas <- rbind(bsI.green[1:3,], bsI.red[7:9,])/((rbind(bsI.green[1:3,], bsI.red[7:9,])) + rbind(bsI.green[4:6,], bsI.red[10:12,]))
  } else { # EPIC
    BSI.betas <- rbind(bsI.green[1:2,], bsI.red[6:7,])/((rbind(bsI.green[1:2,], bsI.red[6:7,])) + rbind(bsI.green[3:4,], bsI.red[ 8:9 ,]))
  }
  
  #calculation of BS con in Type II data
  BSII.betas <- bsII.red/(bsII.red + bsII.green)
  
  apply(rbind(BSI.betas, BSII.betas), 2, median)*100 ## this is the value you are interested in

}



#' Compute internal bisulfite conversion control
#' Compute GCT score for internal bisulfite conversion control. 
#' The higher the GCT score, the more likely the incomplete conversion.
sesame_bisConversionControl <- function(meth, use.median=FALSE) {
  #platform = "EPIC"
  
  platform = meth$platform
  if(platform=="450K") platform = "HM450" #change name for compatibility with sesame
  
  #OBtain base extension from package sesameData
  bool = require(sesameData,quietly=TRUE)
  if(!bool) return(NULL)
  obj = sesameData::sesameDataGet(paste0(platform,".probeInfo"))  
  extC <- obj$typeI.extC ##45427
  extT <- obj$typeI.extT ##15148
  #NOT POSSIBLE TO EXTRACT FROM THIS:  
  #extC2 = with(meth$manifest, probe_id[next_base=="C"]) #NOT THIS 
  #extT2 = with(meth$manifest, probe_id[next_base=="T"]) #NOT THIS 
  
  #hm450.manifest.hg38$nextBaseRef
#  table(meth$manifest$next_base)
  prbs = with(meth$manifest, probe_id[channel=="Red"]) #obtain sites which is used in oob (Green)
  # 89203 for HM450K
  extC <- intersect(prbs, extC) #probes to use
  extT <- intersect(prbs, extT)
  indC = with(meth$manifest, OOBi[probe_id%in%extC & channel=="Red"]) #probe index (NA's removed)
  indT = with(meth$manifest, OOBi[probe_id%in%extT & channel=="Red"]) #probe index (NA's removed)
  #meth$oobG$U[indC,1]
  #meth$oobG$M[indT,1]
  oobG = meth$oobG #obtain out-of-band probe green signals
  nsamples = ncol(oobG$M)
  bscon = rep(NA,nsamples) #values
  for(i in 1:nsamples) {
    extCvals = c(oobG$M[indC,i],oobG$U[indC,i]) #collapse M and U
    extTvals = c(oobG$M[indT,i],oobG$U[indT,i]) #collapse M and U
    if (use.median) {
      bscon[i] = median(extCvals, na.rm=TRUE)/median(extTvals, na.rm=TRUE)
    } else {
      bscon[i] = mean(extCvals, na.rm=TRUE)/mean(extTvals, na.rm=TRUE)
    }
  }
  return(bscon)
}

#Pipeline for quality control (returns list of failing samples/sites)
#plott=TRUE;  bsconScoreThresh=85; minCount=3; minBeadCount=3;bsconOOBsignif=0.1; nbeadsignif=0.1;betaDevThresh=c(0.08,0.6)
doQC = function(meth, plott=TRUE, bsconScoreThresh=85, minCount=3, minBeadCount=3, bsconOOBsignif=0.1, nbeadsignif=0.1, betaDevThresh=c(0.08,0.6)) {
  #minCount number of tolerated samples (used for estNmodes)
  #bsconOOBsignif sifnif level to indicate outlier in bisulfit   
  
  plottDensity = function(x,ind,main="",xlab="") {
	y = rep(0,length(x)) #to plot
     plot(density(x),main=main,xlab=xlab,ylab="pdf");
	points(x,y,cex=0.5,pch=19)
     if(length(ind)>0) {
       points(x[ind],y[ind],col=2,pch=19,cex=0.5)
       text(x[ind],y[ind],ind,adj=1,pos=3,col=2,cex=0.5)
     } 
  }

  #base on robust measure of mean/sd
  getZscoreOutlier = function(x,signif=0.1,plott=TRUE,main="",xlab="") {
    #Zscore = (x - mean(x))/sd(x)
    MED = median(x)
    MAD = median(abs(x-MED))
    scale = 1.4826*MAD
    Zscore = (x - MED)/scale #robust variant
    failed <- which( Zscore > qnorm(1 - signif/length(Zscore)))       
    if(plott) plottDensity(x,failed,main=main,xlab=xlab)
    return(failed)
  }
  
  #Obtain sampleID and siteID-
  sampleIDs = meth$meta$sample_id
  probeIDs = meth$manifest$probe_id

  #PART I: METRICS BASED ON Control-probes
  
  #A) Investigate control probes:
  control_metric = control_metrics(meth)
  for(con in names(control_metric)) {
#    con=names(control_metric)[2]
    vals = control_metric[[con]]
    belowInd = which(vals<attr(vals,"threshold"))
    if(length(belowInd)>0) print(paste0(con,":", paste0(belowInd,collapse="/")))
  }
  if(plott) plotControlMetrics(meth) #plot to pdf for inspection
  failed_control = which(sample_failure(control_metric )) #samples failing to control
  
  #Focus score on bisulfitconversion:
  bsconScore = wm_bscon(meth) #obtain bisulfit score, #should be above 80%
  failed_bscon = which(bsconScore<bsconScoreThresh)
  if(plott) {
    pdf("QC_bscon.pdf")
    plottDensity(bsconScore,failed_bscon,xlab="Score",main="Bisulfit conversion score (wateRmelon)");
    abline(v=bsconScoreThresh,col=2,lty=2);abline(v=80,col=2) #should be above 80%
    dev.off()
  }
  #outliers = getPCAoutlier_control(meth,kk=3,prob=0.9,plott = TRUE)
  
  #ALternative score for bisulfitconversion based on OOB (require sesameData package):
  bsconOOB = sesame_bisConversionControl(meth)
  if(!is.null(bsconOOB)) {
    if(plott) pdf("QC_bsconSesame.pdf")
    failed_bsconOOB <- getZscoreOutlier(bsconOOB,bsconOOBsignif,main="Bisulfit conversion score (SeSAMe)",xlab="bisulfitconversion score")
    if(plott) dev.off()
  } else {
    failed_bsconOOB <- as.character() #Not defined
  }
  
  #PART II: METRICS BASED ON Beta-values
  
  #Indicate bad samples using minfi::showQC based on median
  if(plott) pdf("QC_median.pdf")
  failed_MedianQC = showQC(meth, plott=plott) #use all data
  if(plott) dev.off()


  #INVESTIGATE beta-values (after excluding sites)
  beta = ewastools::dont_normalize(meth) #obtain beta-values (no normalization)
  
  ###########################################
  #REMOVE SITES WHICH ARE CROSS-REACTIVE ETC#
  ###########################################
  
  #REMOVE SITES WITH rs (SNP)
  SNPsites = with(meth,manifest$probe_id[manifest$probe_type=="rs"])
  
  #Obtain list of samples to remove (due to X/Y-chromosome, SNP related, Cross-reactive)
  siteremovalFILE =  paste0( path.package("normtools"),"/SiteRemoval450k.csv") # Get file name in package
  siteRemoval = unlist(read.csv(siteremovalFILE,header=F))
  #sum(siteRemoval450k%in%siteOverview$SITE[siteOverview$CHR%in%c("chrX","chrY")])
  
  sitesWithNA = setdiff(  rownames(beta)[rowSums(is.na(beta))>0],c(siteRemoval,SNPsites))
  
  #following sites contained NA.
  print(paste0("Number of sites with any NA=",length(sitesWithNA)))
  
  siteKeepInd = !rownames(beta)%in%c(siteRemoval,SNPsites,sitesWithNA)
  beta = beta[siteKeepInd,] #update
#  sum(is.na(beta))
  #dim(beta)
  
  #SEARCH MULTIMODAL SITES: beta . tolerate 3 out-side counts
  print("Calculating number of modes for each site....")  
  nmodes = normtools::estNmodes(beta,transposed = FALSE, minCount = minCount) 
  
  site_multiMode = rownames(beta)[nmodes>1] #obtain name of sites being multimodal
  print(paste0("Number of multimodal sites found=",length(site_multiMode)))
  
  #table(nmodes)
  beta = beta[nmodes==1,,drop=FALSE] #update beta again
  
  #FINAL LIST OF SITES TO USE IN CONTROL:
  siteUse = rownames(beta) #this is sites to use in the analysis
  print(paste0("Number of sites to evaluate in QC=",length(siteUse)))
  
  #OVERVIEW OF CHROMOSOME POSITIONS
  #sort(table(with(meth, manifest$chr[manifest$probe_id%in%siteUse])))
  
  #C) Calculate outlyx:
  print("Calculating outlyx, this may take a while...")
  if(plott) pdf("QC_outlyx.pdf")
  outlyx = wm_outlyx(beta,plot = plott)
  if(plott) dev.off()
  failed_outlyx = which(outlyx$outliers)
  #which(outlyx$x > -0.05)
  
  #obtain probe type
  probeType = setNames( ifelse(meth$manifest$channel=="Both",'II','I'),probeIDs)

  #D) Indicate outlying global beta-distributions
  #PLOT beta, highlight failed samples
  print("Checking beta distributions...")
  probeTypes = probeType[rownames(beta)]
  dist = getBetaDistrDev(beta,probeTypes  ) #get distance deviation (squared distance against median)
  meanProbeTypes = exp(colMeans(log(dist))) #get average of each probe type
  
  #dist2 = (dist[,1] + dist[,2])/2 #take average over both probe types
  if(length(betaDevThresh)==1) betaDevThresh = rep(betaDevThresh,2)
  failed_betadistr = which(dist[,1]>betaDevThresh[1] | dist[,2]>betaDevThresh[2])
  
  if(plott) {    #QC_betadist plots are build here
    pdf("QC_betadist.pdf",width=10,height=6)

    #Show density plots of beta
    types = sort(unique(probeTypes))
    for(type in types) {
    #type=types[2]
    	indUse =  type==probeTypes 
       plot(density(beta[indUse,1],from=0,to=1),ylim=c(0,5),xlab="beta-values",main="Density of beta values")
    	mtext(paste0("Probe type ",type))
       for(i in 2:ncol(beta)) lines(density(beta[indUse,i],from=0,to=1))
       for(i in failed_betadistr) lines(density(beta[indUse,i]),col=2,lty=2)  
      } #end foreach probe type  
    
    #plot(checkBeta$dist,xlab="sample index",ylab="Mean Squared Error",main="Deviation against global beta-density");abline(h=0.5,col=2);
    plot(dist[,1],dist[,2],xlab="Type 1",ylab="Type 2",main="Mean Squared Error (per Probe type)")
    abline(v=betaDevThresh[1] ,col=2,lty=2);abline(h=betaDevThresh[2] ,col=2,lty=2);
    if(length(failed_betadistr)>0) text(dist[failed_betadistr,1],dist[failed_betadistr,2],failed_betadistr,col=2,cex=0.5,adj=1,pos = 1)
    
  }
  #Look on PCA (relative (high-dimensional) distance between sample:
  #ONLY FOR CONFIRMING beta_failed (2D)
  pca = prcomp(t(beta[sample(1:nrow(beta),10000),]),center=T,scale=T)
  pca_outl1 =  getPCAoutlier(pca,plott = plott,prob = 0.95);mtext("PCA with 95% outliers")
  pca_outl2 =  getPCAoutlier(pca,plott = plott,prob = 0.7);mtext("PCA with 70% outliers")
  if(plott) dev.off()
  
  rm(beta);gc()
  #PART III: Metrics based on OOB and BeadCounts
  
  #E) Obtain number of counted beads (per probe per sample)
  print("Checking nbeads < 3 distributions...")
  bead_failed = checkBeadCounts(meth, nreq=minBeadCount, sites=siteUse)#, samplerm = failed_betadistr)
  colnames(bead_failed) = meth$meta$sample_id #insert sample names
  nsites_beadfail = colSums(bead_failed) #obtain number of failed sites per sample
  #nsamples_beadfail = rowSums(bead_failed) #sort(,decreasing=T)[1:100]
  
  if(plott) pdf("QC_nbeads.pdf")
  failed_nbeads <- getZscoreOutlier(nsites_beadfail,nbeadsignif,plott=plott,main="Number of failed beads per sample",xlab="#failed beads (per sample)")
  if(plott) dev.off()
  
  #plot(nsites_beadfail)
  #plot(nsamples_beadfail)
  
  #Look on uncertainty in beta-values
  #beta_uncertainty(meth,4,3)
  
  #F) CALC P-VALUES with different methods
  #Calc p-val:
  M = meth$M; U = meth$U;iR = which(meth$manifest$channel=="Red");iG = which(meth$manifest$channel=="Grn");i2 = which(meth$manifest$channel=="Both");
  oobR = meth$oobR;oobG = meth$oobG #obtain out-of-band probe signals
  indNEG =  which(meth$controls$group=="NEGATIVE")
  negR = meth$ctrlR[indNEG,];negG = meth$ctrlG[indNEG,]    #obtain negative control probes
 
  rm(meth);gc()

  print("Checking 4 different pvalue approaches...")
  #CONSIDER DIFFERENT APPROACHES FOR CALC p-values
  #detPoob = ewastools::detectionP(meth)$detP #check that it becomes the same

  failed_pval = list()  #store in list
  mets = c("pval_oob", "pval_neg", "pval_ecdfOOB", "pval_ecdfNEG") #Methods
  for(met in mets) {
#   met=mets[1]	
    print(paste0("Calculate P-values for method ",met))
   if(met=="pval_oob") detP = ewastools_detectionP(M,U,iR,iG,i2) #oob based on summit-closest
   if(met=="pval_neg") detP = ewastools_detectionP.neg(M,U,iR,iG,i2,negR,negG) #conventional way (improved minfi) with negs
   if(met=="pval_ecdfOOB") detP = sesame_detectionPecdf(M,U,iR,iG,i2,oobR,oobG) #ecdf version (oob)
   if(met=="pval_ecdfNEG") detP = sesame_detectionPecdf(M,U,iR,iG,i2,negR,negG) #ecdf version (neg)
   colnames(detP) <- sampleIDs 
   #rownames(detP) <- probeIDs 

   failed_pval[[met]] = getFailedPval(detP)$samples
  }

  #Look closer on P-values
  #indUse = probeIDs%in%siteUse #obtain sites to use
  #plot(density(log10(detPoob[indUse,i]),from=0,to=1))
  #plot(density(log10(detPneg[indUse,i]),from=0,to=1))
  #plot(density(log10(detPecdfOOB[indUse,i]),from=0,to=1))
  #plot(density(log10(detPecdfNEG[indUse,i]),from=0,to=1))    
  #look on failed samples
  #failed_Samples = lapply(failed_pval,function(x) x$samples)
  #failed_Sites = lapply(failed_pval,function(x) x$sites)
  
  
  #WRAP UP and provide output lists
  failed_Samples = failed_pval 
  failed_Samples[["medianQC"]] = sampleIDs[as.integer(names(failed_MedianQC))]
  failed_Samples[["control"]] = sampleIDs[failed_control]
  failed_Samples[["bscon"]] = sampleIDs[failed_bscon]
  failed_Samples[["bsconOOB"]] = sampleIDs[failed_bsconOOB]
  failed_Samples[["outlyx"]] = sampleIDs[failed_outlyx]
  failed_Samples[["betadistr"]] = sampleIDs[failed_betadistr]
  failed_Samples[["nBeads"]]  = names(failed_nbeads)
  
  table = NULL
  for(met in names(failed_Samples)) {
   if( length(failed_Samples[[met]])>0) {
     sampleID = match(failed_Samples[[met]],sampleIDs) #obtain index
     new = cbind(met, sampleID , failed_Samples[[met]]) #new table
     table  = rbind(table , new)
   }
  }
  if(!is.null(table)) colnames(table) = c("Method","ID","Sample")

  return( table )

}

