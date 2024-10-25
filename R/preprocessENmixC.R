
##background correction
preprocessENmixC  <- function(rgSet, bgParaEst="oob", dyeCorr="RELIC",QCinfo=NULL, exQCsample=TRUE,exQCcpg=TRUE, exSample=NULL, exCpG=NULL, nCores=2) {
  
  #This is a wrapper function calling another wrapper function enmix_adj which again calls functions implemented in C
  enmix <- function(meth,bg,bgParaEst,nCores) {
    colnm <- colnames(meth) #keep copy
    meth.o <- foreach(i=1:ncol(meth),.combine=cbind,.export=c("EM_estimateC","new_cm","enmix_adjC")) %dopar% {
      i=i;enmix_adjC(meth[,i],bg[i,],bgParaEst)}
    if(is.matrix(meth.o)){if(sum(is.na(meth.o))>0){stop("Computation ran out of memory, try to set nCores with a smaller value")}}else{
      stop("Computation ran out of memory, try to set nCores with a smaller value")}
    colnames(meth.o)=colnm
    gc(); meth.o
  }
  
#################################################################################
  if(is(rgSet, "rgDataSet")) {
      if(!is.null(QCinfo)){exSample=unique(c(QCinfo$badsample, exSample))}
      exSample=exSample[exSample %in% colnames(rgSet)]
      if(length(exSample)>0){
        rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
        cat(length(exSample), " samples were excluded before ENmix correction\n")
      }
      mdat <- getmeth(rgSet)
    }else if(is(rgSet, "RGChannelSet")){
      if(!is.null(QCinfo)){exSample=unique(c(QCinfo$badsample, exSample))}
      exSample=exSample[exSample %in% colnames(rgSet)]
      if(length(exSample)>0){
        rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
        cat(length(exSample), " samples were excluded before ENmix correction\n")
      }
      mdat <- preprocessRaw(rgSet)
    }else if(is(rgSet, "methDataSet") | is(rgSet, "MethylSet")){
      if(!is.null(QCinfo) & exQCsample){exSample=unique(c(QCinfo$badsample,exSample))}
      exSample=exSample[exSample %in% colnames(rgSet)]
      if(length(exSample)>0){
        rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
        cat(length(exSample), " samples were excluded before ENmix correction\n")
        bgParaEst
      }
      mdat=rgSet; bgParaEst="est"; 
      if(!(dyeCorr=="none")){cat("Warning: Input data need to be a rgDataSet or RGChannelSet to perform dye-bias correction\n");
      cat("Warning: dye-bias correction will not be performed\n")}
      dyeCorr="none"
    }else{stop("Error: object needs to be of class 'RGChannelSet' or 'MethylSet'")}
  
    if(nCores>detectCores()) {
      nCores=detectCores(); 
      cat("Only ",detectCores(), " cores avialable, nCores was reset to ",detectCores(),"\n")
    }
    if(!is.null(QCinfo) & exQCcpg) {exCpG=unique(c(exCpG, QCinfo$badCpG))}
    exCpG=exCpG[exCpG %in% rownames(mdat)]
    if(length(exCpG)>0){
      mdat=mdat[!(rownames(mdat) %in% exCpG),]
      cat(length(exCpG), " CpGs were excluded before ENmix correction\n")
    }
    rm(QCinfo)
    
    if(is(mdat, "methDataSet")){
      probe_type=rowData(mdat)$Infinium_Design_Type
      col=rowData(mdat)$Color_Channel
      probe_type[probe_type %in% c("I","snpI") & col=="Grn"]="IGrn"
      probe_type[probe_type %in% c("I","snpI") & col=="Red"]="IRed"
      probe_type[probe_type %in% c("snpII")]="II"
    }else if(is(mdat, "MethylSet")){
      probe_type <- getProbeType(mdat, withColor=TRUE)
    }
    
    cat("Analysis is running, please wait...!\n")
    
    ##estimate background parameters
    if(bgParaEst == "neg" | bgParaEst == "subtract_neg") {
      if(is(rgSet,"rgDataSet")){
        ctrls <- getCGinfo(rgSet,type="ctrl")
      }else if(is(rgSet,"RGChannelSet")){
      ctrls <- getProbeInfo(rgSet,type="Control")
      }
      ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
      ctrl_address <- as.vector(ctrls$Address[ctrls$Type %in% "NEGATIVE"])

      ctrl_r <- assays(rgSet)$Red[ctrl_address,]
      ctrl_g <- assays(rgSet)$Green[ctrl_address,]
      ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06
      temp <- apply(ctrl_r,2,huber_mus)
      mu <- temp["mu",];sigma <- temp["s",]
      bgRI <- as.data.frame(cbind(mu,sigma))
      temp <- apply(ctrl_g,2,huber_mus);
      mu <- temp["mu",];sigma <- temp["s",]
      bgGI <- as.data.frame(cbind(mu,sigma))
      bgRII <- bgRI;bgGII <- bgGI
    } else if (bgParaEst == "oob" | bgParaEst == "subtract_oob") {
      if(is(rgSet,"rgDataSet")) {
        I_probe <- getCGinfo(rgSet,type="I")
        I_probe=I_probe[I_probe$Color_Channel=="Red",]
      } else if(is(rgSet,"RGChannelSet")) {
        I_probe <- getProbeInfo(rgSet, type="I-Red")
      }
      I_green_bg_M <-  assays(rgSet)$Green[I_probe$AddressB,]
      I_green_bg_U <-  assays(rgSet)$Green[I_probe$AddressA,]
      ctrl_g <- rbind(I_green_bg_M,I_green_bg_U)
      
      if(is(rgSet,"rgDataSet")){
        I_probe <- getCGinfo(rgSet,type="I")
        I_probe=I_probe[I_probe$Color_Channel=="Grn",]
      }else if(is(rgSet,"RGChannelSet")){
        I_probe <- getProbeInfo(rgSet, type="I-Green")
      }
      I_red_bg_M <- assays(rgSet)$Red[I_probe$AddressB,]
      I_red_bg_U <- assays(rgSet)$Red[I_probe$AddressA,]
      ctrl_r <- rbind(I_red_bg_M,I_red_bg_U)
      ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06
      temp <- apply(ctrl_r,2,huber_mus)
      mu <- temp["mu",];sigma=temp["s",]
      bgRI <- as.data.frame(cbind(mu,sigma))
      temp <- apply(ctrl_g,2,huber_mus);
      mu <- temp["mu",];sigma=temp["s",]
      bgGI <- as.data.frame(cbind(mu,sigma))
      bgRII <- bgRI;bgGII=bgGI
      rm(list=c("I_green_bg_M","I_green_bg_U","ctrl_g","I_red_bg_M","I_red_bg_U","ctrl_r"))
    } else if (bgParaEst == "subtract_q5neg") {
      if(is(rgSet,"rgDataSet")){
          ctrls <- getCGinfo(rgSet,type="ctrl")
      }else if(is(rgSet,"RGChannelSet")) {
          ctrls <- getProbeInfo(rgSet,type="Control")
      }
      ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
      ctrl_address <- as.vector(ctrls$Address[ctrls$Type %in% "NEGATIVE"])
      ctrl_r <- assays(rgSet)$Red[ctrl_address,]
      ctrl_g <- assays(rgSet)$Green[ctrl_address,]
      ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06  ## may not need this
      mu <- apply(ctrl_r,2,function(x) quantile(x,probs=0.05,na.rm=TRUE));
      sigma <- apply(ctrl_r,2,function(x)sd(x,na.rm=TRUE));
      bgRI <- as.data.frame(cbind(mu,sigma))
      mu <- apply(ctrl_g,2,function(x) quantile(x,probs=0.05,na.rm=TRUE));
      sigma <- apply(ctrl_g,2,function(x)sd(x,na.rm=TRUE));
      bgGI <- as.data.frame(cbind(mu,sigma))
      bgRII <- bgRI;bgGII=bgGI
    } else if(bgParaEst == "est" | bgParaEst == "subtract_estBG") {
      
        mdat_subset <- mdat[probe_type == "IRed",]
        m_I_red <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
        mdat_subset <- mdat[probe_type == "IGrn",]
        m_I_grn <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
        mdat_subset <- mdat[probe_type == "II",]
        mII <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
        rm(mdat_subset)
        bgRI <- as.data.frame(t(apply(m_I_red,2,estBG)));names(bgRI) <- c("mu","sigma","perc")
        bgGI <- as.data.frame(t(apply(m_I_grn,2,estBG)));names(bgGI) <- c("mu","sigma","perc")
        bgII <- as.data.frame(t(apply(mII,2,estBG)));names(bgII) <- c("mu","sigma","perc")
        
        ##empirically adjusting the estimates
        pp <- apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,max) - apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,min)
        avgp <- apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,mean)
        for(i in 1:nrow(bgGI)){
          if(pp[i]>=0.04){
            bgRI$mu[i] <- quantile(m_I_red[,i],probs=avgp[i])
            bgGI$mu[i] <- quantile(m_I_grn[,i],probs=avgp[i])
            bgII$mu[i] <- quantile(mII[,i],probs=avgp[i])
            bgRI$perc[i] <- avgp[i];bgGI$perc[i]=avgp[i];bgII$perc[i]=avgp[i]
          }
        }
        bgRI <- bgRI[,c("mu","sigma")]
        bgGI <- bgGI[,c("mu","sigma")]
        bgII <- bgII[,c("mu","sigma")]
        A1=sum(bgII$mu)*2/(sum(bgRI$mu)+sum(bgGI$mu))
        A2=sum(bgII$sigma)*2/(sum(bgRI$sigma)+sum(bgGI$sigma))
        bgGII=bgGI;bgGII$mu=bgGII$mu*A1;bgGII$sigma=bgGII$sigma*A2
        bgRII=bgRI;bgRII$mu=bgRII$mu*A1;bgRII$sigma=bgRII$sigma*A2
        rm(list=c("m_I_red","m_I_grn","mII"))
    }
    
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)

    if(dyeCorr == "mean"){
      if(is(rgSet,"rgDataSet")){
          ctrls <- getCGinfo(rgSet,type="ctrl")
      }else if(is(rgSet,"RGChannelSet")){
          ctrls <- getProbeInfo(rgSet,type="Control")
      }
      ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
      ctrl_r <- assays(rgSet)$Red[ctrls$Address,]
      ctrl_g <- assays(rgSet)$Green[ctrls$Address,]
      CG.controls <- ctrls$Type %in% c("NORM_C", "NORM_G")
      AT.controls <- ctrls$Type %in% c("NORM_A", "NORM_T")

      cg_grn=ctrl_g[CG.controls,]
      at_red=ctrl_r[AT.controls,]
      cg_grn <- enmix(cg_grn,bgGI,bgParaEst,nCores)
      at_red <- enmix(at_red,bgRI,bgParaEst,nCores)

      Green.avg <- apply(cg_grn,2,huber_mu)
      Red.avg <- apply(at_red,2,huber_mu)
      ref <- mean(c(Red.avg,Green.avg))
      Grn.factor <- ref/Green.avg
      Red.factor <- ref/Red.avg
    }else if(dyeCorr =="RELIC"){
      if(is(rgSet,"rgDataSet")){
          ctrls <- getCGinfo(rgSet,type="ctrl")
      }else if(is(rgSet,"RGChannelSet")){
        ctrls <- getProbeInfo(rgSet,type="Control")
      }
      ctrls<-ctrls[ctrls$Address %in% rownames(rgSet),]
      ctrl_r<-assays(rgSet)$Red[ctrls$Address,]
      ctrl_g<-assays(rgSet)$Green[ctrls$Address,]
      CG.controls<-ctrls$Type %in% c("NORM_C","NORM_G")
      AT.controls<-ctrls$Type %in% c("NORM_A","NORM_T")
      cg_grn<-ctrl_g[CG.controls,];rownames(cg_grn)=ctrls$ExtendedType[CG.controls]
      at_red<-ctrl_r[AT.controls,];rownames(at_red)=ctrls$ExtendedType[AT.controls]
      cg_grn <- enmix(cg_grn,bgGI,bgParaEst,nCores)
      at_red <- enmix(at_red,bgRI,bgParaEst,nCores)
    }
    rm(rgSet)

    methData <- assays(mdat)$Meth
    N=ceiling(ncol(methData)/(nCores*10))
    parts=rep(1:N,each = ceiling(ncol(methData)/N))[1:ncol(methData)]
    for(i in 1:N){
      id=which(parts==i)
      methD=methData[,id]
      methD[probe_type == "IGrn",] <- enmix(methD[probe_type == "IGrn",],bgGI[id,], bgParaEst, nCores)
      methD[probe_type == "IRed",] <- enmix(methD[probe_type == "IRed",],bgRI[id,], bgParaEst, nCores)
      methD[probe_type == "II",] <- enmix(methD[probe_type == "II",],bgGII[id,], bgParaEst, nCores)
      methData[,id]=methD;
    }
    if(dyeCorr == "mean"){
      methData[probe_type == "IGrn",] <- sweep(methData[probe_type == "IGrn",],2, FUN="*", Grn.factor)
      methData[probe_type == "II",] <- sweep(methData[probe_type == "II",], 2,FUN="*", Grn.factor)
      methData[probe_type == "IRed",] <- sweep(methData[probe_type == "IRed",],2, FUN="*", Red.factor)
    }
    assays(mdat)$Meth <- methData
    rm(methData)

    unmethData <- assays(mdat)$Unmeth
    for(i in 1:N){
      id=which(parts==i)
      unmethD=unmethData[,id]
      unmethD[probe_type == "IGrn",] <- enmix(unmethD[probe_type == "IGrn",],bgGI[id,], bgParaEst, nCores)
      unmethD[probe_type == "IRed",] <- enmix(unmethD[probe_type == "IRed",],bgRI[id,], bgParaEst, nCores)
      unmethD[probe_type == "II",] <- enmix(unmethD[probe_type == "II",],bgRII[id,], bgParaEst, nCores)
      unmethData[,id]=unmethD;
    }
    if(dyeCorr == "mean"){
      unmethData[probe_type == "IGrn",] <- sweep(unmethData[probe_type =="IGrn",], 2, FUN="*", Grn.factor)
      unmethData[probe_type == "IRed",] <- sweep(unmethData[probe_type =="IRed",], 2, FUN="*", Red.factor)
      unmethData[probe_type == "II",] <- sweep(unmethData[probe_type =="II",], 2, FUN="*", Red.factor)
    }
    assays(mdat)$Unmeth <- unmethData
    rm(unmethData)
    stopCluster(c1)
    if(dyeCorr =="RELIC"){mdat=relic(mdat,at_red,cg_grn)}

    annotation=c(paste("Backgroud_corr: ENmix,",bgParaEst,sep=""),paste("dyeBiasCorrection: ",dyeCorr,sep=""))
    if(is(mdat, "methDataSet")){
      metadata(mdat)$preprocessMethod=annotation
    }else if(is(mdat, "MethylSet")){
      mdat@preprocessMethod <- annotation
    }
    mdat
}