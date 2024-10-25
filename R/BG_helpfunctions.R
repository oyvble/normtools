#https://github.com/schalkwyk/wateRmelon/blob/master/R/bgeqot.R

wm_dasen <- function(meth,fudge=100, samplerm=NULL,siteUse=NULL){
   #meth is output from ewastools::read_idats

   # Robust solution
   #https://github.com/schalkwyk/wateRmelon/blob/master/R/dfsfit.R
   dfsfit <- function(mn, onetwo,roco = substring(colnames(mn), regexpr("R0[1-9]C0[1-9]", colnames(mn))), ...){
      
      #https://github.com/schalkwyk/wateRmelon/blob/master/R/dfs2.R
      dfs2 <- function(x, onetwo){  # x is a matrix of intensities
            # onetwo is a character vector 
            # of same order and length 
            # indicating assay I or II 
            one <- density(x[onetwo=='I'], na.rm=T, n = 2^15, from = 0, to = 5000)
            two <- density(x[onetwo=='II'],na.rm=T, n = 2^15, from = 0, to = 5000)
            one$x[which.max(one$y)] - two$x[which.max(two$y)]
         }
      
      mdf<-apply(mn,2,dfs2,onetwo)
      
      if (! is.null(roco) ) {
         scol  <- as.numeric(substr(roco,6,6))
         srow  <- as.numeric(substr(roco,3,3))
         fit   <- try(  lm( mdf ~ srow + scol ), silent=TRUE) 
         if (! inherits (fit, "try-error") ) {mdf   <- fit$fitted.values}
         else { message ('Sentrix position model failed, skipping') }
      }
      otcor <-  matrix(
         rep( mdf, sum(onetwo=='I')),
         byrow=T, 
         nrow=sum(onetwo=='I')
      )
      mn[onetwo=='I',] <- mn[onetwo=='I',] - otcor
      mn
   }
   
   #CONTINUE
   onetwo = ifelse(meth$manifest$channel=="Both",'II','I')
   #table(onetwo)
   mns = meth$M #methylated matrix
   uns = meth$U #unmethylated matrix

   sampleNames = meth$meta[,1][[1]] #obtain sample names
   siteNames = meth$manifest$probe_id #obtain site names
   if(!is.null(samplerm) && length(samplerm)>0) {
      sampleUse = setdiff(1:length(sampleNames),samplerm)
      mns = mns[,sampleUse,drop=F]
      uns = uns[,sampleUse,drop=F]
      sampleNames = sampleNames[sampleUse]
   }
   if(!is.null(siteUse) && length(siteUse)>0) {
      siteKeep = siteNames%in%siteUse
      mns = mns[siteKeep,,drop=F]
      uns = uns[siteKeep,,drop=F]
      onetwo = onetwo[siteKeep]
      siteNames = siteNames[siteKeep]
   }
   
   mnsc <- dfsfit(mns,  onetwo)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   mnsc[onetwo=='I' ,] <- limma::normalizeQuantiles(mnsc[onetwo=='I', ])
   mnsc[onetwo=='II',] <- limma::normalizeQuantiles(mnsc[onetwo=='II',])
   
   unsc[onetwo=='I' ,] <- limma::normalizeQuantiles(unsc[onetwo=='I', ])
   unsc[onetwo=='II',] <- limma::normalizeQuantiles(unsc[onetwo=='II',])
   beta <- mnsc/( mnsc + unsc + fudge )
   colnames(beta) = sampleNames
   rownames(beta) = siteNames #insert site name
   
   return( beta )
}

getBeta_Raw = function (raw, fudge = 100, samplerm=NULL, siteUse=NULL)  {
   if (!all(c("manifest", "M", "U", "meta") %in% names(raw))) stop("Invalid argument")
   with(raw, {
      M[M < 1] = 1
      U[U < 1] = 1
      meth = M/(M + U + fudge)
      rownames(meth) = manifest$probe_id
      colnames(meth) = meta$sample_id
      
      if(!is.null(samplerm) && length(samplerm)>0) {
         sampleUse = setdiff(1:ncol(meth),samplerm)
         meth = meth[,sampleUse,drop=F]
      }
      if(!is.null(siteUse) && length(siteUse)>0) {
         meth = meth[rownames(meth)%in%siteUse,,drop=FALSE]
      }
      return(meth)
   })
}

