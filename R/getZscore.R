#Example:
#n = 1000
#p = 100000
#X = matrix( rnorm(n*p),ncol=n) #sample per column
#Zscore = getZscore(X)

#Transposed version:
#Zscore =  getZscore(t(X),transposed=T) #sample given per-row

#Also includes ROBUST version based on MEDIAN and MAD
getZscore = function(X,mad2sd=1.4826,robust=FALSE,transposed=FALSE) {
 # X (p x n) Data matrix (samples are given per column)
 # transposed whether samples are given per row. See X argument
 if(transposed) X <- t(X) #transpose matrix (keep)

 #EFFICIENT MEDIAN CALCULATIONS:
 if(robust) {
   bool = require(robustbase)
   if(!bool) stop("Install robustbase R-package")
 
   CENTER = robustbase::rowMedians(X) #get Median for each sites
   ABSDIFF_MED = abs( X - CENTER ) #obtain values subtracted with median and taking absolute
   MAD = robustbase::rowMedians(ABSDIFF_MED) #get Median of absolute deviations 
   rm(ABSDIFF_MED);gc()
   SCALE = mad2sd*MAD #robust measure of SD
 } else {
   ntot = ncol(X) #total number of samples
   CENTER = rowSums(X)/ntot #get average per site
   sumsq = rowSums(X^2)
   SCALE = sqrt((sumsq - ntot*CENTER^2)/(ntot-1)) #get scaling
   rm(sumsq);gc()
 }
 
 #Post transformation:
 X = (X-CENTER)/SCALE #center
 if(transposed) X <- t(X) #transpose matrix back to same as input
 return(X)
}