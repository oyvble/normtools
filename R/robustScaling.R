#Example:
#n = 1000
#p = 10000
#X = matrix( rnorm(n*p),ncol=n) #sample per column
#X2 = robustScaling(X) 
#X2 = robustScaling(t(X),transposed=T) #sample given per-row

robustScaling = function(X,mad2sd=1.5,transposed=FALSE,addAttr=TRUE) {
 return(scaling(X,mad2sd=mad2sd,transposed=transposed,addAttr=addAttr,robust=TRUE))
}