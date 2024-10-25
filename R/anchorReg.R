#target_variable=1; gamma=2; lambda = "CV"; alpha=0.5; useIntercept=TRUE
anchorReg <- function (X, anchor, target_variable=1, gamma=2, lambda = "CV", alpha=0.5, useIntercept=TRUE, returnCV=FALSE) {
  require(glmnet,quietly = TRUE)
  if (ncol(X) < 3) {
    print("unsufficient number of columns")
  }
  if(class(X)[1]!="matrix") X <- as.matrix(X)
  if(class(anchor)[1]!="matrix") anchor <- as.matrix(anchor)
  #fit_const <- lm(x ~ 1)
  #fit <- lm(x ~ anchor)
  #cv_data <- fit_const$fitted.values + fit$residuals
  #anchor_data <- fit_const$fitted.values + fit$residuals + sqrt(gamma) * (fit$fitted.values - fit_const$fitted.values)
  
  if(useIntercept) anchor = cbind(1,anchor) #insert intercept
  fit = normtools::margResponseModel(X,anchor, pvalReturn = FALSE,residReturn = TRUE)
  fit_const <- colMeans(X)
  
  cv_data = t( fit_const + t(fit$resid) ) #Const + resid (care must be taken)
  fit$resid=NULL;gc() #remove residuals (no more requireed)
  anchor_data <-  cv_data + sqrt(gamma) *t( t(anchor%*%fit$coef) - fit_const ) #transpose back
  j = target_variable #copy position of target variable
  if(!is.numeric(target_variable)) j <- match(target_variable, colnames(X))
  
  if (lambda == "CV") {
    fit_glmnet_cv <- cv.glmnet(x = cv_data[, -j,drop=FALSE], y=cv_data[, j],alpha=alpha)
    if(returnCV) return(fit_glmnet_cv)
    lambda_cv <- fit_glmnet_cv$lambda.1se
  } else {
    lambda_cv = lambda
  }
  fit_glmnet_anchor <- glmnet(x = anchor_data[, -j,drop=FALSE], y = anchor_data[, j], lambda = lambda_cv,alpha=alpha)
  return(fit_glmnet_anchor)
}
