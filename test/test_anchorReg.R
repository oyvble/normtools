


#CREATE SIM DATASET

#CpG
p = 10
n = 1000
n_train = n*0.1
x <- as.data.frame(matrix(data = rnorm(n*p),nrow = n,ncol = p))
anchor <- as.data.frame(matrix(data = sample(c(0,1), replace=TRUE,size=n ),ncol = 1))
colnames(anchor) <- c('X1')
gamma <- 2 #Anchor penalty
target_variable <- 'V1'#the Y
x[,target_variable] = x[,2] + 3*anchor + rnorm(n)

plot(x[,2],x[,1]);points(x[anchor==1,2],x[anchor==1,1],col=2)

#anchor_regression(x, anchor, gamma, target_variable)

#Store full matrix
x0 = x
anchor0 = anchor
trainind = 1:n_train
x = x[trainind,]
anchor = anchor[trainind,]


#Script approach
fit_const <- lm(x ~ 1); #=colMeans(x)
fit <- lm(x ~ anchor); 
cv_data <- fit_const$fitted.values + fit$residuals

#My approach
designmat = cbind(1,anchor)
fit_const2 <-colMeans(x)
fit2 = normtools::margResponseModel(x,designmat,pvalReturn = FALSE,residReturn = TRUE)
cv_data2 = t(fit_const2 + t(fit2$resid)) #Const + resid (care must be taken)
fit2$resid=NULL;gc() #remove residuals


#max(abs(fit2$coef-fit$coef))
#max(abs(cv_data-cv_data2))
indices <- 1:nrow(cv_data)
j <- match(target_variable, colnames(cv_data))
fit_glmnet_lasso <- cv.glmnet(x = cv_data[indices, -c(j)], cv_data[indices, j])
lambda_cv <- fit_glmnet_lasso$lambda.1se
plot(fit_glmnet_lasso)

#anchor_data <- fit_const$fitted.values + fit$residuals + sqrt(gamma) * (fit$fitted.values - fit_const$fitted.values)
anchor_data <- cv_data + sqrt(gamma) * (fit$fitted.values - fit_const$fitted.values)
anchor_data2 <- cv_data2 + sqrt(gamma) * t( t(designmat%*%fit2$coef)  - fit_const2)
#max(abs(anchor_data-anchor_data2))
