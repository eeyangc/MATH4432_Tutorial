# A linear regression function
# input: n by p design matrix X; n - vector y
# output: coefficient estimates beta; 
#         residual variance estimate sig2; 
#         standard errors of beta; 
#         t-statistics of beta; 
#         p-values of beta
#         Rsquare: 1- RSS/TSS 
#         RSS: residual sum of squares; TSS: Total sum of squares
linReg <- function(X,y){
  n <- nrow(X)
  p <- ncol(X)
  
  X <- cbind(1,X)
  invK <- solve(t(X)%*%X)
  beta <- invK%*%(t(X)%*%y)
  residual <- y-X%*%beta
  
  Rsq <- 1-sum(residual^2)/sum( (y-mean(y))^2 )
  sig2 <- sum(residual^2) / (n-p-1)
  Sig_beta <- sig2*invK
  se <- sqrt(diag(Sig_beta))
  t <- beta/se
  pval <- pt(abs(t),n-p-1,lower.tail = F)*2
  
  return(list(beta=beta,sig2=sig2,se=se,t=t,pval=pval,Rsq = Rsq))
}