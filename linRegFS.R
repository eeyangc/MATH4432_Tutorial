# Forward stagewise regression: Boosting for linear regression
linRegFS <- function(X, y, eps=NULL, center = T, scale = T, tol=1e-8, maxIter=1000){
  
  n  <- nrow(X)
  p  <- ncol(X)
  # ns <- length(lambda)
  
  #Check standardization consistency 
  if(!center & scale){
    warning("Centering is required when 'scale' is TRUE. 'center' is automatically set to be TRUE.")
    center <- TRUE
  }
  
  #Centering, Scaling or do nothing
  if(center){
    # beta <- matrix(0,p,ns)
    ym <- mean(y)
    y  <- y-ym
    Xm <- colMeans(X)
    X <- scale(X,center = T,scale = F)
    if(scale){
      Xsd <- sqrt(colMeans(X^2))
      X <- apply(X,MARGIN = 1,function(X,sd) X/sd,sd=Xsd)
      X <- t(X)
    }
  }
  
  # #default step size
  # if(is.null(eps)) {
  #   XX <- t(X)%*%X
  #   #pre-calculate eigenvalues and eigenvectors
  #   eigenXX <- eigen(XX)
  #   D <- eigenXX$values
  #   U <- eigenXX$vectors
  #   eps <- max(D)/1000
  # }
  
  beta <- matrix(0,p,maxIter)    # solution path for beta
  beta0 <- rep(0,maxIter)        # solution path for beta0/intercept
  betat <- rep(0,p)         # initialize beta by zero
  yhat <- rep(0,n)
  
  for(iter in 1:maxIter){
    r <- y - yhat
    g <- t(X)%*%r
    
    if(is.null(eps)) eps <- 0.01*max(abs(g))/sqrt(sum(r^2))
    idx <- which.max(abs(g))
    betat[idx] <- betat[idx] + eps*sign(g[idx])
    yhat <- yhat + eps*sign(g[idx])*X[,idx]
    beta[,iter] <- betat
    
    #corr <- max(abs(g))/sqrt(sum(r^2)*n)
    corr <- max(cor(X,r))
    cat(iter,"-th iteration;\t updating",idx,"-th coefficient;\t eps:",eps,";\t max corr:",corr,"\n")
    if(corr<tol){
      beta <- beta[,1:iter]
      break
    }
  }
  
  #Recover the intercept
  if(scale){
    beta0 <- ym - t(Xm/Xsd)%*%beta
    beta <- beta/Xsd
  } else if(center){
    beta0 <- ym - Xm%*%beta
  } else{
    beta0 <- ym
  }
  
  beta <- list(beta=beta,beta0=beta0,eps=eps,iteration=iter)
  attr(beta, "class") <- "linRegFS"
  return(beta)
}
