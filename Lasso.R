# Ridge regression. Both X and y are not scaled/centerized
Lasso <- function(X, y, lambda=NULL, s=NULL, center = T, scale = T, tol=1e-8, maxIter=1000, method="ACD") {
  
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
  
  #default lambda setting
  if(is.null(lambda)) {
    XX <- t(X)%*%X
    #pre-calculate eigenvalues and eigenvectors
    eigenXX <- eigen(XX)
    D <- eigenXX$values
    U <- eigenXX$vectors
    loglam <- seq(5*log10(max(D)),-0.5*log10(max(D)),length.out = 100)
    lambda <- exp(loglam)
  }
  
  ns <- length(lambda)
  beta <- matrix(0,p,ns)    # solution path for beta
  beta0 <- rep(0,ns)        # solution path for beta0/intercept
  betat <- rep(0,p)         # initialize beta by zero
  obj <- vector("list",ns)
  
  # using gradient descent
  if(method=="GD" | method=="AGD"){
    if(is.null(s)){
      XX <- t(X)%*%X
      eigenXX <- eigen(XX)
      s <- 1/max(eigenXX$values)
    }
    for(i_lam in 1:ns){
      if(method=="GD"){
        fit <- LassoSolver3(X,y,lambda[i_lam],s,betat,tol,maxIter) 
      } else if(method=="AGD"){
        fit <- LassoSolver4(X,y,lambda[i_lam],s,betat,tol,maxIter) 
      }
      betat <- fit$betat
      obj[[i_lam]] <- fit$obj
      
      #Recover the intercept
      if(scale){
        beta0[i_lam] <- ym - t(Xm/Xsd)%*%fit$betat
        beta[,i_lam] <- fit$betat/Xsd
      } else if(center){
        beta0[i_lam] <- ym - Xm%*%fit$betat
        beta[,i_lam] <- fit$betat
      } else{
        beta0[i_lam] <- ym
        beta[,i_lam] <- fit$betat
      }
    }
    
  }
  
  # using coordinate descent
  if(method=="ACD" | method=="CD"){
    Xvar <- colSums(X^2)
    
    if(method=="ACD"){
      mm <- rep(0,p)    # record whether the variable in the active set
      m <- rep(0,p)     # record active variables
      nin <- 0          # initialize the number of all zeros
    }
    
    #Fit the path
    for(i_lam in 1:ns){
      if(method=="ACD"){
        fit <- LassoSolver1(X,y,lambda[i_lam],Xvar,betat,mm,m,nin,tol,maxIter)
        mm <- fit$mm
        m <- fit$m
        nin <- fit$nin
      } else if(method=="CD"){
        fit <- LassoSolver2(X,y,lambda[i_lam],Xvar,betat,tol,maxIter)
      }
      
      betat <- fit$betat
      obj[[i_lam]] <- fit$obj
      
      #Recover the intercept
      if(scale){
        beta0[i_lam] <- ym - t(Xm/Xsd)%*%fit$betat
        beta[,i_lam] <- fit$betat/Xsd
      } else if(center){
        beta0[i_lam] <- ym - Xm%*%fit$betat
        beta[,i_lam] <- fit$betat
      } else{
        beta0[i_lam] <- ym
        beta[,i_lam] <- fit$betat
      }
    }
  }
  
  beta <- list(beta=beta,beta0=beta0,obj=obj,lambda=lambda)
  attr(beta, "class") <- "Lasso"
  return(beta)
}





###################################################### Lasso solver using CD with active set ######################################################
LassoSolver1 <- function(X,y,lam,Xvar,beta_old,mm,m,nin,tol,maxIter){
  n  <- nrow(X)
  p  <- ncol(X)
  
  betat <- beta_old
  yhat <- X%*%betat
  obj_outer <- rep(0,maxIter)
  obj_outer[1] <-  evalObj(y,yhat,betat,lam)
  total_inner <- 0
  
  for(iter_outer in 2:maxIter){
    for(j in 1:p){
      betat[j] <- t(X[,j]) %*% (y-yhat) + Xvar[j]*beta_old[j]
      betat[j] <- softThresh(betat[j],lam) / Xvar[j]
      diff_beta <- beta_old[j] - betat[j]
      
      if(abs(diff_beta)>0){
        yhat <- yhat - X[,j]*diff_beta
        beta_old[j] <- betat[j]
        if(mm[j]==0){
          nin <- nin+1
          mm[j] <- nin
          m[nin] <- j
        }
      }
    }
    obj_outer[iter_outer] <- evalObj(y,yhat,betat,lam)
    if(abs(obj_outer[iter_outer]-obj_outer[iter_outer-1])<tol){
      obj_outer <- obj_outer[1:iter_outer]
      break
    }
    
    # inner loop for active set
    obj_inner <- rep(0,maxIter)
    obj_inner[1] <-  evalObj(y,yhat,betat,lam)
    for(iter_inner in 2:maxIter){
      total_inner <- total_inner + 1
      for(k in 1:nin){
        j <- m[k]
        
        betat[j] <- t(X[,j]) %*% (y-yhat) + Xvar[j]*beta_old[j]
        betat[j] <- softThresh(betat[j],lam) / Xvar[j]
        diff_beta <- beta_old[j] - betat[j]
        if(abs(diff_beta)>0){
          yhat <- yhat - X[,j]*diff_beta
          beta_old[j] <- betat[j]
        }
      }
      obj_inner[iter_inner] <- evalObj(y,yhat,betat,lam)
      if(abs(obj_inner[iter_inner]-obj_inner[iter_inner-1])<tol) break
    }
  }
  cat(" Finished!\t",iter_outer-1," outer loops;\t",total_inner," inner loops\n")
  out <- list(betat = betat, mm=mm, m=m, nin=nin, obj=obj_outer)
}



###################################################### Lasso solver using CD without active set ######################################################
LassoSolver2 <- function(X,y,lam,Xvar,beta_old,tol,maxIter){
  n  <- nrow(X)
  p  <- ncol(X)
  
  betat <- beta_old
  yhat <- X%*%betat
  obj <- rep(0,maxIter)
  obj[1] <-  evalObj(y,yhat,betat,lam)
  
  for(iter in 2:maxIter){
    for(j in 1:p){
      betat[j] <- t(X[,j]) %*% (y-yhat) + Xvar[j]*beta_old[j]
      betat[j] <- softThresh(betat[j],lam) / Xvar[j]
      diff_beta <- beta_old[j] - betat[j]
      
      if(abs(diff_beta)>0){
        yhat <- yhat - X[,j]*diff_beta
        beta_old[j] <- betat[j]
      }
    }
    obj[iter] <- evalObj(y,yhat,betat,lam)
    if(abs(obj[iter]-obj[iter-1])<tol){
      obj <- obj[1:iter]
      break
    }
  }
  cat(" Finished!\t",iter-1," loops\n")
  out <- list(betat = betat, obj=obj)
}


###################################################### Lasso solver using GD ######################################################
LassoSolver3 <- function(X,y,lam,s,betat,tol,maxIter){
  n  <- nrow(X)
  p  <- ncol(X)
  
  yhat <- X%*%betat
  obj <- rep(0,maxIter)
  obj[1] <-  evalObj(y,yhat,betat,lam)
  
  for(iter in 2:maxIter){
    betat <- betat - s * t(X) %*% (yhat-y) 
    betat <- sapply(betat,softThresh,lam=s*lam)#softThresh(betat,s*lam)
    yhat <- X%*%betat
    
    obj[iter] <- evalObj(y,yhat,betat,lam)
    # cat(iter-1," loops: obj = ",obj[iter],"\n")
    if(abs(obj[iter]-obj[iter-1])<tol){
      obj <- obj[1:iter]
      break
    }
  }
  cat(" Finished!\t",iter-1," loops\n")
  out <- list(betat = betat, obj=obj)
}


###################################################### Lasso solver using accelerated GD with momentum ######################################################
LassoSolver4 <- function(X,y,lam,s,beta_old,tol,maxIter){
  n  <- nrow(X)
  p  <- ncol(X)
  
  thetat <- betat <- beta_old
  yhat <- X%*%betat
  obj <- rep(0,maxIter)
  obj[1] <-  evalObj(y,yhat,betat,lam)
  
  for(iter in 2:maxIter){
    betat <- thetat - s * t(X) %*% (yhat-y) 
    betat <- sapply(betat,softThresh,lam=s*lam)#softThresh(betat,s*lam)
    thetat <- betat + (iter-1)/(iter+2) * (betat-beta_old)
    yhat <- X%*%betat
    beta_old <- betat
    
    obj[iter] <- evalObj(y,yhat,betat,lam)
    # cat(iter-1," loops: obj = ",obj[iter],"\n")
    if(abs(obj[iter]-obj[iter-1])<tol){
      obj <- obj[1:iter]
      break
    }
  }
  cat(" Finished!\t",iter-1," loops\n")
  out <- list(betat = betat, obj=obj)
}

###################################################### Soft Threshold function ######################################################
softThresh <- function(b,lam){
  out <- max(0,b-lam) - max(0,-b-lam)
  out
}

###################################################### evaluate objective function ######################################################
evalObj <- function(y,yhat,betat,lam){
  # 1 / (2 * nobs) RSS + lambda * penalty
  n <- length(y)
  out <- 0.5 * sum((y-yhat)^2) + lam*sum(abs(betat))
  out
}