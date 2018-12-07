
G_sl = function(UA, UB, C, D){
  return(UA %*% (( t(UA) %*% C %*% UB)*D) %*% t(UB))
}



S_sl = function(A, lambda, p){
  for(i in 1:(p-1))
    for(j in (i+1):p){
      x = abs(A[i,j])
      A[i,j] = ifelse(x <= lambda, 0, sign(A[i,j])*(x - lambda))
      A[j,i] = A[i,j]
    }
  return(A)
}

R_sl = function(A, p, eps = 1e-8){
  decomp = eigen(A, only.values = FALSE)
  D = diag(p)
  diag(D) = pmax(decomp$values, eps)
  U = decomp$vectors
  ret = tcrossprod(U %*% D, U)
  return((ret+t(ret))/2)
}

cond = function(X, lastX){
  return(norm(X-lastX, type = "F") / max(c(1, norm(X,type = 'F'), norm(lastX, type = 'F'))))
}



cdtrace = function(FsigmaF, lambda, rho){
  p = nrow(FsigmaF)
  F = diag(p) - matrix(1,p,p)/p
  I = diag(p)
  
  #G solusion
  eigF = eigen(F)
  eigFsigmaF = eigen(FsigmaF)
  eA = eigF$values
  UA = eigF$vectors
  eB = eigFsigmaF$values
  UB = eigFsigmaF$vectors
  D1 = 1/(matrix(eA, ncol=1) %*% matrix(eB, nrow=1)+6*rho)
  D2 = 1/(matrix(eB, ncol=1) %*% matrix(eA, nrow=1)+6*rho)
  
  Theta1 = Theta2 = Theta3 = Theta4 =solve(FsigmaF+5*I)
  Lambda1 = Lambda2 = Lambda3 = Lambda4 = Lambda5 = matrix(0,p,p)
  
  lastTheta1 = lastTheta2 = lastTheta3 = lastTheta4 = Theta1
  
  niter = 0
  while(niter <= 100000){
    Theta1 = G_sl(UA, UB, 2*rho*(Theta3 +Theta2+Theta4)+F+2*(Lambda1-Lambda3+Lambda4),D1)
    Theta2 = G_sl(UB, UA, 2*rho*(Theta3 +Theta1+Theta4)+F+2*(Lambda3-Lambda2-Lambda5),D2)
    Theta3 = S_sl((rho*(Theta1+Theta2)-Lambda1+Lambda2)/(2*rho), lambda/(2*rho), p)
    Theta4 = R_sl((rho*(Theta1+Theta2)-Lambda4+Lambda5)/(2*rho),p)
    
    Lambda1 = Lambda1+rho*(Theta3-Theta1)
    Lambda2 = Lambda2+rho*(Theta2-Theta3)
    Lambda3 = Lambda3+rho*(Theta1-Theta2)
    Lambda4 = Lambda4+rho*(Theta4-Theta1)
    Lambda5 = Lambda5+rho*(Theta2-Theta4)
    
    niter = niter+1
    
    if(all(c(cond(Theta1,lastTheta1),cond(Theta2,lastTheta2),
             cond(Theta3,lastTheta3),cond(Theta4,lastTheta4))<0.001)){
      break
    }
    lastTheta1 = Theta1
    lastTheta2 = Theta2
    lastTheta3 = Theta3
    lastTheta4 = Theta4
  }
  Theta3
}

CompDtrace = function(clr.data, lambda = NULL, rho = 0.5,
                      nlambda = 20, lambda.min.ratio = 0.004, 
                      lambda.max.ratio = 10){
  n <- nrow(clr.data)
  p <- ncol(clr.data)
  FsigmaF <- cov(clr.data)*(1-1/n)
  F = diag(p) - matrix(1,p,p)/p
  
  if(!is.null(lambda)){ nlambda <- length(lambda) }
  if(is.null(lambda)){
    lambda.max <- lambda.max.ratio*max(abs(FsigmaF))
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.max),log(lambda.min),length=nlambda))
  }
  
  ret = list()
  path = list()
  bic = c()
  for(i in 1:nlambda){
    ret[[i]] <- cdtrace(FsigmaF, lambda = lambda[i], rho)
    path[[i]] <- 1*(ret[[i]]!=0)
    bic = c(bic, n*norm((F%*%ret[[i]]%*%FsigmaF+FsigmaF%*%ret[[i]]%*%F)/2-F,"1")+
              log(n)*sum(path[[i]][upper.tri(path[[i]])]))
    print(i)
  }
  
  return(list(Theta=ret, path = path, lambda = lambda, nlambda = nlambda, bic = bic))
}




# set.seed(666)
# p=5
# Omega = 2*diag(p)
# Omega[1,2] = Omega[2,1] = 1
# Omega[2,3] = Omega[3,2] = 1
# Omega[3,4] = Omega[4,3] = 1
# Omega[4,5] = Omega[5,4] = 1
# Omega[1,5] = Omega[5,1] = 1
# 
# require(MASS)
# Sigma<-solve(Omega)
# n = 200
# p=5
# Y<-exp(mvrnorm(n,rep(0,p),Sigma))
# Y <- Y/rowSums(Y)
# F = diag(p) - matrix(1,p,p)/p
# x.f = Y
# FsigmaF = F%*%cov(log(Y))%*%F*(1-1/n)
# cov(clr(Y))/ cov((log(Y)%*%F))
# 
# log(2, base = exp(1))
