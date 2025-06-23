#################
###packages######
#################
library(MASS)
library(glmnet)


##################
###functions######
##################
muf<-function(x){# baseline function of X
  y=sin(rowSums(x^2))
  return(y)
}

hinge<-function(x){# hinge loss function
  y=pmax(1-x,0)
  return(y)
}

ind<-function(x){# indicated function
  return(x>0)
}

sg<-function(x,s){#sigmoid function
  y=1/(1+exp(-x/s))
  y[is.na(y)]=0
  return(y)
}

snp<-function(x,s){#first derivative of the sigmoid function sn=1/(1+exp(-u/s))
  y=exp(-x/s)/(s*(1+exp(-x/s))^2)
  y[is.na(y)]=0
  return(y)
}

lsg<-function(x,s,x0){
  y=sg(x0,s)+snp(x0,s)*(x-x0)
}

snpp<-function(x,s){#second derivative of the sigmoid function sn=1/(1+exp(-u/s))
  y=(exp(-2*x/s)-exp(-x/s))/(s^2*(1+exp(-x/s))^3)
  y[is.na(y)]=0
  return(y)
}

varf<-function(x,s){
  y=snp(x,s)^2/snpp(x,s)^2
  return(y)
}

loss_ind<-function(data,beta){# negative rank correlation function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(ind(data$y-data$y[i]))%*%ind((data$x-xim)%*%beta)
  }
  l=-1*l/(n*(n-1))
  return(l)
}

loss_cal<-function(data,beta){#negative concordance-assistant learning objective function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(data$y-data$y[i])%*%ind((data$x-xim)%*%beta)
  }
  l=-1*l/(n*(n-1))
  return(l)
}

loss_linear<-function(data,beta){# negative rank correlation function with indicate function of X^T beta replaced by a linear function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(ind(data$y-data$y[i]))%*%((data$x-xim)%*%beta)
  }
  l=-1*l/(n*(n-1))
  return(l)
}

loss_sg<-function(data,beta,s){# negative rank correlation function with indicate function of X^T beta replaced by sigmoid function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(ind(data$y-data$y[i]))%*%sg((data$x-xim)%*%beta,s)
  }
  l=-1*l/(n*(n-1))
  return(l)
}

loss_psg<-function(data,beta,s,lambda,penalty="ALasso",betat=NA,gam=1){#loss function for penalized mrc with sigmoid function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(ind(data$y-data$y[i]))%*%sg((data$x-xim)%*%beta,s)
  }
  if(penalty=="Lasso"){
    l=-1*l/(n*(n-1))+lambda*sum(abs(beta))
  }
  if(penalty=="ALasso"){
    l=-1*l/(n*(n-1))+lambda*sum(abs(beta)/abs(betat)^gam)
  }
  return(l)
}

loss_prsg<-function(data,beta,s,lambda,penalty="ALasso",betat=NA,gam=1){#loss function for penalized maximum rank estimator with sigmoid function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(data$y-data$y[i])%*%sg((data$x-xim)%*%beta,s)
  }
  if(penalty=="Lasso"){
    l=-1*l/(n*(n-1))+lambda*sum(abs(beta))
  }
  if(penalty=="ALasso"){
    l=-1*l/(n*(n-1))+lambda*sum(abs(beta)/abs(betat)^gam)
  }
  return(l)
}

loss_scal<-function(data,beta,lambda){#loss function for sparse concordance-assistant learning
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(data$y-data$y[i])%*%hinge((data$x-xim)%*%beta)
  }
  l=l/(n*(n-1))+lambda*sum(abs(beta))
  return(l)
}

score_pcl<-function(data,beta,s){
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t((data$y[i]>data$y)*snp((xim-data$x)%*%beta,s))%*%(xim-data$x)[,-1]
  }
  l=-1*l/(n*(n-1))
  return(l)
}

score_sg<-function(data,beta,s){#score function for the MRC function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t((data$y[i]>data$y)*snp((xim-data$x)%*%beta,s))%*%(xim-data$x)[,-1]
  }
  l=-1*l/(n*(n-1))
  return(l)
}

score_cal<-function(data,beta,s){#score function for the CAL function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t((data$y[i]-data$y)*snp((xim-data$x)%*%beta,s))%*%(xim-data$x)[,-1]
  }
  l=-1*l/(n*(n-1))
  return(l)
}

scorep_sg<-function(data,beta,s){
  n=length(data$y)
  l=0
  for(i in 1:n){
    for(j in 1:n){
      l=l+ind(data$y[j]-data$y[i])*as.numeric(snpp((data$x[j,]-data$x[i,])%*%beta,s))*t(t((data$x[j,]-data$x[i,])[-1]))%*%((data$x[j,]-data$x[i,])[-1])
    }
  }
  l=-1*l/(n*(n-1))
  return(l)
}

loss_sg_penal<-function(data,beta,s,lambda){#penalized MRC function
  n=length(data$y)
  l=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    l=l+t(ind(data$y-data$y[i]))%*%sg((data$x-xim)%*%beta,s)
  }
  l=-1*l/(n*(n-1))+lambda*sum(abs(beta[-1]))
  return(l)
}

ev<-function(data,betahat){#empirical value function
  n=length(data$y)
  dhat=sign(data$x%*%betahat)
  ev=2*sum(data$y[dhat==data$a])/n
  return(ev)
}

evc<-function(data,betahat,c){#empirical value function with c
  n=length(data$y)
  dhat=sign(data$x%*%betahat-c)
  ev=2*sum(data$y[dhat==data$a])/n
  return(ev)
}

cv_mrc<-function(data,nfold=5,lambda,s,betai,stepsize,sh,method="optim",criteria="value",penalty="ALasso",betat=NA,gam=1){
  #cross validation with empirical value function
  n=length(data$y)
  cv=0
  ntest=n/nfold
  for(i in 1:nfold){
    data_train=data
    data_test=data
    data_train$x=data$x[-((1+(i-1)*ntest):(i*ntest)),]
    data_train$y=data$y[-((1+(i-1)*ntest):(i*ntest))]
    data_train$a=data$a[-((1+(i-1)*ntest):(i*ntest))]
    data_test$x=data$x[((1+(i-1)*ntest):(i*ntest)),]
    data_test$y=data$y[((1+(i-1)*ntest):(i*ntest))]
    data_test$a=data$a[((1+(i-1)*ntest):(i*ntest))]
    if(method=="gd"){
      ytil=data_train$y*data_train$a
      data_train_r=data_train
      data_train_r$y=ytil
      betahat=as.numeric(rmrc(data=data_train_r,lambda = lambda,s=s,betai=betai,stepsize = stepsize,sh=sh,steps = 1000))
    }
    if(method=="scal"){
      ytil=data_train$y*data_train$a
      data_train_r=data_train
      data_train_r$y=ytil
      betahat=as.numeric(cal_subgradient(data=datar,lambda=lambda,betai=betai,stepsize=stepsize,sh=0.01,steps=steps))
    }
    if(method=="optim_mrc"){
      ytil=data_train$y*data_train$a
      data_train_r=data_train
      data_train_r$y=ytil
      betahat=c(1,optim(betai[-1],function(x) loss_psg(data=data_train_r,beta=c(1,x),s=s,lambda=lambda,betat=betat,gam=gam),method = "BFGS")$par)
      betahat[abs(betahat)<sh]=0
    }
    if(method=="optim_mre"){
      ytil=data_train$y*data_train$a
      data_train_r=data_train
      data_train_r$y=ytil
      betahat=c(1,optim(betai[-1],function(x) loss_prsg(data=data_train_r,beta=c(1,x),s=s,lambda=lambda,betat=betat,gam=gam),method = "BFGS")$par)
      betahat[abs(betahat)<sh]=0
    }
    if(method=="optim_scal"){
      ytil=data_train$y*data_train$a
      data_train_r=data_train
      data_train_r$y=ytil
      betahat=c(1,optim(betai[-1],function(x) loss_scal(data=data_train_r,beta=c(1,x),lambda=lambda))$par)
      betahat[abs(betahat)<sh]=0
    }
    if(criteria=="value"){
      cv=cv+ev(data=data_test,betahat=betahat)
    }
    if(criteria=="loss_ind"){
      ytil=data_test$y*data_test$a
      data_test_r=data_test
      data_test_r$y=ytil
      cvk=loss_ind(data=data_test_r,beta=betahat)
      print(cvk)
      cv=cv+cvk
    }
    if(criteria=="loss_cal"){
      ytil=data_test$y*data_test$a
      data_test_r=data_test
      data_test_r$y=ytil
      cvk=loss_cal(data=data_test_r,beta=betahat)
      print(cvk)
      cv=cv+cvk
    }
  }
  cv=cv/nfold
  return(cv)
}

##################
###Optimization###
##################

lse<-function(data){#least square estimator
  n=length(data$y)
  xtil=cbind(rep(1,n),data$x)
  betahat=solve(t(xtil)%*%xtil)%*%t(xtil)%*%data$y
  return(betahat[-1,])
}



rmrc<-function(data,lambda,s,betai,stepsize,sh=0.01,steps=100){#maximum rank correlation estimation
  betai=matrix(betai,nrow=length(betai))
  n=length(data$y)
  for(i in 1:steps){
    stepi=stepsize/i
    gd=0
    for(j in 1:n){
      xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
      gd=gd+(-1)/(n*(n-1))*t((data$y[j]>data$y)*snp((xjm-data$x)%*%betai,s))%*%(xjm-data$x)
    }
    betai1=betai
    betai1[1,]=0
    gd=t(gd)+lambda*sign(betai1)
    gd[1]=0
    betahat=betai-gd*stepsize
    betahat[abs(betahat)<sh]=0
    print(as.numeric(betahat))
    if(norm(betahat-betai,type = "1")<0.0001){
      break
    }
    betai=betahat
  }
  return(betahat)
}

cal<-function(data,lambda,s,betai,stepsize,sh=0.01,steps=100){#maximum rank correlation estimation
  betai=matrix(betai,nrow=length(betai))
  n=length(data$y)
  for(i in 1:steps){
    stepi=stepsize/i
    gd=0
    for(j in 1:n){
      xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
      gd=gd+(-1)/(n*(n-1))*t((data$y[j]-data$y)*snp((xjm-data$x)%*%betai,s))%*%(xjm-data$x)
    }
    betai1=betai
    betai1[1,]=0
    gd=t(gd)+lambda*sign(betai1)
    gd[1]=0
    betahat=betai-gd*stepsize
    betahat[abs(betahat)<sh]=0
    print(as.numeric(betahat))
    if(norm(betahat-betai,type = "1")<0.0001){
      break
    }
    betai=betahat
  }
  return(betahat)
}

rmrc_constraint<-function(data,lambda,s,betai,xnew,stepsize,sh=0.01,steps=100){#maximum rank correlation estimation
  betai=matrix(betai,nrow=length(betai))
  n=length(data$y)
  for(i in 1:steps){
    stepi=stepsize/i
    gd=0
    for(j in 1:n){
      xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
      gd=gd+(-1)/(n*(n-1))*t((data$y[j]>data$y)*snp((xjm-data$x)%*%betai,s))%*%(xjm-data$x)
    }
    betai1=betai
    betai1[1,]=0
    gd=t(gd)+lambda*sign(betai1)
    gd[1]=0
    betahat=betai-gd*stepsize
    betahat[abs(betahat)<sh]=0
    xnew1=xnew
    xnew1[1]=0
    if(as.numeric(betahat)%*%xnew>0){
      betahat=betahat-as.numeric((xnew%*%betahat)/(xnew[-1]%*%xnew[-1]))*xnew1
    }
    print(as.numeric(betahat))
    #print(as.numeric(betahat)%*%xnew)
    if(norm(betahat-betai,type = "1")<0.0001){
      break
    }
    betai=betahat
  }
  return(betahat)
}

rmrc_subgradient<-function(data,lambda,s,betai,stepsize,sh=0.01,steps=100){#penalized maximum rank correlation estimation
  betai=matrix(betai,nrow=length(betai))
  n=length(data$y)
  for(i in 1:steps){
    stepi=stepsize/i
    gd=0
    for(j in 1:n){
      xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
      gd=gd+(-1)/(n*(n-1))*t((data$y[j]>data$y)*snp((xjm-data$x)%*%betai,s))%*%(xjm-data$x)
    }
    betai1=betai
    betai1[1,]=0
    gd=t(gd)+lambda*sign(betai1)
    betahat=betai-gd*stepi
    betahat=betahat/abs(betahat[1,])
    betahat[abs(betahat)<sh]=0
    #print(as.numeric(betahat))
    if(norm(betahat-betai,type = "1")<0.01){
      break
    }
    betai=betahat
  }
  return(betahat)
}

cal_subgradient<-function(data,lambda,betai,stepsize=1,sh=0.01,steps=100){#sparse concordance-assisted learning
  betai=matrix(betai,nrow=length(betai))
  n=length(data$y)
  for(i in 1:steps){
    gd=0
    for(j in 1:n){
      xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
      gd=gd+1/(n*(n-1))*t((data$y[j]-data$y)*dhinge((xjm-data$x)%*%betai))%*%(xjm-data$x)
    }
    betai1=betai
    betai1[1,]=0
    gd=t(gd)+lambda*sign(betai1)
    betahat=betai-gd*stepsize
    betahat=betahat/abs(betahat[1,])
    betahat[abs(betahat)<sh]=0
    print(as.numeric(betahat))
    print(i)
    if(norm(betahat-betai,type = "1")<0.01){
      break
    }
    betai=betahat
  }
  return(betahat)
}

rmrc_la<-function(data){#maximum rank correlation estimation
  n=length(data$y)
  la=0
  for(j in 1:n){
    xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
    la=la+(-1)/(n*(n-1))*t(data$y[j]>data$y)%*%(xjm-data$x)
  }
  betahat=-t(la)/norm(la,type = "2")
  return(betahat)
}

evf<-function(data,pihat,betahat,c){#empirical value function
  evf=sum(data$y*(data$a==sign(data$x%*%betahat-c))/pihat)
  return(evf)
}
#####################
###hypothesis test###
#####################
#empirical likelihood ratio test function
emplik = function(  z,  # matrix with one data vector per row, a column vector is ok when d=1
                    mu,  # hypothesized mean, default (0 ... 0) in R^d
                    lam,  # starting lambda, default (0 ... 0)
                    eps,  # lower cutoff for -log( ), default 1/nrow(z)
                    M,  # upper cutoff for -log( ), default Inf
                    thresh=1e-30,  # convergence threshold for log likelihood (default is aggressive)
                    itermax=100,  # upper bound on number of Newton steps (seems ample)
                    verbose=FALSE)  # controls printed output
{ 
  # empirical likelihood test for whether
  # mean of (rows of) z is mu
  
  
  # Internal function mllog, modified -log( ) with derivatives
  
  mllog = function( x, eps, M, der=0 ){
    # minus log and its first der derivatives, on  eps < x < M
    # 4th order Taylor approx to left of eps and right of M
    # der = 0 or 1 or 2
    # 4th order is lowest that gives self concordance
    
    if( missing(M) )
      M = Inf
    if( eps>M )
      stop("Thresholds out of order")
    
    lo = x < eps
    hi = x > M
    md = (!lo) & (!hi)
    
    # Coefficients for 4th order Taylor approx below eps
    coefs      = rep(0,5)
    coefs[1]   = -log(eps)
    coefs[2:5] = (-eps)^-(1:4)/(1:4)
    
    # Coefficients for 4th order Taylor approx above M
    Coefs      = rep(0,5)
    Coefs[1]   = -log(M)
    Coefs[2:5] = (-M)^-(1:4)/(1:4)
    
    # degree 4 polynomial approx to log
    h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
      # sum c[t+1] y^t
      tee = 1:4
      ans = y*0
      ans = ans + cvals[1]
      for( j in tee )
        ans = ans + y^j*cvals[j+1]
      ans
    }
    
    # first derivative of h at y, from approx at pt
    hp = function(y,pt){
      tee = 0:3
      ans = y*0
      for( j in tee )
        ans = ans + (-y/pt)^j
      ans = ans * (-pt)^-1
      ans
    }
    
    # second derivative of h at y, from approx at pt
    hpp = function(y,pt){
      tee = 0:2
      ans = y*0
      for( j in tee )
        ans = ans + (j+1) * (-y/pt)^j
      ans = ans *(-pt)^-2 
      ans
    }
    
    # function value
    f      = x*0
    f[lo]  = h( x[lo]-eps, coefs )
    f[hi]  = h( x[hi]-M,   Coefs )
    f[md]  = -log(x[md])
    
    if( der<1 )return(cbind(f))
    
    # first derivative
    fp     = x*0
    fp[lo] = hp( x[lo]-eps, eps )
    fp[hi] = hp( x[hi]-M, M )
    fp[md] = -1/x[md]
    
    if( der<2 )return(cbind(f,fp))
    
    # second derivative
    fpp     = x*0
    fpp[lo] = hpp( x[lo]-eps, eps )
    fpp[hi] = hpp( x[hi]-M, M )
    fpp[md] = 1/x[md]^2
    
    return( cbind(f,fp,fpp) )
    # End of mllog()
  } 
  
  # Internal function to do a linear model via SVD
  # Empirical likelihood's Newton steps are of
  # least squares type.
  
  svdlm = function(X,y){
    # Linear model regression coefficient via SVD
    
    # Tolerances for generalized inverse via SVD
    RELTOL = 1e-9     
    ABSTOL = 1e-100
    
    # Get Xplus = generalized inverse of X
    # If svd algorithm failures are encountered
    # it sometimes helps to try svd(t(X)) and
    # translate back. First check to ensure that
    # X does not contain NaN or Inf or -Inf.
    svdX     = svd(X)
    d        = svdX$d
    lo       = d < (RELTOL * max(d) + ABSTOL)
    dinv     = 1/d
    dinv[lo] = 0
    Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
    # taking care with diag when dinv is 1x1
    # to avoid getting the identity matrix of
    # size floor(dinv)
    
    # Return X^+ y
    Xplus %*% matrix(y,ncol=1)
  }
  # end of svdlm
  
  
  # Backtracking line search parameters [Tweak only with extreme caution.] 
  # See Boyd and Vandenberghe, pp464-466.
  ALPHA = 0.3  # seems better than 0.01 on some 2d test data (sometimes fewer iters)
  BETA  = 0.8
  # We need  0 < ALPHA < 1/2   and 0 < BETA < 1
  
  # Backtrack threshold: you can miss by this much.
  BACKEPS = 0
  # Consider replacing 0 by 1e-10 if backtracking seems to be
  # failing due to round off.
  
  if( is.vector(z) )
    z = matrix(z,ncol=1)
  
  n = nrow(z)
  d = ncol(z)
  
  if( missing(mu)  )
    mu  = rep(0,d)
  z = t( t(z)-mu ) # subtract mu from each z[i,]
  
  if( missing(eps) )eps = 1/n
  if( missing(M)   )M   = Inf
  
  #
  # Use lam = 0 or initial lam, whichever is best
  #
  
  init0 = mllog( rep(1,n), eps=eps, M=M, der=2 ) # i.e. lam = 0
  
  if( missing(lam) ){
    init = init0
    lam  = rep(0,d)
  }else{
    init = mllog( 1+z%*%lam, eps=eps, M=M, der=2 )
    if( sum(init0[,1]) < sum(init[,1]) ){
      lam = rep(0,d)
      init = init0
    }
  }
  
  # Initial f, g
  fold = sum(init[,1])
  gold = apply( z * init[,2],2,sum )
  
  converged = FALSE
  iter      = 0
  oldvals   = init
  while( !converged ){
    iter = iter + 1
    
    # Get Newton Step
    rootllpp = sqrt(oldvals[,3])  # sqrt 2nd deriv of -llog lik
    zt = z
    for( j in 1:d )
      zt[,j] = zt[,j] * rootllpp
    yt   = oldvals[,2] / rootllpp
    step = -svdlm(zt,yt)  #  more reliable than step = -lm( yt~zt-1 )$coef
    
    backtrack = FALSE 
    s = 1   # usually called t, but R uses t for transpose
    while( !backtrack ){
      newvals = mllog( 1+z%*%(lam+s*step),eps=eps,M=M,der=2 )
      fnew    = sum(newvals[,1])
      targ    = fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
      if(  fnew <= targ ){
        # backtracking has converged
        backtrack = TRUE
        oldvals   = newvals
        fold      = fnew
        gold      = apply( z * oldvals[,2],2,sum )
        # take the step
        lam       = lam + s*step
      }else{
        s = s * BETA
      }
    }
    
    # Newton decrement and gradient norm
    ndec     = sqrt( sum( (step*gold)^2 ) )
    gradnorm = sqrt( sum(gold^2))
    
    if(verbose)print(c(fold,gradnorm,ndec,lam))
    
    converged = ( ndec^2 <= thresh)
    if( iter > itermax )break
  }
  
  wts    = (1/n)/(1+z%*%lam)
  logelr = sum( mllog(1+z%*%lam,eps=eps,M=M,der=0) )
  
  list(logelr=logelr,lam=lam, wts=wts,
       converged=converged,iter=iter,ndec=ndec,gradnorm=gradnorm)
}

#empirical likelihood ratio test
elt<-function(data,s,betahat){
  n=length(data$y)
  eta=matrix(nrow=(n*(n-1)/2),ncol=length(betahat))
  for(i in 2:n){
    for(j in 1:(i-1))
    eta[(sum(1:(i-1))-i+1+j),]=t((data$y[i]>data$y[j])*snp((data$x[i]-data$x[j])%*%betahat,s))%*%(data$x[i]-data$x[j])+t((data$y[i]<data$y[j])*snp((data$x[j]-data$x[i])%*%betahat,s))%*%(data$x[j]-data$x[i])
  }
  llr=emplik(eta)$logelr
  ts=-2*llr
  pvalue=1-pchisq(ts,nx)
  return(pvalue)
}


#jackknife empirical likelihood
jeltb<-function(data,s,betahat){
  n=length(data$y)
  u=score_sg(data=data,beta=betahat,s=s)
  v=matrix(nrow = n,ncol = (length(betahat)-1))
  for(i in 1:n){
    datai=data
    datai$y=datai$y[-i]
    datai$x=datai$x[-i,]
    v[i,]=n*u-(n-1)*score_sg(data=datai,beta = betahat,s=s)
  }
  llr=emplik(v)$logelr
  ts=-2*llr
  pvalue=1-pchisq(ts,length(betahat))
  return(pvalue)
}

vhat<-function(data,s,i,beta){
  n=length(data$y)
  u=score_sg(data=data,beta=betahat,s=s)
  datai=data
  datai$y=datai$y[-i]
  datai$x=datai$x[-i,]
  vi=n*u-(n-1)*score_sg(data=datai,beta = beta,s=s)
  return(vi)
}

logstar<-function(x,n){
  if(x>1/n){
    y=log(x)
  }
  if(x<=(1/n)){
    y=2*n*x-(n*x)^2/2+log(1/n)-1.5
  }
  return(y)
}

jel<-function(data,s,betas,method="mrc"){
  n=length(data$y)
  beta=c(1,betas)
  if(method=="mrc"){
    u=score_sg(data=data,beta=beta,s=s)
  }
  if(method=="mre"){
    u=score_cal(data=data,beta=beta,s=s)
  }
  
  v=matrix(nrow = n,ncol = (length(beta)-1))
  for(i in 1:n){
    datai=data
    datai$y=datai$y[-i]
    datai$x=datai$x[-i,]
    if(method=="mrc"){
      v[i,]=n*u-(n-1)*score_sg(data=datai,beta = beta,s=s)
    }
    if(method=="mre"){
      v[i,]=n*u-(n-1)*score_cal(data=datai,beta = beta,s=s)
    }
  }
  jel=emplik(v)$logelr
  return(jel)
}

jel_constraint<-function(data,s,xnew,betas){
  n=length(data$y)
  beta2=c(1,0,betas)%*%xnew/xnew[2]*(-1)
  beta=c(1,beta2,betas)
  u=score_sg(data=data,beta=beta,s=s)
  v=matrix(nrow = n,ncol = (length(beta)-1))
  for(i in 1:n){
    datai=data
    datai$y=datai$y[-i]
    datai$x=datai$x[-i,]
    v[i,]=n*u-(n-1)*score_sg(data=datai,beta = beta,s=s)
  }
  jel=emplik(v)$logelr
  return(jel)
}
#quantile
quantile<-function(data,betahat,s,xnew){
  n=length(data$y)
  u=score_sg(data=data,beta=beta,s=s)
  up=scorep_sg(data=data,beta=beta,s=s)
  v=matrix(nrow = n,ncol = (length(beta)-1))
  vp=list(NULL)
  for(i in 1:n){
    datai=data
    datai$y=datai$y[-i]
    datai$x=datai$x[-i,]
    v[i,]=n*u-(n-1)*score_sg(data=datai,beta = beta,s=s)
    vp[[i]]=n*up-(n-1)*scorep_sg(data=datai,beta = beta,s=s)
  }
  A=0
  B=0
  for(i in 1:n){
    xim=matrix(rep(data$x[i,],n),nrow=n,byrow = T)
    A=A+vp[[i]]
    B=B+t(t(v[i,]))%*%v[i,]
  }
  A=A/n
  B=B/n
  Sigma0=xnew[-1]%*%t(A)%*%solve(B)%*%A%*%xnew[-1]
}
#single beta test
bitest<-function(data,betahat,index,lambda,s){
  n=length(data$y)
  nx=length(betahat)
  d1=0
  for(i in 1:n){
    for(j in 1:n){
      d1=d1+(-1)/(n*(n-1))*(data$y[i]>data$y[j])*(snp((data$x[i,]-data$x[j,])%*%betahat,s)*(data$x[i,index]-data$x[j,index]))*(snp((data$x[i,]-data$x[j,])%*%betahat,s)*(data$x[i,index]-data$x[j,index]))
    }
  }
  var=-d1
  
}

xbtest<-function(xnew,betahat,data,lambda,s,method="mrc",penalty="ALasso",betat=NA,gam=0.1,c=0,less=T){#
  n=length(data$y)
  nx=length(betahat)
  bindex=which(betahat!=0)
  bindex=bindex[-1]
  d1=0
  d2=0
  for(i in 1:n){
    d1i=0
    if(method=="mrc"){
      for(j in 1:n){
      #d1=d1+(data$y[i]>data$y[j])*t(snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex]))%*%(snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex]))
      d1i=d1i+(data$y[i]>data$y[j])*snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex])+(data$y[j]>data$y[i])*snp((data$x[j,]-data$x[i,])%*%betahat,s)%*%((data$x[j,]-data$x[i,])[bindex])
      d2=d2+((data$y[i]>data$y[j])*snpp((data$x[i,]-data$x[j,])%*%betahat,s))[1,1]*(((data$x[i,]-data$x[j,])[bindex])%*%t((data$x[i,]-data$x[j,])[bindex]))+((data$y[j]>data$y[i])*snpp((data$x[j,]-data$x[i,])%*%betahat,s))[1,1]*(((data$x[j,]-data$x[i,])[bindex])%*%t((data$x[j,]-data$x[i,])[bindex]))
    }
    }
    if(method=="mre"){
      for(j in 1:n){
        #d1=d1+(data$y[i]>data$y[j])*t(snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex]))%*%(snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex]))
        d1i=d1i+(data$y[i]-data$y[j])*snp((data$x[i,]-data$x[j,])%*%betahat,s)%*%((data$x[i,]-data$x[j,])[bindex])+(data$y[j]-data$y[i])*snp((data$x[j,]-data$x[i,])%*%betahat,s)%*%((data$x[j,]-data$x[i,])[bindex])
        d2=d2+((data$y[i]-data$y[j])*snpp((data$x[i,]-data$x[j,])%*%betahat,s))[1,1]*(((data$x[i,]-data$x[j,])[bindex])%*%t((data$x[i,]-data$x[j,])[bindex]))+((data$y[j]-data$y[i])*snpp((data$x[j,]-data$x[i,])%*%betahat,s))[1,1]*(((data$x[j,]-data$x[i,])[bindex])%*%t((data$x[j,]-data$x[i,])[bindex]))
      }
    }
    d1=d1+t(d1i)%*%d1i/(n-1)^2
  }
  B=d1/n
  A=d2/((-2)*n*(n-1))
  D=lambda*sign(betahat[bindex])/abs(betat[bindex])^gam
  Ainv=solve(A)
  ndebias=xnew%*%betahat+xnew[bindex]%*%Ainv%*%D-c
  nsd=sqrt(xnew[bindex]%*%Ainv%*%B%*%Ainv%*%xnew[bindex]/n)#sd of xnew*beta
  if(less==T){
    pvalue=1-pnorm(ndebias,sd=nsd)#test H0: xnew*beta<0
  }
  if(less==F){
    pvalue=pnorm(ndebias,sd=nsd)#test H0: xnew*beta>0
  }
  dt=Ainv%*%D
  result=list(pvalue,nsd,ndebias,dt)
  return(result)
}

xbtest_bootstrap_reweight<-function(xnew,betahat,data,lambda,s,stepsize,sh=0.01,steps=100){#
  n=length(data$y)
  nx=length(betahat)
  bindex=which(betahat!=0)
  betai=betahat[bindex]
  ts=rep(0,100)
  for(k in 1:100){
    w=10*rbeta(n,0.125,1.125)
    for(i in 1:steps){
      gd=0
      for(j in 1:n){
        xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
        gd=gd+(-1)/(n*(n-1))*((rep(w[j],n)+w)*t((data$y[j]>data$y)*snp((xjm[,bindex]-data$x[,bindex])%*%betai,s)))%*%(xjm[,bindex]-data$x[,bindex])
      }
      betahati=betai-t(gd)*stepsize
      betahati=betahati/norm(betahati,type = "2")
      betai=betahati
    }
    ts[k]=xnew[bindex]%*%betahati
  }
  nsd=sd(ts)
  d2=0
  for(i in 1:n){
    for(j in 1:n){
      d2=d2+((data$y[i]>data$y[j])*snpp((data$x[i,]-data$x[j,])%*%betahat,s))[1,1]*(((data$x[i,]-data$x[j,])[bindex])%*%t((data$x[i,]-data$x[j,])[bindex]))+((data$y[j]>data$y[i])*snpp((data$x[j,]-data$x[i,])%*%betahat,s))[1,1]*(((data$x[j,]-data$x[i,])[bindex])%*%t((data$x[j,]-data$x[i,])[bindex]))
    }
  }
  A=d2/((-2)*n*(n-1))
  D=lambda*sign(betahat[bindex])
  Ainv=solve(A)
  ndebias=xnew%*%betahat-(1/sqrt(n))*xnew[bindex]%*%Ainv%*%D
  #nsd=sqrt(xnew[bindex]%*%Ainv%*%B%*%Ainv%*%xnew[bindex]/n)
  pvalue=1-pnorm(ndebias,sd=nsd)#test xnew*beta>0
  result=list(pvalue,nsd)
  return(result)
}

xbtest_bootstrap<-function(xnew,betahat,chat,pihat,betai,data,lambda,s,stepsize,sh=0.01,steps=100){#
  n=length(data$y)
  nx=length(betahat)
  bindex=which(betahat!=0)
  ts=rep(0,100)
  for(k in 1:100){
    dataindex=sample(1:n,n,replace=TRUE)
    xb=data$x[dataindex,]
    yb=data$y[dataindex,]
    datab=list(xb=xb,yb=yb)
    for(i in 1:steps){
      gd=0
      for(j in 1:n){
        xjm=matrix(rep(datab$xb[j,],n),nrow=n,byrow = T)
        gd=gd+(-1)/(n*(n-1))*t((datab$yb[j]>datab$yb)*snp((xjm-datab$xb)%*%betai,s))%*%(xjm-datab$xb)
      }
      betahati=betai-t(gd)*stepsize
      betahati=betahati/abs(betahati[1,])
      betai=betahati
    }
    cv=seq(-2,2,0.01)
    ev=sapply(cv, function(x) evf(data=data,pihat=pihat,betahat=betahat,c=x))
    chat=cv[ev==max(ev)]
    if(length(chat)>1){
      chat=mean(chat)
    }
    ts[k]=xnew%*%betahati-chat
  }
  nsd=sd(ts)
  d2=0
  for(i in 1:n){
    for(j in 1:n){
      d2=d2+((data$y[i]>data$y[j])*snpp((data$x[i,]-data$x[j,])%*%betahat,s))[1,1]*(((data$x[i,]-data$x[j,])[bindex])%*%t((data$x[i,]-data$x[j,])[bindex]))+((data$y[j]>data$y[i])*snpp((data$x[j,]-data$x[i,])%*%betahat,s))[1,1]*(((data$x[j,]-data$x[i,])[bindex])%*%t((data$x[j,]-data$x[i,])[bindex]))
    }
  }
  A=d2/((-2)*n*(n-1))
  D=lambda*sign(betahat[bindex])
  Ainv=solve(A)
  ndebias=xnew%*%betahat-(1/sqrt(n))*xnew[bindex]%*%Ainv%*%D
  #nsd=sqrt(xnew[bindex]%*%Ainv%*%B%*%Ainv%*%xnew[bindex]/n)
  pvalue=1-pnorm(ndebias-chat,sd=nsd)#test xnew*beta>0
  result=list(pvalue,nsd)
  return(result)
}
################################
###tuning parameter selection###
################################
vsnp<-function(data,betahat,s){
  n=length(data$y)
  xbd=matrix(nrow = n,ncol = n)
  for(i in 1:n){
    xbd[i,]=data$x%*%betahat-(data$x[i,]%*%betahat)[1,1]
  }
  ehat=rep(0,1000)
  for(i in 1:1000){
    samp=rnorm(n=(n^2/2),mean=0,sd=sd(xbd))
    ehat[i]=mean(snp(samp,s=s))
  }
  result=list(mv=mean(ehat),sdv=sd(ehat))
  return(result)
}

vsnpp<-function(data,betahat,s){
  n=length(data$y)
  xbd=matrix(nrow = n,ncol = n)
  for(i in 1:n){
    xbd[i,]=data$x%*%betahat-(data$x[i,]%*%betahat)[1,1]
  }
  ehat=rep(0,1000)
  for(i in 1:1000){
    samp=rnorm(n=(n^2/2),mean=0,sd=sd(xbd))
    ehat[i]=mean(snpp(samp,s=s))
  }
  result=list(mv=mean(ehat),sdv=sd(ehat))
  return(result)
}

vsnpp<-function(data,betahat,s){
  
  ehat=rep(0,100)
  for(i in 1:100){
    data=generate_data_mrc(n=100,beta=beta0,sigma_x=1,sigma=0.5)
    n=length(data$y)
    xbd=matrix(nrow = n,ncol = n)
    for(j in 1:n){
      xbd[j,]=(data$x%*%betahat-(data$x[j,]%*%betahat)[1,1])*(data$y>data$y[j])
    }
    ehat[i]=mean(snpp(c(xbd),s=s))
  }
  result=list(mv=mean(ehat),sdv=sd(ehat))
  return(result)
}

xbetas<-function(data,s,betahat){
  n=length(data$y)
  xbetam=matrix(nrow = n,ncol = n)
  for(j in 1:n){
    xjm=matrix(rep(data$x[j,],n),nrow=n,byrow = T)
    xbetam[j,]=abs((xjm-data$x)%*%betahat)/s
  }
  rate=sum(as.numeric(xbetam)>5)/(n^2-n)
  return(rate)
}
#####################
###data generation###
#####################
generate_data_mrc<-function(n,beta,Sigma_x,sigma,error_type="normal",method="pnorm"){#data generation for MRC
  nx=length(beta)
  x=mvrnorm(n=n,mu=rep(0,nx),Sigma=Sigma_x)
  if(error_type=="normal"){
    epsi=rnorm(n,0,sigma)
  }
  if(error_type=="mix"){
    epsi=rnorm(n,0,sigma)+rbinom(n,size=1,prob=0.05)*rnorm(n,-10,sigma)+rbinom(n,size=1,prob=0.05)*rnorm(n,10,sigma)
  }
  if(error_type=="t2"){
    epsi=rt(n,df=2)
  }
  if(method=="pnorm"){
    y=psif(x%*%beta)+epsi
  }
  if(method=="linear"){
    y=(x%*%beta)+epsi
  }
  data=list(x=x,y=y)
  return(data)
}

generate_data<-function(n,alpha,beta,Sigma_x,sigma){
  #data generation function
  nx=length(beta)
  x=mvrnorm(n=n,mu=rep(0,nx),Sigma=Sigma_x)
  pia=exp(x%*%alpha[-1]+alpha[1])/(1+exp(x%*%alpha[-1]+alpha[1]))
  a=rbinom(n=n,size=1,prob=pia)
  a[a==0]=-1
  epsi=rnorm(n,0,sigma)
  y=muf(x)+psif(x%*%beta)*a+epsi
  data=list(x=x,a=a,y=y)
  return(data)
}

generate_data_fixpi<-function(n,pi,beta,Sigma_x,sigma,psi="linear",noise="Gaussian"){
  #data generation function with fixed propensity score
  nx=length(beta)
  x=mvrnorm(n=n,mu=rep(0,nx),Sigma=Sigma_x)
  av=c(rep(1,pi*n),rep(-1,(1-pi)*n))
  a=sample(av)
  if(noise=="Gaussian"){
    epsi=rnorm(n,0,sigma)
  }
  if(noise=="GaussianC1"){
    epsi=rnorm(n,0,sigma)
    epsi[sample(1:n,(n/200),replace = F)]=runif(n=(n/200),min=-50,max=50)
  }
  if(noise=="GaussianC2"){
    epsi=rnorm(n,0,sigma)
    epsi[sample(1:n,(n/200),replace = F)]=runif(n=(n/200),min=-20,max=20)
  }
  if(noise=="GaussianC3"){
    epsi=rnorm(n,0,sigma)
    epsi[sample(1:n,(n/50),replace = F)]=runif(n=(n/50),min=-20,max=20)
  }
  if(noise=="t"){
    epsi=rt(n,df=2)
  }
  if(noise=="t2"){
    epsi=rt(n,df=2)/2
  }
  if(psi=="linear"){
    y=muf(x)+(x%*%beta)*a+epsi
  }
  if(psi=="linear0"){
    y=(x%*%beta)*a+epsi
  }
  if(psi=="2linear"){
    y=muf(x)+(2*x%*%beta)*a+epsi
  }
  if(psi=="4linear"){
    y=muf(x)+(4*x%*%beta)*a+epsi
  }
  if(psi=="5linear"){
    y=muf(x)+(5*x%*%beta)*a+epsi
  }
  if(psi=="10linear"){
    y=muf(x)+(10*x%*%beta)*a+epsi
  }
  if(psi=="cubic"){
    y=muf(x)+(x%*%beta)^3*a+epsi
  }
  if(psi=="cubic8"){
    y=muf(x)+(x%*%beta)^3/8*a+epsi
  }
  if(psi=="cubic2"){
    y=muf(x)+(x%*%beta)^3/2*a+epsi
  }
  if(psi=="exp"){
    y=muf(x)+(exp(x%*%beta)-1)*a+epsi
  }
  if(psi=="exp5"){
    y=muf(x)+(exp(x%*%beta)-1)/5*a+epsi
  }
  if(psi=="exp2"){
    y=muf(x)+(exp(x%*%beta)-1)/2*a+epsi
  }
  if(psi=="01"){
    y=muf(x)+(x%*%beta>1)*a+epsi
  }
  if(psi=="normalcdf"){
    y=muf(x)+(pnorm(x%*%beta)-0.5)*a+epsi
  }
  data=list(x=x,a=a,y=y)
  return(data)
}

generate_data_rct<-function(n,pi,beta,Sigma_x,sigma,psi="linear",noise="Gaussian"){
  nx=length(beta)
  x=mvrnorm(n=n,mu=rep(0,nx),Sigma=Sigma_x)
  a=rbinom(n=n,size=1,prob=pi)
  a[a==0]=-1
  if(noise=="Gaussian"){
    epsi=rnorm(n,0,sigma)
  }
  if(noise=="t"){
    epsi=rt(n,df=4)
  }
  if(psi=="linear"){
    y=muf(x)+5*(x%*%beta)*a+epsi
  }
  if(psi=="cubic"){
    y=muf(x)+(x%*%beta)^3*a+epsi
  }
  if(psi=="exp"){
    y=muf(x)+exp(x%*%beta)*a+epsi
  }
  if(psi=="01"){
    y=muf(x)+(x%*%beta>1)*a+epsi
  }
  if(psi=="normalcdf"){
    y=muf(x)+(pnorm(x%*%beta)-0.5)*a+epsi
  }
  data=list(x=x,a=a,y=y)
  return(data)
}
##################
####Algorithms####
##################
S<-function(a,b){#soft-thresholding operator
  S=sign(a)*max(abs(a)-b,0)
  return(S)
}

pgd<-function(data,beta0,s,lambda,betat,gam=1,stepsize=1,method="mrc"){#proximal gradient descent 
  betai=beta0
  nx=length(beta0)
  for(i in 1:10000){
    betai0=betai
    if(method=="mrc"){
      betai[-1]=betai[-1]-stepsize*score_sg(data,beta=betai,s)
    }
    if(method=="cal"){
      betai[-1]=betai[-1]-stepsize*score_cal(data,beta=betai,s)
    }
    for(j in 2:nx){
      betai[j]=S(betai[j],stepsize*lambda/abs(betat[j])^gam)
    }
    print(betai)
    if(norm(betai-betai0,type = "2")<0.001){
      break
    }
  }
  return(betai)
}