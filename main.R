source("/geode2/home/u020/cyishan/Quartz/Project2/scripts/Functions.R")
###### setting
n=1000#sample size
nx=5#number of covariates
xnew = c(1,2,2,1,-1)#new patients covariates for test
beta0=c(1,-1,0.5,0,0)#true beta
pi0=0.5#propensity score
betai=beta0+c(0,0.5,0.5,0.5,0.5)#initial guess for beta
sigma=1#error term variance
Sigma_x=diag(rep(1,length(beta0)))
gam=1#gamma for adaptive Lasso
sv=c(1/2,1,2,4)/sqrt(n)#tuning parameter for sigmoid function for tets
s1=1/sqrt(n)#tuning parameter for sigmoid function for beta estimation
ns=length(sv)
sh=0.01# hard threshold for beta estimation if proximal gradient descent is not used

###### data generation
data=generate_data_fixpi(n=n,pi=pi0,beta=beta0,Sigma_x=Sigma_x,sigma=sigma,psi = "2linear")

###### initial process
ytil=data$y*data$a
datar=data
datar$y=ytil

###### tuning parameter lambda selection
betahatm1=matrix(rep(0,ns*nx),ncol = nx)
betahatm2=matrix(rep(0,ns*nx),ncol = nx)
p1=rep(0,times=ns)
p2=rep(0,times=ns)

lambdav1=c(0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001,0)/10
lambdav2=c(0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001,0)
nlam=length(lambdav1)
cvv1=rep(0,nlam)
cvv2=rep(0,nlam)

###### betat estimation for adaptive Lasso
betat1=c(1,optim(betai[-1],function(x) loss_psg(data=datar,beta=c(1,x),s=s1,lambda=0,betat=betai),method = "BFGS")$par)
betat2=c(1,optim(betai[-1],function(x) loss_prsg(data=datar,beta=c(1,x),s=s1,lambda=0,betat=betai),method = "BFGS")$par)

###### tuning parameter selection with 5-folder CV
for(j in 1:nlam){
  lambda1=lambdav1[j]
  lambda2=lambdav2[j]
  cv1=cv_mrc(data=data,nfold=5,lambda=lambda1,s=s1,betai=betai,stepsize=stepsize,sh=sh,criteria="loss_ind",method = "optim_mrc",betat=betat1,gam=gam)
  cvv1[j]=cv1
  cv2=cv_mrc(data=data,nfold=5,lambda=lambda2,s=s1,betai=betai,stepsize=stepsize,sh=sh,criteria="loss_cal",method = "optim_mre",betat=betat2,gam=gam)
  cvv2[j]=cv2
}
lam1=mean(lambdav1[cvv1==min(cvv1)])
lam2=mean(lambdav2[cvv2==min(cvv2)])
###### beta estimation
betahat1=c(1,optim(betai[-1],function(x) loss_psg(data=datar,beta=c(1,x),s=s1,lambda=lam1,betat=betat1,gam=gam),method = "BFGS")$par)
betahat1[abs(betahat1)<sh]=0
betahat2=c(1,optim(betai[-1],function(x) loss_prsg(data=datar,beta=c(1,x),s=s1,lambda=lam2,betat=betat2,gam=gam),method = "BFGS")$par)
betahat2[abs(betahat2)<sh]=0

###### beta estimation with proximal gradient descent
#betahat1=pgd(data=datar,beta0=betai,s=s1,lambda=lam1,betat=betat1,gam=gam,stepsize=1,method = "mrc")
#betahat2=pgd(data=datar,beta0=betai,s=s1,lambda=lam2,betat=betat2,gam=gam,stepsize=1,method = "cal")
for(i in 1:ns){
  ###### s value
  s2=sv[i]
  ###### de-biased Wald test
  testhat1=xbtest(xnew = xnew,betahat = betahat1,data = datar,lambda = lam1,s=s2,method="mrc",betat=betat1,gam=gam)
  pvalue1=testhat1[[1]][1,1]
  testhat2=xbtest(xnew = xnew,betahat = betahat2,data = datar,lambda = lam2,s=s2,method="mre",betat=betat2,gam=gam)
  pvalue2=testhat2[[1]][1,1]
  ###### pvalues
  p1[i]=pvalue1
  p2[i]=pvalue2
}
