library(speff2trial)
data(ACTG175)

#compare ZDV+ddI vs ZDV+zalcitabine
y=ACTG175$cd420[ACTG175$arms%in%c(1,2)]
a=ACTG175$arms[ACTG175$arms%in%c(1,2)]
a[a==2]=-1
x1=as.matrix(ACTG175[ACTG175$arms%in%c(1,2),c("age","wtkg","karnof","cd40","cd80","hemo","homo","drugs","race","gender","str2","symptom")])

#Standardize X
x1[,1:5]=scale(x1[,1:5])
data1<-list(x=x1,y=y,a=a)

#pre-process
ytil=data1$y*data1$a
data1r=data1
data1r$y=ytil

mean(data1$y[data1$a==1])
mean(data1$y[data1$a==-1])

mod<-lm(data1r$y~data1r$x)
betai=as.numeric(mod$coefficients[-1])
betai=betai/betai[1]
print(betai)
sh=0.0001
n=length(data1r$y)
gam=1
s=1/sqrt(n)

#beta initial estimation without penalty
betahat1=c(1,optim(betai[-1],function(x) loss_psg(data=data1r,beta=c(1,x),s=s,lambda=0,betat=betai),method = "BFGS")$par)
betahat2=c(1,optim(betai[-1],function(x) loss_prsg(data=data1r,beta=c(1,x),s=s,lambda=0,betat=betai),method = "BFGS")$par)


#lambda selection

lambdav1=c(0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0)/10
lambdav2=c(0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0)
nlam=length(lambdav1)
cvv12=rep(0,nlam)
cvv22=rep(0,nlam)

for(j in 1:nlam){
  lambda1=lambdav1[j]
  lambda2=lambdav2[j]
  cv12=cv_mrc(data=data1,nfold=5,lambda=lambda1,s=s,betai=betahat1,stepsize=stepsize,sh=sh,criteria="loss_ind",method = "optim_mrc",betat=betahat1,gam=gam)
  cvv12[j]=cv12
  cv22=cv_mrc(data=data1,nfold=5,lambda=lambda2,s=s,betai=betahat2,stepsize=stepsize,sh=sh,criteria="loss_cal",method = "optim_mre",betat=betahat2,gam=gam)
  cvv22[j]=cv22
}
lam12=mean(lambdav1[cvv12==min(cvv12)])
lam22=mean(lambdav2[cvv22==min(cvv22)])

#beta estimation with penalty
betahat12=c(1,optim(betahat1[-1],function(x) loss_psg(data=data1r,beta=c(1,x),s=s,lambda=lam12,betat=betahat1,gam=gam),method = "BFGS")$par)
betahat12[abs(betahat12)<sh]=0
betahat22=c(1,optim(betahat2[-1],function(x) loss_prsg(data=data1r,beta=c(1,x),s=s,lambda=lam22,betat=betahat2,gam=gam),method = "BFGS")$par)
betahat22[abs(betahat22)<sh]=0

#c estimation
cvec=seq(-50,50,by=0.2)
evv12=sapply(cvec,function(x) evc(data=data1,betahat=betahat12,c=x))
evv22=sapply(cvec,function(x) evc(data=data1,betahat=betahat22,c=x))
maxindex12=which.max(evv12)
c12=mean(cvec[maxindex12])
maxindex22=which.max(evv22)
c22=mean(cvec[maxindex22])


#Predict ITRs for the data
ITRs12=sign(x1%*%betahat12-c12)
ITRs22=sign(x1%*%betahat22-c22)

#Wald test
testhat12=xbtest(xnew = xnew,betahat = betahat12,data = data1r,lambda = lam12,s=1,method="mrc",betat=betahat1,gam=1,c=c12)
pvalue12=testhat12[[1]][1,1]
esd12=testhat12[[2]][1,1]
testhat22=xbtest(xnew = xnew,betahat = betahat22,data = data1r,lambda = lam22,s=1,method="mre",betat=betahat2,gam=1,c=c22)
pvalue22=testhat22[[1]][1,1]
esd22=testhat22[[2]][1,1]
