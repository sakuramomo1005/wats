## Prof. Ying Liu (Columbia Univsersity)'s codes for OWL
library(glmnet)
library(kernlab)
library(ggplot2)
library(MASS)


wsvm<-function(X,A,wR,kernel='linear',sigma=0.05,C=1,e=0.00001, ridge.reg =FALSE){
  wAR=A*wR
  if (kernel=='linear'){
    K=X%*%t(X)
  }
  else if (kernel=='rbf'){
    rbf=rbfdot(sigma=sigma)
    K=kernelMatrix(rbf,X)
  } else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))

  H=K*(wAR%*%t(wAR))
  if(ridge.reg) H <-  H  + diag(10^{-4}, nrow(H))
  n=length(A)
  solution=ipop(-abs(wR),H,t(A*wR),0,numeric(n),rep(C,n),0,maxiter=100)
  alpha=primal(solution)
  alpha1=alpha*wR*A
  if (kernel=='linear'){
    w=t(X)%*%alpha1 #parameter for linear
    fitted=X%*%w
    rm=sign(wR)*A-fitted
  } else if (kernel=='rbf'){
    #there is no coefficient estimates for gaussian kernel
    #but there is fitted value, first we compute the fitted value without adjusting for bias
    fitted=K%*%alpha1
    rm=sign(wR)*A-fitted
  }
  Imid =(alpha < C-e) & (alpha > e)
  rmid=rm[Imid==1];
  if (sum(Imid)>0){
    bias=mean(rmid)
  } else {
    Iup=((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
    Ilow=((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
    rup=rm[Iup]
    rlow=rm[Ilow]
    bias=(min(rup)+max(rlow))/2}
  fit=bias+fitted
  if (kernel=='linear') {
    model=list(alpha1=alpha1,bias=bias,fit=fit,beta=w)
    class(model)<-'linearcl'
  } else if (kernel=='rbf') {
    model=list(alpha1=alpha1,bias=bias,fit=fit,sigma=sigma,X=X)
    class(model)<-'rbfcl'}
  return (model)
}



Qlearning_Single<-function(H,A,R,pentype='lasso',m=4){
  n=length(A)
  X=cbind(H,A,diag(A)%*%H)
  if (pentype=='lasso'){
    cvfit=cv.glmnet(X,R,nfolds=m)
    co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
  }else if (pentype=='LSE'){
    co=coef(lm(R~X))
  }else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',
it can also be 'LSE' without penalization"))

  XX1=cbind(rep(1,n),H,rep(1,n),diag(n)%*%H)
  XX2=cbind(rep(1,n),H,rep(-1,n),-diag(n)%*%H)
  Q1=XX1%*%co
  Q2=XX2%*%co
  Q=apply(cbind(Q1,Q2),1,max)
  Qsingle=list(co=co,Q=Q)
  class(Qsingle)='qlearn'
  Qsingle
}

Qlearning<-function (X,AA,RR,K,pentype='lasso',m=4) {
  R_future=0
  coef=list()
  models=list()
  if (is.matrix(X)){
    for (j in K:1){
      R=RR[[j]]+R_future
      if (min(R)!=max(R)){
        models[[j]]=Qlearning_Single(X,AA[[j]],R,pentype=pentype)
        R_future=models[[j]]$Q}
      else {
        models[[j]]=list(co=rep(0,2+2*dim(X)[2]),Q=R)
        R_future=R
      }}
  }

  if (is.list(X)){
    for (j in K:1){
      R=RR[[j]]+R_future
      if (min(R)!=max(R)){
        models[[j]]=Qlearning_Single(X[[j]],AA[[j]],R,pentype=pentype,m=4)
        R_future=models[[j]]$Q}
      else {
        models[[j]]=list(co=rep(0,2+2*dim(X[[j]])[2]),Q=R)
        R_future=R }}
  }
  models}

predict.linearcl<-function(object,x,...){
  predict=sign(object$bias+x%*%object$beta)
}

predict.rbfcl<-function(object,x,...){
  rbf=rbfdot(sigma=object$sigma)
  n=dim(object$X)[1]
  if (is.matrix(x)) xm=dim(x)[1]
  else if (is.vector(x)) xm=1
  else stop('x must be vector or matrix')
  if (xm==1){ K <- apply(object$X,1,rbf,y=x)
  }else{   K<- matrix(0, xm, n)
  for (i in 1:xm) {K[i,]=apply(object$X,1,rbf,y=x[i,]) }}
  predict=sign(object$bias+K%*%object$alpha1)
}

# predict.rbfcl<-function(object,x){
#   rbf=rbfdot(sigma=object$sigma)
#   n=dim(object$X)[1]
#   if (is.matrix(x)) xm=dim(x)[1]
#   if (is.vector(x)) xm=1
#   K<- matrix(0, xm, n)
#   for (i in 1:n) {
#   K[,i]=apply(x,1,rbf,y=object$X[i,]) }
#   predict=sign(object$bias+K%*%object$alpha1)
# }

predict.qlearn<-function(object,x,...){
  p=dim(x)[2]
  predict=sign(object$co[p+2]+x%*%object$co[(p+3):(2*p+2)])
}




#this draw the coefficient plots for single stage coefficients
plot.linearcl<-function(x,index=NULL,names=NULL,ylab='std coefficients',xlab='',col='gray',...){
  p=length(x$beta)
  if (is.null(index)) index=1:p
  if (is.null(names)) names=paste('V',index,sep='')
  names=factor(names,levels=names)
  co=x$beta[index]/sqrt(sum(x$beta[index]^2))
  method=rep('Olearn',p)
  data1=data.frame(co,method,names)
  f<-ggplot(data1, aes(names, co))+
    ylab(ylab) +
    xlab(xlab) +
    #scale_x_discrete(labels=names1)
    geom_bar(stat="identity",width=0.5,fill=col) +
    facet_wrap(~method)+
    coord_flip() +
    #scale_fill_grey() +
    theme_bw(base_size=11)
  suppressWarnings(print(f))
}

plot.qlearn<-function(x,index=NULL,names=NULL,ylab='std coefficients',xlab='',col='gray',...){
  p=length(x$co)
  if (is.null(index)) index=(p/2+1):p
  if (is.null(names)) names=paste('V',index-p/2,sep='')
  names=as.factor(names)
  co=x$co[index]/sqrt(sum(x$co[index]^2))
  method=rep('Qlearning',p)
  data1=data.frame(co,method,names)
  f<-ggplot(data1, aes(names, co))+
    ylab(ylab) +
    xlab(xlab) +
    #scale_x_discrete(labels=names1)
    geom_bar(stat="identity",width=0.5,fill=col) +
    facet_wrap(~method)+
    coord_flip() +
    #scale_fill_grey() +
    theme_bw(base_size=11)
  suppressWarnings(print(f))
}




Plearning<-function(X,AA,RR,n,K,pi,pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=0.00001){
  select=matrix(1,n,1)
  QL=matrix(0,n,K)
  M=matrix(1,n,K)
  C=matrix(1,n,K)
  models=list()
  prob=matrix(1,n,K)
  QLproj=matrix(0,n,K+1)
  Qspecify=matrix(0,n,K)
  QR_future=0
  Rsum=0
  if (is.matrix(X)){
    for (k in K:1){
      A=AA[[k]]
      output_Q=Qlearning_Single(X,A,RR[[k]]+QR_future,pentype=pentype,m=m)
      QR_future=output_Q$Q
      #subsititute the outcome by expected outcome of best treatment
      QL[,k]=output_Q$Q
      if(k<K) R_p=Rsum*select/prob[,K]+apply(QLproj[,(k+1):K]%*%Qspecify[,(k+1):K],2,sum)
      if(k==K) R_p=Rsum*select/prob[,K]
      R=(RR[[k]]+R_p)
      if (kernel=='linear'){
        models[[k]]=Olearning_Single(X,A,R,pi[[k]],pentype=pentype,clinear=clinear,e=e,m=m)
      }else if (kernel=='rbf'){
        models[[k]]=Olearning_Single(X,A,R,pi[[k]],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
      }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))

      right=(sign(models[[k]]$fit)==A)
      #update fo next stage
      M[,k:K]=M[,k:K]*(right%*%rep(1,K-k+1))
      if (k>1) C[,k:K]=M[,k-1:K-1]-M[,k:K]
      if (k==1){
        C[,2:K]=M[,1:(K-1)]-M[,2:K]
        C[,1]=rep(1,n)-M[,1]
      }

      select=select*right
      prob[,k:K]=prob[,k:K]*(pi[[k]]*rep(1,K-k+1))
      Rsum=rep(1,n)
      for (j in k:K){
        if (j>1) {QLproj[,j]=(C[,j]-(1-pi[[j]])*M[,j-1])/prob[,j]
        } else QLproj[,1]=(C[,j]-(1-pi[[j]]))/prob[,j]

        Qspecify[,j]=QL[,j]+Rsum
        Rsum=Rsum+RR[[j]]
      }}}

  if (is.list(X)){
    for (k in K:1){
      A=AA[[k]]
      output_Q=Qlearning_Single(X[[k]],A,RR[[k]]+QR_future,pentype=pentype)
      QR_future=output_Q$Q
      #subsititute the outcome by expected outcome of best treatment
      QL[,k]=output_Q$Q
      if(k<K) R_p=Rsum*select/prob[,K]+apply(QLproj[,(k+1):K]%*%Qspecify[,(k+1):K],2,sum)
      if(k==K) R_p=Rsum*select/prob[,K]
      R=(RR[[k]]+R_p)
      if (kernel=='linear'){
        models[[k]]=Olearning_Single(X[[k]],A,R,pi[[k]],pentype=pentype)
      }else if (kernel=='rbf'){
        models[[k]]=Olearning_Single(X[[k]],A,R,pi[[k]],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
      }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))

      right=(sign(models[[k]]$fit)==A)
      #update fo next stage
      M[,k:K]=M[,k:K]*(right%*%rep(1,K-k+1))
      if (k>1) C[,k:K]=M[,k-1:K-1]-M[,k:K]
      if (k==1){
        C[,2:K]=M[,1:(K-1)]-M[,2:K]
        C[,1]=rep(1,n)-M[,1]
      }

      select=select*right
      prob[,k:K]=prob[,k:K]*(pi[[k]]*rep(1,K-k+1))
      Rsum=rep(1,n)
      for (j in k:K){
        if (j>1) {QLproj[,j]=(C[,j]-(1-pi[[j]])*M[,j-1])/prob[,j]
        } else QLproj[,1]=(C[,j]-(1-pi[[j]]))/prob[,j]

        Qspecify[,j]=QL[,j]+Rsum
        Rsum=Rsum+RR[[j]]
      }}}
  models
}


Olearning_Single<-function(H,A,R2,pi=rep(1,n),pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=1e-5,
                           ridge.reg =TRUE)
{
  npar=length(clinear)
  n=length(A)
  p=dim(H)[2]

  if (max(R2)!=min(R2)){
    if (pentype=='lasso'){
      cvfit=cv.glmnet(H,R2,nfolds=m)
      co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
    }else if (pentype=='LSE'){
      co=coef(lm(R2~H))
    }else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',
it can also be 'LSE' without penalization"))
    r=R2-cbind(rep(1,n),H)%*%co
  } else r=R2
  rand=sample(m,n,replace=TRUE)

  r=r/pi
  if (kernel=='linear'){
    V=matrix(0,m,npar)
    for (i in 1:m){
      this=(rand!=i)
      X=H[this,]
      Y=A[this]
      R=r[this]
      Xt=H[!this,]
      Yt=A[!this]
      Rt=r[!this]
      for (j in 1:npar){
        model=wsvm(X,Y,R,C=clinear[j],e=e, ridge.reg=ridge.reg)
        YP=predict(model,Xt)
        V[i,j]=sum(Rt*(YP==Yt))/sum(YP==Yt)
      }}
    mimi=colMeans(V)
    best=which.max(mimi)
    bestC= cbest=clinear[best]
    bestSig = NULL
    model=wsvm(H,A,r,C=cbest,e=e)}

  if (kernel=='rbf'){
    nsig=length(sigma)
    V=array(0,c(npar,nsig,m))
    for (i in 1:m){
      this=(rand!=i)
      X=H[this,]
      Y=A[this]
      R=r[this]
      Xt=H[!this,]
      Yt=A[!this]
      Rt=r[!this]
      for (j in 1:npar){
        for (s in 1:nsig){
          model=wsvm(X,Y,R,'rbf',sigma=sigma[s],C=clinear[j],e=e, ridge.reg=ridge.reg)
          YP=predict(model,Xt)
          V[j,s,i]=sum(Rt*(YP==Yt))/sum(YP==Yt)
        }}}
    mimi=apply(V,c(1,2),mean)
    best=which(mimi==max(mimi),arr.ind=TRUE)
    bestC=clinear[best[1]]
    bestSig=sigma[best[2]]
    print(bestC)
    print(bestSig)
    model=wsvm(H,A,r,'rbf', bestSig, C=bestC, e=e, ridge.reg=ridge.reg)}

  # HP modified this
  model$bestC <- bestC
  model$bestSig <- bestSig
  model
}

Olearning<-function(X,AA,RR,n,K,pi,pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=0.00001){
  select=rep(TRUE,n)
  R_future=0
  prob=rep(1,n)
  models=list()
  if (is.matrix(X)){
    for (j in K:1){
      R=(RR[[j]]+R_future)/(prob*pi[[j]])
      if (kernel=='linear'){
        models[[j]]=Olearning_Single(X[select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,clinear=clinear,e=e,m=m)
      }else if (kernel=='rbf'){
        models[[j]]=Olearning_Single(X[select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
      }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))



      select[which(select==1)]=(sign(models[[j]]$fit)==AA[[j]][select])
      R_future=R_future+RR[[j]]
      prob=prob*pi[[j]]
    }
  }

  if (is.list(X)){
    for (j in K:1){
      R=(RR[[j]]+R_future)/(prob*pi[[j]])
      if (kernel=='linear'){
        models[[j]]=Olearning_Single(X[[j]][select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,clinear=clinear,e=e,m=m)
      }else if (kernel=='rbf'){
        models[[j]]=Olearning_Single(X[[j]][select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
      }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))
      select[which(select==1)]=(sign(models[[j]]$fit)==AA[[j]][select])
      R_future=R_future+RR[[j]]
      prob=prob*pi[[j]]
    }
  }

  class(models)=paste('Olearning',pentype,sep='')
  # sumR=RR[[1]]
  # if (K>=2){
  # for (j in 2:K) sumR=sumR+RR[[j]]
  # }
  # value=sum(sumR*select)/sum(select)
  models
}



#this function makes simulation data set for single stage
make_classification<-function(n_cluster,pinfo,pnoise,n_sample,centroids=0){
  mu=numeric(pinfo)
  Sigma=diag(pinfo)
  if(is.vector(centroids)){
    centroids=5*mvrnorm(n_cluster,mu,Sigma)}

  X=matrix(0,n_sample,pinfo+pnoise)
  y=numeric(n_sample)
  z=numeric(n_sample)

  stopp=1
  A=2*rbinom(n_sample,1,0.5)-1


  for (k in 1:n_cluster){
    start=stopp
    stopp=stopp+n_sample/n_cluster
    y[start:(stopp-1)]=2*(k%%2)-1
    X[start:(stopp-1),1:pinfo]=matrix(1,n_sample/n_cluster,1)%*%centroids[k,]+mvrnorm(n_sample/n_cluster,mu,Sigma)
    z[start:(stopp-1)]=1.5*(2*(k%%2)-1)*A[start:(stopp-1)]+rnorm(n_sample/n_cluster)
  }

  mun=numeric(pnoise)
  Sigman=diag(pnoise)
  X[,pinfo+(1:pnoise)]=mvrnorm(n_sample,mun,Sigman)
  example=list(X=X,A=A,y=y,R=z,centroids=centroids)
  class(example)<-'examplesingle'
  return(example)
}



make_2classification<-function(n_cluster,pinfo,pnoise,n_sample,centroids=0){
  mu=numeric(pinfo)
  Sigma=diag(pinfo)
  if(is.vector(centroids)){
    centroids=5*mvrnorm(n_cluster,mu,Sigma)}

  X=matrix(0,n_sample,pinfo+pnoise)
  y=list()
  z=list()
  A=list()
  for (i in 1:2){
    y[[i]]=rep(1,n_sample)
    z[[i]]=rep(1,n_sample)
  }

  stopp=1
  A[[1]] = 2*rbinom(n_sample,1,0.5)-1
  A[[2]] = 2*rbinom(n_sample,1,0.5)-1



  for (k in 1:n_cluster){
    start=stopp
    stopp=stopp+n_sample/n_cluster
    y[[1]][start:(stopp-1)]=2*(k%%2)-1
    y[[2]][start:(stopp-1)]=2*(floor(k/2)%%2)-1
    X[start:(stopp-1),1:pinfo]=matrix(1,n_sample/n_cluster,1)%*%centroids[k,]+mvrnorm(n_sample/n_cluster,mu,Sigma)
    z[[2]][start:(stopp-1)]=y[[1]][start:(stopp-1)]*A[[1]][start:(stopp-1)]+y[[2]][start:(stopp-1)]*A[[2]][start:(stopp-1)]+rnorm(n_sample/n_cluster)
  }
  z[[1]]=rep(0,n_sample)
  mun=numeric(pnoise)
  Sigman=diag(pnoise)
  X[,pinfo+(1:pnoise)]=mvrnorm(n_sample,mun,Sigman)
  example=list(X=X,A=A,y=y,R=z,centroids=centroids)
  class(example)<-'example2'
  return(example)
}


