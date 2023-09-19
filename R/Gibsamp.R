#Y: observed group status (length: sample size)
#Yi: observed group status (length: group number)
#X: covariates in matrix form
#c: test time
#grid: can define grid and find corresponding baseline survival function
#theta0, Sigma0, m0, v0, a0, b0: For priors
#alpha: Sensitivity
#beta: Specificity
#order: order for knots (usually 2 or 3)
#knots: If NULL default, can set up own knots
#maxiter: maximum interation number
#burn.in: burn-in number
#' @export

Gibbs=function(Y=NULL,Yi=NULL,X=X,c=c,grid=NULL,groupID,theta,theta0=numeric(ncol(X)),Sigma0=length(Y)*solve(t(X)%*%X),
               m0=0,v0=0.1,a0=1,b0=1,alpha,beta,order=2,knots=NULL,maxiter=6000,burn.in=2000){
  #data type Y or Yi consider both.

  #Ispline function
  Ispline<-function(x,order,knots){
    # get I spline matrix with order
    # x is a row vector
    # k is the order of I spline
    # knots are a sequence of increasing points
    # the number of free parameters in M spline is the length of knots plus 1.


    m=length(knots)
    ### get Ispline bases ###
    k=order+1
    n=m-2+k # number of parameters
    t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

    yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
    for (l in k:n){
      yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
    }

    yytem1=yy1
    for (ii in 1:order){
      yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
      for (i in (k-ii):n){
        yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
      }
      yytem1=yytem2
    }


    index=rep(0,length(x))
    for (i in 1:length(x)){
      index[i]=sum(t<=x[i])
    }

    ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

    if (order==1){
      for (i in 2:n){
        ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
      }
    }else{
      for (j in 1:length(x)){
        for (i in 2:n){
          if (i<(index[j]-order+1)){
            ibases[i-1,j]=1
          }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
          }else{
            ibases[i-1,j]=0
          }
        }
      }
    }
    return(ibases)
  }

  #package
  require(truncnorm);require(mvtnorm);require(dplyr);

  if(is.null(knots)) {k=10; knots=c(0,quantile(c,seq(0.1,0.9,length.out=k-order)),max(c)+0.01)}
  if(is.null(grid)) {grid=knots}
  knots1=seq(0,max(grid),length=length(knots))
  count=table(groupID)
  if(is.null(Y)) {Y=rep(Yi,count[,2])}

  N=nrow(X); p=ncol(X); gamma=alpha+beta-1; k=length(knots)-2+order
  I=Ispline(c,order,knots=knots); dij=numeric(N); zij=dij; THETA.sample=matrix(0,maxiter,k+p+1)
  I1=Ispline(grid,order,knots=knots1);ST=matrix(0,maxiter,length(grid))
  #ST1=ST;ST2=ST;ST3=ST # erase

  #prior
  eta=1; gam=rep(0.1,k); gam0=-3

  for(iter in 1:maxiter){
    u=t(gam0+gam%*%I)+X%*%theta
    Yi=aggregate(Y,list(groupID),FUN=mean)[,2];Di=Yi

    #Calculate probability
    m1=aggregate(1-pnorm(u),list(groupID),FUN=prod)[,2]
    prY1=alpha-gamma*m1; prD1=1-m1; prY0=1-prY1; prD0=1-prD1
    prob1=alpha*prD1/prY1; prob2=1-beta*prD0/prY0

    #Derive True group status
    Di[Yi==1]=sapply(prob1[Yi==1],function(x){rbinom(1,1,x)})
    Di[Yi==0]=sapply(prob2[Yi==0],function(x){rbinom(1,1,x)})
    groupsize=aggregate(Y,list(groupID),FUN=length)[,2]
    D=rep(Di,times=groupsize)

    #Derive True individual status
    pij=pnorm(u)/rep(1-m1,times=groupsize); pij[pij>1]=1
    dij=numeric(N); dij[D==1]=sapply(pij[D==1],function(x){rbinom(1,1,x)})

    index=rep(aggregate(dij[D==1],list(groupID[D==1]),FUN=sum)[,2]==0,times=groupsize[Di==1])
    while(sum(index)>0){
      dij[D==1][index]=sapply(pij[D==1][index],function(x){rbinom(1,1,x)})
      index=rep(aggregate(dij[D==1],list(groupID[D==1]),FUN=sum)[,2]==0,times=groupsize[Di==1])}

    zij[dij==0]=sapply(u[dij==0],function(x){rtruncnorm(1,-Inf,0,x,1)})
    zij[dij==1]=sapply(u[dij==1],function(x){rtruncnorm(1,0,Inf,x,1)})

    #gam0 posterior
    W0=v0+N; E0=(v0*m0+sum(unlist(zij)-t(I)%*%gam-X%*%theta))/W0
    gam0=rnorm(1,E0,sqrt(1/W0))

    #gam posterior
    W=rowSums(I^2)
    for(l in 1:length(W)){
      if(W[l]==0){
        gam[l]=rexp(1,eta)
      } else{
        El=(sum(I[l,]*(zij-gam0-t(I[-l,])%*%gam[-l]-X%*%theta))-eta)/W[l]
        gam[l]=rtruncnorm(1,0,Inf,El,sqrt(1/W[l]))
      }
    }

    #theta posterior
    Sigmatilde=solve(solve(Sigma0)+t(X)%*%X)
    thetatilde=c(Sigmatilde%*%(solve(Sigma0)%*%theta0+t(X)%*%(zij-gam0-t(gam%*%I))))
    theta=c(rmvnorm(1,thetatilde,Sigmatilde))

    #eta posterior
    eta=rgamma(1,1+k,1+sum(gam))

    #Save Gibbs result
    THETA.sample[iter,1:p]=theta
    THETA.sample[iter,p+1]=gam0
    THETA.sample[iter,-(1:(p+1))]=gam

    ST[iter,]=1-pnorm(t(gam0+gam%*%I1))
    #ST1[iter,]=1-pnorm(c(theta%*%c(1,0))+t(gam0+gam%*%I1))
    #ST2[iter,]=1-pnorm(c(theta%*%c(0,1))+t(gam0+gam%*%I1))
    #ST3[iter,]=1-pnorm(c(theta%*%c(1,1))+t(gam0+gam%*%I1))

    if( iter %% 100 == 0 ) {cat(paste("iteration", iter, "complete\n"));print(theta)}

  }

  #Throw away burn-in
  THETA.sample=THETA.sample[-(1:burn.in),]# HDI=apply(THETA.sample[,1:p],2,hdi)

  #Results to show
  return(list(THETA.sample=THETA.sample, mean=colMeans(THETA.sample[,1:p]),
              q.025=apply(THETA.sample[,1:p],2,function(x){quantile(x,0.25)}),
              q.975=apply(THETA.sample[,1:p],2,function(x){quantile(x,0.975)}),
              grid=grid, St=colMeans(ST),# St1=colMeans(ST1), St2=colMeans(ST2), St3=colMeans(ST3),
              plot=plot(grid, colMeans(ST), xlab="t", ylab=expression(S[0](t)),lwd=2, col=2, type="l")))

  #summary (list: mean/ quantiles / baseline cdf /grid and S(t)) # plot grid(t) and S(t)

}

GibbsInd=function(Y=NULL,X=X,c=c,grid=NULL,theta,theta0=numeric(ncol(X)),Sigma0=length(Y)*solve(t(X)%*%X),
               m0=0,v0=0.1,a0=1,b0=1,alpha,beta,order=2,knots=NULL,maxiter=6000,burn.in=2000){
  #data type Y or Yi consider both.

  #Ispline function
  Ispline<-function(x,order,knots){
    # get I spline matrix with order
    # x is a row vector
    # k is the order of I spline
    # knots are a sequence of increasing points
    # the number of free parameters in M spline is the length of knots plus 1.


    m=length(knots)
    ### get Ispline bases ###
    k=order+1
    n=m-2+k # number of parameters
    t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

    yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
    for (l in k:n){
      yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
    }

    yytem1=yy1
    for (ii in 1:order){
      yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
      for (i in (k-ii):n){
        yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
      }
      yytem1=yytem2
    }


    index=rep(0,length(x))
    for (i in 1:length(x)){
      index[i]=sum(t<=x[i])
    }

    ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

    if (order==1){
      for (i in 2:n){
        ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
      }
    }else{
      for (j in 1:length(x)){
        for (i in 2:n){
          if (i<(index[j]-order+1)){
            ibases[i-1,j]=1
          }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
          }else{
            ibases[i-1,j]=0
          }
        }
      }
    }
    return(ibases)
  }

  #package
  require(truncnorm);require(mvtnorm);require(dplyr);

  if(is.null(knots)) {k=10; knots=c(0,quantile(c,seq(0.1,0.9,length.out=k-order)),max(c)+0.01)}
  if(is.null(grid)) {grid=knots}

  N=nrow(X); p=ncol(X); gamma=alpha+beta-1; k=length(knots)-2+order
  I=Ispline(c,order,knots=knots); dij=numeric(N); zij=dij; THETA.sample=matrix(0,maxiter,k+p+1)
  I1=Ispline(grid,order,knots=knots);ST=matrix(0,maxiter,length(grid))

  #prior
  eta=1; gam=rep(0.1,k); gam0=-3

  for(iter in 1:maxiter){
    u=t(gam0+gam%*%I)+X%*%theta

    #Calculate probability
    prob1=alpha*pnorm(u)/(alpha-gamma*(1-pnorm(u)))
    prob2=1-beta*(1-pnorm(u))/(1-alpha+gamma*(1-pnorm(u)))

    dij[Y==1]=sapply(prob1[Y==1],function(x){rbinom(1,1,x)})
    dij[Y==0]=sapply(prob2[Y==0],function(x){rbinom(1,1,x)})

    zij[dij==0]=sapply(u[dij==0],function(x){rtruncnorm(1,-Inf,0,x,1)})
    zij[dij==1]=sapply(u[dij==1],function(x){rtruncnorm(1,0,Inf,x,1)})

    zij=unlist(zij)

    #gam0 posterior
    W0=v0+N; E0=(v0*m0+sum(unlist(zij)-t(I)%*%gam-X%*%theta))/W0
    gam0=rnorm(1,E0,sqrt(1/W0))

    #gam posterior
    W=rowSums(I^2)
    for(l in 1:length(W)){
      if(W[l]==0){
        gam[l]=rexp(1,eta)
      } else{
        El=(sum(I[l,]*(zij-gam0-t(I[-l,])%*%gam[-l]-X%*%theta))-eta)/W[l]
        gam[l]=rtruncnorm(1,0,Inf,El,sqrt(1/W[l]))
      }
    }

    #theta posterior
    Sigmatilde=solve(solve(Sigma0)+t(X)%*%X)
    thetatilde=c(Sigmatilde%*%(solve(Sigma0)%*%theta0+t(X)%*%(zij-gam0-t(gam%*%I))))
    theta=c(rmvnorm(1,thetatilde,Sigmatilde))

    #eta posterior
    eta=rgamma(1,1+k,1+sum(gam))

    #Save Gibbs result
    THETA.sample[iter,1:p]=theta
    THETA.sample[iter,p+1]=gam0
    THETA.sample[iter,-(1:(p+1))]=gam

    ST[iter,]=1-pnorm(t(gam0+gam%*%I1))

    if( iter %% 100 == 0 ) {cat(paste("iteration", iter, "complete\n"));print(theta)}

  }

  #Throw away burn-in
  THETA.sample=THETA.sample[-(1:burn.in),1:p]# HDI=apply(THETA.sample[,1:p],2,hdi)

  #Results to show
  return(list(THETA.sample=THETA.sample, mean=colMeans(THETA.sample[,1:p]),
              q.025=apply(THETA.sample[,1:p],2,function(x){quantile(x,0.25)}),
              q.975=apply(THETA.sample[,1:p],2,function(x){quantile(x,0.975)}),
              grid=grid, St=colMeans(ST),
              plot=plot(grid, colMeans(ST), xlab="t", ylab=expression(S[0](t)),lwd=2, col=2, type="l")))

  #summary (list: mean/ quantiles / baseline cdf /grid and S(t)) # plot grid(t) and S(t)

}
