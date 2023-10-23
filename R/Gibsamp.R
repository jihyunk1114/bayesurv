#Y: individual status by group test result (length: sample size)
#Yi: group status (length: number of groups)
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
#make a note that k will be redefined when knots are specified.
#warning if both Y and Yi not provided
#' @export
Gibbs = function(Y = NULL, Yi = NULL, X, c, groupID, grid = NULL, n.grid = 100, k = 10,
                 init.theta = numeric(ncol(X)), eta = 1, gam0 = -3, gam = NULL,
                 theta0 = NULL, Sigma0 = NULL, m0 = 0, v0 = 0.1, a0 = 1, b0 = 1,
                 alpha = 1, beta = 1, order = 2, knots = NULL, quantile = TRUE,
                 maxiter = 6000, burn.in = 2000) {

  #make a note that k will be re-defined when knots are specified.
  if(is.null(knots)==F) {print("k will be redefined when knots are specified")}
  if(is.null(Y) && is.null(Yi)) stop("Warning: Both Y and Yi are not provided.")

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
  require(truncnorm); require(mvtnorm); require(dplyr)

  # when user do not specify knots or grid:
  if (is.null(knots)) {
    if (quantile == TRUE) {knots = c(0, quantile(c, seq(0.1, 0.9, length.out = k - order)), max(c) + 0.01)
    } else{knots = seq(0, max(c)+ 0.01, length.out = k + 2 - order)}
  } else{k = length(knots) - 2 + order}
  if (is.null(grid)) {grid = seq(0, max(c)+0.01, length.out = n.grid)}

  n.knots = length(knots); n.grid = length(grid)
  knots1 = seq(0, max(grid), length = n.knots)

  # if Yi, group test result with length as number of group, is present
  count = table(groupID); if (is.null(Y)) {Y = rep(Yi, count)}

  N = nrow(X); p = ncol(X); gamma = alpha + beta - 1
  I = Ispline(c, order, knots = knots); I1 = Ispline(grid, order, knots = knots1)
  dij = numeric(N); zij = dij

  #prior
  if(is.null(theta0)) {theta0 = numeric(ncol(X))}
  if(is.null(Sigma0)) {Sigma0 = nrow(X)*solve(t(X)%*%X)}

  #init
  theta = init.theta; gam = rep(0.1,k)

  # To save chains
  theta.mat = matrix(0, maxiter, k + p + 1)
  g.mat = matrix(0, maxiter, n.grid)

  for (iter in 1:maxiter) {
    u = t(gam0 + gam %*% I) + X %*% theta
    Yi = aggregate(Y, list(groupID), FUN = mean)[, 2]
    Di = Yi

    #Calculate probability
    m1 = aggregate(1 - pnorm(u), list(groupID), FUN = prod)[, 2]
    prY1 = alpha - gamma * m1; prD1 = 1 - m1
    prY0 = 1 - prY1; prD0 = 1 - prD1
    prob1 = alpha * prD1 / prY1; prob2 = 1 - beta * prD0 / prY0

    #Derive true group status
    Di[Yi == 1] = sapply(prob1[Yi == 1], function(x) {rbinom(1, 1, x)})
    Di[Yi == 0] = sapply(prob2[Yi == 0], function(x) {rbinom(1, 1, x)})
    groupsize = aggregate(Y, list(groupID), FUN = length)[, 2]
    D = rep(Di, times = groupsize)

    #Derive True individual status
    pij = pnorm(u) / rep(1 - m1, times = groupsize); pij[pij > 1] = 1
    dij = numeric(N); dij[D == 1] = sapply(pij[D == 1], function(x) {rbinom(1, 1, x)})

    index = rep(aggregate(dij[D == 1], list(groupID[D == 1]), FUN = sum)[, 2] == 0, times = groupsize[Di == 1])
    while (sum(index) > 0) {
      dij[D == 1][index] = sapply(pij[D == 1][index], function(x) {rbinom(1, 1, x)})
      index = rep(aggregate(dij[D == 1], list(groupID[D == 1]), FUN = sum)[, 2] == 0, times = groupsize[Di == 1])
    }

    zij[dij == 0] = sapply(u[dij == 0], function(x) {rtruncnorm(1,-Inf, 0, x, 1)})
    zij[dij == 1] = sapply(u[dij == 1], function(x) {rtruncnorm(1, 0, Inf, x, 1)})

    #gam0 posterior
    W0 = v0 + N; E0 = (v0 * m0 + sum(unlist(zij) - t(I) %*% gam - X %*% theta)) / W0
    gam0 = rnorm(1, E0, sqrt(1 / W0))

    #gam posterior
    W = rowSums(I ^ 2)
    for (l in 1:length(W)) {
      if (W[l] == 0) {
        gam[l] = rexp(1, eta)
      } else{
        El = (sum(I[l,] * (zij - gam0 - t(I[-l,]) %*% gam[-l] - X %*% theta)) - eta) / W[l]
        gam[l] = rtruncnorm(1, 0, Inf, El, sqrt(1 / W[l]))
      }
    }

    #theta posterior
    Sigmatilde = solve(solve(Sigma0) + t(X) %*% X)
    thetatilde = c(Sigmatilde %*% (solve(Sigma0) %*% theta0 + t(X) %*% (zij - gam0 - t(gam %*% I))))
    theta = c(rmvnorm(1, thetatilde, Sigmatilde))

    #eta posterior
    eta = rgamma(1, 1 + k, 1 + sum(gam))

    #Save Gibbs result
    theta.mat[iter, 1:p] = theta
    theta.mat[iter, p + 1] = gam0
    theta.mat[iter,-(1:(p + 1))] = gam

    g.mat[iter,] = t(gam0 + gam %*% I1)

    if (iter %% 100 == 0) {cat(paste("iteration", iter, "complete\n"))}

  }

  #Baseline survival function Plot
  St = colMeans(1 - pnorm(g.mat[-(1:burn.in),]))
  plot(grid, St, xlab = "t", ylab = expression(S[0](t)), lwd = 2, col = 2, type = "l")

  #Summary table
  mean = colMeans(theta.mat[-(1:burn.in), 1:p])
  q.025 = apply(theta.mat[, 1:p], 2, function(x) {quantile(x, 0.25)})
  q.975 = apply(theta.mat[, 1:p], 2, function(x) {quantile(x, 0.975)})
  rname = NULL
  for (ii in 1:p) {rname = c(rname, paste0("theta", ii))}
  summary.theta = as.data.frame(cbind(mean, q.025, q.975), row.names = rname)

  #Results to show
  return(list(theta.mat = theta.mat, theta.mean = mean, theta.q = c(q.025, q.975), alpha = alpha, beta = beta,
      g.mat = g.mat, grid = grid, St = St, summary.theta = summary.theta))

  #summary (list: mean/ quantiles / baseline cdf /grid and S(t)) # plot grid(t) and S(t)

}
