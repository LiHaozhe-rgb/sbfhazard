                  
SBF.MH.LC<-function(data,b.grid,x.grid=NULL,n.grid.additional=0, x.min=NULL, x.max=NULL, integral.approx='midd',it=100,kern=function(u){return(0.75*(1-u^2)*(abs(u)<1))},initial=NULL,convergence_tol=0.001,verbose=TRUE,phi=NULL)
{                               
  smooth.alpha<-function(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
  {
    alpha.smooth.i<-array(dim=c(d,n))
    
    alpha.smooth.i[1,] <-rep( 1,n)
    
    alpha.smooth.i.0<-(K.b/k.b)%*%(alpha[[1]]*dx[[1]])
    
    
    for (k in 2:d){
      for (i in 1:n)
      {
        alpha.smooth.i[k,i] <- as.numeric((K.X.b[[k]][i,]/k.X.b[[k]][i])%*%(alpha[[k]]*dx[[k]]))
      }}
    return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))
  }
  get.new.alpha<-function(data,alpha,K.X.b,k.X.b,K.b,k.b,Y,k,x.grid,dx,n.grid,d,n)
  {
    
    alpha.smooth.i<-smooth.alpha(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
    
    if (k==1) alpha.minusk.smooth<-numeric(n) else  alpha.minusk.smooth<-array(dim=c(n,n.grid[1])) 
    
    for (i in 1:n)
    {
      if (k==1) alpha.minusk.smooth[i] <-  prod(alpha.smooth.i$alpha.smooth.i[-1,i]) else
      {alpha.minusk.smooth[i,] <-  prod(alpha.smooth.i$alpha.smooth.i[-k,i])*alpha.smooth.i$alpha.smooth.i.0
      }
    } 
    
    if (k==1)
    {
      ###########################
      D<- rowSums(sapply(1:n, function (i) { return((dx[[1]]*(alpha.minusk.smooth[i]*Y[i,]))%*%(K.b/k.b))}))
    }else D<- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(alpha.minusk.smooth[i,])))*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    O<- rowSums(sapply(1:n, function (i) { return((data$status[i])*(K.X.b[[k]][i,]/k.X.b[[k]][i]))}))
    
    # print(paste(O))
    # print(paste(D))
    return(O/D)
  }
  ##### Note: The adjusted kernel,  K.b[k,b,,]/k.b[k,b,], has rowSums equal one.
  
  # K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
  # k.b<-array(0,dim=c(n.grid[1]))
  d<-ncol(data)-1
  n<-nrow(data)
  
  if(is.null(x.grid)) x.grid<-lapply(1:d,function(k) data[order(data[,k]),k])
 
  
  if(is.null(x.min)) x.min<-sapply(x.grid,head,1)
  if(is.null(x.max)) x.max<-sapply(x.grid,tail,1)
  
   x.grid2<-lapply(1:d, function(k) seq(x.min[k],x.max[k],length=n.grid.additional))
   x.grid<-lapply(1:d,function(k) sort(c(x.grid[[k]],x.grid2[[k]])))
  
   
   n.grid<-sapply(x.grid, length)
   
 
  
  
   if (integral.approx=='midd'){
                dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
                              c(  (ddx[1]/2)+(x.grid[[k]][1]-x.min[k])  ,(ddx[-(n.grid[k]-1)]+ddx[-1])/2,  
                                  (ddx[n.grid[k]-1]/2)+(x.max[k]-x.grid[[k]][n.grid[k]])  
                                )
                             }
                           )
  }
  
  if (integral.approx=='left'){
    dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
                                 c(  ddx,  x.max[k]-x.grid[[k]][n.grid[k]])  
                              }
             )
  }
  
  
  if (integral.approx=='right'){
    dx<-lapply(1:d,function(k){  ddx<-diff(x.grid[[k]])
                                  c(  (x.grid[[k]][1]-x.min[k])  ,ddx)
                                }
               )
  }
  
  x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 
   u<-x.grid.array-t(x.grid.array) # u is x_j-u_j for all x,u on the grid considered
  
  
  K.b<-apply(u/b.grid[1],1:2,kern)/(b.grid[1])
  k.b<-colSums(dx[[1]]*apply(u/b.grid[1],1:2,kern)/(b.grid[1]))
  # ## change
  # k.b <- as.vector(K.b %*% dx[[1]])
  
  
  X<-t(data[,-(d+1)])
  
  Y<-t(sapply(1:n,function(i) { temp<-numeric(n.grid[1])
                                 for (l in 1:n.grid[1]) {temp[l]<-as.numeric((x.grid[[1]][l]<=data$time[i]))
                                                          } 
                               return(temp) 
                              }
        ))
  

 
  K.X.b<-k.X.b<-list()
  for( k in 1:d){
    K.X.b[[k]]<-array(0,dim=c(n,n.grid[k]))
    k.X.b[[k]]<-numeric(n)
    for (i in 1:n)
    {
      u<-x.grid[[k]]
      u<-(X[k,i]-u)
      K.X.b[[k]][i,]<-sapply(u/b.grid[k],kern)/(b.grid[k])
      k.X.b[[k]][i]<-sum(dx[[k]]*sapply(u/b.grid[k],kern)/(b.grid[k]))
    } }
  
  
  alpha_backfit<-list()
  
  if (is.null(initial)){
  for(k in 1:d){
    alpha_backfit[[k]]<-rep(1, n.grid[k])
  }
  } else  alpha_backfit<-initial
 
   final_delta <- NA_real_
   converged <- FALSE
   for (l in 2:it)
  {  
    alpha_backfit_old<-alpha_backfit
    for(k in 1:d)
    {
  
      
      
      alpha_backfit[[k]]<- get.new.alpha(data,alpha_backfit,K.X.b,k.X.b,K.b,k.b,Y,k,x.grid,dx,n.grid,d,n)
                                       
      ###########
     # for (j in 1:(d-1))
     # {   alpha_backfit[[j]]<- 1*(alpha_backfit[[j]]/ alpha_backfit[[j]][1]) }
      
        
      alpha_backfit[[k]][is.nan(alpha_backfit[[k]])]<-0
      alpha_backfit[[k]][alpha_backfit[[k]]==Inf]<-0
      
    }
   
    final_delta <- max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE)
    if (is.finite(final_delta) && final_delta <= convergence_tol) {
      converged <- TRUE
      break
    }
    if (isTRUE(verbose)) print(c(l,final_delta))
  }

   if (d>1 && !is.null(phi)){  CorrSum<-1
     for (k in 1:(d-1))
     {  
       corr <-log(alpha_backfit[[k+1]][n.grid[k+1]/2])-phi[[k]](x.grid[[k+1]][n.grid[k+1]/2])
       alpha_backfit[[k+1]] <- alpha_backfit[[k+1]]/exp(corr) 
       alpha_backfit[[1]] <- alpha_backfit[[1]]*  exp(corr)     
       } 

   }

  eps <- 1e-12
  fit.diagnostics <- .sbf_kernel_norm_diagnostics(
    kernel_norm = list(k.b, k.X.b),
    floor_value = eps
  )
  k.b.safe <- pmax(as.numeric(k.b), eps)
  k.X.b.safe <- lapply(k.X.b, function(x) pmax(as.numeric(x), eps))
  alpha.smooth <- smooth.alpha(alpha_backfit,K.X.b,k.X.b.safe,K.b,k.b.safe,x.grid,dx,n.grid,d,n)
  alpha.smooth.obs <- alpha.smooth$alpha.smooth.i
  alpha.smooth.time <- as.numeric(alpha.smooth$alpha.smooth.i.0)

  dx1 <- as.numeric(dx[[1]])
  baseline.multiplier <- if (d > 1) {
    apply(alpha.smooth.obs[-1, , drop = FALSE], 2, prod)
  } else {
    rep(1, n)
  }
  baseline.exposure <- as.numeric(crossprod(baseline.multiplier, Y)) * dx1
  time.exposure <- as.numeric(Y %*% (alpha.smooth.time * dx1))

  covariate.scalar <- vector("list", d)
  if (d >= 2) {
    for (k in 2:d) {
      prod.minusk <- apply(alpha.smooth.obs[-k, , drop = FALSE], 2, prod)
      covariate.scalar[[k]] <- prod.minusk * time.exposure
    }
  }

  multiplicative.predict.state <- list(
    bandwidth = as.numeric(b.grid),
    kernel = kern,
    train_x = as.matrix(data[, seq_len(d), drop = FALSE]),
    dx = dx,
    k.b = k.b.safe,
    k.X.b = k.X.b.safe,
    status = as.numeric(data$status),
    baseline.exposure = baseline.exposure,
    covariate.scalar = covariate.scalar
  )

  return(list(
    alpha_backfit=alpha_backfit,
    l=l,
    x.grid=x.grid,
    data=data,
    fit_diagnostics=fit.diagnostics,
    multiplicative_predict_state=multiplicative.predict.state,
    converged=converged,
    iterations_used=l,
    iterations_max=it,
    final_delta=final_delta,
    convergence_tol=convergence_tol
  ))
}

do.sim<-function(n,b.grid,model=1,rho=0,d, distr,it,kern,seed.n)
{ 
  require(MASS)
  ptm <- proc.time()
  set.seed(seed.n)
  
  
  
  if (distr=='normal')
  {
    stddev<-rep(1,d-1)
    corMat<-matrix(rho, nrow=d-1,ncol=d-1)
    corMat[col(corMat)==row(corMat)]<-1
    covMat<-stddev %*% t(stddev) * corMat
    Z<-mvrnorm(n=n,  mu=rep(0,d-1), Sigma = covMat, empirical = FALSE)
    Z<-2.5*atan(Z)/pi
  } 
  
  if (distr=='uniform')
  {
    Z<-runif(n*(d-1),-1,1)
    Z<-matrix(Z,ncol=(d-1),nrow=n) 
  }
  
  phi<-list()  
  if (model==3)
  {
    for(k in 1:(d-1)){
      if ((k%%2)==1)  {phi[[k]]<- function(z) 0.2*log(z+2.5)} else {phi[[k]]<-function(z) 2*z}
    }
    
  }
  
  if (model==5)
  {  
    for(k in 1:(d-1)){
      phi[[k]]<-function(z) 2*z}
  }
  
  if (model==1)
  {
    for(k in 1:(d-1)){
      if ((k%%2)==1)  {phi[[k]]<- function(z) -z} else {phi[[k]]<-function(z) 2*z}
    }
  }
  
  if (model==4)
  {  
    for(k in 1:(d-1)){
      if ((k%%2)==1)  {phi[[k]]<- function(z) z^2} else {phi[[k]]<-function(z) 2*z}
    }
  }
  
  if (model==2)
  {  
    for(k in 1:(d-1)){
      if ((k%%2)==1)  {phi[[k]]<- function(z) 2*sin(pi*z)} else {phi[[k]]<-function(z) 2*z}
    }
  }
  
  top<-0
  for (k in 1:(d-1))
  {
    top<-  top+phi[[k]](Z[,k])
  }
  Time<-rexp(n,exp(top))
  C<-rexp(n,exp(top)/1.75)
  TT<-Time*(Time<=C)+C*(Time>C)
  status<-(Time<=C)*1  ## censoring indicator
  
  
  
  
  ################# Data generation 2: Comaprison with Lin.
  data<-data.frame(cbind(TT,Z))
  data$status<-status
  names(data)[1]<-'time'
  name<-names(data)
  
  ptm <- proc.time()
  alpha_backfit<-SBF.MH.LC(data,b.grid,it=it,phi=phi)
  comp.time.sbf<-(proc.time() - ptm)[[3]]
  
  l.sbf<-alpha_backfit$l
  x.grid<-alpha_backfit$x.grid
  n.grid<-sapply(x.grid, length)
  alpha_backfit<-alpha_backfit$alpha_backfit
  
  return(list(SBF=alpha_backfit,l.sbf=l.sbf,comp.time.sbf=comp.time.sbf,x.grid=x.grid,n.grid=n.grid, phi=phi))
}

# 
# d<-3
# rho<-0
# n<-200
# model<-2
# b.grid<-numeric(d)
# b.grid[1]<-30
# b.grid[2:d]<-rep(0.3,d-1)
# distr<-'uniform'
# it=10
# 
# seed.n=round(runif(1,0,999999999))
# 
# result<-do.sim(n, b.grid, model,rho,d, distr, it, seed.n=seed.n) 
# 
# k=1
# plot(result$x.grid[[k+1]],result$phi[[k]](result$x.grid[[k+1]]),lty=3,col=1,lwd=2) 
# lines(result$x.grid[[k+1]], log(result$SBF[[k+1]]),col='red',lwd=2)
# 
# 
# k=2
# plot(result$x.grid[[k+1]],result$phi[[k]](result$x.grid[[k+1]]),lty=3,col=1,lwd=2) 
# lines(result$x.grid[[k+1]], log(result$SBF[[k+1]]),col='red',lwd=2)
