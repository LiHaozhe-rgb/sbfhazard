
SBF.MH.CLL<-function(formula,
                     data,bandwidth,
                     weight='sw',
                     x.grid=NULL, n.grid.additional=0,
                     x.min=NULL, x.max=NULL,
                     integral.approx='right',
                     it=100,
                     kern=function(u){return(0.75*(1-u^2)*(abs(u)<1))},
                     initial=NULL,
                     convergence_tol=0.001,
                     kcorr=kcorr,
                     LC,
                     wrong,
                     classic.backfit=FALSE, print=FALSE){  
  

# 
# formula <- formula(formula)
# if (class(formula) != "formula") {
#   stop("Error: Invalid formula.")
# }
# data.selected <- as.list(attr(terms(frmla), "variables"))[-1]
#   
"div" <- function(x,y) ifelse(y==0&x==0,0,base:::"/"(x,y))

Terms <- terms(x=formula,data=data)
mm <- na.omit(get_all_vars(formula(Terms),data=data))
if (NROW(mm) == 0) stop("No (non-missing) observations")

response <- model.response(model.frame(update(formula,".~1"),data=mm))
X       <- prodlim::model.design(Terms,
                                 data=mm,
                                 maxOrder=1,
                                 dropIntercept=TRUE)[[1]]


time <- as.vector(response[,"time"])
status <- as.vector(response[,"status"])
X <- cbind( time,X)


## Grid construction

d <- ncol(X) 
n <- nrow(X)

if(is.null(x.grid)) x.grid<-lapply(1:d,function(k) X[order(X[,k]),k])
if(is.null(x.min)) x.min<-sapply(x.grid,head,1)
if(is.null(x.max)) x.max<-sapply(x.grid,tail,1)



x.grid.additional <- lapply(1:d, function(k) seq(x.min[k],x.max[k], length=n.grid.additional))
x.grid <- lapply(1:d, function(k) sort(c(x.grid[[k]], x.grid.additional[[k]])))
n.grid <- sapply(x.grid, length)

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



## Kernel weighting and exposure process

K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
k.b<-array(0,dim=c(n.grid[1]))



x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 


u<- x.grid.array-t(x.grid.array) # u[t,s]= t-s for all t,s on the grid considered
K.b<-apply(u/bandwidth[1],1:2,kern)/(bandwidth[1])
k.b<-colSums(dx[[1]]*K.b)                      # small k for normailzing kernel, sincs grid points are symmetric normalization can be used for both row and column
if (kcorr==FALSE) {k.b<-rep(1,n.grid[1])}                     


dX0.b<- u#/bandwidth[1]     # dX0.b[t,s]= t-s for all t,s on the grid considered
K.b <- K.b %*% diag(1/k.b)  ### row-wise division --> colSums(dx[[1]]*K.b)=1
K.b[is.na(K.b)] <-0
### define exposure=Y[i,s]

Y<-t(sapply(1:n,function(i) { temp<-numeric(n.grid[1])
for (l in 1:n.grid[1]) {temp[l]<-as.numeric((x.grid[[1]][l]<=time[i]))
} 
return(temp) 
}
))



dX.b<-K.X.b<-k.X.b<-list()
for( k in 1:d){
  k.X.b[[k]]<-numeric(n)
  
  x.grid.array<-matrix(-X[,k],nrow =n.grid[k], ncol =n ,byrow=TRUE)   # 
  u<- x.grid.array+x.grid[[k]]      ####  u[,i]=x_k - X_{ik}
  K.X.b[[k]] <- apply(u/bandwidth[k],1:2,kern)/(bandwidth[k])
  k.X.b[[k]]<-colSums(dx[[k]]*K.X.b[[k]])   
  if (kcorr==FALSE) k.X.b[[k]] <- rep(1,n)
  dX.b[[k]]<- u               #### dX.b[[k]] [,i] = x_k-X_{ik} 
  # K.X.bC[[k]] <- K.X.b[[k]]/ k.X.b[[k]]
  K.X.b[[k]]<- K.X.b[[k]] %*%  diag(1/k.X.b[[k]])  ### row-wise division ---> col integration=1
  #(dx[[k]]*K.X.b[[k]]) =1
  K.X.b[[k]][is.na(K.X.b[[k]])] <-0
  
}


## Kernel-smoothed occurrences and other variables

get.D<-function(k,dx,K.X.b,K.b,Y){
  # V_{0,0}^j, the denominator or eq(31)
  if (k==1)
  {
    D <-     colSums( Y    %*%  (t(K.b) *dx[[1]])  ) 
    
  }else 
    D <-   K.X.b[[k]]  %*%  Y   %*%   as.matrix(dx[[1]]) 
  
  return(D)
}

get.O1<-function(k,status,K.X.b){
  # U_0^j in the numerator of eq(31), the first term
  O1 <- as.numeric(  K.X.b[[k]] %*% status )
  
  return(O1)
}


get.O3<-function(k,dx,K.X.b,K.b,Y){
  # V_{j,0}^j in the numerator of eq(31), the second term
  if (k==1)
  {
    O3 <- colSums( Y    %*%  (t(dX0.b*K.b) *dx[[1]])  )  
    
  }else O3 <- (( dX.b[[k]]*K.X.b[[k]] ) %*%  Y   %*%  as.matrix(dx[[1]]) ) 
  
  return(O3) 
}

get.D.1<-function(k,dx,K.X.b,K.b,dX.b,dX0.b,Y){
  
  if (k==1)
  {
    D <-     colSums( Y    %*%  (t(dX0.b^2*K.b) *dx[[1]])  ) 
    
  }else
    D <-   (dX.b[[k]]^2*K.X.b[[k]])  %*%  Y   %*%  as.matrix(dx[[1]] )
  return(D)
}


get.O1.1<-function(k,status,K.X.b,dX.b){
  O1<- as.numeric((dX.b[[k]]*K.X.b[[k]]) %*% status)
  return(O1)
}


get.O3.1<-function(k,dx,K.X.b,K.b,dX.b,dX0.b,Y){
  
  if (k==1)
  {
    O3 <- colSums( Y    %*%  (t(dX0.b*K.b)*dx[[1]])  )  
  }else O3 <- ( dX.b[[k]]*K.X.b[[k]] ) %*%  Y    %*%  as.matrix(dx[[1]]) 
  
  return(O3)
}







taylor.alpha<-function(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
  # the last two terms without taking the sum in the numerator of eq(31)
{
  
  taylor.alpha.i <- array(dim=c(n,d))
  
  taylor.alpha.i[ ,1] <- rep(0,n)
  
  # the wrong variable indicates whether the difference is the same order as alpha DL
  if(wrong==TRUE) taylor.alpha.i.0 <-   as.numeric( t(as.matrix(dx[[1]]))   %*%   (( alpha[[1]]+(t(dX0.b)*alpha.1[[1]] ))*K.b) )
  
  if(wrong==FALSE) taylor.alpha.i.0 <-   as.numeric( t(as.matrix(dx[[1]]))   %*%   (( alpha[[1]]+(dX0.b*alpha.1[[1]] ))*K.b) )
  
 
  for (k in 2:d){
    
    taylor.alpha.i[,k] <-    as.numeric( t(as.matrix(dx[[k]])) %*% ((alpha[[k]]+alpha.1[[k]]*dX.b[[k]])* K.X.b[[k]])  )
    
  }
  
  
  return(list( taylor.alpha.i= taylor.alpha.i,taylor.alpha.i.0=taylor.alpha.i.0))
}

get.alpha.new<-function(O1,O3,D,alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
  # taking the sum in the numerator of eq(31)
{

  if(classic.backfit==FALSE){
  
    taylor.alpha <- taylor.alpha(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    
    if (k==1) taylor.alpha.minusk <-rowSums(cbind(taylor.alpha$taylor.alpha.i[,-1],0)) else  
    {taylor.alpha.minusk <-  rowSums(cbind(taylor.alpha$taylor.alpha.i[,-k],0)) + matrix(taylor.alpha$taylor.alpha.i.0, ncol=n.grid[[1]], nrow=n, byrow = TRUE)
    }
  }
  
  
  if(classic.backfit==TRUE){
    temp <- array(dim=c(n,d)) 
    for (j in 1:d){
     index <- sapply(1:length(X[,j]), function(i) which.min(abs(X[i,j]-x.grid[[j]])))
     temp[,j] <- alpha[[j]][index]
    }
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
   
    if (k==1) taylor.alpha.minusk <-    rowSums(cbind(temp[,-1],0))  else  
    {taylor.alpha.minusk <-  rowSums(cbind(temp[,-k],0)) + matrix(alpha[[1]], ncol=n.grid[[1]], nrow=n, byrow = TRUE)
    }
  }
    
    
    

  
    
    if (k==1)
    { # the sum of last two terms
      O2 <- colSums( (taylor.alpha.minusk*Y)    %*%  (t(K.b) *dx[[1]])  )  ### transpose k or not? old: not

    }else 
      O2 <-  K.X.b[[k]]  %*%  (taylor.alpha.minusk*Y)    %*%  as.matrix(dx[[1]]) 
    
    
    

    O <- O1[[k]] - alpha.1[[k]]*O3[[k]]-O2


  
  
  return(as.numeric(div(O,D[[k]])))
} 

get.alpha.1.new<-function(O1.1,O3.1,D.1,alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)
{    
  if(classic.backfit==FALSE){


    taylor.alpha <- taylor.alpha(alpha,alpha.1,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    
    
    
    
    if (k==1) taylor.alpha.minusk <-  rowSums(cbind(taylor.alpha$taylor.alpha.i[,-1],0)) else  
    {taylor.alpha.minusk <- rowSums(cbind(taylor.alpha$taylor.alpha.i[,-k],0)) + matrix(taylor.alpha$taylor.alpha.i.0, ncol=n.grid[[1]], nrow=n, byrow = TRUE)
    }
    
  }

    if(classic.backfit==TRUE){
      temp <- array(dim=c(n,d)) 
      for (j in 1:d){
        index <- sapply(1:length(X[,j]), function(i) which.min(abs(X[i,j]-x.grid[[j]])))
        temp[,j] <- alpha[[j]][index]
      }
      
      if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
      
      if (k==1) taylor.alpha.minusk <-    rowSums(cbind(temp[,-1],0))  else  
      {taylor.alpha.minusk <-  rowSums(cbind(temp[,-k],0)) + matrix(alpha[[1]], ncol=n.grid[[1]], nrow=n, byrow = TRUE)
      }
    }
    
    
    if (k==1)
    { 
      O2.1 <- colSums( (taylor.alpha.minusk*Y)    %*%  (t(dX0.b*K.b) *dx[[1]])  )  
    }else 
      O2.1 <- ( dX.b[[k]]*K.X.b[[k]] ) %*%  (taylor.alpha.minusk*Y)    %*%  as.matrix(dx[[1]])
    


    
    O.1 <- O1.1[[k]]-alpha[[k]]*O3.1[[k]]-O2.1
    
  
  
    return(as.numeric(div(O.1,D.1[[k]])))
} 


D<-O1<-O3<-D1.1<-O1.1<-O3.1<-list()

D <- lapply(1:d, get.D, dx,K.X.b,K.b,Y)
O1 <- lapply(1:d, get.O1, status,K.X.b)
O3 <- lapply(1:d, get.O3, dx,K.X.b,K.b,Y)

D.1 <- lapply(1:d, get.D.1, dx,K.X.b,K.b,dX.b,dX0.b,Y)
O1.1 <- lapply(1:d, get.O1.1,status,K.X.b,dX.b)
O3.1 <- lapply(1:d, get.O3.1,dx,K.X.b,K.b,dX.b,dX0.b,Y)

alpha_backfit<-list()
alpha.1_backfit<-list()

if (is.null(initial)){
  for(k in 1:d){
    alpha_backfit[[k]]<-rep(0, n.grid[k])
  }
} else  alpha_backfit<-initial

for(k in 1:d){
  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
  
}


count<-rep(1,d)
final_delta <- NA_real_
converged <- FALSE
for (l in 2:it)
{  
  
  
  
  alpha_backfit_old<-alpha_backfit
  alpha.1_backfit_old<-alpha.1_backfit
  for(k in 1:d)
  {
  
  
    alpha_backfit[[k]]<- get.alpha.new(O1,O3,D,alpha_backfit,alpha.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n,l)
    
  
     #alpha_backfit[[k]][is.nan(alpha_backfit[[k]])]<-0
    #alpha_backfit[[k]][is.na(alpha_backfit[[k]])]<-0
     # alpha_backfit[[k]][alpha_backfit[[k]]==Inf]<-0
     # alpha_backfit[[k]][alpha_backfit[[k]]==-Inf]<-0
  # alpha_backfit[[k]][alpha_backfit[[k]]<0]<-0.0001
    #  }
    #  for(k in 1:d)
    #   {
    #    if (k!=1) {
    
    #     }
    if (LC==TRUE) alpha.1_backfit[[k]]<-rep(0,n.grid[k]) else{
      alpha.1_backfit[[k]]<- get.alpha.1.new(O1.1,O3.1,D.1,alpha_backfit,alpha.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,Y,k,x.grid,n.grid,d,n)}
    alpha.1_backfit[[k]][is.nan(alpha.1_backfit[[k]])]<-0
    alpha.1_backfit[[k]][is.na(alpha.1_backfit[[k]])]<-0
    # alpha.1_backfit[[k]][alpha.1_backfit[[k]]>30]<-30
    #  alpha.1_backfit[[k]][alpha.1_backfit[[k]]<=(-30)]<--30
    #   print(!(prod(alpha.1_backfit_old[[k]] == 0) ))
    
    # 
    # if ((prod(alpha.1_backfit_old[[k]] == 0) ))  {### if old vaues all zero
    #   if (max(abs(alpha.1_backfit_old[[k]]- alpha.1_backfit[[k]])) >(100))    {alpha.1_backfit[[k]] <-rep(0, n.grid[k])
    #   #  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
    #   count[k]<-1
    #   }}
    # 
    # if (!(prod(alpha.1_backfit_old[[k]] == 0) )){ ### if old vaues NOT all zero
    #   if (max(abs(alpha.1_backfit_old[[k]]- alpha.1_backfit[[k]])) >(100/log(count[k])))    alpha.1_backfit[[k]] <-rep(0, n.grid[k])
    #   #  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
    # }
   if(k!=1){
      alpha_backfit[[1]]<-alpha_backfit[[1]]+ as.numeric((t(alpha_backfit[[k]])%*%dx[[k]])/sum(dx[[k]]))
      alpha_backfit[[k]]<-alpha_backfit[[k]]- as.numeric((t(alpha_backfit[[k]])%*%dx[[k]])/sum(dx[[k]]))
   }
    
    count[k]<-count[k]+1
    # plot(x.grid[[k]],alpha.1_backfit[[k]],lty=3,col=1,lwd=2) 
    

  }
  #}
  
  
  
  #plot(x.grid[[2]],phi[[1]](x.grid[[2]]),lty=3,col=1,lwd=2) 
  #  if (l==2) plot(x.grid[[2]],log(alpha_backfit[[2]]),ylim=c(-2,2) ) else lines(x.grid[[2]],log(alpha_backfit[[2]] ) )
  #plot(x.grid[[1]], log(alpha_backfit[[1]]),col='blue',lwd=2) 
  
  
    
    
  final_delta <- max(max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE),max(abs(unlist(alpha.1_backfit_old)-unlist(alpha.1_backfit)),na.rm=TRUE))
  if (is.finite(final_delta) && final_delta <= convergence_tol) {
    converged <- TRUE
    break
  }
  for(k in 1:d){
   if (print==TRUE) print(c(l,max(abs(unlist((alpha_backfit_old[[k]]))-unlist((alpha_backfit[[k]]))),na.rm=TRUE),max(abs(unlist((alpha.1_backfit_old[[k]]))-unlist((alpha.1_backfit[[k]]
    ))),na.rm=TRUE)))
  }
}



for(k in 2:d){
  alpha_backfit[[1]]<-alpha_backfit[[1]]+ as.numeric((t(alpha_backfit[[k]])%*%dx[[k]])/sum(dx[[k]]))
  alpha_backfit[[k]]<-alpha_backfit[[k]]- as.numeric((t(alpha_backfit[[k]])%*%dx[[k]])/sum(dx[[k]]))
}

Y.colsum <- colSums(Y)
Y.dx.1 <- as.numeric(Y %*% as.matrix(dx[[1]]))
taylor.terms <- taylor.alpha(alpha_backfit,alpha.1_backfit,dx,K.X.b,K.b,dX.b,dX0.b,x.grid,n.grid,d,n)
row.total <- rowSums(taylor.terms$taylor.alpha.i)
nuisance.baseline <- as.numeric(crossprod(row.total, Y))
nuisance.covariate <- vector("list", d)
if (d >= 2) {
  weighted.Y <- sweep(Y, 2, dx[[1]], "*")
  time.nuisance <- as.numeric(weighted.Y %*% taylor.terms$taylor.alpha.i.0)
  for (k in 2:d) {
    nuisance.covariate[[k]] <- (row.total - taylor.terms$taylor.alpha.i[, k]) * Y.dx.1 + time.nuisance
  }
}

fit.data <- data.frame(time = time, X[, -1, drop = FALSE], status = status, check.names = FALSE)
fit.diagnostics <- .sbf_kernel_norm_diagnostics(
  kernel_norm = list(k.b, k.X.b),
  floor_value = 1e-12
)
additive.predict.state <- list(
  bandwidth = as.numeric(bandwidth),
  kernel = kern,
  kernel_correction = isTRUE(kcorr),
  local_constant = isTRUE(LC),
  X = X,
  dx = dx,
  k.b = pmax(as.numeric(k.b), 1e-12),
  k.X.b = lapply(k.X.b, function(x) pmax(as.numeric(x), 1e-12)),
  Y.colsum = Y.colsum,
  Y.dx.1 = Y.dx.1,
  status = status,
  nuisance.baseline = nuisance.baseline,
  nuisance.covariate = nuisance.covariate
)

return(list(
  alpha_backfit=alpha_backfit,
  alpha.1_backfit=alpha.1_backfit,
  l=l,
  x.grid=x.grid,
  data=fit.data,
  O1=O1,
  D=D,
  fit_diagnostics=fit.diagnostics,
  additive_predict_state=additive.predict.state,
  converged=converged,
  iterations_used=l,
  iterations_max=it,
  final_delta=final_delta,
  convergence_tol=convergence_tol
))
}
