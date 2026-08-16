# Additive SBF engine 
.sbf_fit_additive_engine <- function(data,bandwidth,
                     x.grid=NULL, n.grid.additional=0L,
                     x.min=NULL, x.max=NULL,
                     integral.approx='right',
                     it=100L,
                     kernel="epanechnikov",
                     initial=NULL,
                     convergence_tol=0.001,
                     kcorr=TRUE,
                     LC,
                     wrong,
                     classic.backfit=FALSE, print=FALSE){  

kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
kernel_fun <- kernel_spec$kernel
kernel_name <- kernel_spec$name

"div" <- function(x, y) ifelse(y == 0 & x == 0, 0, x / y)

mm <- data

time <- as.numeric(mm$time)
status <- as.numeric(mm$status)
X <- as.matrix(mm[, seq_len(ncol(mm) - 1L), drop = FALSE])


## Grid construction

d <- ncol(X) 
n <- nrow(X)

grid_spec <- .sbf_prepare_grid_spec(
  observed_grid = lapply(seq_len(d), function(k) X[order(X[, k]), k]),
  x_grid = x.grid,
  n_grid_additional = n.grid.additional,
  x_min = x.min,
  x_max = x.max,
  integral_approx = integral.approx
)
x.grid <- grid_spec$x_grid
n.grid <- grid_spec$n_grid
dx <- grid_spec$dx

dx.sum <- vapply(dx, sum, numeric(1))
dx.1 <- dx[[1]]
dx.1.matrix <- matrix(dx.1, ncol = 1L)



## Kernel weighting and exposure process

K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
k.b<-array(0,dim=c(n.grid[1]))



x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 


u<- x.grid.array-t(x.grid.array) # u[t,s]= t-s for all t,s on the grid considered
K.b<-.sbf_eval_kernel(u, bandwidth[1], kernel_fun)
k.b<-colSums(dx[[1]]*K.b)                      # small k for normailzing kernel, sincs grid points are symmetric normalization can be used for both row and column
if (kcorr==FALSE) {k.b<-rep(1,n.grid[1])}                     


dX0.b<- u#/bandwidth[1]     # dX0.b[t,s]= t-s for all t,s on the grid considered
K.b <- sweep(K.b, 2, k.b, "/")  ### row-wise division --> colSums(dx[[1]]*K.b)=1
K.b[is.na(K.b)] <-0
K.b.dx <- sweep(K.b, 1, dx.1, "*")
K.b.t.dx <- sweep(t(K.b), 1, dx.1, "*")
dX0.K.b.dx <- dX0.b * K.b.dx
dX0.t.K.b.dx <- t(dX0.b) * K.b.dx
dX0.K.b.t.dx <- sweep(t(dX0.b * K.b), 1, dx.1, "*")
dX0.sq.K.b.t.dx <- sweep(t((dX0.b^2) * K.b), 1, dx.1, "*")
### define exposure=Y[i,s]

Y <- (outer(time, x.grid[[1]], FUN = ">=") * 1) #Ymatrix
Y.colsum <- colSums(Y) #Yi(s)
Y.dx.1 <- as.numeric(Y %*% dx.1.matrix)



dX.b<-K.X.b<-k.X.b<-list()
K.X.dx <- vector("list", d)
dX.K.X <- vector("list", d)
dX.K.X.dx <- vector("list", d)
dX.sq.K.X <- vector("list", d)
for( k in 1:d){
  k.X.b[[k]]<-numeric(n)
  
  x.grid.array<-matrix(-X[,k],nrow =n.grid[k], ncol =n ,byrow=TRUE)   # 
  u<- x.grid.array+x.grid[[k]]      ####  u[,i]=x_k - X_{ik}
  K.X.b[[k]] <- .sbf_eval_kernel(u, bandwidth[k], kernel_fun)
  k.X.b[[k]]<-colSums(dx[[k]]*K.X.b[[k]])   
  if (kcorr==FALSE) k.X.b[[k]] <- rep(1,n)
  dX.b[[k]]<- u               #### dX.b[[k]] [,i] = x_k-X_{ik} 
  # K.X.bC[[k]] <- K.X.b[[k]]/ k.X.b[[k]]
  K.X.b[[k]]<- sweep(K.X.b[[k]], 2, k.X.b[[k]], "/")  ### row-wise division ---> col integration=1
  #(dx[[k]]*K.X.b[[k]]) =1
  K.X.b[[k]][is.na(K.X.b[[k]])] <-0
  dX.K.X[[k]] <- dX.b[[k]] * K.X.b[[k]]
  dX.sq.K.X[[k]] <- dX.b[[k]] * dX.K.X[[k]]
  if (k > 1L) {
    K.X.dx[[k]] <- sweep(K.X.b[[k]], 1, dx[[k]], "*")
    dX.K.X.dx[[k]] <- dX.b[[k]] * K.X.dx[[k]]
  }
  
}


## Kernel-smoothed occurrences and other variables

get.D<-function(k){
  # V_{0,0}^j, the denominator or eq(31)
  if (k==1)
  {
    D <- as.numeric(crossprod(Y.colsum, K.b.t.dx))
    
  }else 
    D <- as.numeric(K.X.b[[k]] %*% Y.dx.1)
  
  return(D)
}

get.O1<-function(k){
  # U_0^j in the numerator of eq(31), the first term
  O1 <- as.numeric(K.X.b[[k]] %*% status)
  
  return(O1)
}


get.O3<-function(k){
  # V_{j,0}^j in the numerator of eq(31), the second term
  if (k==1)
  {
    O3 <- as.numeric(crossprod(Y.colsum, dX0.K.b.t.dx))
    
  }else O3 <- as.numeric(dX.K.X[[k]] %*% Y.dx.1)
  
  return(O3) 
}

get.D.1<-function(k){
  
  if (k==1)
  {
    D <- as.numeric(crossprod(Y.colsum, dX0.sq.K.b.t.dx))
    
  }else
    D <- as.numeric(dX.sq.K.X[[k]] %*% Y.dx.1)
  return(D)
}


get.O1.1<-function(k){
  O1<- as.numeric(dX.K.X[[k]] %*% status)
  return(O1)
}

taylor.alpha<-function(alpha,alpha.1)
  # the last two terms without taking the sum in the numerator of eq(31)
{
  
  taylor.alpha.i <- array(dim=c(n,d))
  
  taylor.alpha.i[ ,1] <- rep(0,n)
  
  # the wrong variable indicates whether the difference is the same order as alpha DL
  if(wrong==TRUE) {
    taylor.alpha.i.0 <- as.numeric(
      crossprod(alpha[[1]], K.b.dx) + crossprod(alpha.1[[1]], dX0.t.K.b.dx)
    )
  }
  
  if(wrong==FALSE) {
    taylor.alpha.i.0 <- as.numeric(
      crossprod(alpha[[1]], K.b.dx) + crossprod(alpha.1[[1]], dX0.K.b.dx)
    )
  }
  
 
  for (k in 2:d){
    
    taylor.alpha.i[,k] <- as.numeric(
      crossprod(alpha[[k]], K.X.dx[[k]]) + crossprod(alpha.1[[k]], dX.K.X.dx[[k]])
    )
    
  }
  
  
  return(list(
    taylor.alpha.i= taylor.alpha.i,
    taylor.alpha.i.0=taylor.alpha.i.0,
    taylor.alpha.i.0.matrix = matrix(taylor.alpha.i.0, ncol = n.grid[1], nrow = n, byrow = TRUE),
    row.total = rowSums(taylor.alpha.i)
  ))
}

get.alpha.new<-function(O1,O3,D,alpha,alpha.1,Y,k)
  # taking the sum in the numerator of eq(31)
{

  if(classic.backfit==FALSE){
  
    taylor.alpha <- taylor.alpha(alpha,alpha.1)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    
    if (k==1) {
      taylor.alpha.minusk <- taylor.alpha$row.total
    } else {
      taylor.alpha.minusk <- taylor.alpha$row.total - taylor.alpha$taylor.alpha.i[, k] + taylor.alpha$taylor.alpha.i.0.matrix
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
      O2 <- as.numeric(crossprod(colSums(taylor.alpha.minusk * Y), K.b.t.dx))

    }else 
      O2 <- as.numeric(K.X.b[[k]] %*% ((taylor.alpha.minusk * Y) %*% dx.1.matrix))
    
    
    

    O <- O1[[k]] - alpha.1[[k]]*O3[[k]]-O2


  
  
  return(as.numeric(div(O,D[[k]])))
} 

get.alpha.1.new<-function(O1.1,O3.1,D.1,alpha,alpha.1,Y,k)
{    
  if(classic.backfit==FALSE){


    taylor.alpha <- taylor.alpha(alpha,alpha.1)
    
    if (k==1)  taylor.alpha.minusk<-numeric(n) else   taylor.alpha.minusk<-array(dim=c(n,n.grid[1])) 
    
    
    
    
    if (k==1) {
      taylor.alpha.minusk <- taylor.alpha$row.total
    } else {
      taylor.alpha.minusk <- taylor.alpha$row.total - taylor.alpha$taylor.alpha.i[, k] + taylor.alpha$taylor.alpha.i.0.matrix
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
      O2.1 <- as.numeric(crossprod(colSums(taylor.alpha.minusk * Y), dX0.K.b.t.dx))
    }else 
      O2.1 <- as.numeric(dX.K.X[[k]] %*% ((taylor.alpha.minusk * Y) %*% dx.1.matrix))
    


    
    O.1 <- O1.1[[k]]-alpha[[k]]*O3.1[[k]]-O2.1
    
  
  
    return(as.numeric(div(O.1,D.1[[k]])))
} 


D<-O1<-O3<-D1.1<-O1.1<-O3.1<-list()

D <- lapply(seq_len(d), get.D)
O1 <- lapply(seq_len(d), get.O1)
O3 <- lapply(seq_len(d), get.O3)

D.1 <- lapply(seq_len(d), get.D.1)
O1.1 <- lapply(seq_len(d), get.O1.1)
# The LL value and slope equations share this first-order exposure cross moment.
O3.1 <- O3
  
alpha_backfit<-list()
alpha.1_backfit<-list()

alpha_backfit <- if (is.null(initial)) {
  lapply(n.grid, function(n) rep(0, n))
} else {
  initial
}

for(k in 1:d){
  alpha.1_backfit[[k]]<-rep(0, n.grid[k])
  
}
  
iterations_used <- 0L
final_delta <- NA_real_
converged <- FALSE
for (iter in seq_len(it))
{  
  
  
  
  alpha_backfit_old<-alpha_backfit
  alpha.1_backfit_old<-alpha.1_backfit
  for(k in 1:d)
  {
  
  
    alpha_backfit[[k]]<- get.alpha.new(O1,O3,D,alpha_backfit,alpha.1_backfit,Y,k)
    
    if (LC==TRUE) alpha.1_backfit[[k]]<-rep(0,n.grid[k]) else{
      alpha.1_backfit[[k]]<- get.alpha.1.new(O1.1,O3.1,D.1,alpha_backfit,alpha.1_backfit,Y,k)}
    alpha.1_backfit[[k]][is.nan(alpha.1_backfit[[k]])]<-0
    alpha.1_backfit[[k]][is.na(alpha.1_backfit[[k]])]<-0
   
   if(k!=1){
      center.shift <- as.numeric(crossprod(alpha_backfit[[k]], dx[[k]]) / dx.sum[k])
      alpha_backfit[[1]]<-alpha_backfit[[1]] + center.shift
      alpha_backfit[[k]]<-alpha_backfit[[k]] - center.shift
   }
  }
  
  
    
  final_delta <- .sbf_fit_delta(
    c(alpha_backfit_old, alpha.1_backfit_old),
    c(alpha_backfit, alpha.1_backfit)
  )
  iterations_used <- iter
  if (is.finite(final_delta) && final_delta <= convergence_tol) {
    converged <- TRUE
    break
  }
  for(k in 1:d){
   if (print==TRUE) print(c(iter,max(abs(unlist((alpha_backfit_old[[k]]))-unlist((alpha_backfit[[k]]))),na.rm=TRUE),max(abs(unlist((alpha.1_backfit_old[[k]]))-unlist((alpha.1_backfit[[k]]
    ))),na.rm=TRUE)))
  }
}



for(k in 2:d){
  center.shift <- as.numeric(crossprod(alpha_backfit[[k]], dx[[k]]) / dx.sum[k])
  alpha_backfit[[1]]<-alpha_backfit[[1]] + center.shift
  alpha_backfit[[k]]<-alpha_backfit[[k]] - center.shift
}

taylor_terms <- taylor.alpha(alpha_backfit, alpha.1_backfit)
row_total <- taylor_terms$row.total
nuisance.baseline <- as.numeric(crossprod(row_total, Y))
nuisance.covariate <- vector("list", d)
if (d >= 2L) {
  weighted_Y <- sweep(Y, 2, dx.1, "*")
  time_nuisance <- as.numeric(weighted_Y %*% taylor_terms$taylor.alpha.i.0)
  for (k in 2:d) {
    nuisance.covariate[[k]] <- (row_total - taylor_terms$taylor.alpha.i[, k]) * Y.dx.1 + time_nuisance
  }
}

fit_diagnostics <- .sbf_kernel_norm_diagnostics(
  kernel_norm = list(k.b, k.X.b),
  floor_value = 1e-12
)

additive_predict_state <- list(
  bandwidth = as.numeric(bandwidth),
  kernel = kernel_fun,
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
  x.grid=x.grid,
  data=mm,
  fit_diagnostics=fit_diagnostics,
  additive_predict_state=additive_predict_state,
  converged=converged,
  iterations_used=iterations_used,
  iterations_max=it,
  final_delta=final_delta,
  convergence_tol=convergence_tol
))
}


# Internal additive SBF fit adapter used by the public API wrapper.
.sbf_add_fit <- function(data,
                               bandwidth,
                               formula = NULL,
                               x_grid = NULL,
                               n_grid_additional = 0L,
                               x_min = NULL,
                               x_max = NULL,
                               integral_approx = "midd",
                               iterations = 200L,
                               kernel = "epanechnikov",
                               initial = NULL,
                               convergence_tol = 0.001,
                               kernel_correction = TRUE,
                               local_constant = FALSE,
                               wrong_taylor = FALSE,
                               classic_backfit = FALSE,
                               print_trace = FALSE,
                               warn_nonconvergence = TRUE) {
  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  feature_names <- setdiff(names(data), c("time", "status"))

  fit <- .sbf_fit_additive_engine(
    data = data,
    bandwidth = bandwidth,
    x.grid = x_grid,
    n.grid.additional = n_grid_additional,
    x.min = x_min,
    x.max = x_max,
    integral.approx = integral_approx,
    it = iterations,
    kernel = kernel,
    initial = initial,
    convergence_tol = convergence_tol,
    kcorr = kernel_correction,
    LC = local_constant,
    wrong = wrong_taylor,
    classic.backfit = classic_backfit,
    print = print_trace
  )

  fit$converged <- isTRUE(fit$converged)
  fit$model_type <- "additive"
  fit["formula"] <- list(formula)
  fit$feature_names <- feature_names
  fit$component_info <- .sbf_component_info_from_grid(
    x_grid = fit$x.grid,
    component_names = c("time", feature_names),
    roles = c("time", rep("covariate", length(feature_names))),
    support_min = x_min,
    support_max = x_max
  )
  fit$fit_settings <- list(
    bandwidth = bandwidth,
    integral_approx = integral_approx,
    iterations = iterations,
    convergence_tol = convergence_tol,
    kernel = kernel_spec$kernel,
    kernel_name = kernel_spec$name,
    kernel_correction = kernel_correction,
    local_constant = local_constant,
    wrong_taylor = wrong_taylor,
    classic_backfit = classic_backfit,
    print_trace = print_trace
  )
  class(fit) <- unique(c("sbf_additive_fit", class(fit)))
  if (warn_nonconvergence && !isTRUE(fit$converged)) {
    .sbf_emit_nonconvergence_warning(
      model_type = fit$model_type,
      iterations_used = fit$iterations_used,
      iterations_max = fit$iterations_max,
      final_delta = fit$final_delta,
      convergence_tol = fit$convergence_tol
    )
  }
  fit
}
