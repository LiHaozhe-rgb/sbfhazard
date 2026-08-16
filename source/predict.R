####################### Prediction #######################

# predict alpha for a new data set by the same formula in the paper
predict.DM <- function(result, data_o, data_n){
  #### INPUTS: 
  # result:  list, result of SBF_MH_LC
  # data_o:  matrix, original data set 
  # data_n:  matrix, new data set
  
  #### OUTPUT:
  # estimated alpha_j 
  
  
  data_n <- as.matrix(data_n)
  if (ncol(data_n)==1) data_n <- t(data_n)                                    # check the format of new data 
  
  # browser()
  
  alpha <- result$SBF
  x.grid <- result$x.grid
  dx <- result$dx
  K.X.b <- result$K.X.b
  k.X.b <- result$k.X.b
  K.b <- result$K.b
  k.b <- result$k.b
  Y <- result$Y
  kern <- result$kern
  b.grid <- result$b.grid
  
  m <- nrow(data_n)                                                           # the number of individuals of new data set                   
  n <- length(x.grid[[1]])                                                    # the number of individuals of the old data set
  d <- length(alpha)                                                          # the dimensions of data set, including Time
  # 
  # # order each covariate of the new data set to make sure the kernel function is coincide with the data
  # data_n <- data.frame(lapply(1:d,function(k) data_n[order(data_n[,k]),k]))   # sort T and Z (each column is ordered increasingly)
  
  hazard <- matrix(0, nrow = m, ncol = (d+1))                                     # store alpha_j for new data and the product
  
  
  # the kernel function for all covariates
  x.grid.array1 <- matrix(rep(data_n[,1], times = n), nrow = m, ncol = n, byrow=FALSE)        # matrix, m*n, [T, ..., T]
  x.grid.array2 <- matrix(rep(x.grid[[1]], times = m),nrow = n, ncol = m, byrow=FALSE)        # matrix, n*m, [T, ..., T]
  u <-  x.grid.array2 - t(x.grid.array1)                                      # T_j-T_k, matrix, n*m
  K.b2 <- apply(u/b.grid[1],1:2, kern)/(b.grid[1])                            # the kernel result for T, matrix, n*m
  
  x.grid.array3 <- matrix(rep(x.grid[[1]], times = n),nrow = n, ncol = n, byrow=FALSE) 
  u_0 <- x.grid.array3 - t(x.grid.array3)
  k.b2 <- colSums(dx[[1]]*apply(u_0/b.grid[1],1:2,kern)/(b.grid[1]))          # the kernel boundary correction, a vector times a matrix return a matrix whose column is the the previous column times with vector, 1*n
  
  K.X.b2 <- k.X.b2 <- list()
  for(k in 1:d){
    K.X.b2[[k]]<-array(0,dim=c(n,m))
    k.X.b2[[k]]<-numeric(n)
    for (i in 1:n){
      # the kernel function at a certain point
      u <-  data_o[i,k] - data_n[,k]                                          # X_ik(s) - x, vector, (m), a scalar minus each ele. of a vector
      K.X.b2[[k]][i,] <- sapply(u/b.grid[k],kern)/(b.grid[k])                 # K(X_ik(s) - x), vector, (m)
      
      u_k <- data_o[i,k] - x.grid[[k]]                                        # X_ik(s) - u, vector, (n+n.add), a scalar minus each ele. of a vector
      k.X.b2[[k]][i] <- sum(dx[[k]]*sapply(u_k/b.grid[k],kern)/(b.grid[k]))   # the integration of k(X_ik(s) - u)
    }
  }
  
  
  # alpha estimation of the new data set
  alpha.smooth.i <- array(dim=c(d,n))
  alpha.smooth.i[1,] <- rep(1,n)                                              # alpha_star in eq.(4), vector, (n), always set it to be 1
  alpha.smooth.i.0 <- (K.b/k.b)%*%(alpha[[1]]*dx[[1]]) 
  
  for(k in 2:d){
    for (i in 1:n){
      alpha.smooth.i[k,i] <- as.numeric((K.X.b[[k]][i,]/k.X.b[[k]][i])%*%(alpha[[k]]*dx[[k]]))
    }
  }
  
  # browser()
  for(k in 1:d){
    if (k==1) alpha.minusk.smooth <- numeric(n) else  alpha.minusk.smooth <- array(dim=c(n,n))
    
    for (i in 1:n){
      if (k==1) alpha.minusk.smooth[i] <- prod(alpha.smooth.i[-1,i]) else       # scalar, exclude first column since it is always 1
      {alpha.minusk.smooth[i,] <- prod(alpha.smooth.i[-k,i])*alpha.smooth.i.0}  # vector, (n+n.add), the last term is adjusted alpha_0 
    }
    
    if (k==1) {
      D <- rowSums(sapply(1:n, function (i) { return((dx[[1]]*(alpha.minusk.smooth[i]*Y[i,]))%*%(K.b2/k.b2))})) # vector, (m)
    }else D <- rowSums(sapply(1:n, function (i) { return(as.numeric(dx[[1]]%*%(Y[i,]*(alpha.minusk.smooth[i,])))*(K.X.b2[[k]][i,]/k.X.b2[[k]][i]))})) # vector, (m)
    O <- rowSums(sapply(1:n, function (i) { return(data_o$status[i]*(K.X.b2[[k]][i,]/k.X.b2[[k]][i]))}))
    
    hazard[,k] <- O/D 
  }
  
  for (i in 1:m){
    hazard[i,d+1] <- prod(hazard[i,-(d+1)])
  }
  
  return(hazard)
}