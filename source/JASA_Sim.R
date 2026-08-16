################## RUN BIG SIMULATION##########################
#library(Rmpi)
#library(snow)

########################## comments by DM ######################################
# This file should produce the same plots as the paper
# The main questions are: 1) why no weight when compute alpha.smooth? 2) the final results were normalized follows what theory?
# answer from Stefan: 1) the weights are constant and the same, so they get canceled in the computation (check more from papers)
graphics.off()
library(parallel)                                                             # comments from DM

d.par<-c(3)                                                                   # the dimension of data, including time
rho.par<-c(0.8)                                                               # the correlation between Z tilde
n.par<-c(200)                                                                 # the number of observations of data
model.par<-c(2)                                                               # the model
n.sim.par<-c(200)                                                             # the number of simulations/ repetitions
n.cores<-c(7)                                                                # number of cores can be used in parallel
r<-1:n.sim.par                                                                # a list associated with number of simulations

#d.par<-c(100)
#rho.par<-c(0)
#model.par<-c(1)
#n.par<-c(1000)
#n.sim.par<-c(1)
#n.cores<-c(2)
#r<-1:n.sim.par
set.seed(23)
jobs<-expand.grid(d.par,rho.par,model.par,n.par,n.sim.par,n.cores,r)          # create a data frame of parameters for running them in parallel
jobs<-lapply(1:nrow(jobs), function(x){list(jobs=jobs[x,],seeds=round(runif(n.sim.par,0,999999999)))})  # assign different seeds to jobs

myparallelfunction<-function(x){
  library("MASS")                                                             # to use the mvnorm function

  smooth.alpha<-function(alpha,K.X.b,k.X.b,K.b,k.b,b,n.grid,d,n,nb,get.baseline)
  {
    alpha.smooth.i<-array(dim=c(d,n))
    alpha.smooth.i.0<-array(dim=c(n,n.grid[1]))
    for (k in 1:d)
    { if (k==1){ for (i in 1:n)     ## (k==1) is outside loop in "SBF_multiplicative_hazard_StSp.R"
      { alpha.smooth.i[k,i] <- 1                                              # equivalent to what he did in SBF_MH_LC.R, for alpha_star 
        alpha.smooth.i.0[i,]<-(K.b/k.b)%*%alpha[[k]]  }                       # no weight dx because they are the same and constant, can be cancelled (TODO)
      } else { for (i in 1:n)
        { alpha.smooth.i[k,i] <- (K.X.b[[k]][i,]/k.X.b[[k]][i])%*%alpha[[k]] } 
      }
    }
    # alpha.smooth.i[2:d,]<-matrix(unlist(mclapply(2:d, function(k){(K.X.b[[k]][b,,]/k.X.b[[k]][b,])%*%(alpha[[k]])})),nrow=(d-1),byrow=TRUE)
    if (get.baseline==FALSE) return(alpha.smooth.i) else return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))                # estimate alpha_0 or not
  }

  get.denominator<-function(data,alpha,K.X.b,k.X.b,K.b,k.b,Y,kk,bb,x.grid,n.grid,d,n,nb,get.baseline)
    ### comments by DM: maybe could be simplified, too much nested ifs????
    ### the output is actually the alpha, not just the denominator
  {
    alpha.smooth.i<-smooth.alpha(alpha,K.X.b,k.X.b,K.b,k.b,bb,n.grid,d,n,nb,get.baseline)    
    if (get.baseline==FALSE)    ## in "SBF_multiplicative_hazard_StSp.R" FALSE, but -1 instead of -kk 
    { alpha.minusk.smooth<-numeric(n)
      for (i in 1:n) { alpha.minusk.smooth[i] <-  prod(alpha.smooth.i[-kk,i]) }                                                                      # alpha.smooth.i[k,i] is 1 when kk=1
    }
    if (get.baseline==TRUE){  ## this is for estimating alpha_0, function of t neq ==1
      alpha.minusk.smooth<-array(dim=c(n,n.grid[1])) 
      for (i in 1:n){				
         if (kk==1){ alpha.minusk.smooth[i,] <- rep(prod(alpha.smooth.i$alpha.smooth.i[-1,i]), n.grid[1])                                            # for update alpha_0, repeat the prod, why bother to do it?
         } else { alpha.minusk.smooth[i,] <- prod(alpha.smooth.i$alpha.smooth.i[-kk,i])*alpha.smooth.i$alpha.smooth.i.0[i,] }                        # equivalent to SBF_MH_LC
    } } 
    if (get.baseline==TRUE){  ## if alpha_0 not constant
      if (kk==1){
        D<- rowSums(sapply(1:n, function(i){ return((
          c((x.grid[[1]][-1]-x.grid[[1]][-n.grid[1]]), x.grid[[1]][n.grid[1]]-x.grid[[1]][(n.grid[1]-1)])
          *(alpha.minusk.smooth[i,]*Y[i,]))
          %*%(K.b/k.b)
          )
          }))  # equivalent, the first part is computing the integration weight dx
        }else D<- rowSums(sapply(1:n, function(i){ return(
          as.numeric(
            c((x.grid[[1]][-1]-x.grid[[1]][-n.grid[1]]), x.grid[[1]][n.grid[1]]-x.grid[[1]][(n.grid[1]-1)])
            %*%(Y[i,]
            *(alpha.minusk.smooth[i,])))
            *(K.X.b[[kk]][i,]/k.X.b[[kk]][i])
          )}))                                              # equivalent to SBF_MH_LC
     }else {
      if (kk==1){
        D<- rowSums(sapply(1:n, function(i){ return((c((x.grid[[1]][-1]-x.grid[[1]][-n.grid[1]]), x.grid[[1]][n.grid[1]]-x.grid[[1]][(n.grid[1]-1)])*(alpha.minusk.smooth[i]*Y[i,]))%*%(K.b/k.b))}))   # the same as when get.baseline==TRUE and kk=1
        }else D<- rowSums(sapply(1:n, function(i){ return(as.numeric(c((x.grid[[1]][-1]-x.grid[[1]][-n.grid[1]]), x.grid[[1]][n.grid[1]]-x.grid[[1]][(n.grid[1]-1)])%*%(Y[i,]*(alpha.minusk.smooth[i])))*(K.X.b[[kk]][i,]/k.X.b[[kk]][i]))}))
    }
    # O <- rowSums(sapply(1:n, function(i){ return((alpha.minusk.smooth[kk,bb,i]*data$status[i])*(K.X.b[[kk]][bb,i,]/k.X.b[[kk]][bb,i]))}))
    O <- rowSums(sapply(1:n, function(i){ return((data$status[i])*(K.X.b[[kk]][i,]/k.X.b[[kk]][i]))}))                                               # equivalent to SBF_MH_LC 
    return(O/D)
  } 


  sbf<-function(data,name,x.grid,b.grid,d,n,nb,n.grid,it=100,kern,phi,initial=matrix(1,nrow=d,ncol=n.grid),get.baseline)
  {
    #### compute the estimated alpha (DM)
    ##### Note: The adjusted kernel,  K.b[k,b,,]/k.b[k,b,], has rowSums equal one.
    
    K.b<-array(0,dim=c(n.grid[1],n.grid[1]))
    k.b<-array(0,dim=c(n.grid[1]))
    
    x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 
    u<-x.grid.array-t(x.grid.array) # u is x_j-u_j for all x,u on the grid considered
    K.b<-apply(u/b.grid[1,1],1:2,kern)/(b.grid[1,1])
    k.b<-rowSums(apply(u/b.grid[1,1],1:2,kern)/(b.grid[1,1]))
     
    ##### get covariances 
    X<-array(dim=(c(d,n)))
    ##### in this example
    for(k in 1:d){
      X[k,]<-as.numeric(data[name[k]][,1])
    }
        
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
        K.X.b[[k]][i,]<-sapply(u/b.grid[k,1],kern)/(b.grid[k,1])
        k.X.b[[k]][i]<-sum(sapply(u/b.grid[k,1],kern)/(b.grid[k,1]))
      } }
            
    alpha_backfit<-list()        
    for(k in 1:d){
      alpha_backfit[[k]]<-rep(1, n.grid[k])
    }
    
    alpha.temp<-list()
    for (l in 2:it)
    {  alpha_backfit_old<-alpha_backfit
    for(k in 1:d)
    {
      # if(k==1) alpha.temp<-lapply(1:d,function(z){alpha_backfit[[z]][b,]}) else 
      #   { alpha.temp[1:(k-1)]<-lapply(1:(k-1),function(z){alpha_backfit[[z]][b,]})
      #     #alpha.temp[k:d]<-lapply(k:d,function(z){alpha_backfit_old[[z]][b,]})
      #     }
            
      alpha.temp<-lapply(1:d,function(z){alpha_backfit[[z]]})      
      for(m in 1:d)
      {
        alpha.temp[[m]][is.nan(alpha.temp[[m]])]<-0
        alpha.temp[[m]][alpha.temp[[m]]==Inf]<-0
      }
      for (j in (1:(d-1))[-k])
      { if (j==1) {alpha.temp[[1]]<-alpha.temp[[1]]*(1/alpha.temp[[1]][1])} else                                       # for time, normalized by the estimation of first observation
        alpha.temp[[j]]<-alpha.temp[[j]]*(exp(phi[[j-1]](x.grid[[j]][n.grid[j]/2]))/alpha.temp[[j]][n.grid[j]/2])      # equivalent to Sim_est normalization 
      }
      if (get.baseline==FALSE)  alpha.temp[[1]]<-rep(1,n.grid[1])
      
      alpha_backfit[[k]]<- get.denominator(data,alpha.temp,K.X.b,k.X.b,K.b,k.b,Y,k,b,x.grid,n.grid,d,n,nb,get.baseline)
      #alpha_backfit[l,k,b,]<- O.k.b[k,b,]/get.denominator( alpha_backfit[l-1,,b,],K.X.b,k.X.b,Y,k,b)
    }
    
    if (max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE)<=0.001) break
    print(c(l,max(abs(unlist(alpha_backfit_old)-unlist(alpha_backfit)),na.rm=TRUE)))
    }
    return(list(alpha_backfit=alpha_backfit,l=l))
  }
  
  MY<-function(Z,Y,status,initial,h0,itera,rrr,phi,it)
    ################### Comments from DM #################
  # the MY... function is the implementation of Lin et al.'s paper: A global partial likelihood estimation
  {
    # Z is a matrix with the covariates values (consider p covariates, for any p)  
    # Y is a vector with the observed times (n=length(Y) is the sample size)
    # status is the censoring indicator
    # initial are initial values for the covariate effects (matrix same dimension as Z)
    # suggestion in the paper is to use Huang (1999)'s estimator but below they use the true values
    # h0 is the bandwidth parameter (scalar, the same bandwidth for all covariates)
    # itera is tolerance criterion (either it is achieved or the number of iterations is > 150)
    # rrr is a cut point of Z-values. The functions eta(z) are estimated at each z in the grid
    #   the number rrr just define how we vary the z's. 
    #   First we estimate at the rrr-th z-point, afterwards at the bigger z's and finally at the smaller ones.        
    n<-length(Y)
    p<-dim(Z)[[2]]
    time<-sort(unique(Y[status==1]))  ## sort the times so the risk set is easier to compute below
    d<-length(time)
    risk<-(outer(Y,rep(1,d))>=outer(rep(1,n),time))  ## the risk set for each time                            # matrix n*d
    # (number of individuals at risk at time Yi, i=1,..,n)
    STATUS<-status*(outer(Y,rep(1,d))==outer(rep(1,n),time))
    score<-rep(1,2)*0  ## for newton-raphson (first derivative)                                               # why not define it as 0 directly?
    imat<-outer(score,score)*0  ## for newton-raphson (hessian)                                               # outer(A,B) return an array with dim c(dim(A), dim(B))
    GBETA<-initial                                                                                            # matrix n*p
    GBETA1<-GBETA+10  ## initialize GBETA1 to be just different from GBETA, 10 or other value is fine
    alpha0<-rep(1,2)*0  ## alpha is the vector of parameters (eta(z),eta'(z)),                                # local linear estimator
    # we are interested in the first component but both are derived solving the equation (2.7)                
    h<-h0
    r<-1   ## number of iteration for the algorith below
    
    ## Now equation (2.7) is solved using Newton-Raphson (tolerance=itera and max iterations = 150)
    while((max(abs(GBETA1-GBETA),na.rm=TRUE)>itera)&(r<=it))
    {  
      #   if (is.na(sum(GBETA))==TRUE){
      #   break
      #}
      GBETA1<-GBETA
      for (k in 1:p)
      {
        z0<-sort(unique(round(Z[,k],2)))  ## the evaluation points are sorted
        alpha0[1]<-initial[round(Z[,k],2)==z0[rrr],k][1]                                                      # obtain the corresponding initial value eta(z)
        alpha1<-alpha0                                                                                        # not updating eta'(z) ????
        j<-rrr
        while (j<=length(z0))
        {
          wei<-3/4*(1-(Z[,k]-z0[j])^2/(h^2))*(abs(Z[,k]-z0[j])<=h)/h  # Epanechnikov evaluations
          Y0<-Y[wei>0]                                                                                        # filter the Y of neighborhood
          Z0<-Z[wei>0,k]                                                                                   
          wei0<-wei[wei>0]                                                                                 
          n0<-length(wei0)                                                                                   
          risk0<-risk[wei>0,]
          STATUS0<- STATUS[wei>0,]
          
          alpha1[1]<-initial[round(Z[,k],2)==z0[j],k][1]
          alpha<-alpha1+1                                                                                     # why plus one for both eta(z) and eta'(z)????
          fenmu<-c(rep(1,n)%*%(risk*exp(c(GBETA%*%rep(1,p)))))                                                # vector 1*d = (1*n)%*%((n*d)*exp(n*1))
          while(sum(abs(alpha-alpha1),na.rm=TRUE)>=0.00001)
          {
            alpha<-alpha1
            wexpxa<-wei0*exp(alpha[1]+alpha[2]*(Z0-z0[j]))                                                    # local linear estimator for alpha, not eta?
            for (jj in 1:p)
            {
              if (jj!=k)
                wexpxa<- wexpxa*exp(GBETA[wei>0,jj])                                                          # updating 
            }
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa))
            score[1]<-sum(c(wei0%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            score[2]<-sum(c(((Z0-z0[j])*wei0)%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa))
            imat[1,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            imat[1,2]<- imat[2,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])^2))
            imat[2,2]<- -sum(c(STATUS%*%c(fz/fenmu)))
            alpha1<-tryCatch(alpha-c(solve(imat)%*%score), error=function(bla) 
            {print(paste(bla,'error in l. 390'))
              return(alpha)}
            )
            if (is.na(sum(alpha1))) alpha1<-GBETA[round(Z[,k],2)==z0[j],k]
            #if (is.na(sum(alpha1))==TRUE){
            #    break
            #  }
            #print(c(sum(abs(alpha-alpha1)),j))
          }
          GBETA[round(Z[,k],2)==z0[j],k]<-alpha1[1]
          j<-j+1
        }
        alpha0[1]<-initial[round(Z[,k],2)==z0[rrr-1],k][1]
        alpha1<-alpha0
        j<-rrr-1
        while (j>=1)
        {
          #if (is.na(sum(alpha1))==TRUE){
          #  break
          # }
          wei<-3/4*(1-(Z[,k]-z0[j])^2/(h^2))*(abs(Z[,k]-z0[j])<=h)/h
          Y0<-Y[wei>0]
          Z0<-Z[wei>0,k]
          wei0<-wei[wei>0]
          n0<-length(wei0)
          risk0<-risk[wei>0,]
          STATUS0<- STATUS[wei>0,]
          
          alpha1[1]<-initial[round(Z[,k],2)==z0[j],k][1]
          alpha<-alpha1+1
          fenmu<-c(rep(1,n)%*%(risk*exp(c(GBETA%*%rep(1,p)))))
          while(sum(abs(alpha-alpha1),na.rm=TRUE)>=0.00001)
          {
            alpha<-alpha1
            wexpxa<-wei0*exp(alpha[1]+alpha[2]*(Z0-z0[j]))
            for (jj in 1:p)
            {
              if (jj!=k)
                wexpxa<- wexpxa*exp(GBETA[wei>0,jj])
            }
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa))
            score[1]<-sum(c(wei0%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            score[2]<-sum(c(((Z0-z0[j])*wei0)%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa))
            imat[1,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            imat[1,2]<- imat[2,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])^2))
            imat[2,2]<- -sum(c(STATUS%*%c(fz/fenmu)))
            alpha1<-tryCatch(alpha-c(solve(imat)%*%score), error=function(bla) 
            {print(paste(bla,'error in l. 439'))
              return(alpha)}
            )
            if (is.na(sum(alpha1))) alpha1<- GBETA[round(Z[,k],2)==z0[j],k]
            #   if (is.na(sum(alpha1))==TRUE){
            #     break
            #   }
            #print(c(sum(abs(alpha-alpha1)),j))
          }
          GBETA[round(Z[,k],2)==z0[j],k]<-alpha1[1]
          j<-j-1
        }
        #GBETA[,k]<- GBETA[,k]-GBETA[round(Z[,k],2)==z0[rrr],k][1]+initial[round(Z[,k],2)==z0[rrr],k][1]        
        GBETA[,k]<- GBETA[,k]-GBETA[,k][order(abs(Z[,k]-0))][1]+phi[[k]](Z[,k][order(abs(Z[,k]-0))][1])        # normalized by the estimation and true value of minimum
      }
      r<-r+1
      print(c(r,max(abs(GBETA1-GBETA))))
    }
    
    #  if(max(abs(GBETA1-GBETA))<=itera)
    return(list(l.lin=r,GBETA=GBETA))
    # if((max(abs(GBETA1-GBETA))>itera))
    #  return(GBETA*0)
  }

  MY.orig<-function(Z,Y,status,initial,h0,itera,rrr,it)
  {
    # Z is a matrix with the covariates values (consider p covariates, for any p)  
    # Y is a vector with the observed times (n=length(Y) is the sample size)
    # status is the censoring indicator
    # initial are initial values for the covariate effects (matrix same dimension as Z)
    # suggestion in the paper is to use Huang (1999)'s estimator but below they use the true values
    # h0 is the bandwidth parameter (scalar, the same bandwidth for all covariates)
    # itera is tolerance criterion (either it is achieved or the number of iterations is > 150)
    # rrr is a cut point of Z-values. The functions eta(z) are estimated at each z in the grid
    #   the number rrr just define how we vary the z's. 
    #   First we estimate at the rrr-th z-point, afterwards at the bigger z's and finally at the smaller ones.
    n<-length(Y)
    p<-dim(Z)[[2]]
    time<-sort(unique(Y[status==1]))  ## sort the times so the risk set is easier to compute below
    d<-length(time)
    risk<-(outer(Y,rep(1,d))>=outer(rep(1,n),time))  ## the risk set for each time 
    # (number of individuals at risk at time Yi, i=1,..,n)
    STATUS<-status*(outer(Y,rep(1,d))==outer(rep(1,n),time))
    score<-rep(1,2)*0  ## for newton-raphson (first derivative)
    imat<-outer(score,score)*0  ## for newton-raphson (hessian)
    GBETA<-initial
    GBETA1<-GBETA+10  ## initialize GBETA1 to be just different from GBETA, 10 or other value is fine
    alpha0<-rep(1,2)*0  ## alpha is the vector of parameters (eta(z),eta'(z)), 
    # we are interested in the first component but both are derived solving the equation (2.7)
    h<-h0
    r<-1   ## number of iteration for the algorith below
    
    ## Now equation (2.7) is solved using Newton-Raphson (tolerance=itera and max iterations = 150)
    while((max(abs(GBETA1-GBETA))>itera)&(r<=it))
    {
      GBETA1<-GBETA
      for (k in 1:p)
      {
        z0<-sort(unique(round(Z[,k],2)))  ## the evaluation points are sorted
        alpha0[1]<-initial[round(Z[,k],2)==z0[rrr],k][1]
        alpha1<-alpha0
        j<-rrr
        while (j<=length(z0))
        {
          wei<-3/4*(1-(Z[,k]-z0[j])^2/(h^2))*(abs(Z[,k]-z0[j])<=h)/h  # Epanechnikov evaluations
          Y0<-Y[wei>0]
          Z0<-Z[wei>0,k]
          wei0<-wei[wei>0]
          n0<-length(wei0)
          risk0<-risk[wei>0,]
          STATUS0<- STATUS[wei>0,]
          
          alpha1[1]<-initial[round(Z[,k],2)==z0[j],k][1]
          alpha<-alpha1+1
          fenmu<-c(rep(1,n)%*%(risk*exp(c(GBETA%*%rep(1,p)))))
          while(sum(abs(alpha-alpha1))>=0.00001)
          {
            alpha<-alpha1
            wexpxa<-wei0*exp(alpha[1]+alpha[2]*(Z0-z0[j]))
            for (jj in 1:p)
            {
              if (jj!=k)
                wexpxa<- wexpxa*exp(GBETA[wei>0,jj])
            }
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa))
            score[1]<-sum(c(wei0%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            score[2]<-sum(c(((Z0-z0[j])*wei0)%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa))
            imat[1,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            imat[1,2]<- imat[2,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])^2))
            imat[2,2]<- -sum(c(STATUS%*%c(fz/fenmu)))
            alpha1<-tryCatch(alpha-c(solve(imat)%*%score), error=function(bla) 
            {print(paste(bla,'error in l. 390'))
              return(NA)}
            )
            if (is.na(sum(alpha1))==TRUE) break
            #print(c(sum(abs(alpha-alpha1)),j))
          }
          if (is.na(sum(alpha1))==TRUE) break
          GBETA[round(Z[,k],2)==z0[j],k]<-alpha1[1]
          j<-j+1
        }
        if (is.na(sum(alpha1))==TRUE) break
        alpha0[1]<-initial[round(Z[,k],2)==z0[rrr-1],k][1]
        alpha1<-alpha0
        j<-rrr-1
        while (j>=1)
        {
          wei<-3/4*(1-(Z[,k]-z0[j])^2/(h^2))*(abs(Z[,k]-z0[j])<=h)/h
          Y0<-Y[wei>0]
          Z0<-Z[wei>0,k]
          wei0<-wei[wei>0]
          n0<-length(wei0)
          risk0<-risk[wei>0,]
          STATUS0<- STATUS[wei>0,]
          
          alpha1[1]<-initial[round(Z[,k],2)==z0[j],k][1]
          alpha<-alpha1+1
          fenmu<-c(rep(1,n)%*%(risk*exp(c(GBETA%*%rep(1,p)))))
          while(sum(abs(alpha-alpha1))>=0.00001)
          {
            alpha<-alpha1
            wexpxa<-wei0*exp(alpha[1]+alpha[2]*(Z0-z0[j]))
            for (jj in 1:p)
            {
              if (jj!=k)
                wexpxa<- wexpxa*exp(GBETA[wei>0,jj])
            }
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa))
            score[1]<-sum(c(wei0%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fenzi<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            score[2]<-sum(c(((Z0-z0[j])*wei0)%*%STATUS0))-sum(c(STATUS%*%c(fenzi/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa))
            imat[1,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])))
            imat[1,2]<- imat[2,1]<- -sum(c(STATUS%*%c(fz/fenmu)))
            fz<-c(rep(1,n0)%*%(risk0*wexpxa*(Z0-z0[j])^2))
            imat[2,2]<- -sum(c(STATUS%*%c(fz/fenmu)))
            alpha1<-tryCatch(alpha-c(solve(imat)%*%score), error=function(bla) 
            {print(paste(bla,'error in l. 444'))
              return(NA)}
            )
            if (is.na(sum(alpha1))==TRUE) break
            
            #print(c(sum(abs(alpha-alpha1)),j))
          }
          if (is.na(sum(alpha1))==TRUE) break
          GBETA[round(Z[,k],2)==z0[j],k]<-alpha1[1]
          j<-j-1
        }
        if (is.na(sum(alpha1))==TRUE) break
        GBETA[,k]<- GBETA[,k]-GBETA[round(Z[,k],2)==z0[rrr],k][1]+initial[round(Z[,k],2)==z0[rrr],k][1]
      }
      if (is.na(sum(alpha1))==TRUE) break
      r<-r+1
      print(c(r,max(abs(GBETA1-GBETA))))
    }
    
  #  if(max(abs(GBETA1-GBETA))<=itera)
  #  if((max(abs(GBETA1-GBETA))>itera))
  #    return(list(l.lin=r,GBETA=GBETA*0))
    if (is.na(sum(alpha1))==TRUE)  return(list(l.lin=r,GBETA=GBETA,broke=TRUE))
    return(list(l.lin=r,GBETA=GBETA,broke=FALSE))    
  }
  
  do.sim<-function(n.cores=2,n.sim=1,n,r,nb, n.grid, b.grid,x.grid,  model=1,rho=0,d, distr, initial=matrix(0, ncol=(d-1),nrow=n),it,kern,seed.n,get.baseline,d1,d2)
  { 
    ptm <- proc.time()
    # LIN<-SBF<-SBF2<-list()    
    # n.grid<-list
    # x.grid<-list()    
    set.seed(seed.n[r])
                  
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
      # 
      xx.grid<-list()
      nn.grid<-numeric(d)
      for (k in 1:d)
      {
        xx.grid[[k]]<-data[order(data[,k]),k]
        # xx.grid[[k]]<-data[,k]
        #xx.grid2<-seq(k.min[k],k.max(k),length=100)
        #xx.grid[[k]]<-c(xx.grid2,xx.grid[[k]])
        #xx.grid[[k]]<- unique(xx.grid[[k]][(order(xx.grid[[k]]))])
        
        nn.grid[k]<-length(xx.grid[[k]])        
      }
      xx.grid2<-seq(k.min[1],max(TT),length=n.grid[1])
      # nn.grid[1]<-n.grid[1]
      xx.grid[[1]]<-c(xx.grid2,xx.grid[[1]])
      xx.grid[[1]]<- xx.grid[[1]][(order(xx.grid[[1]]))]
      nn.grid[1]<-length(xx.grid[[1]])
            
      #  alpha_backfit<-sbf(data,name,xx.grid,b.grid,d,n,nb,nn.grid,it,kern,phi,get.baseline=FALSE)       
      #    alpha_backfit.exp<-sbf.exp(data,name,x.grid,b.grid,d,n,nb,n.grid,it,kern,phi)
      if(get.baseline==TRUE)
        ptm <- proc.time()
      alpha_backfit.with.baseline<-sbf(data,name,xx.grid,b.grid,d,n,nb,nn.grid,it,kern,phi,get.baseline=TRUE)
      comp.time.sbf<-(proc.time() - ptm)[[3]]
      l.sbf<-alpha_backfit.with.baseline$l
      alpha_backfit.with.baseline<-alpha_backfit.with.baseline$alpha_backfit
      
      #l.exp<-alpha_backfit.exp$l
      # alpha_backfit.exp<-alpha_backfit.exp$alpha_backfit
                  
      #  if (d>2){
      #    for (k in 1:(d-1))
      #     {  
      #      corr <-log(alpha_backfit[[k+1]][2,nn.grid[k+1]/2])-phi[[k]](xx.grid[[k+1]][nn.grid[k+1]/2])
      #      alpha_backfit[[k+1]][2,] <- alpha_backfit[[k+1]][2,]/ exp(corr)      
      #      #   corr <-alpha_backfit.exp[[k+1]][2,n.grid[k+1]/2]-phi[[k]](x.grid[[k+1]][n.grid[k+1]/2])
      #      #   alpha_backfit.exp[[k+1]][2,] <- alpha_backfit.exp[[k+1]][2,]- corr
      #      #  alpha_backfit[l,k+2,2,]<-alpha_backfit[l,k+2,2,] / exp(corr)
      #     }
      #  }      
      # 
      
      if(get.baseline==TRUE)
      {
        if (d>2){
          for (k in 1:(d-1))
          {  
            corr <-log(alpha_backfit.with.baseline[[k+1]][nn.grid[k+1]/2])-phi[[k]](xx.grid[[k+1]][nn.grid[k+1]/2])
            alpha_backfit.with.baseline[[k+1]] <- alpha_backfit.with.baseline[[k+1]]/ exp(corr)                                  # the same normalization as SBF_MH_LC
            
            #   corr <-alpha_backfit.exp[[k+1]][2,n.grid[k+1]/2]-phi[[k]](x.grid[[k+1]][n.grid[k+1]/2])
            #   alpha_backfit.exp[[k+1]][2,] <- alpha_backfit.exp[[k+1]][2,]- corr
            #  alpha_backfit[l,k+2,2,]<-alpha_backfit[l,k+2,2,] / exp(corr)
          }
        }
      }
      
      ptm <- proc.time()
      #re<-MY(Z,TT,status,initial,b.grid[2,1],0.001,100,phi,it) ## tolerance is 0.001 and 100 is an arbitrary choice
      re<-MY.orig(Z,TT,status,initial,b.grid[2,1],0.001,100,it) ## tolerance is 0.001 and 100 is an arbitrary choice

      comp.time.lin<-(proc.time() - ptm)[[3]]
      ## eta will be estimated at each grid point z0. z0[100] would divide the points into two groups
      ## first it estimates at the right side of it and then at the left side. 
      l.lin<-re$l.lin
      broke<-re$broke
      re<-re$GBETA
      LIN1<-list()

      # same as (1:200)%*%t(Z0), dim(ETA1)=c(200,201)

      for (k in 1:(d-1))
      { LIN1[[k]]<-numeric(nn.grid[k+1])
      for(j in 1:nn.grid[k+1])
      {
        LIN1[[k]][j]<-re[,k][order(abs(Z[,k]-xx.grid[[k+1]][j]))][1]
      }
      }

      if (d>2){
        for (k in 1:(d-1))
        {
          LIN1[[k]]<-LIN1[[k]]-LIN1[[k]][nn.grid[k+1]/2]+phi[[k]](xx.grid[[k+1]][nn.grid[k+1]/2])
          #LIN1[[k+1]]<-LIN1[[k+1]]+LIN1[[k]][n.grid[k]/2]-phi[[k]](xx.grid[[k]][n.grid[k]/2])
        }
      }
      #win.graph()

      #  LIN1[r,]=ETA[r,1,]
      #  LIN2[r,]=ETA[r,2,]
      # SBF1[r,]=alpha_backfit[l,2,2,]
      # SBF2[r,]=alpha_backfit[l,3,2,]
      # dput(A1,"ETA1.r")
      # dput(A2,"ETA2.r")
  result<-list(LIN=LIN1,SBF=alpha_backfit.with.baseline,l.lin=l.lin,l.sbf=l.sbf,comp.time.lin=comp.time.lin,comp.time.sbf=comp.time.sbf,x.grid=xx.grid,n.grid=nn.grid, phi=phi,broke=broke)

    
    myfilename<-paste('n.sim=',n.sim,'Model',model, 'Dim=',d,'distribution=', distr,'n=',n, 'rho=',rho)
    
   # pdf(file=paste('plots',myfilename,'.pdf'),height = 4*d1, width = 15, onefile=TRUE)
    
      #   quartz()
    # par(mfrow=(c(d1,d2)))
      
    #  for (k in 1:(d-1))
    #  {
    #    plot(result$x.grid[[k+1]],result$LIN[[k]],col=4,type="l",lwd=2,xlab=expression('z'['k']),ylab="covariate effect",ylim=c(-2.3,4),main=paste('SIM',r,' ','Model',model,' ','dim=',d,' ','n=',n,' ', 'rho=',rho,' ','k=',k))
    #    lines(result$x.grid[[k+1]],result$phi[[k]](result$x.grid[[k+1]]),lty=3,col=1,lwd=2) 
     #  lines(result$x.grid[[k+1]], log(result$SBF[[k+1]][2,]),col='green')
     #    if(get.baseline==TRUE) lines(result$x.grid[[k+1]], log(result$SBF[[k+1]]),col='red',lwd=2)
     #   lines(x.grid[[k+1]], (alpha_backfit.exp[[k+1]][2,]),col='purple')
     #  legend("topleft",c("Lin He Huang (2016)","SBF","True function"),col=c(4,'red',1),lty=c(1,1,3),lwd=2,bt="n")
     #  }
     #  if(get.baseline==TRUE) plot(result$x.grid[[1]], ylim=c(-2.3,4),log(result$SBF2[[1]][2,]),col='green',lwd=2,type='l')
    
############## quality check and ISE calculation (see errors[2*k+2])     
     errors<-numeric(2*(d))
     errors[1]<-r
     errors[2]<-seed.n[r]
     for (k in 1:(d-1))
     {    
       errors[2*k+1]<-(1/length(xx.grid[[k+1]]))*sum((phi[[k]](xx.grid[[k+1]]) -LIN1[[k]])^2)
       # errors[[k]][s,2]<-sum((aa[[s]]$phi[[k]](aa[[s]]$x.grid[[k+1]]) - log(aa[[s]]$SBF2[[k+1]][2,]))^2)
       errors[2*k+2]<-(1/length(xx.grid[[k+1]]))*sum((phi[[k]](xx.grid[[k+1]]) - log(alpha_backfit.with.baseline[[k+1]]))^2)
      }
     
   write(errors,file=paste('ERRORS',myfilename,'.txt'),append=TRUE,ncolumns=(2*d))
     
   write(seed.n[r],file=paste('SEED',myfilename,'.txt'),append=TRUE,ncolumns=1)
 
   write(broke,file=paste('BROKE DOWN',myfilename,'.txt'),append=TRUE,ncolumns=1)
  
  comp.time<-c(r,seed.n[r],l.lin,l.sbf,comp.time.lin,comp.time.sbf)
  write(comp.time,file=paste('Comptime',myfilename,'.txt'),append=TRUE,ncolumns=6)
     
  bias<-numeric(2*d)
  bias[1]<-r
  bias[2]<-seed.n[r]
  for (k in 1:(d-1))
  {    
    bias[2*k+1]<-sum((phi[[k]](xx.grid[[k+1]]) - LIN1[[k]]))
    bias[2*k+2]<-sum((phi[[k]](xx.grid[[k+1]]) - log(alpha_backfit.with.baseline[[k+1]])))
  }
  write(bias,file=paste('BIAS',myfilename,'.txt'),append='TRUE',ncolumns=(2*d))
  
  LIN1<-matrix(unlist(LIN1),ncol=(d-1))
  SBF<-matrix(log(unlist(alpha_backfit.with.baseline)[-(1:nn.grid[1])]),ncol=(d-1))                                  # remove the estimation for time
  write.table(LIN1,file=paste('LIN_hazard',myfilename,'.txt'),append='TRUE',row.names=FALSE,col.names=FALSE)
  write.table(SBF,file=paste('SBF_hazard',myfilename,'.txt'),append='TRUE',row.names=FALSE,col.names=FALSE)
  
  hu<-length(read.delim(paste('SEED',myfilename,'.txt'),header = FALSE)[,1])
  me<-paste('Sim',hu,'of',n.sim,'Model',model, 'Dim=',d,'n=',n, 'rho=',rho)
  
  #if(hu==n.sim) dev.off()
  write(me, file='howfar.txt',append=TRUE)
  
    return(result)
  }
   
  seed.n<-x$seeds
  jobs<-x$jobs
  
  d<-jobs[[1]]
  rho<-jobs[[2]]
  model<-jobs[[3]]
  n<-jobs[[4]]
  n.sim<-jobs[[5]]
  n.cores<-jobs[[6]]
  r<-jobs[[7]]
    
  n.grid=rep(1,d)
  nb=2
  b.grid<-matrix(nrow=d,ncol=nb)
  b.grid[1,]<-c(5,5)
  b.grid[2:d,]<-matrix(c(0.2,0.2), byrow = TRUE, nrow=(d-1), ncol=2)
  k.min<-c(0,rep(-1,d-1))
  k.max<-c(28, rep(1,d-1))
  x.grid<-lapply(1:d, function (z) {seq(k.min[z],k.max[z],length=n.grid[z])})
  distr<-'normal'
  it=150
  kern<-function(u){return(0.75*(1-u^2)*(abs(u)<1))}
  
  d2<-3    ### plot columns
  d1<-ceiling((d-1)/3)###plot rows
  
  do.sim(n.cores=n.cores,n.sim=n.sim,n=n,r=r,nb=nb,n.grid=n.grid, b.grid=b.grid, x.grid=x.grid,model=model,rho=rho,d=d, it=it,kern=kern,distr=distr,seed.n=seed.n, get.baseline=TRUE,d1=d1,d2=d2) 
}

dots <- function(res, i, args = NULL){
  cat(".")
  if (i == length(res)){
    cat("\n")
  }
}
# Run code parallel
#numslaves<-detectCores()
#numslaves<-31
#cl <- makeCluster(numslaves)
#clusterExport(cl,'jobs')
start.time <- Sys.time() 
result<- mclapply(jobs, myparallelfunction,mc.cores=7)
end.time <-Sys.time()
time.taken <- end.time - start.time
time.taken  ## 7 cores 45 mins
#stopCluster(cl)
#mpi.exit()
for(ss in 1:length(jobs)){

distr<-'normal'
jobss<-jobs[[ss]]$jobs
d<-jobss[[1]]
rho<-jobss[[2]]
model<-jobss[[3]]
n<-jobss[[4]]
n.sim<-jobss[[5]]
n.cores<-jobss[[6]]
r<-jobss[[7]]
d1<-ceiling((d-1)/3)

myfilename<-paste('n.sim=',n.sim,'Model',model, 'Dim=',d,'distribution=', distr,'n=',n, 'rho=',rho)

if (ss==1) pdf(file=paste('plots',myfilename,'.pdf'),height = 4*d1, width = 15, onefile=TRUE)
#   quartz()
par(mfrow=(c(ceiling((d-1)/3),3)))
if (result[[ss]]$broke==TRUE){
for (k in 1:(d-1))
{
  y.min<-min(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))-1.5
  y.max<-max(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))+1.5
  plot(result[[ss]]$x.grid[[k+1]],result[[ss]]$LIN[[k]],col=4,type="l",lwd=2,xlab=expression('z'['k']),ylab="covariate effect",ylim=c(y.min,y.max),main=paste('SIM',r,' ','Model',model,' ','dim=',d,' ','n=',n,' ', 'rho=',rho,' ','k=',k))
  lines(result[[ss]]$x.grid[[k+1]],result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]),lty=3,col=1,lwd=2) 
  #lines(result$x.grid[[k+1]], log(result$SBF[[k+1]][2,]),col='green')
  lines(result[[ss]]$x.grid[[k+1]], log(result[[ss]]$SBF[[k+1]]),col='red',lwd=2)
  #   lines(x.grid[[k+1]], (alpha_backfit.exp[[k+1]][2,]),col='purple')
  legend("topleft",c("Lin He Huang (2016) (ITERATION BROKE DOWN)","SBF","True function"),col=c(4,'red',1),lty=c(1,1,3),lwd=2,bt="n")
}
}else{
  for (k in 1:(d-1))
  {
    y.min<-min(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))-1.5
    y.max<-max(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))+1.5
    plot(result[[ss]]$x.grid[[k+1]],result[[ss]]$LIN[[k]],col=4,type="l",lwd=2,xlab=expression('z'['k']),ylab="covariate effect",ylim=c(y.min,y.max),main=paste('SIM',r,' ','Model',model,' ','dim=',d,' ','n=',n,' ', 'rho=',rho,' ','k=',k))
    lines(result[[ss]]$x.grid[[k+1]],result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]),lty=3,col=1,lwd=2) 
    #lines(result$x.grid[[k+1]], log(result$SBF[[k+1]][2,]),col='green')
    lines(result[[ss]]$x.grid[[k+1]], log(result[[ss]]$SBF[[k+1]]),col='red',lwd=2)
    #   lines(x.grid[[k+1]], (alpha_backfit.exp[[k+1]][2,]),col='purple')
    legend("topleft",c("Lin He Huang (2016)","SBF","True function"),col=c(4,'red',1),lty=c(1,1,3),lwd=2,bt="n")
  }
}

#if(get.baseline==TRUE) plot(result$x.grid[[1]], ylim=c(-2.3,4),log(result$SBF2[[1]][2,]),col='green',lwd=2,type='l')

}
dev.off()

for(ss in 1:length(jobs))
  {  
  distr<-'normal'
  jobss<-jobs[[ss]]$jobs
  d<-jobss[[1]]
  rho<-jobss[[2]]
  model<-jobss[[3]]
  n<-jobss[[4]]
  n.sim<-jobss[[5]]
  n.cores<-jobss[[6]]
  r<-jobss[[7]]
  d1<-ceiling((1)/3)
  
  myfilename<-paste('n.sim=',n.sim,'Model',model, 'Dim=',d,'distribution=', distr,'n=',n, 'rho=',rho)
  
  if (ss==1) pdf(file=paste('ALLplots',myfilename,'.pdf'),height = 15, width = 15, onefile=TRUE)
  #   quartz()
  if(ss==1)par(mfrow=(c(ceiling((1)/3),1)))
  if (result[[ss]]$broke==TRUE) next
  
  x.0<-x.0.index<-numeric(d)
  x.0<-sapply(1:d, function(j) min(abs(result[[ss]]$x.grid[[j]]+0.5)))
  x.0.index<-sapply(1:d, function(j) which.min(abs(result[[ss]]$x.grid[[j]]+0.5)))                # why use the mins of absolute values of plus 0.5 and -0.5???? 
  
  x.1<-x.1.index<-numeric(d)
  x.1<-sapply(1:d, function(j) min(abs(result[[ss]]$x.grid[[j]]-0.5)))
  x.1.index<-sapply(1:d, function(j) which.min(abs(result[[ss]]$x.grid[[j]]-0.5)))
    
  if (d>=2){
    for (k in 1:(d-1))
    {  
      corr0 <-log(result[[ss]]$SBF[[k+1]][x.0.index[k+1]])-result[[ss]]$phi[[k]]( x.0[k+1])
      corr1 <-log(result[[ss]]$SBF[[k+1]][x.1.index[k+1]])-result[[ss]]$phi[[k]]( x.1[k+1])
      corr<-mean(c(corr1,corr0))                                                                 # use the mean of two correlation
      result[[ss]]$SBF[[k+1]] <- result[[ss]]$SBF[[k+1]]/ exp(corr)
      result[[ss]]$SBF[[1]] <- result[[ss]]$SBF[[1]]*  exp(corr)           
      #   corr <-alpha_backfit.exp[[k+1]][2,n.grid[k+1]/2]-phi[[k]](x.grid[[k+1]][n.grid[k+1]/2])
      #   alpha_backfit.exp[[k+1]][2,] <- alpha_backfit.exp[[k+1]][2,]- corr
      #  alpha_backfit[l,k+2,2,]<-alpha_backfit[l,k+2,2,] / exp(corr)
    }
  }    
  
  if (d>=2){
    for (k in 1:(1))
    { corr0 <-result[[ss]]$LIN[[k]][x.0.index[k+1]]-result[[ss]]$phi[[k]](x.0[k+1])
      corr1 <-result[[ss]]$LIN[[k]][x.1.index[k+1]]-result[[ss]]$phi[[k]](x.1[k+1])
      corr<-mean(c(corr1,corr0))      
      result[[ss]]$LIN[[k]]<-result[[ss]]$LIN[[k]]-corr
      #LIN1[[k+1]]<-LIN1[[k+1]]+LIN1[[k]][n.grid[k]/2]-phi[[k]](xx.grid[[k]][n.grid[k]/2])
    }
  }
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442")
  
  col1 <-    adjustcolor( cbbPalette [1], alpha.f = 0.4)
  col2 <-    adjustcolor( cbbPalette [2], alpha.f = 0.4)
  col3 <-    adjustcolor( cbbPalette [3], alpha.f = 0.4) 
    
    for (k in 1:(1))
    {
      y.min<-min(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))-1.5
      y.max<-max(result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]))+1.5
    if(ss==1)  plot(result[[ss]]$x.grid[[k+1]],result[[ss]]$LIN[[k]],col=col2,type="l",xlab=expression('z'['k']),ylab="covariate effect",ylim=c(y.min,y.max),main=paste('Model',model,' ','dim=',d,' ','n=',n,' ', 'rho=',rho,' ','k=',k)) else {
      lines(result[[ss]]$x.grid[[k+1]],result[[ss]]$LIN[[k]],col=col2,type="l",xlab=expression('z'['k']),ylab="covariate effect",ylim=c(y.min,y.max),main=paste('SIM',r,' ','Model',model,' ','d=',d-1,' ','n=',n,' ', 'rho=',rho,' ','k=',k))
    }
      lines(result[[ss]]$x.grid[[k+1]],result[[ss]]$phi[[k]](result[[ss]]$x.grid[[k+1]]),col=1,lwd=3) 
      #lines(result$x.grid[[k+1]], log(result$SBF[[k+1]][2,]),col='green')
      lines(result[[ss]]$x.grid[[k+1]], log(result[[ss]]$SBF[[k+1]]),col=col3)
      #   lines(x.grid[[k+1]], (alpha_backfit.exp[[k+1]][2,]),col='purple')
    if(ss==1)  legend("topleft",c("Lin He Huang (2016)","SBF","True function"),col=c(col2,col3,1),lty=c(1,1,1),lwd=2,bt="n")
    }
  
  }
  
dev.off()
save(result,file='result.RData')
########## END ##################
