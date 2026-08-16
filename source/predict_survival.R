############ survival prediction ############
# This function indicates that: for each combination of covariates Z, we should consider different time values
# that is why the dimension of survival.prob is either m*m or m*(n-1)
# when the time vector is provided, it could be not ordered
# only compute the prediction for time in this function

predict.survival<-function(result, data, times=NULL){
## INPUTS:
  # result:                  list, the result of function SBF.MH_LC
  # data:        m*(d-1)     matrix, the data set for prediction, without status and time
  # times:       m           vector, the time of data, if null, using the x.grid[[1]] instead
## OUTPUTS:
  # list(survival probability, time)
  
  data<-as.matrix(data)                                                       # new data set, without time, m*(d-1)
  x.grid<-result$x.grid
  alpha_sbf<-result$alpha_backfit
  if (is.null(times)) surv.times<-x.grid[[1]][-1] else surv.times<-times      # vector, ordered time, (n-1)
  surv.prob <- matrix(nrow=nrow(data), ncol=length(surv.times))               # empty matrix , m*(n-1)
  hazard.values<-predict.hazard(result,cbind(rep(1,nrow(data)),data))         # prediction with nearest points, m*(d+1), time =1 for all i
  
  for(i in 1:nrow(data)){
    # find.index<- error<-numeric(length(alpha_sbf))
    # for (j in 2:length(alpha_sbf)){
    #   error[j] <- min(abs(as.numeric(as.numeric(data[i,j-1]))-x.grid[[j]]))
    #   find.index[j]<- which.min(abs(as.numeric(data[i,j-1])-x.grid[[j]]))
    # }
    
    par<-1
    for (j in 2:length(alpha_sbf)){                                           # exclude time
      par <- par * hazard.values[i,j]                                         # the product of predicted alpha without time, scalar
    }
    
    if (is.null(times)){
      # why taking exp and using the alpha_2*...*alpha_d ??? cumulative prob?
      # the cumulative sum is affected by the normalization method, here, he did not use any normalization
      # minus sign in the begining to make the plot on the right hand
      surv.prob[i,]<-exp(-par*cumsum(alpha_sbf[[1]][-1]*diff(x.grid[[1]])))   # cumulative sum of alpha*dx, only using the computed to approximate integration, vector, (n-1)
    } else{ 
      for (j in 1:length(times)) {                                            # if time is provided, using the provided time
        index <- which.min(abs(result$x.grid[[1]]-times[j]))                  # index of nearest point
        error <-  result$x.grid[[1]][index] - times[j]                        # difference between the nearest point and the one we wanted to predict
        
        if (error==0)  surv.prob[i,j] <- exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index]))) else{   if (error>0&index==2) surv.prob[i,j] <- exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index])))
        if ((error>0&index>2) | (error<0& index>1)) {
          if (error>0) {index2<-index-1} 
          if (error<0)  {index2<-index+1}
          error2<-   result$x.grid[[1]][index2] - times[j]          
          temp1<-  exp(-par*sum(alpha_sbf[[1]][2:index]*diff(x.grid[[1]][1:index])))
          temp2<-  exp(-par*sum(alpha_sbf[[1]][2:index2]*diff(x.grid[[1]][1:index2])))
          surv.prob[i,j] <- (abs(error)*temp1 + abs(error2)*temp2)/(abs(error)+abs(error2))}
        }
        if (index==1) surv.prob[i,j]  <- 1
      }}
    
  }
  return(list(surv.prob=surv.prob, surv.times= surv.times))
}
