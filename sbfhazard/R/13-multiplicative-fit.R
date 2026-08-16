# Multiplicative SBF engine ---------------------------------------------------

                  
.sbf_fit_multiplicative_engine <- function(data,
                                           b.grid,
                                           x.grid = NULL,
                                           n.grid.additional = 0L,
                                           x.min = NULL,
                                           x.max = NULL,
                                           integral.approx = "midd",
                                           it = 100L,
                                           kernel = "epanechnikov",
                                           initial = NULL,
                                           convergence_tol = 0.001,
                                           identification = "sample_mean",
                                           truth_functions = NULL,
                                           verbose = FALSE)
{                               
  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  kernel_fun <- kernel_spec$kernel
  kernel_name <- kernel_spec$name
  min_component <- 1e-8

  smooth.alpha<-function(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
  {
    alpha.smooth.i<-array(dim=c(d,n))
    
    alpha.smooth.i[1,] <-rep( 1,n)
    
    alpha.smooth.i.0 <- as.numeric((K.b/k.b) %*% (alpha[[1]] * dx[[1]]))
    
    
    for (k in 2:d){
      K.norm <- K.X.b[[k]] / k.X.b[[k]]
      alpha.smooth.i[k,] <- as.numeric(K.norm %*% (alpha[[k]] * dx[[k]]))
    }
    return(list(alpha.smooth.i=alpha.smooth.i,alpha.smooth.i.0=alpha.smooth.i.0))
  }
  get.new.alpha<-function(data,alpha,K.X.b,k.X.b,K.b.norm,Y,k,x.grid,dx,n.grid,d,n)
  {
    
    alpha.smooth.i<-smooth.alpha(alpha,K.X.b,k.X.b,K.b,k.b,x.grid,dx,n.grid,d,n)
    
    if (k==1) alpha.minusk.smooth<-numeric(n) else  alpha.minusk.smooth<-array(dim=c(n,n.grid[1])) 
    
    if (k==1) {
      alpha.minusk.smooth <- apply(alpha.smooth.i$alpha.smooth.i[-1, , drop = FALSE], 2, prod)
    } else {
      prod.minusk <- apply(alpha.smooth.i$alpha.smooth.i[-k, , drop = FALSE], 2, prod)
      alpha.minusk.smooth <- outer(prod.minusk, alpha.smooth.i$alpha.smooth.i.0)
    }
    
    if (k==1)
    {
      # A <- sweep(Y, 1, alpha.minusk.smooth, "*")
      # A <- sweep(A, 2, dx[[1]], "*")
      # D <- colSums(A %*% K.b.norm)
      weighted_risk <- as.numeric(crossprod(alpha.minusk.smooth, Y)) * dx[[1]]
      D <- as.numeric(weighted_risk %*% K.b.norm)
    }else {
      # scalar <- rowSums((Y * alpha.minusk.smooth) * matrix(dx[[1]], nrow = n, ncol = n.grid[1], byrow = TRUE))
      scalar <- as.numeric((Y * alpha.minusk.smooth) %*% dx[[1]])
      K.norm <- K.X.b[[k]] / k.X.b[[k]]
      # D <- colSums(K.norm * scalar)
      D <- as.numeric(crossprod(K.norm, scalar))
    }
    
    K.norm <- K.X.b[[k]] / k.X.b[[k]]
    # O <- colSums(K.norm * status)
    O <- as.numeric(crossprod(K.norm, status))
    
    # print(paste(O))
    # print(paste(D))
    return(O/D)
  }
  d<-ncol(data)-1
  n<-nrow(data)
  status <- as.numeric(data$status)
  
  grid_spec <- .sbf_prepare_grid_spec(
    observed_grid = lapply(seq_len(d), function(k) data[order(data[, k]), k]),
    x_grid = x.grid,
    n_grid_additional = n.grid.additional,
    x_min = x.min,
    x_max = x.max,
    integral_approx = integral.approx
  )
  x.grid <- grid_spec$x_grid
  n.grid <- grid_spec$n_grid
  dx <- grid_spec$dx
  
  x.grid.array<-matrix(rep(x.grid[[1]], times=n.grid[1]),nrow = n.grid[1], ncol = n.grid[1],byrow=FALSE)   # 
   u<-x.grid.array-t(x.grid.array) # u is x_j-u_j for all x,u on the grid considered
  
  
  K.b<-.sbf_eval_kernel(u, b.grid[1], kernel_fun)
  k.b<-colSums(dx[[1]]*K.b)
  K.b.norm <- K.b / k.b
  # ## change
  # k.b <- as.vector(K.b %*% dx[[1]])
  
  
  X <- as.matrix(data[, seq_len(d), drop = FALSE])
  
  Y <- (outer(data$time, x.grid[[1]], FUN = ">=") * 1)
  

 
  K.X.b<-k.X.b<-list()
  for( k in 1:d){
    K.X.b[[k]]<-array(0,dim=c(n,n.grid[k]))
    k.X.b[[k]]<-numeric(n)
    u <- outer(X[, k], x.grid[[k]], FUN = "-")
    K.X.b[[k]] <- .sbf_eval_kernel(u, b.grid[k], kernel_fun)
    # k.X.b[[k]] <- rowSums(K.X.b[[k]] * matrix(dx[[k]], nrow = n, ncol = n.grid[k], byrow = TRUE))
    k.X.b[[k]] <- as.numeric(K.X.b[[k]] %*% dx[[k]])
  }

  normalize_alpha_components <- function(alpha) {
    for (j in 2:d) {
      scale.j <- .sbf_multiplicative_identification_scale(
        alpha_component = alpha[[j]],
        grid = x.grid[[j]],
        widths = dx[[j]],
        identification = identification,
        min_component = min_component,
        truth_functions = truth_functions,
        covariate_index = j - 1L
      )
      alpha[[j]] <- alpha[[j]] / scale.j
      alpha[[1]] <- alpha[[1]] * scale.j
    }

    alpha
  }
  
  
  alpha_backfit<-list()
  
  alpha_backfit <- if (is.null(initial)) {
    lapply(n.grid, function(n) rep(1, n))
  } else {
    initial
  }

  iterations_used <- 0L
  final_delta <- NA_real_
  converged <- FALSE
  component_floor_count <- 0L
 
   for (iter in seq_len(it))
  {  
    alpha_backfit_old<-alpha_backfit
    for(k in 1:d)
    {
  
      
      
      updated <- get.new.alpha(data,alpha_backfit,K.X.b,k.X.b,K.b.norm,Y,k,x.grid,dx,n.grid,d,n)
                                       
      ###########
     # for (j in 1:(d-1))
     # {   alpha_backfit[[j]]<- 1*(alpha_backfit[[j]]/ alpha_backfit[[j]][1]) }
      
        
      floor_idx <- !is.finite(updated) | updated < min_component
      component_floor_count <- component_floor_count + sum(floor_idx)
      updated[floor_idx] <- min_component
      alpha_backfit[[k]] <- updated
      
    }

    alpha_backfit <- normalize_alpha_components(alpha_backfit)

    final_delta <- .sbf_fit_delta(alpha_backfit_old, alpha_backfit)
    iterations_used <- iter
    if (is.finite(final_delta) && final_delta <= convergence_tol) {
      converged <- TRUE
      break
    }
    if (verbose) print(c(iter, final_delta))
  }

  eps <- 1e-12
  fit_diagnostics <- .sbf_kernel_norm_diagnostics(
    kernel_norm = list(k.b, k.X.b),
    floor_value = eps
  )
  fit_diagnostics$component_floor_count <- as.integer(component_floor_count)
  k.b.safe <- pmax(as.numeric(k.b), eps)
  k.X.b.safe <- lapply(k.X.b, function(x) pmax(as.numeric(x), eps))
  alpha.smooth <- smooth.alpha(alpha_backfit, K.X.b, k.X.b.safe, K.b, k.b.safe, x.grid, dx, n.grid, d, n)
  alpha.smooth.obs <- alpha.smooth$alpha.smooth.i
  alpha.smooth.time <- alpha.smooth$alpha.smooth.i.0

  dx1 <- as.numeric(dx[[1]])
  baseline.multiplier <- if (d > 1L) {
    apply(alpha.smooth.obs[-1, , drop = FALSE], 2, prod)
  } else {
    rep(1, n)
  }
  baseline.exposure <- as.numeric(crossprod(baseline.multiplier, Y)) * dx1
  time.exposure <- as.numeric(Y %*% (alpha.smooth.time * dx1))

  covariate.scalar <- vector("list", d)
  if (d >= 2L) {
    for (k in 2:d) {
      prod.minusk <- apply(alpha.smooth.obs[-k, , drop = FALSE], 2, prod)
      covariate.scalar[[k]] <- prod.minusk * time.exposure
    }
  }

  multiplicative_predict_state <- list(
    bandwidth = as.numeric(b.grid),
    kernel = kernel_fun,
    train_x = X,
    dx = dx,
    k.b = k.b.safe,
    k.X.b = k.X.b.safe,
    status = status,
    baseline.exposure = baseline.exposure,
    covariate.scalar = covariate.scalar
  )

  return(list(
    alpha_backfit=alpha_backfit,
    x.grid=x.grid,
    data=data,
    converged=converged,
    iterations_used=iterations_used,
    iterations_max=it,
    final_delta=final_delta,
    convergence_tol=convergence_tol,
    fit_diagnostics=fit_diagnostics,
    multiplicative_predict_state=multiplicative_predict_state
  ))
}


# Internal multiplicative SBF fit adapter used by the public API wrapper.
.sbf_mult_fit <- function(data,
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
                                identification = "sample_mean",
                                truth_functions = NULL,
                                verbose = FALSE,
                                warn_nonconvergence = TRUE) {
  kernel_spec <- .sbf_kernel_spec(kernel, arg_name = "kernel")
  feature_names <- setdiff(names(data), c("time", "status"))

  fit <- .sbf_fit_multiplicative_engine(
    data = data,
    b.grid = bandwidth,
    x.grid = x_grid,
    n.grid.additional = n_grid_additional,
    x.min = x_min,
    x.max = x_max,
    integral.approx = integral_approx,
    it = iterations,
    kernel = kernel,
    initial = initial,
    convergence_tol = convergence_tol,
    identification = identification,
    truth_functions = truth_functions,
    verbose = verbose
  )

  fit$converged <- isTRUE(fit$converged)
  fit$model_type <- "multiplicative"
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
    min_component = 1e-8,
    identification = identification
  )
  class(fit) <- unique(c("sbf_multiplicative_fit", class(fit)))
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
