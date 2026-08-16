MY.orig <- function(Z, Y, status, initial, h0, itera = 0.001, rrr = 100, it = 150) {
  # Z: n x p covariate matrix
  # Y: length n observed times
  # status: length n event indicator (1=event)
  # initial: n x p initial values for component functions (often 0 or "true" in sims)
  # h0: bandwidth (scalar, same for all covariates)
  # itera: outer-loop tolerance
  # rrr: starting grid index (original code uses 100)
  # it: max outer iterations (original code uses 150)
  
  Z <- as.matrix(Z)
  n <- length(Y)
  p <- ncol(Z)
  
  # Unique event times
  time <- sort(unique(Y[status == 1]))
  d_time <- length(time)
  
  # Risk set indicator matrix: n x d_time
  risk <- (outer(Y, rep(1, d_time)) >= outer(rep(1, n), time))
  
  # STATUS matrix: n x d_time, 1 at (i, t_j) if i fails at time t_j
  STATUS <- status * (outer(Y, rep(1, d_time)) == outer(rep(1, n), time))
  
  # Newton-Raphson objects
  score <- numeric(2)
  imat <- matrix(0, 2, 2)
  
  GBETA <- initial
  GBETA1 <- GBETA + 10  # just to start outer loop
  alpha0 <- c(0, 0)
  h <- h0
  r <- 1
  
  broke <- FALSE
  
  # Outer loop
  while ((max(abs(GBETA1 - GBETA), na.rm = TRUE) > itera) && (r <= it)) {
    GBETA1 <- GBETA
    
    for (k in 1:p) {
      # evaluation points for covariate k
      z0 <- sort(unique(round(Z[, k], 2)))
      
      # guard rrr
      if (length(z0) < 2) {
        broke <- TRUE
        return(list(l.lin = r, GBETA = GBETA, broke = broke))
      }
      if (rrr >= length(z0)) rrr_eff <- length(z0) - 1 else if (rrr <= 1) rrr_eff <- 2 else rrr_eff <- rrr
      
      # initial alpha at z0[rrr]
      alpha0[1] <- initial[round(Z[, k], 2) == z0[rrr_eff], k][1]
      alpha1 <- alpha0
      
      ## ---------- Forward sweep: j = rrr .. length(z0) ----------
      j <- rrr_eff
      while (j <= length(z0)) {
        # Epanechnikov kernel weights (same as original)
        wei <- 3/4 * (1 - (Z[, k] - z0[j])^2 / (h^2)) * (abs(Z[, k] - z0[j]) <= h) / h
        
        # local subset
        idx <- which(wei > 0)
        if (length(idx) == 0) {
          j <- j + 1
          next
        }
        
        Y0 <- Y[idx]
        Z0 <- Z[idx, k]
        wei0 <- wei[idx]
        n0 <- length(wei0)
        risk0 <- risk[idx, , drop = FALSE]
        STATUS0 <- STATUS[idx, , drop = FALSE]
        
        alpha1[1] <- initial[round(Z[, k], 2) == z0[j], k][1]
        alpha <- alpha1 + 1
        
        # fenmu: 1 x d_time (risk set denominator at each event time)
        # original: rep(1,n) %*% (risk * exp( GBETA %*% rep(1,p) ))
        linpred <- exp(as.numeric(GBETA %*% rep(1, p)))
        fenmu <- as.numeric(rep(1, n) %*% (risk * linpred))
        
        # inner NR loop
        while (sum(abs(alpha - alpha1), na.rm = TRUE) >= 1e-5) {
          alpha <- alpha1
          
          # wexpxa = kernel weights * exp(local linear) * exp(other components)
          wexpxa <- wei0 * exp(alpha[1] + alpha[2] * (Z0 - z0[j]))
          if (p > 1) {
            for (jj in 1:p) {
              if (jj != k) wexpxa <- wexpxa * exp(GBETA[idx, jj])
            }
          }
          
          # score[1]
          fenzi <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa))
          score[1] <- sum(as.numeric(wei0 %*% STATUS0)) - sum(as.numeric(STATUS %*% (fenzi / fenmu)))
          
          # score[2]
          fenzi2 <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])))
          score[2] <- sum(as.numeric(((Z0 - z0[j]) * wei0) %*% STATUS0)) - sum(as.numeric(STATUS %*% (fenzi2 / fenmu)))
          
          # Hessian imat
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa))
          imat[1, 1] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])))
          imat[1, 2] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          imat[2, 1] <- imat[1, 2]
          
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])^2))
          imat[2, 2] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          
          # NR update
          alpha1 <- tryCatch(
            alpha - as.numeric(solve(imat) %*% score),
            error = function(e) NA
          )
          
          if (is.na(sum(alpha1))) break
        }
        
        if (is.na(sum(alpha1))) break
        
        # update function values at points that equal z0[j] (rounded)
        GBETA[round(Z[, k], 2) == z0[j], k] <- alpha1[1]
        j <- j + 1
      }
      
      if (is.na(sum(alpha1))) {
        broke <- TRUE
        return(list(l.lin = r, GBETA = GBETA, broke = broke))
      }
      
      ## ---------- Backward sweep: j = rrr-1 .. 1 ----------
      alpha0[1] <- initial[round(Z[, k], 2) == z0[rrr_eff - 1], k][1]
      alpha1 <- alpha0
      
      j <- rrr_eff - 1
      while (j >= 1) {
        wei <- 3/4 * (1 - (Z[, k] - z0[j])^2 / (h^2)) * (abs(Z[, k] - z0[j]) <= h) / h
        
        idx <- which(wei > 0)
        if (length(idx) == 0) {
          j <- j - 1
          next
        }
        
        Y0 <- Y[idx]
        Z0 <- Z[idx, k]
        wei0 <- wei[idx]
        n0 <- length(wei0)
        risk0 <- risk[idx, , drop = FALSE]
        STATUS0 <- STATUS[idx, , drop = FALSE]
        
        alpha1[1] <- initial[round(Z[, k], 2) == z0[j], k][1]
        alpha <- alpha1 + 1
        
        linpred <- exp(as.numeric(GBETA %*% rep(1, p)))
        fenmu <- as.numeric(rep(1, n) %*% (risk * linpred))
        
        while (sum(abs(alpha - alpha1), na.rm = TRUE) >= 1e-5) {
          alpha <- alpha1
          
          wexpxa <- wei0 * exp(alpha[1] + alpha[2] * (Z0 - z0[j]))
          if (p > 1) {
            for (jj in 1:p) {
              if (jj != k) wexpxa <- wexpxa * exp(GBETA[idx, jj])
            }
          }
          
          fenzi <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa))
          score[1] <- sum(as.numeric(wei0 %*% STATUS0)) - sum(as.numeric(STATUS %*% (fenzi / fenmu)))
          
          fenzi2 <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])))
          score[2] <- sum(as.numeric(((Z0 - z0[j]) * wei0) %*% STATUS0)) - sum(as.numeric(STATUS %*% (fenzi2 / fenmu)))
          
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa))
          imat[1, 1] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])))
          imat[1, 2] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          imat[2, 1] <- imat[1, 2]
          
          fz <- as.numeric(rep(1, n0) %*% (risk0 * wexpxa * (Z0 - z0[j])^2))
          imat[2, 2] <- -sum(as.numeric(STATUS %*% (fz / fenmu)))
          
          alpha1 <- tryCatch(
            alpha - as.numeric(solve(imat) %*% score),
            error = function(e) NA
          )
          if (is.na(sum(alpha1))) break
        }
        
        if (is.na(sum(alpha1))) break
        
        GBETA[round(Z[, k], 2) == z0[j], k] <- alpha1[1]
        j <- j - 1
      }
      
      if (is.na(sum(alpha1))) {
        broke <- TRUE
        return(list(l.lin = r, GBETA = GBETA, broke = broke))
      }
      
      # Anchoring / normalization (original MY.orig does this, using "rrr point" anchoring in some variants;
      # your pasted MY.orig uses: subtract at z0[rrr] and add initial at z0[rrr].
      # Here I keep the exact line you had in the pasted MY.orig:
      GBETA[, k] <- GBETA[, k] - GBETA[round(Z[, k], 2) == z0[rrr_eff], k][1] + initial[round(Z[, k], 2) == z0[rrr_eff], k][1]
    }
    
    r <- r + 1
    # print progress like original
    # print(c(r, max(abs(GBETA1 - GBETA), na.rm = TRUE)))
  }
  
  return(list(l.lin = r, GBETA = GBETA, broke = broke))
}