# sandwich variance estimator for linear case
sandwich_normal <- function(data, 
                            model = NULL, # requires either GLM model (warm-start) or X, Y, w
                            beta_hat = NULL, # beta estimate
                            X = NULL,
                            Y = NULL,
                            w = NULL,
                            id_name = "id"
                          ){
  
  if(is.null(model) & (is.null(X) | is.null(Y) | is.null(w) ) )  stop("Must provide either lm() fit or X, Y, w")
  if(is.null(model) & (is.null(beta_hat)))  stop("Must provide either lm() fit or beta_hat")
  setDF(model$dat)
  
  n <- length(unique(data$id))
  m <- function(x, beta){ (x %*% beta)[,1] }
  id_vec <- as.vector(model$dat[,id_name]) # rename id
  ids <- unique(id_vec)
  
  if(is.null(beta_hat)){
    beta_hat <- as.numeric(coef(model$model))
  }
  
  if(!is.null(model)){
    X <- model$model$x # design matrix
    Y <- model$model$y # outcome vector -- ** assumes outcome vector is always first vector in model$model object
    w <- model$model$weights # weights
  }
  
  ### standard error calculations
  phi <- function(id, beta, X, Y, w, id_vec) {
    idx <- which(id_vec == id)
    x_sub <- X[idx,]
    ww_sub <- w[idx]
    YY_sub <- Y[idx]
    
    m_res <- m(x_sub, beta)
    return ( colSums(x_sub * ww_sub * (YY_sub - m_res)) )
    # cnst <- ww_sub * (YY_sub - m_res)
    # return ( colSums( apply(x_sub, 2, function(x) x * cnst) ) )
  }
  
  phi_grad <- function(id, beta, X, Y, w, id_vec) {
    idx <- which(id_vec == id)
    x_sub <- X[idx,]
    ww_sub <- w[idx]
    YY_sub <- Y[idx]
    
    m_res <- m(x_sub, beta)
    # out_prods <- t(x_sub) %*% diag(-ww_sub) %*% x_sub # slow
    out_prods <- crossprod(x_sub, -ww_sub * x_sub )  
    
    # lapply(1:nrow(x_sub), 
      #                   function(i){ x_sub[i,] %*% t(x_sub[i,])}  )  # tcrossprod(x_sub[i,]) }
    #scls <- -ww_sub
    
    return( out_prods )
      # Reduce('+', 
      #           lapply(1:nrow(x_sub), function(i) { out_prods[[i]] * scls[i] }) ) )

  }
  
  ## estimating functions themselves
  phi_mat <- t( sapply(ids, 
                      FUN = function(id){ phi(id, beta_hat, X, Y, w, id_vec) }, 
                      simplify = 0) )
  
  ## Jacobians of estimating functions
  phi_grad_list <- lapply(ids, FUN = function(id){ phi_grad(id, beta_hat, X, Y, w, id_vec) })
  
  A_hat <- crossprod(phi_mat) / n
  B_hat <- Reduce("+", phi_grad_list) / n
  
  B_inv <- solve(B_hat)
  V_hat <- B_inv %*% A_hat %*% B_inv
  
  SEs <- sqrt(diag(V_hat) / n) ## hand-coded standard errors
  
  # convergence_check <- rowSums(phi_mat) # should all be close to 0
  # if( any(convergence_check > convg_tol) ){
  #   message("Convergence not reached!")
  #   print(convergence_check)
  # }
  
  # confidence intervals
  CIs <- cbind( beta_hat - qnorm(0.975) * SEs,
                beta_hat + qnorm(0.975) * SEs) ## hand-coded CIs
  
  rownames(CIs) <- names(beta_hat)
  names(SEs) <- names(beta_hat)
  
  return( list(SEs = SEs, CIs = CIs, cov = V_hat/n) )
}
