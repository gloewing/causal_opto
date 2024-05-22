# binary SEs
binary_SE <- function(data, 
                      beta_hat, # beta estimate
                      model = NULL, # requires either GLM model (warm-start) or X, Y, w
                      X = NULL,
                      Y = NULL,
                      w = NULL,
                      id_name = "id",
                      id = NULL
                      ){
  
  if(is.null(model) & (is.null(X) | is.null(Y) | is.null(w) ) ){
    message("Must provide either GLM fit or X, Y, w")
  }
  
  setDF(data)
  n <- length(unique(data$id))
  m <- function(x, beta){ plogis(x %*% beta)[,1] }
  if(is.null(id)){
    id_vec <- as.vector(model$dat[,id_name]) # rename id
    ids <- unique(id_vec)
  }else{
    id_vec <- as.vector(id) # rename id
    ids <- unique(id_vec)
  }

  
  if(!is.null(model)){
    X <- as.matrix(model.matrix(model)) # design matrix
    Y <- model$y # outcome vector
    w <- model$prior.weights # weights
  }
  
  ### standard error calculations
  phi <- function(id, beta, X, Y, w, id_vec) {
    idx = which(id_vec == id)
    x_sub <- X[idx,]
    ww_sub <- w[idx]
    YY_sub <- Y[idx]
    
    m_res <- m(x_sub, beta)
    return ( colSums(x_sub * m_res * (1 - m_res) * ww_sub * (YY_sub - m_res)) )
  }
  
  phi_grad <- function(id, beta, X, Y, w, id_vec) {
    idx = which(id_vec == id)
    x_sub <- X[idx,]
    ww_sub <- w[idx]
    YY_sub <- Y[idx]
    
    m_res <- m(x_sub, beta)
    scls <- ww_sub * m_res * (1 - m_res) * (
      (1 - 2*m_res) * (YY_sub - m_res) - m_res * (1 - m_res) )
    #out_prods <- t(x_sub) %*% diag(scls) %*% x_sub # could be faster: crossprod(x_sub, diag(-ww_sub)) %*% x_sub  #t(x_sub) %*% diag(-ww_sub) %*% x_sub 
    out_prods <- crossprod(x_sub, x_sub * scls)
    # out_prods <- lapply(1:nrow(x_sub), 
    #                     function(i){ x_sub[i,] %*% t(x_sub[i,]) } )  # tcrossprod(x_sub[i,])

    # Reduce('+', 
    #        lapply(1:nrow(x_sub), function(i) { out_prods[[i]] * scls[i]})
    #        )
  }
  
  ## estimating functions themselves
  phi_mat <- t( sapply(ids, 
                    FUN = function(id){ phi(id, beta_hat, X, Y, w, id_vec) }, 
                    simplify = 0) )
  
  ## Jacobians of estimating functions
  phi_grad_list <- lapply(ids, 
                          FUN = function(id){phi_grad(id, beta_hat, X, Y, w, id_vec)})
  
  A_hat <- crossprod(phi_mat) / n #( t(phi_mat) %*% phi_mat ) / n 
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
  
  return( list(SEs = SEs, CIs = CIs, cov = V_hat/n ) )
  
}

