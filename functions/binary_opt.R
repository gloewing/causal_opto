# binary optimization
bin_msm_splines <- function(model,
                    data, # original dataset GLM warm start applied to
                    opt = "root", # optimization applied "gee", "loss"
                    method = "L-BFGS-B", # optimization method applied if using "loss"
                    id, # ID variable if GEE applied
                    convg_tol = 1e-3 # check convergence
                    ){
  
  source("~/Desktop/NIMH Research/Causal/Code/Functions/binary_SEs.R")
  
  beta0 <- coef(model) # warm start fit
  X <- model.matrix(model) # design matrix
  Y <- model$y # outcome vector
  w <- model$prior.weights # weights
  n <- length(unique(id)) # number of clusters/subjects

  ## fitted value function
  m <- function(x, beta){ plogis(x %*% beta)[,1] }
  
  ## estimating function at  beta
  est_function <- function(beta) {
    fit.vals <- m(X, beta)
    colSums( X * (fit.vals * (1 - fit.vals) * w * (Y - fit.vals)) )
  }
  
  estim.precis <- NA
  
  if(opt == "root"){
    # root solve
    beta_hat <- rootSolve::multiroot(f = est_function, start = beta0)
    estim.precis <- beta_hat$f.root # estimating equation evaluated at optimum
    beta_hat <- beta_hat$root
    
    if( any( abs(estim.precis) > convg_tol ) ){
      message("Convergence not reached!")
      print(estim.precis)
    }
  }else if(opt == "loss"){
    # use loss function
    obj_function <- function(beta){ sum(ww*(Y - m(X, beta))^2) }
    obj_gr <- function(beta){ -2 * est_function(beta) }

    beta_hat <- optimx::optimr(par = beta0, 
                               fn = obj_function, 
                               gr = obj_gr,
                               method = method)$par
  }else if(opt == "gee"){
    mod <- geese(YY ~ d_sum * preCount,
                 mean.link = "logit",
                 variance = "gaussian",
                 scale.fix = TRUE,
                 scale.value = 1,
                 weights = w,
                 id = id,
                 corstr = "independence",
                 data = data)
    beta_hat <- mod$beta
  }
  
  return(list(beta_hat = beta_hat,
              estim.precis = estim.precis,
              opt = opt) )
}