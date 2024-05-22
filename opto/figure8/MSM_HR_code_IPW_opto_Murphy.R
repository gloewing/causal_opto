# IPW for optoDA delta = 1 conditional on X = 1

msm_HR_optoIPW_Murphy <- function(data,
                   formula,
                   family = binomial,
                   delta = 1,
                   A_name,
                   Y_name,
                   X_name,
                   W_name,
                   id_name,
                   trial_name,
                   session_name = NA,
                   h = NA # numerator weights - to update later
){
  #######################################################################
  
  ids <- sort(unique(data[,id_name]))
  Wt <- data[,W_name]
  W <- rep(NA, length(Wt))
  A <- data[,A_name] # save original treatment and time-varying confounder A, X
  X <- data[,X_name] # save original treatment and time-varying confounder A, X
  data$A <- A
  data$X <- X
  data$I0 <- 1 * I(A == 0) 
  data$I1 <- 1 * I(A == X)
  data$ids <- data[,id_name]
  n <- nrow(data)
  data$weight <- rep(NA, n)
  if(is.na(session_name)){
    data$session <- 1 # only sessio n
  }else{
    data$session <- data[,session_name]
  }
  
  ####################
  # shift by 1 
  ####################
  # calculate shift on the biggest k to start with so do not have to do this repeatedly below
  setDT(data)
  start <- 1 # 0 means start summing at current trial, 1 starts at next trial and 2 skips the subsequent trial
  k <- 1
  lagged_vars <- c("A_lag", "X_lag")
  data[, c(lagged_vars) := (shift( .(A, X), type = "lag", fill = NA)), keyby = .(ids, session) ] # goes back 1
  data <- data[data$X_lag == 1,] # remove 0s
  
  #######################################################################################################
  # calculate denominator weights as just inverse prob weights because of conditioning on X_{t-1}
  #######################################################################################################
  message("denominator weights")
  
  # numerator/stabilizing weights
  if(is.na(h))   h <- 1
  data$weight <- h / data$weights
 
  # # define variables inside dataset
  data$ww <- data$weight #* regime_ind
  # data$trial <- data[, trial_name]
  data$YY <- data[,..Y_name]
  data <- data[data$ww != 0,] # remove 0s
  
  ####################
  # MSM
  ####################
  # make more general in mean specification with products of d's using string in `formula` R syntax
  # data$d_sum <- as.factor(rowSums(data[,c("d_1", "d_2", "d_3", "d_4", "d_5")]))
  
  # in general should be a glm
  # data[,Y_name]
  if(family == "gaussian"){
    mod <- stats::lm(YY ~ A_lag * preCount, # splines::bs(trial, df = 4)
                      weights = ww,
                      data = data)
  }else{
    data$YY <- as.integer(data$YY) # ensure converted
    
    mod <- stats::glm(YY ~ A_lag * preCount, # splines::bs(trial, df = 4)
                     family = family,
                     weights = ww,
                     data = data)
  }
  
  return(list(dat = data,
              model1 = mod))
  
}

