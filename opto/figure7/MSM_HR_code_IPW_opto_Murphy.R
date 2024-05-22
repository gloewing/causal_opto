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
  # saved in /home/folder/msm_hr
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
  
  # indicator matrix for treatment regime 
  # indicator_matrix <- expand.grid( rep(list(0:1), delta))
  # colnames(indicator_matrix) <- paste0("d_", 1:delta)
  # q <- nrow(indicator_matrix)
  # data_orig <- data
  # data <- data_orig[rep(seq_len(nrow(data_orig)), each = q),] # copy
  # data <- cbind(data, indicator_matrix[rep(1:q, n),]) # copies each row 2^delta
  # N <- nrow(data)
  # regime_ind <- rep(NA, N)
  # 
  # for(i in 1:length(ids)){
  #   message(paste0("regime calculation for subject: ", i))
  #   
  #   idx <- which(data[id_name] == ids[i]) # copied data indices
  #   idx_orig <- which(data_orig[id_name] == ids[i]) # original dataset indices
  #   for(jj in 1:q){
  #     # iterate through indicator matrix to find specific regime (a row of indicator_matrix)
  #     # use those regimes (+1) as indicators to subset I0 and I1 columns of original data
  #     # then keep track of staggering between data and data_orig
  #     regime <- as.numeric(indicator_matrix[jj,] + 1) # add one because we use as indices below
  #     idx_copy <- idx[seq(jj, length(idx), by = q)] # find rows of copied matrix corresponding to the current regime
  #     data_subset <- data_orig[idx_orig,] # use ids of current animal on original (non-copied) dataset
  #     regime_ind[idx_copy] <- c( rep(0, delta - 1), 
  #                                sapply( 1:(length(idx_orig)-delta+1),  function(dd){
  #                                 prod(data_subset[dd:(dd+delta-1), c("I0", "I1")][cbind(1:delta, regime)]) #
  #                                } 
  #                               ))
  #   }
  # }
  # 
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
    mod <- stats::lm(YY ~ A_lag, 
                      weights = ww,
                      data = data)
  }else{
    data$YY <- as.integer(data$YY) # ensure converted
    
    mod <- stats::glm( YY ~ A_lag, 
                     family = family,
                     weights = ww,
                     data = data)
  }
  
  setDF(data)

  return(list(dat = data,
              model1 = mod))
  
}

