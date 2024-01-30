# IPW for optoDA delta = 5

msm_HR_optoIPW <- function(data,
                   formula,
                   family = binomial,
                   delta = 2,
                   A_name,
                   Y_name,
                   X_name,
                   W_name,
                   id_name,
                   trial_name,
                   session_name = NA,
                   h = NA # numerator weights - to update later
){
  # saved in /home/loewingergc/msm_hr
  #######################################################################
  
  ids <- sort(unique(data[,id_name]))
  Wt <- data[,W_name]
  W <- rep(NA, length(Wt))
  A <- data[,A_name] # save original treatment and time-varying confounder A, X
  X <- data[,X_name] # save original treatment and time-varying confounder A, X
  data$I0 <- 1 * I(A == 0) 
  data$I1 <- 1 * I(A == X)
  n <- nrow(data)
  data$weight <- rep(NA, n)
  
  ###################################################################################
  # calculate denominator weights as product of moving window of delta time-points
  ###################################################################################
  message("denominator weights")
  for(i in 1:length(ids)){
    idx <- which(data[id_name] == ids[i])
    WWtt <- Wt[idx]
    W[idx] <- c( rep(Inf, delta - 1), # use Inf so it will make weights 0 when invert
                 sapply( 1:(length(idx)-delta+1),  function(dd) prod(WWtt[dd:(dd+delta-1)])  ) 
                )
  }
  
  # numerator/stabilizing weights
  if(is.na(h))   h <- 1
  data$weight <- h / W
  
  # indicator matrix for treatment regime 
  indicator_matrix <- expand.grid( rep(list(0:1), delta))
  colnames(indicator_matrix) <- paste0("d_", 1:delta)
  q <- nrow(indicator_matrix)
  data_orig <- data
  data <- data_orig[rep(seq_len(nrow(data_orig)), each = q),] # copy
  data <- cbind(data, indicator_matrix[rep(1:q, n),]) # copies each row 2^delta
  N <- nrow(data)
  regime_ind <- rep(NA, N)
  
  for(i in 1:length(ids)){
    message(paste0("regime calculation for subject: ", i))
    
    idx <- which(data[id_name] == ids[i]) # copied data indices
    idx_orig <- which(data_orig[id_name] == ids[i]) # original dataset indices
    for(jj in 1:q){
      # iterate through indicator matrix to find specific regime (a row of indicator_matrix)
      # use those regimes (+1) as indicators to subset I0 and I1 columns of original data
      # then keep track of staggering between data and data_orig
      regime <- as.numeric(indicator_matrix[jj,] + 1) # add one because we use as indices below
      idx_copy <- idx[seq(jj, length(idx), by = q)] # find rows of copied matrix corresponding to the current regime
      data_subset <- data_orig[idx_orig,] # use ids of current animal on original (non-copied) dataset
      regime_ind[idx_copy] <- c( rep(0, delta - 1), 
                                 sapply( 1:(length(idx_orig)-delta+1),  function(dd){
                                  prod(data_subset[dd:(dd+delta-1), c("I0", "I1")][cbind(1:delta, regime)]) #
                                 } 
                                ))
    }
  }
  
  # define variables inside dataset
  data$ww <- data$weight * regime_ind
  data$trial <- data[, trial_name]
  data$YY <- data[,Y_name]
  data <- data[data$ww != 0,] # remove 0s
  
  ####################
  # MSM
  ####################
  # make more general in mean specification with products of d's using string in `formula` R syntax
  data$d_sum <- as.factor(rowSums(data[,c("d_1", "d_2", "d_3", "d_4", "d_5")]))
  
  # in general should be a glm
  # data[,Y_name]
  if(family == "gaussian"){
    mod <- stats::lm(YY ~ d_1*d_2, # splines::bs(trial, df = 4)
                      weights = ww,
                      data = data)
  }else{
    data$YY <- as.integer(data$YY) # ensure converted
    
    mod <- stats::glm(YY ~ d_1*d_2*d_3*d_4*d_5, # splines::bs(trial, df = 4)
                     family = family,
                     weights = ww,
                     data = data)
  }
  
  if(family == "gaussian"){
    mod2 <- stats::lm(YY ~ d_sum, # splines::bs(trial, df = 4)
                     weights = ww,
                     data = data)
  }else{
    data$YY <- as.integer(data$YY) # ensure converted
    
    mod2 <- stats::glm(YY ~ d_sum, # splines::bs(trial, df = 4)
                      family = family,
                      weights = ww,
                      data = data)
  }

  return(list(dat = data,
              model1 = mod,
              model2 = mod2,
              regime = regime_ind))
  
}

