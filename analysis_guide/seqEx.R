seqEx <- function(data,
                   A_name,
                   Y_name,
                   X_name,
                   W_name,
                   id_name,
                   trial_name,
                   seq.mat,
                   session_name = NA,
                   h = NA 
){
  
  library(data.table)
  delta <- ncol(seq.mat)
  ids <- sort(unique(data[,id_name]))
  Wt <- data[,W_name]
  data$A <- data[,A_name]
  data$X <- data[,X_name]
  data$trial <- data[,trial_name]
  W <- rep(NA, nrow(data))
  A <- data[,A_name] # save original treatment and time-varying confounder A, X
  X <- data[,X_name] # save original treatment and time-varying confounder A, X
  data$I0 <- as.integer(1 * I(A == 0))
  data$I1 <- as.integer(1 * I(A == X))
  n <- nrow(data)
  data$weight <- rep(NA, n)
  if(is.na(session_name)){
    data$session <- 1 # treat it as only session
  }else{
    data$session <- data[,session_name]
  }
  
  # numerator/stabilizing weights
  if(is.na(h))   data$h <- 1
  data <- data[order(data$id, data$session, data$trial),] # sort rows
  
  ###################################################################################
  # calculate denominator weights as product of moving window of delta time-points
  ###################################################################################
  for(i in 1:length(ids)){
    sess <- unique( data$session[ data[id_name] == ids[i] ] )
    for(ss in sess){
      idx <- which(data[id_name] == ids[i] & data$session == ss)
      WWtt <- Wt[idx]
      W[idx] <- c( rep(Inf, delta - 1), # use Inf so it will make weights 0 when invert
                   sapply( 1:(length(idx)-delta+1),  
                           function(dd) prod(WWtt[dd:(dd+delta-1)]) ) )
    }
  }
  
  # # numerator/stabilizing weights
  data$weight <- data$h / W
  
  # indicator matrix for treatment regime 
  message("Formatting regimes")
  indicator_matrix <- seq.mat
  rm(seq.mat)
  colnames(indicator_matrix) <- paste0("d_", 1:delta)
  q <- nrow(indicator_matrix)
  data <- data[rep(seq_len(nrow(data)), each = q),] # copy
  data <- cbind(data, indicator_matrix[rep(1:q, n),]) # copies each row 2^delta
  indicator_matrix <- as.matrix(indicator_matrix) # use below as a matrix
  if(delta == 1)   colnames(data)[ncol(data)] <- "d_1" # colnames of indicator matrix do not work when delta = 1
  
  ############################
  ## construct regimes
  ############################
  data.table::setDT(data)
  data[,keyby = .(id, session), regime_number := rep(1:q, .N / q)] # create regime number column
  mycol <- paste0("r", 1:delta)
  data[, keyby = .(id, session, regime_number), c(mycol) := (data.table::shift(I1, n = (delta-1):0, type = "lag", fill = 0))] 
  data.table::setkey(data, regime_number)
  data[order(id, session, trial, regime_number)]
  for(jj in 1:q){
    idx <- which(data$regime_number == jj)
    data.table::set(data, idx, "regime_p1", (1*as.matrix(data[idx,..mycol])) %*% indicator_matrix[jj,])
  } 
  
  # I0
  indicator_matrix <- 1 - indicator_matrix # switch which columns we choose (original 0s correspond to the I0 columns)
  data[, keyby = .(id, session, regime_number), c(mycol) := (data.table::shift(I0, n = (delta-1):0, type = "lag", fill = 0))] 
  for(jj in 1:q){
    idx <- which(data$regime_number == jj)
    data.table::set(data, idx, "regime_p0", (1*as.matrix(data[idx,..mycol])) %*% indicator_matrix[jj,])
  } 

  # define variables 
  data[, c(W_name) := NULL ] # remove initial weight variable to avoid issues
  data[, c(mycol) := NULL] # remove regime numbers
  data[, regime_p := regime_p0 + regime_p1]
  data[, w := 1 * (regime_p == delta)]
  regime_ind <- data$w 
  data[, w := w * weight]
  data <- data[ w != 0 ]

  # remove columns
  remove.nm <- c("I0", "I1", "regime_number", "regime_p1", "regime_p0", "regime_p", 
                 "weight", "h") # "indicator_matrix[rep(1:q, n), ]", 
  data[, (remove.nm) := NULL]
  
  return(data)
}

