msm_HR <- function(data,
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
  
  library(data.table)
  #######################################################################
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
    data$session <- 1 # only sessio n
  }else{
    data$session <- data[,session_name]
  }
  
  # numerator/stabilizing weights
  if(is.na(h))   data$h <- 1
  
  data <- data[order(data$id, data$session, data$trial),] # sort rows
  
  ###################################################################################
  # calculate denominator weights as product of moving window of delta time-points
  ###################################################################################
  message("Calculate weights")
  # setDT(data)
  # data[, keyby = .(id, session), W := sapply(shift(Wt, n = (delta-1):0, type = "lag", fill = Inf), prod)] #[,regime_prod1 := Reduce(`*`, .SD), .SDcols = idx_regI1] #[,regime_prod1]
  
  for(i in 1:length(ids)){
    # idx <- which(data[id_name] == ids[i])
    sess <- unique( data$session[ data[id_name] == ids[i] ] )
    for(ss in sess){
      idx <- which(data[id_name] == ids[i] & data$session == ss)
      WWtt <- Wt[idx]
      W[idx] <- c( rep(Inf, delta - 1), # use Inf so it will make weights 0 when invert
                   sapply( 1:(length(idx)-delta+1),  function(dd) prod(WWtt[dd:(dd+delta-1)])  )
      )
    }
  }
  
  # 
  # # data[, keyby = .(id, session), W := prod(shift(Wt, n = (delta-1):0, type = "lag", fill = Inf))] #[,regime_prod1 := Reduce(`*`, .SD), .SDcols = idx_regI1] #[,regime_prod1]
  # 
  # # numerator/stabilizing weights
  data$weight <- data$h / W
  
  # indicator matrix for treatment regime 
  message("Formatting regimes")
  indicator_matrix <- expand.grid( rep(list(0:1), delta))
  colnames(indicator_matrix) <- paste0("d_", 1:delta)
  q <- nrow(indicator_matrix)
  data <- data[rep(seq_len(nrow(data)), each = q),] # copy
  data <- cbind(data, indicator_matrix[rep(1:q, n),]) # copies each row 2^delta
  indicator_matrix <- as.matrix(indicator_matrix) # use below as a matrix
  ################
  # also works byt seems slower and more memory intesive
  # setDT(data)
  # setDT(ind_matrix)
  # data <- data[rep(seq_len(nrow(data)), each = 2^delta)]
  # ind_matrix <- ind_matrix[rep(1:q, times = nrow(data))]
  # data <- cbind(data, ind_matrix)
  
  ############################
  ## data.table approach
  ############################
  setDT(data)
  data[,keyby = .(id, session), regime_number := rep(1:q, .N / q)] # create regime number column
  mycol <- paste0("r", 1:delta)
  
  data[, keyby = .(id, session, regime_number), c(mycol) := (shift(I1, n = (delta-1):0, type = "lag", fill = 0))] #[,regime_prod1 := Reduce(`*`, .SD), .SDcols = idx_regI1] #[,regime_prod1]
  
  setkey(data, regime_number)
  data[order(id, session, trial, regime_number)]
  for(jj in 1:q){
    idx <- which(data$regime_number == jj)
    set(data, idx, "regime_p1", (1*as.matrix(data[idx,..mycol])) %*% indicator_matrix[jj,])
  } 
  
  # for I0
  indicator_matrix <- 1 - indicator_matrix # switch which columns we choose (original 0s correspond to the I0 columns)
  data[, keyby = .(id, session, regime_number), c(mycol) := (shift(I0, n = (delta-1):0, type = "lag", fill = 0))] #[,regime_prod1 := Reduce(`*`, .SD), .SDcols = idx_regI1] #[,regime_prod1]
  for(jj in 1:q){
    idx <- which(data$regime_number == jj)
    set(data, idx, "regime_p0", (1*as.matrix(data[idx,..mycol])) %*% indicator_matrix[jj,])
  } 
  
  # delete columns
  data[, c(mycol) := NULL]
  
  # define variables inside dataset
  message(paste0("delta", delta))
  data[, regime_p := regime_p0 + regime_p1]
  data[, ww := 1 * (regime_p == delta)]
  regime_ind <- data$ww 
  data[, ww := ww * weight]
  data$YY <- data[,..Y_name]
  data <- data[ww!=0]
  
  ####################
  # MSM
  ####################
  message("Fitting MSM")
  
  # MSM is a glm
  if(family == "gaussian"){
    mod <- stats::lm(YY ~ d_1*d_2, # splines::bs(trial, df = 4)
                     weights = ww,
                     model = FALSE,
                     x = TRUE,
                     y = TRUE,
                     data = data)
  }else{
    data$YY <- as.integer(data$YY) # ensure converted
    
    mod <- stats::glm(YY ~ d_1*d_2, # splines::bs(trial, df = 4)
                      family = family,
                      weights = ww,
                      data = data)
  }
 
  return(list(dat = data,
              model = mod))

  
}

