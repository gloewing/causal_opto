dtr_sim <- function(n,
                    timepoints = 10,
                    sigmaX_base = c(1,1), # sd of X0 and X1
                    sigmaX = 1,
                    A_prob = 1/2, # IF WE ADJUST REQUIRES RE-DOING ALGEBRA AND RECODING
                    beta_div = 1 # scale all elements of treatment-related component of beta
                    ){
    
  
  # draw V - baseline covariates
  #V <- runif(n, min = 0, max = 1)
  
  # draw X0 and X1 (assumes delta = 2), and A0
  # X1 <- rbinom(n, size = 1, prob = 0.5) #rnorm(n, mean = V, sd = sigmaX_base[1] )
  # X2 <- #rnorm(n, mean = V, sd = sigmaX_base[2] )
  # A1 <- rbinom(n, size = 1, prob = I(X1 > 0) * A_prob) # A0 and A1 (assumes delta = 2)

  # draw beta
  beta <- c(0, 0.25, 2, 0.175, 0.5) / beta_div # make time-points further back have smaller beta magnitudes
  # X_j A_j X_{j-1} A_{j-1}
  # c(-1, 0.25, 1, 0.175, 0.5)
  # beta2 <- 2
    
  dat_lst <- list(length = n)
  for(i in 1:n){
    dat <- data.frame(matrix(NA, nrow = timepoints, ncol = 6) )
    colnames(dat) <- c("id", "trial", "X", "A", "Y", "weights")
    X1 <- rbinom(1, size = 1, prob = 0.5)
    Aj <- rbinom(1, size = 1, prob = X1 * A_prob)
    w <- ifelse(X1 == 0, 1, Aj*A_prob + (1-Aj)*(1-A_prob) ) # if X_j == 0, no laser (then weights=1) otherwise use bernoulli likelihood
    dat[1,] <- c(i, 1, X1, Aj, NA, w) # no model for outcome at time-point 1 so leave Y = NA
    
    for(j in 2:timepoints){
      dat[j, "X"] <- Xj <- rbinom(1, size = 1, prob = 0.3*(1 - dat$A[j-1]) + 0.7*dat$A[j-1] + 0.1  )  # addition of 0.1 ensures E[X_j] = 0.5
      dat[j, "A"] <- Aj <- rbinom(1, size = 1, prob = Xj * A_prob)
      # draw outcome using linear model
      Y_j <- rnorm(1, 
                  mean = beta[1] + as.numeric(t(dat[j:(j-1), c("X","A")])) %*% beta[-1], # + beta2 * j  # additive model: 1st component is laser-related outcome model and  2nd component is time-related (laser independent)
                  sd = sigmaX)
      w <- ifelse(Xj == 0, 1, Aj*A_prob + (1-Aj)*(1-A_prob) ) # if X_j == 0, no laser (then weights=1) otherwise use bernoulli likelihood
      dat[j,] <- c(i, j, Xj, Aj, Y_j, w)
    }
    dat_lst[[i]] <-  dat
  }
  
  return(list(data = do.call(rbind, dat_lst),
              beta = beta,
              A_prob = A_prob)
          )
}
