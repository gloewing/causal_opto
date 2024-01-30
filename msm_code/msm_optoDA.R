# HR-MSM for optoDA
# Gabe Loewinger and Alex Levis 1-25-24

curium <- TRUE
if(curium){
  wd <- "/home/loewingergc/optoDA_preprocess/"
  code_wd <- "/home/loewingergc/msm_hr/"
  save_wd <- "/home/loewingergc/optoDA_results/"

  # learners
  toml.list <- configr::read.config(file = "/lscratch/SpontaneousBehaviour/optoda_intermediate_results/closed_loop_learners.toml")
  learner_ids <- toml.list[[1]]$learners
  
}else{
  wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/preprocessed_data/"
  code_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/msm_hr/"
  save_wd <- "/Users/loewingergc/Desktop/Research/optoDA_results/"
  # learners
  toml.list <- configr::read.config(file = "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoda_intermediate_results/closed_loop_learners.toml")
  learner_ids <- toml.list[[1]]$learners
  }

Y_nm <- "Y5" # use k = 5 for outcome
delta <- 5
cov_type <- "HC0" # "HC3" is too computationally intensive
library(data.table)
dat <- data.table::fread(paste0(wd, "optoDA.csv"))
source(paste0(code_wd, "MSM_HR_code_IPW_opto_d5.R"))
source(paste0(code_wd, "MSM_HR_code_IPW_opto_Murphy.R"))
# define outcomes based on k ahead

#####################################
# calculate outcomes with different k
#####################################
# calculate shift on the biggest k to start with so do not have to do this repeatedly below
start <- 2 # 0 means start summing at current trial, 1 starts at next trial and 2 skips the subsequent trial
kVec <- c(2,5,10)
k <- max(kVec)
mycol <- paste0("V", 1:k)
dat[, c(mycol) := (shift(available, n = start:(start+k), type = "lead", fill = NA)), keyby = .(mouse_id, session) ]

# iterate through and save for each
Ymat <- matrix(NA, nrow = nrow(dat), ncol = length(kVec))
colnames(Ymat) <- paste0("Y", kVec)
d <- dat
for(kk in 1:length(kVec)){
  k <- kVec[kk]
  mycol <- paste0("V", 1:k)
  Ymat[,kk] <- apply(dat[,..mycol], 1, max) # calculate function of variables
}

dat <- as.data.frame(cbind(d, Ymat))
# remove columns
dat <- dat[,!(colnames(dat) %in% mycol)]

#############################################
# calculate IPW weights for current trial
#############################################
treat_prob <- 0.75 # P(A | X = 1) # experimental treatment probability
dat$weights <- 1 # initialize for all available = 0 trials
dat$weights <- ifelse(dat$available == 1 & dat$A == 1, treat_prob, dat$weights) # if treated and available
dat$weights <- ifelse(dat$available == 1 & dat$A == 0, 1 - treat_prob, dat$weights) # if not treated but available

######################################
# fit models for different targets
######################################
targets <- unique(dat$target)
dat$seshID <- 1 # add fake session variable since we are not consider session here

# change id name for sandwich package
id_col <- which(colnames(dat) == "mouse_id")
colnames(dat)[id_col] <- "id"

for(cond in c("condit", "uncondit")){
  # "condit is conditional on A_{t-1} == 1 (i.e., Murphy approach)
  for(lrn in c(TRUE, FALSE, NA)){
    # only learners = TRUE, ALL animals = FALSE
    
    for(Y_nm in c("Y2", "Y5", "Y10")){
      # iterate through definitions of outcome
      
      # remove trials at start or end of sessions (indicated by NAs in outcome -- this is how I constructed Ys)
      d <- dat[ !is.na(dat[Y_nm]), ] # Y_nm is the specific Y I specified above
      
      if(cond == "condit"){
        # full flexible model
        coef_dim <- 2 # number of coefs for full factorial model
        seMat <- matrix(NA, ncol = coef_dim, nrow = length(targets))
        betaMat <- matrix(NA, ncol = coef_dim, nrow = length(targets))
        colnames(seMat) <- paste0("se_var_", 1:coef_dim)
        colnames(betaMat) <- paste0("beta_var_", 1:coef_dim)
      }else if (cond == "uncondit"){
      # full flexible model
      coef_dim <- 32 # number of coefs for full factorial model
      seMat <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      betaMat <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      colnames(seMat) <- paste0("se_var_", 1:coef_dim)
      colnames(betaMat) <- paste0("beta_var_", 1:coef_dim)
      
      # sum models
      coef_dim <- 6 # number of coefs for full factorial model
      seMat2 <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      betaMat2 <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      colnames(seMat2) <- paste0("se2_var_", 1:coef_dim)
      colnames(betaMat2) <- paste0("beta2_var_", 1:coef_dim)
      }
      
      for(tt in 1:length(targets)){
        trgt <- targets[tt]
        set.seed(1)
        
        if(is.na(lrn)){
          # only include non-learner mice
          rows_incl <- which(d$target == trgt & !d$id %in% learner_ids)
        }else if(lrn){
          # include only learners
          rows_incl <- which(d$target == trgt & d$id %in% learner_ids)
        }else{
          # all mice
          rows_incl <- which(d$target == trgt)
        }
        
        
        if(cond == "uncondit"){
          # our approach
          mod <- msm_HR_optoIPW(data = d[rows_incl,],
                                formula = NA,
                                family = "binomial",
                                delta = delta,
                                A_name = "A",
                                Y_name = Y_nm,
                                X_name = "available",
                                W_name = "weights",
                                id_name = "id",
                                trial_name = "trial",
                                session_name = "seshID",
                                h = NA)
          
          # standard errors
          seMat[tt,] <- as.numeric(sqrt(diag(sandwich::vcovCL(mod$model1, type = cov_type, cluster = ~ id))))
          seMat2[tt,] <- as.numeric(sqrt(diag(sandwich::vcovCL(mod$model2, type = cov_type, cluster = ~ id))))
          betaMat[tt,] <- as.numeric(coef(mod$model1))
          betaMat2[tt,] <- as.numeric(coef(mod$model2))
          
        }else if(cond == "condit"){
          # Murphy approach: conditional on A_{t-1} == 1
          mod <- msm_HR_optoIPW_Murphy(data = d[rows_incl,],
                                      formula = NA,
                                      family = "binomial",
                                      delta = 1,
                                      A_name = "A",
                                      Y_name = Y_nm,
                                      X_name = "available",
                                      W_name = "weights",
                                      id_name = "id",
                                      trial_name = "trial",
                                      session_name = "seshID",
                                      h = NA)
          
          # standard errors
          seMat[tt,] <- as.numeric(sqrt(diag(sandwich::vcovCL(mod$model1, type = cov_type, cluster = ~ id))))
          betaMat[tt,] <- as.numeric(coef(mod$model1))
        }
        
      }
      
      # save results for learner/all and outcome (Y_nm)
      if(cond == "condit"){
        flNm <- paste0(save_wd, "optoDA_lrn_", lrn, "_Y_",  Y_nm, "_conditional") # filename with specifications with Murphy approach
        res <- cbind(betaMat, seMat) # only one model
      }else if(cond == "uncondit"){
        flNm <- paste0(save_wd, "optoDA_lrn_", lrn, "_Y_",  Y_nm) # filename with specifications 
        res <- cbind(betaMat, betaMat2, seMat, seMat2)
      }  
      
      # write results file
      write.csv(res, file = flNm, row.names = FALSE)
  
    }
  }
}







# mean(dat$Y10,na.rm=TRUE)
# mean(dat$Y5,na.rm=TRUE)
# mean(dat$Y2,na.rm=TRUE)
# 
# dat %>% as_tibble() %>%
#   group_by(mouse_id, session) %>%
#   summarise(trials = n())
# 
# dat[,
#    lapply(.SD, max),
#    .SDcols = ..mycol]
# 
# dat[, myMax := dat[,..mycol] ]
# dat[, myMax := lapply(..mycol, max), .SDcols = ..c(mycol)]

# dat[, Y2 := Reduce(`+`, shift(pose_length, n = start:(start+k), type = "lead", fill = 0)), keyby = .(mouse_id, session) ]

#as.matrix(dat[, shift(pose_length, n = start:(start+k), type = "lead", fill = 0), keyby = .(mouse_id, session) ][,-c("mouse_id", "session")])
#Rfast::rowMaxs()


# 
# # k = 5
# k <- 5
# dat[, Y5 := Reduce(`+`, shift(pose_length, n = start:(start+k), type = "lead")), keyby = .(mouse_id, session) ]
# 
