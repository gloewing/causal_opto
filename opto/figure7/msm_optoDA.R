# Dose-Response Only (pooled across stim sessions)
# HR-MSM for optoDA
# Gabe Loewinger and Alex Levis 1-25-24

curium <- FALSE
if(curium){
  wd <- "/home/loewingergc/optoDA_preprocess/"
  code_wd <- "/home/loewingergc/msm_hr/"
  save_wd <- "/home/loewingergc/optoDA_results/"

  # learners
  toml.list <- configr::read.config(file = "/lscratch/SpontaneousBehaviour/optoda_intermediate_results/closed_loop_learners.toml")
  learner_ids <- toml.list[[1]]$learners
  
}else{
  wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Data_processed/"
  code_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_MSM/Code/"
  save_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/Final/Dose/"
  # learners
  learner_ids <- read.csv(paste0(wd, "learners.csv"))$mouse_id 
  }

delta <- 5
cov_type <- "HC0" # "HC3" is too computationally intensive
library(data.table)
library(dplyr)
control_ids <- read.csv(paste0(wd, "ctrl_ids.csv"))$mouse_id
dat <- data.table::fread(paste0(wd, "optoDA_data.csv")) %>%
            table.express::filter(!mouse_id %in% control_ids) # treatment group only
source(paste0(code_wd, "MSM_HR_code_IPW_opto_Murphy.R"))
source(paste0(code_wd, "MSM_HR_code_IPW_opto_DT.R"))

setDT(dat)

#####################################
# calculate outcomes with different k
#####################################
# calculate shift on the biggest k to start with so do not have to do this repeatedly below
start <- 2 # 0 means start summing at current trial, 1 starts at next trial and 2 skips the subsequent trial
kVec <- c(1,5,10)
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
setDF(d)

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
# dat$seshID <- 1 # add fake session variable since we are not consider session here
targets <- targets[targets > 0]

# change id name for sandwich package
id_col <- which(colnames(dat) == "mouse_id")
colnames(dat)[id_col] <- "id"

for(cond in c("condit", "uncondit")){ # 
  # "condit is conditional on A_{t-1} == 1 (i.e., Murphy approach)
  for(lrn in c(NA)){
    # only learners = TRUE, non-learns = FALSE, ALL animals = NA
    
    for(kk in c(1,5,10)){ #  
      # iterate through definitions of outcome
      Y_nm <- paste0("Y", kk)
      # remove trials at start or end of sessions (indicated by NAs in outcome -- this is how I constructed Ys)
      d <- dat[ !is.na(dat[Y_nm]), ] # Y_nm is the specific Y I specified above
      d$YY_fail <- as.integer(kk - d[,Y_nm]) # failures for GLM
      
      if(cond == "condit"){
        # full flexible model
        coef_dim <- 4 # number of coefs for full factorial model
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
      coef_dim <- 4 # number of coefs for full factorial model
      seMat2 <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      betaMat2 <- matrix(NA, ncol = coef_dim, nrow = length(targets))
      colnames(seMat2) <- paste0("se2_var_", 1:coef_dim)
      colnames(betaMat2) <- paste0("beta2_var_", 1:coef_dim)
      }
      
      for(tt in 1:length(targets)){
        trgt <- targets[tt]
        set.seed(1)
        
        if(is.na(lrn)){
          # all mice
          rows_incl <- which(d$target == trgt)
        }else if(lrn){
          # include only learners
          rows_incl <- which(d$target == trgt & d$id %in% learner_ids)
        }else if(!lrn){
          # only include non-learner mice
          rows_incl <- which(d$target == trgt & !d$id %in% learner_ids)
        }
        
        
        if(cond == "uncondit"){
          # our approach
          mod <- msm_HR_optoIPW_DT(data = d[rows_incl,],
                                formula = NA,
                                family = "binomial",
                                delta = delta,
                                A_name = "A",
                                Y_name = Y_nm,
                                X_name = "available",
                                W_name = "weights",
                                id_name = "id",
                                remove_ind = c(4, 7:8, 12:16, 20, 23:32), # seq(1,32)[-c(22, 21, 17, 1)], 
                                trial_name = "trial",
                                session_name = "session_repeat",
                                h = NA)
          
          # standard errors
          seMat2[tt,] <- mod$SEs
          betaMat2[tt,] <- mod$final_model$beta_hat
          
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
                                      session_name = "session_repeat",
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
        flNm <- paste0(save_wd, "optoDA_Int_lrn_", lrn, "_Y_",  Y_nm) # filename with specifications 
        res <- cbind(betaMat, betaMat2, seMat, seMat2)
      }  
      
      # write results file
      write.csv(res, file = flNm, row.names = FALSE)
  
    }
  }
}
