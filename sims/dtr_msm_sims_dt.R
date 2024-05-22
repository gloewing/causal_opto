# DTR MSM simulations
# 11/24/23

library(sandwich)
 
cluserInd <- TRUE
if(cluserInd){

  # paths
  wd <- "insert"
  data_path <-"insert" # path to original Matlab files for data
  save.folder <- "insert" # msm_sims_dt4" #
  
  args <- commandArgs(TRUE)
  seed <- as.integer( Sys.getenv('SLURM_ARRAY_TASK_ID') ) # seed index from array id
  n <- as.integer( as.numeric(args[1]) ) # sample size
  time_points <- as.integer( as.numeric(args[2]) ) # time points per animal
  
}else{
  wd <- "insert"
  n <- 100
  time_points <- 500
  seed <- 7
}

# load code
source(paste0(wd,"dtr_sim.R"))
source(paste0(wd,"MSM_HR_code_DT.R"))
source(paste0(wd,"sandwich_normal.R"))
source(paste0(wd,"saveFn.R"))

# simulation parameters
boots <- 5000
coef_dim <- 4 # dimension of true parameter vector
var_array  <- array(NA, dim = c(4, coef_dim, coef_dim) )
numSims <- 1000 # simulation replicates

# wait time for cluster
wait_min <- 10

# file name
fileNm <- paste0("dtr_msm",
                 "_n_", n,
                 "_tm_", time_points,
                 "_totSims_", numSims,
                 "_boots_", boots)

##################################
# simulate data
##################################
set.seed(seed)

dtr_data <- dtr_sim(n = n,
                    timepoints = time_points,
                    sigmaX_base = c(1,1), # sd of X0 and X1
                    sigmaX = 1,
                    beta_div = 1,
                    A_prob = 1/2)

# calculate true parameters
theta_truth <- rep(NA, coef_dim)
theta_truth[1] <- dtr_data$beta[1] + 0.5 * (dtr_data$beta[4]) + 0.4 * dtr_data$beta[2]
theta_truth[2] <- 0.5 * dtr_data$beta[5] + 0.2 * dtr_data$beta[2]
theta_truth[3] <- 0.4 * dtr_data$beta[3]
theta_truth[4] <- 0.2 * dtr_data$beta[3]

##################################
# fit model
##################################
set.seed(seed)
dtr_data$data$session <- dtr_data$data$seshID <- 1 # add fake session variable

mod <- msm_HR(data = dtr_data$data,
              formula = NA,
              family = "gaussian",
              delta = 2,
              A_name = "A",
              Y_name = "Y",
              X_name = "X",
              W_name = "weights",
              id_name = "id",
              trial_name = "trial",
              session_name = "session",
              h = NA)

#########################
# variance estimates
#########################
num_CIs <- 14
seMat <- matrix(NA, ncol = coef_dim, nrow = num_CIs)
resMat <- matrix(NA, ncol = (2+num_CIs) * coef_dim, nrow = numSims)
sMat <- powerMat <- matrix(NA, ncol = num_CIs * coef_dim, nrow = numSims)
colnames(resMat) <- c(paste0(paste0("res_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs)), 
                      paste0("bias_",1:coef_dim), 
                      paste0("est_",1:coef_dim))
colnames(powerMat) <- paste0(paste0("power_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs))
colnames(sMat) <- paste0(paste0("SE_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs))

# HC0-HC3
HC_vec <- paste0("HC", c(0,2,3))
for(j in 1:3)  var_array[j,,] <- sandwich::vcovCL(mod$model, type = HC_vec[j], cluster = ~ id)
var_array[4,,] <- sandwich_normal(data = dtr_data$data, model = mod)$cov # our estimator

seMat[1,] <- sqrt(diag(var_array[1,,])) # HC0
seMat[2,] <- sqrt(diag(var_array[4,,])) # our estimator
seMat[3,] <- sqrt(diag(var_array[2,,])) # "HC2"
seMat[4,] <- sqrt(diag(var_array[3,,])) # "HC3"

# Clustered Panel
seMat[5,] <- sqrt(diag(sandwich::vcovPL(mod$model, cluster = ~ id, order.by = ~ trial, lag = "max", adjust = FALSE ) ))
seMat[6,] <- sqrt(diag(sandwich::vcovPL(mod$model, cluster = ~ id, order.by = ~ trial, lag = "NW1994" ) ))

# bootstrap 1-way
b_methods <- c("wild-rademacher", "wild-mammen", "wild-webb", "wild-norm") # "residual" -- is not working
for(meth in 1:length(b_methods)){
  set.seed(seed)
  seMat[6+meth,] <- sqrt(diag(sandwich::vcovBS(mod$model, 
                                               cluster = ~ id, 
                                               type = b_methods[meth],
                                               R = boots) ))
}

# bootstrap 2-way
b_methods <- c("wild-rademacher", "wild-mammen", "wild-webb", "wild-norm") # "residual" -- is not working
for(meth in 1:length(b_methods)){
  set.seed(seed)
  seMat[10+meth,] <- sqrt(diag(sandwich::vcovBS(mod$model, 
                                               cluster = ~ id + trial, 
                                               type = b_methods[meth],
                                               R = boots) ))
}

sMat[seed,] <- as.vector(seMat) # save SEs in matrix

#########################
# calculate summary stats
#########################
cnt <- 1
beta <- coef(mod$model) # estimate of coefficients
# iterate through covariance estimators
for(varNum in 1:num_CIs){
  # iterate through coefficient dimensions
  for(coefDim in 1:coef_dim){
    se_est <- seMat[varNum, coefDim] # estimate for coefficient #<coefDim> and variance estimate #<varNum>
    upper <- beta[coefDim] + 1.96 * se_est
    lower <- beta[coefDim] - 1.96 * se_est
    resMat[seed, cnt] <- 1 * (theta_truth[coefDim] <= upper & theta_truth[coefDim] >= lower) # coverage
    powerMat[seed, cnt] <- 1 * (0 < upper & 0 < lower) + 1 * (0 > upper & 0 > lower) # power
    cnt <- cnt + 1 # update column counter
  }
}

# bias and raw coefficient estimates
colIdx <- (num_CIs * coef_dim + 1):( (1+num_CIs) * coef_dim) # indices for bias
colIdx2 <- (max(colIdx)+1):ncol(resMat) # indices for coefficient estimates
resMat[seed, colIdx] <- (coef(mod$model) - theta_truth) / theta_truth # % bias
resMat[seed, colIdx2] <- coef(mod$model) # estimates

resMat1 <- cbind(resMat, powerMat, sMat)
rm(resMat, powerMat, sMat)

#############################################################################
##############################################
# specific contrasts / excursion effects
##############################################
c_mat <- rbind(c(0,0,1,0), 
               c(0,-1,1,0),
               c(1,0,0,0), # dose 0
               c(1,0.5,0.5,0), # dose 1
               rep(1,4) ) # dose 2

theta_truth <- c_mat %*% theta_truth # update based on contrasts of true values
beta <- as.numeric( c_mat %*% coef(mod$model) ) # estimate of coefficients

num_CIs <- 4 # new number of CIs to construct for just ours and HCs (not doing bootstrap)
coef_dim <- nrow(c_mat) # number of new constrasts (no longer just number of coefficients)
seMat <- matrix(NA, ncol = coef_dim, nrow = num_CIs)
resMat <- matrix(NA, ncol = (2+num_CIs) * coef_dim, nrow = numSims)
sMat <- powerMat <- matrix(NA, ncol = num_CIs * coef_dim, nrow = numSims)
colnames(resMat) <- c(paste0(paste0("res_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs)), 
                      paste0("bias_",1:coef_dim), 
                      paste0("est_",1:coef_dim))
colnames(powerMat) <- paste0(paste0("power_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs))
colnames(sMat) <- paste0(paste0("SE_var_", rep(1:num_CIs, each = coef_dim), "_coef_"), rep(1:coef_dim, num_CIs))

# SEs - # our estimator is j = 4
for(j in 1:4)   seMat[j,] <- sqrt(diag(c_mat %*% var_array[j,,] %*% t(c_mat))) 

cnt <- 1
# iterate through covariance estimators
for(varNum in 1:num_CIs){
  # iterate through coefficient dimensions
  for(coefDim in 1:coef_dim){
    se_est <- seMat[varNum, coefDim] # estimate for coefficient #<coefDim> and variance estimate #<varNum>
    upper <- beta[coefDim] + 1.96 * se_est
    lower <- beta[coefDim] - 1.96 * se_est
    resMat[seed, cnt] <- 1 * (theta_truth[coefDim] <= upper & theta_truth[coefDim] >= lower) # coverage
    powerMat[seed, cnt] <- 1 * (0 < upper & 0 < lower) + 1 * (0 > upper & 0 > lower) # power
    cnt <- cnt + 1 # update column counter
  }
}

# bias and raw coefficient estimates
colIdx <- (num_CIs * coef_dim + 1):( (1+num_CIs) * coef_dim) # indices for bias
colIdx2 <- (max(colIdx)+1):ncol(resMat) # indices for coefficient estimates
resMat[seed, colIdx] <- (beta - theta_truth) / theta_truth # % bias
resMat[seed, colIdx2] <- beta # estimates

resMat2 <- cbind(resMat, powerMat, sMat)
rm(resMat, powerMat, sMat)

#############################################################################
# wait time for cluster
wait_time <- round(runif(1, wait_min * 60, (wait_min + 3) * 60)) # random draw to prevent saving issues
sys_sleep(wait_time, "s") # wait time in seconds

########################
# save results
########################
# print("setWD to save file")

print(paste("Save 1:", seed))

saveFn_Indiv(file = resMat1, 
             fileNm = fileNm, 
             iterNum = seed, 
             iters = numSims,
             save.folder = save.folder)

print(paste("Save 2:", seed))
# contrasts
#saveFn_new
saveFn_Indiv(file = resMat2, 
             fileNm = paste0(fileNm, "_c"), 
             iterNum = seed, 
             iters = numSims,
             save.folder = save.folder)

# count to see if ready
print(paste("count:", seed))
#file_cnt_total <- file_cntr(fileNm = fileNm, save.folder = save.folder)

# see if all completed
# if(file_cnt_total == numSims){
#   process_files(fileNm = fileNm, save.folder = save.folder, colnm = NULL, itrs = numSims)
# }

print(paste("Save 3:", seed))
# also save aggregate
# saveFn(file = resMat1,
#        fileNm = paste0(fileNm, "_agg"),
#        iterNum = seed,
#        save.folder = save.folder)
# 
# print(paste("Save 4:", seed))
# saveFn(file = resMat2,
#        fileNm = paste0(fileNm, "_c_agg2"),
#        iterNum = seed,
#        save.folder = save.folder)