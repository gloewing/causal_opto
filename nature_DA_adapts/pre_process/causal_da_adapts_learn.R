library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(refund)
library(dplyr)

# read data and extract variable names
setwd("~/Desktop/NIMH Research/Causal/DA_adapts_learn_rate")
matdat <- R.matlab::readMat("seshMerge.mat")
matdat <- matdat$seshMerge


# indices in dataset of variables
photo_idx <- 7
opt_idx <- 4
session_idx <- 2
trial_idx <- 1
sess_idx <- c(1:6)
iti_idx <- 15

# indices in photometry data of time points (e.g., cue onset/offset)
# cued trials (0 and 2) -- not sure about cued trials==1
cue_onset <- 151
cue_offset <- 200
post_cue_end <- 250 # analyze window from cue_onset:post_cue_end according to Luke
reward_onset <- 301
reward_off <- 500

# animal group ids
control <-  c(1, 2, 4, 6, 9, 11, 15, 16, 19)
slPlus <- c(7, 8, 12, 13, 17) 
slMinus <- c(3, 5, 10, 14, 18, 20) 
sPlusLPlus <- c(21, 22, 23, 24, 25)
n <- 24 # sample size


id_mat <- data.frame( ids = 1:n, group = rep(NA, n) )
id_mat[control,2] <- "control"
id_mat[slPlus,2] <- "slPlus"
id_mat[slMinus,2] <- "slMinus"
id_mat[sPlusLPlus,2] <- "sPlusLPlus"

# session data
nms <- rownames(matdat[,1,])
cue_nms <- nms[sess_idx] # names with same dimension about trial type
dat <- do.call(cbind, lapply(sess_idx, function(x) do.call(c, matdat[x,,])  ) ) # seq_along(cue_nms)
colnames(dat) <- cue_nms
ids <- do.call(c, lapply(seq_along(1:n), function(x) rep(x, length( matdat[1,,][[x]])  ) ))
trial_num <- do.call(c, sapply(seq_along(1:n), function(x)  1:length( ids[ids == x] ) )   ) # trial number (ignores session structure) # assumes in order within animal ID
ids <- data.frame(ids = ids, trial_num = trial_num)

# treatment groups
ids <- left_join(ids, id_mat, by = "ids")

# ITIs (not provided in all sessions)
ITIs <- vector(length = nrow(ids) )
iti_vec <- do.call(c, matdat[iti_idx,,])  # ITIs 
ITIs[1:length(iti_vec)] <- iti_vec # assumes (correctly) ITIs happen for first part of vector and then ends

# photometry
photo <- do.call(rbind, matdat[photo_idx,,] )
colnames(photo) <- paste0("photometry.", 1:ncol(photo) )

# summary of photometry for outcomes
### In cued trials (trialID==0 or 2), cues start at data point 151 in the trials and end at data point 200 (i.e., 0.5 sec cue that ends 1 s before reward delivery)
## "Analysis windows were chosen to capture the extent of mean phasic activations following each kind of stimulus. For NAc–DA and VTA–DA, reward responses were quantified from 0 to 2 s after reward delivery and 
## cue responses were quantified from 0 to 1 s after cue delivery. DS–DA exhibited much faster kinetics, and thus reward and cue responses were quantified from 0 to 0.75 s after delivery."
# 2) Quantify DA (subtract out mean baseline): 
#       -mean baseline DA[1:150]
#       -mean cue DA[151:250] # includes
#       -mean reward DA[301:500]

base_DA <- rowMeans( photo[ , 1:(cue_onset-1)] ) # baseline DA before cue onset
cue_DA <- rowMeans( photo[ , cue_onset:post_cue_end] ) # cue DA 
reward_DA <- rowMeans( photo[ , reward_onset:reward_off] ) # reward DA 

# DA during periods relative to pre-cue baseline
DA_mat <- data.frame( cue_DA =    cue_DA - base_DA,
                      reward_DA = cue_DA - base_DA)

# join data together
dat <- as.data.frame( cbind(ids, ITIs, dat, DA_mat) )
# dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(iti_vec, ids, ITIs, DA_mat, photo, trial_num)




# information from Luke, paper and README

# summary of photometry for outcomes
### In cued trials (trialID==0 or 2), cues start at data point 151 in the trials and end at data point 200 (i.e., 0.5 sec cue that ends 1 s before reward delivery)
## "Analysis windows were chosen to capture the extent of mean phasic activations following each kind of stimulus. For NAc–DA and VTA–DA, reward responses were quantified from 0 to 2 s after reward delivery and 
## cue responses were quantified from 0 to 1 s after cue delivery. DS–DA exhibited much faster kinetics, and thus reward and cue responses were quantified from 0 to 0.75 s after delivery."


# Mouse identity
# Individual mice belonged to different experimental groups, as listed below:
#   Control (no exogenous stimulation): 1, 2, 4, 6, 9, 11, 15, 16, 19
# stimLick+ (DA stim at reward on lick+ trials): 7, 8, 12, 13, 17
# stimLick- (DA stim at reward on lick- trials): 3, 5, 10, 14, 18, 20
# stim+Lick+ (large, uncalibrated DA stim at reward on lick+ trials; as in Fig 6): 21, 22, 23, 24, 25
# 
# Definitions for each field (column) in the “seshMerge” data structure:
#   %first are some trialwise identifiers
# trialID: Type of trial delivered. 0 == cued, 1==uncued, 2==omitted.
# seshID: Which session this trial came from.
# lickState: Whether there was licking during the delay just preceding reward. 0==lick-, 1==lick+
#   stimState: Whether exogenous VTA stimulation was delivered on this trial. 0==no stim, 1==stim
# 
# %Raw data for 7 seconds across each trial: from 3 sec before reward to 4 sec after reward
# *** notes on this data formatting: Sampling rate is 100 Hz, in other words, 100 data points per second. On cued trials, the 0.5 cue starts 1.5 s before reward. Thus when cues happened on a trial they occur during data points 151-200. Reward or reward omission occurs at data point 301 in each trial. Fully preprocessed data that in the case of jRCaMP1b data has been movement and bleaching corrected to produce a dF/F measurement, and then z-scored to enable comparisons between different mice.

## to Do:
# 1) make treatment group indicator variable(s) for Control  (1, 2, 4, 6, 9, 11, 15, 16, 19), stimLick+ (7, 8, 12, 13, 17), stimLick- (3, 5, 10, 14, 18, 20) and stim+Lick+ (21, 22, 23, 24, 25)
# 2) Quantify DA (subtract out mean baseline): 
#       -mean baseline DA[1:150]
#       -mean cue DA[151:250] # includes the .5second cue and 0.5 second after the cue
#       -mean reward DA[301:500]
# 3) ***identify which trials/sessions to use for figure 5 (and which sessions): identify what right estimand is
# 4) ***confirm with Luke that "cue responses were quantified from 0 to 1 s after cue delivery " means "0 seconds post cue onset to 1 second post cue onset" (i.e., from time points 151:250)
