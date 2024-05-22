# dataset preprocess of Spontaneous DA on Curium for Both treatment and control animals

## ***sign on***
# ssh loewingergc@curium.nimh.nih.gov

## ***load R***
# R

# whether to discard data with timestamp switches
# discard_data <- TRUE
# sess_length <- 1800 # min length for sessions
# min_sess_length <- 10 # make 10 but could probably be as low as 4 sec because most seem to be 3 sec 
min_trial_frames <- 3 # 3 minimum consecutive frames
Hz <- 1 / 30
post_threshold <- FALSE # if true then cut off based on <min_trial_frames>
# **********post_threshold = FALSE yielded estimates identical to what is in the   learning_timecourse_binsize-30.parquet and TRUE was NOT consistent (my estimates end up being too low as expected since we don't count the short ones) despite what is written in paper **********
# read dataset
library(arrow)
library(configr)
library(data.table)
library(table.express)
source("/home/loewingergc/optoDA_preprocess_new/data_preprocess.R")
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")
# dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/learning_timecourse_binsize-30.parquet")

# learners
toml.list <- configr::read.config(file = "/lscratch/SpontaneousBehaviour/optoda_intermediate_results/closed_loop_learners.toml")
learner_ids <- toml.list[[1]]$learners

sess_rm <- c("Habituation 2", "habituation 3", 
             "Habituation 3", "habituation 4", 
             "Habituation 4",
             "Habituation 1", "habituation 1")

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

# ---------------------------------------------------------------------
###################
# filter variables
###################
# experiment_type == "reinforcement", session 1 & 2, and feedback_status > 0, supported by this file below:
# https://github.com/dattalab/dopamine-reinforces-spontaneous-behavior/blob/14ed714c4f49be1f4728a080eb65c4ee42f74c23/notebooks_panels/fig_3j%2C%203k-transition%20changes.ipynb#L121
# subset to specific experiment type
setDT(dat)
sess_opt <- c(1, 2) # equivalently gives sessions 1-2:  unique(dat$session_number[dat$feedback_status == 1]) # sessions with stimulation (feedback_status == 1 is stim)
dat <- dat %>% 
  table.express::filter(experiment_type %in% c("reinforcement", "reinforcement_photometry"), # reinforcement seems to be the only experiment type with all 9 "learner" mice
                        session_number %in% sess_opt, # only sessions where there is sometimes stimulation
                        # area == "snc (axon)", # "snc (axon)" -- opto animals & "ctrl" -- control animals
                        !dat$SessionName %in% sess_rm ) %>% # remove rows that can't be used
  table.express::mutate(learner = 1*I(mouse_id %in% learner_ids), # label who is a learner
                        available = 1*I(target_syllable == predicted_syllable & feedback_status >= 0)) %>%
  table.express::mutate(A = 1*I( available == 1 & feedback_status == 1)) # treatment/stim

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

# **************************************************
idx_final <- idx # save final file to add below
# **************************************************

# write data
if(post_threshold){
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_data_thresh.csv") # write file
}else{
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_data.csv") # write file
}  
rm(dat, d, idx, A_idx, avail_idx)
