# dataset preprocess of Spontaneous DA on Curium

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
                            area == "snc (axon)", # "snc (axon)" -- opto animals & "ctrl" -- control animals
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
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_thresh.csv") # write file
}else{
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA.csv") # write file
}  
rm(dat, d, idx, A_idx, avail_idx)

# ---------------------------------------------------------------------
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
# experiment_type == "reinforcement", session 1 & 2, and feedback_status > 0, supported by this file below:
# https://github.com/dattalab/dopamine-reinforces-spontaneous-behavior/blob/14ed714c4f49be1f4728a080eb65c4ee42f74c23/notebooks_panels/fig_3j%2C%203k-transition%20changes.ipynb#L121
# subset to specific experiment type
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- c(1, 2) # equivalently gives sessions 1-2:  unique(dat$session_number[dat$feedback_status == 1]) # sessions with stimulation (feedback_status == 1 is stim)
dat <- dat[dat$session_number %in% sess_opt,] # only sessions where there is sometimes stimulation

# filter
dat <- dat[dat$area == "snc (axon)",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# label who is a learner
dat$learner <- 1*I(dat$mouse_id %in% learner_ids)

# availability indicator
dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable &
                         dat$feedback_status >= 0) # make all target poses where feedback_status = -1 into equivalent of wrong pose

# treatment (stim)
dat$A <- 1*I( dat$available == 1 & dat$feedback_status == 1) # 

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

# **************************************************
idx_final <- idx # save final file to add below
# **************************************************

# write data
if(post_threshold){
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_thresh.csv") # write file
}else{
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA.csv") # write file
}  

rm(dat, d, idx, A_idx, avail_idx)




# **************************************************
# find counts during opto session
# **************************************************
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")
# dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/learning_timecourse_binsize-30.parquet")

# learners
toml.list <- configr::read.config(file = "/lscratch/SpontaneousBehaviour/optoda_intermediate_results/closed_loop_learners.toml")
learner_ids <- toml.list[[1]]$learners

###################
# filter variables
###################
# experiment_type == "reinforcement", session 1 & 2, and feedback_status > 0, supported by this file below:
# https://github.com/dattalab/dopamine-reinforces-spontaneous-behavior/blob/14ed714c4f49be1f4728a080eb65c4ee42f74c23/notebooks_panels/fig_3j%2C%203k-transition%20changes.ipynb#L121
# subset to specific experiment type
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- c(1, 2) # equivalently gives sessions 1-2:  unique(dat$session_number[dat$feedback_status == 1]) # sessions with stimulation (feedback_status == 1 is stim)
dat <- dat[dat$session_number %in% sess_opt,] # only sessions where there is sometimes stimulation

# filter
dat <- dat[dat$area == "snc (axon)",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# label who is a learner
dat$learner <- 1*I(dat$mouse_id %in% learner_ids)

# availability indicator
dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable &
                        dat$feedback_status >= 0) # make all target poses where feedback_status = -1 into equivalent of wrong pose

# treatment (stim)
dat$A <- 1*I( dat$available == 1 & dat$feedback_status == 1) # 
dat$XX <- 1*I( dat$target_syllable == dat$predicted_syllable) # correct pose WITHOUT feedback_status threshold

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

idx2 <- idx %>% 
  table.express::group_by(mouse_id, target, session_number) %>%
  table.express::summarise(#total_counts = sum(available),
                           total_counts = sum(XX)) %>% # calculate proportions of being in right pose (available is indicator that met threshold of confidence)
  tidyr::pivot_wider(names_from = session_number, values_from = total_counts ) # pivot wider to make columns

colnames(idx2)[3:4] <- paste0("Sess_", c(1,2))
setDF(idx2)

idx2$optoCount <- rowMeans(idx2[,c("Sess_1", "Sess_2")]) # pre stim days

#data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_opto_counts.csv") # write file
# write data
if(post_threshold){
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_opto_counts_thresh.csv") # write file
}else{
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_opto_counts.csv") # write file
}  

rm(dat, d, idx, idx2, A_idx, avail_idx)




# **************************************************
# find baseline rates and re-process
# **************************************************
# for each id and target:
# -- sessions: c(-1, 0, 3, 4)
# -- preprocess as above
# -- find proportion of trials that have target
# -- make table with id, target, proportions-before, proportion-after
# -- left_join according to id and target 

dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- c(-1,0,3,4) # -1:0 are baselines and 3:4 are post (followup)
dat <- dat[dat$session_number %in% sess_opt,] # only baseline sessions 
# filter
dat <- dat[dat$area == "snc (axon)",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# availability/treatment indicator -- ****Note: all feedback_status = 0 for baseline experiments
dat$A <- dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable) # Note: no feedback_status in baseline experiments

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

idx2 <- idx %>% 
          table.express::group_by(mouse_id, target, session_number) %>%
          table.express::summarise(#prop = mean(available),
                                   total_counts = sum(available)) %>% # calculate proportions of being in right pose (available is indicator that met threshold of confidence)
          #tidyr::pivot_wider(names_from = session_number, values_from = prop ) %>% # pivot wider to make columns
          tidyr::pivot_wider(names_from = session_number, values_from = total_counts ) # pivot wider to make columns
colnames(idx2)[3:6] <- paste0("Sess_", c(-1,0,3,4))
setDF(idx2)
setDF(idx_final)

# join and calculate averages of pre-post proportions
idx_mrg <- left_join(idx_final, idx2, c("mouse_id", "target")) # join by target (because same for all stim days)
# idx_mrg$preProp <- rowMeans(idx_mrg[,c("Sess_0", "Sess_-1")]) # pre stim days
# idx_mrg$postProp <- rowMeans(idx_mrg[,c("Sess_3", "Sess_4")]) # post stim days
idx_mrg$preCount <- rowMeans(idx_mrg[,c("Sess_0", "Sess_-1")]) # pre stim days
idx_mrg$postCount <- rowMeans(idx_mrg[,c("Sess_3", "Sess_4")]) # post stim days


# write data
if(post_threshold){
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_prePost_thresh.csv", row.names = FALSE)
}else{
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_prePost.csv", row.names = FALSE)
}  



# **************************************************
# find control rates and re-process
# **************************************************
# for each id and target:
# -- sessions: c(-1, 0, 3, 4)
# -- preprocess as above
# -- find proportion of trials that have target
# -- make table with id, target, proportions-before, proportion-after
# -- left_join according to id and target 

dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- -1:4 # -1:0 are baselines and 3:4 are post (followup)
dat <- dat[dat$session_number %in% sess_opt,] # only baseline sessions 
# filter
dat <- dat[dat$opsin == "ctrl",] #"ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# label who is a learner
dat$learner <- 1*I(dat$mouse_id %in% learner_ids)

# **** Note that control animals *DO* have feedback_status for stim days
# availability/treatment indicator -- ****Note: all feedback_status = 0 for baseline experiments
dat$A <- dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable) # Note: no feedback_status in baseline experiments

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

idx2 <- idx %>% 
  table.express::group_by(mouse_id, target, session_number) %>%
  table.express::summarise(total_counts = sum(available)) %>% # calculate proportions of being in right pose (available is indicator that met threshold of confidence)
  tidyr::pivot_wider(names_from = session_number, values_from = total_counts ) # pivot wider to make columns
colnames(idx2)[3:8] <- paste0("Sess_", -1:4)

setDF(idx2)
setDF(idx_final)

idx2$optoCount <- rowMeans(idx2[,c("Sess_1", "Sess_2")]) # pre stim days
idx2$preCount <- rowMeans(idx2[,c("Sess_0", "Sess_-1")]) # pre stim days
idx2$postCount <- rowMeans(idx2[,c("Sess_3", "Sess_4")]) # post stim days


# write data
# data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_prePost_ctrl.csv", row.names = FALSE)
if(post_threshold){
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_prePost_ctrl_thresh.csv", row.names = FALSE)
}else{
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_prePost_ctrl.csv", row.names = FALSE)
}  



# ---------------------------------------------------------------------
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
# experiment_type == "reinforcement", session 1 & 2, and feedback_status > 0, supported by this file below:
# https://github.com/dattalab/dopamine-reinforces-spontaneous-behavior/blob/14ed714c4f49be1f4728a080eb65c4ee42f74c23/notebooks_panels/fig_3j%2C%203k-transition%20changes.ipynb#L121
# subset to specific experiment type
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- c(1, 2) # equivalently gives sessions 1-2:  unique(dat$session_number[dat$feedback_status == 1]) # sessions with stimulation (feedback_status == 1 is stim)
dat <- dat[dat$session_number %in% sess_opt,] # only sessions where there is sometimes stimulation

# filter
dat <- dat[dat$area == "ctrl",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# label who is a learner
dat$learner <- 1*I(dat$mouse_id %in% learner_ids)

# availability indicator
dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable &
                        dat$feedback_status >= 0) # make all target poses where feedback_status = -1 into equivalent of wrong pose

# treatment (stim)
dat$A <- 1*I( dat$available == 1 & dat$feedback_status == 1) # 

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

# **************************************************
idx_final <- idx # save final file to add below
# **************************************************
if(post_threshold){
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_control_thresh.csv") # write file
  
}else{
  data.table::fwrite(idx, "/home/loewingergc/optoDA_preprocess_new/optoDA_control.csv") # write file
}  

rm(dat, d, idx, A_idx, avail_idx)


#--------------------------------------------------------
# moving sum - active
#--------------------------------------------------------
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- -1:4#c(-1,0,3,4) # -1:0 are baselines and 3:4 are post (followup)
dat <- dat[dat$session_number %in% sess_opt,] # only baseline sessions 
# filter
dat <- dat[dat$area == "snc (axon)",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# availability/treatment indicator -- ****Note: all feedback_status = 0 for baseline experiments
dat$A <- dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable) # Note: no feedback_status in baseline experiments

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

# preprocess
idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

# calculate counts
idx_counts <- data_bin_conts(idx, bin_size = 30, Hz = 1/30, max_len = 1800)

# idx_trials <- data_bin_trials(idx, bin_size = 30, Hz = 1/30, max_len = 1800)
# d <- idx_trials %>% table.express::filter(target == 17,
#                                 mouse_id == 3440,
#                                 session_number == 2) ###%>%
                    ###table.express::select(count, bin_start, bin_end, "__index_level_0__", mouse_id)

  
  
idx2 <- idx_counts %>% 
  tidyr::pivot_wider(names_from = session_number, values_from = bin_counts ) #%>%# pivot wider to make columns
  
colnames(idx2)[6:11] <- paste0("S", 0:5)

idx_mrg = idx2 %>%
            dplyr::group_by(mouse_id, target) %>%
            dplyr::summarise(mouse_id = unique(mouse_id), target = unique(target),
                               "S_pre" = mean(c(S0, S1), na.rm=T),
                               "S_opt" = mean(c(S2, S3), na.rm=T),
                               "S_post" = mean(c(S4, S5), na.rm=T))
 
# write data
if(post_threshold){
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_summary_thresh.csv", row.names = FALSE)
}else{
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_summary.csv", row.names = FALSE)
}  


#--------------------------------------------------------
# moving sum - ctrl
#--------------------------------------------------------
dat <- arrow::read_parquet(file = "/lscratch/SpontaneousBehaviour/optoda_raw_data/closed_loop_behavior_transfer.parquet")

###################
# filter variables
###################
dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
sess_opt <- -1:4#c(-1,0,3,4) # -1:0 are baselines and 3:4 are post (followup)
dat <- dat[dat$session_number %in% sess_opt,] # only baseline sessions 
# filter
dat <- dat[dat$area == "ctrl",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used

# availability/treatment indicator -- ****Note: all feedback_status = 0 for baseline experiments
dat$A <- dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable) # Note: no feedback_status in baseline experiments

####################
# columns to remove
####################
col_rm <- which(grepl( "likes" , colnames( dat ) ))
col_rm <- c(col_rm, which(grepl( "pc" , colnames( dat ) )) )
dat <- dat[,-col_rm] # remove

index_col <- which(colnames(dat) == "__index_level_0__")
colnames(dat)[index_col] <- "index"

# preprocess
idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials

# calculate counts
idx_counts <- data_bin_conts(idx, bin_size = 30, Hz = 1/30, max_len = 1800)

idx2 <- idx_counts %>% 
  tidyr::pivot_wider(names_from = session_number, values_from = bin_counts ) #%>%# pivot wider to make columns

colnames(idx2)[6:11] <- paste0("S", 0:5)

idx_mrg = idx2 %>%
  dplyr::group_by(mouse_id, target) %>%
  dplyr::summarise(mouse_id = unique(mouse_id), target = unique(target),
                   "S_pre" = mean(c(S0, S1), na.rm=T),
                   "S_opt" = mean(c(S2, S3), na.rm=T),
                   "S_post" = mean(c(S4, S5), na.rm=T))

# write data
if(post_threshold){
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_summary_thresh_ctrl.csv", row.names = FALSE)
}else{
  data.table::fwrite(idx_mrg, "/home/loewingergc/optoDA_preprocess_new/optoDA_summary_ctrl.csv", row.names = FALSE)
}  




