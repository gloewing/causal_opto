# dataset preprocess of Spontaneous DA on Curium
# pre-process data for replicating original analysis

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


# **************************************************
# find baseline rates and re-process
# **************************************************
# for each id and target:
# -- sessions: c(-1, 0)
# -- preprocess as above
# -- find proportion of trials that have target
# -- make table with id, target, proportions-before, proportion-after
# -- left_join according to id and target 

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
sess_opt <- 1:2 # -1:0 are baselines and 3:4 are post (followup)
# dat <- dat[dat$experiment_type %in% c("reinforcement", "reinforcement_photometry"),] # reinforcement seems to be the only experiment type with all 9 "learner" mice
# dat <- dat[dat$session_number %in% sess_opt,] # only baseline sessions 
# # filter
# dat <- dat[dat$area == "snc (axon)",] # "snc (axon)" -- opto animals & "ctrl" -- control animals
# dat <- dat[ !dat$SessionName %in% sess_rm, ] # remove rows that can't be used
# 
# # availability/treatment indicator -- ****Note: all feedback_status = 0 for baseline experiments
# dat$A <- dat$available <- 1*I( dat$target_syllable == dat$predicted_syllable) # Note: no feedback_status in baseline experiments
dat <- dat %>% 
  table.express::filter(experiment_type %in% c("reinforcement", "reinforcement_photometry"), # reinforcement seems to be the only experiment type with all 9 "learner" mice
                        session_number %in% sess_opt, # only sessions where there is sometimes stimulation
                        # area == "snc (axon)", # "snc (axon)" -- opto animals & "ctrl" -- control animals
                        !dat$SessionName %in% sess_rm ) %>% # remove rows that can't be used
  table.express::mutate(learner = 1*I(mouse_id %in% learner_ids), # label who is a learner
                        available = 1*I(target_syllable == predicted_syllable),
                        A = 1*I(target_syllable == predicted_syllable)) # A is just dummy variable here

idx <- data_preprocess(dat, post_threshold = post_threshold) # pre-process data to turn into trials



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
  tidyr::pivot_wider(names_from = session_number, values_from = total_counts ) # pivot wider to make columns
colnames(idx2)[3:4] <- paste0("Sess_", sess_opt)
setDF(idx2)

# join and calculate averages of pre-post proportions
# idx_mrg <- left_join(idx_final, idx2, c("mouse_id", "target")) # join by target (because same for all stim days)
idx2$optoCount <- rowMeans(idx2[,c("Sess_1", "Sess_2")]) # pre stim days
# idx_mrg$postCount <- rowMeans(idx_mrg[,c("Sess_3", "Sess_4")]) # post stim days


# write data
if(post_threshold){
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_active_thresh.csv", row.names = FALSE)
}else{
  data.table::fwrite(idx2, "/home/loewingergc/optoDA_preprocess_new/optoDA_active.csv", row.names = FALSE)
}  
