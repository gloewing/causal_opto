# data pre-process code
data_preprocess <- function(dat, post_threshold = FALSE){
  
  # identify total number of unique sessions and convert dates to session # for each animal
  dat <- as.data.table(dat)
  dat <- dat[, by = mouse_id, session := as.numeric(factor(date))]
  
  d <- dat[order(mouse_id, session, date, target_syllable, index, timestamp)] # order rows
  ids <- unique(d$mouse_id)
  n <- length(ids)
  
  d <- d[, keyby = .(mouse_id, session), t := seq(1, .N)]
  A_idx <- d[which(d$A == 1), c("mouse_id", "session", "t", "A", "available", "predicted_syllable", "target_syllable")]
  avail_idx <- d[which(d$available == 1), c("mouse_id", "session", "t", "A", "available", "predicted_syllable", "target_syllable")]
  
  # find places where there are transitions
  idx <- d[, keyby = .(mouse_id, session),
           .(pose_diff = diff(predicted_syllable),
             pose = predicted_syllable[-.N],
             A = A[-.N],
             t = t[-.N],
             session_number = session_number[-.N],
             session_repeat = session_repeat[-.N],
             target = target_syllable[-.N],
             available = available[-.N])]
  
  idx <- idx[, keyby = .(mouse_id, session),
             .(pose_diff_idx = which(pose_diff != 0),
               pose_diff = pose_diff[which(pose_diff != 0)],
               pose = pose[which(pose_diff != 0)],
               t = t[which(pose_diff != 0)],
               A = A[which(pose_diff != 0)],
               session_repeat = session_repeat[which(pose_diff != 0)],
               session_number = session_number[which(pose_diff != 0)],
               target = target[which(pose_diff != 0)],
               available = available[which(pose_diff != 0)]) ]
  
  idx[,keyby = .(mouse_id, session), 
      c("t1", "t2") := .(c(1, t[-.N] + 1),
                         t)]
  
  # find length of poses
  idx <- idx[, keyby = .(mouse_id, session),
             pose_length := diff(c(0, pose_diff_idx))]#,
  
  
  for(i in ids){
    # sessions for this animal
    index <- which(d$mouse_id == i) # rows for this animal
    sessions <- unique(d$session[index]) # unique sessions
    
    for(s in sessions){
      # iterate over sessions
      # indices for different datasets
      index2 <- which(idx$mouse_id == i        &      idx$session == s) # working dataset: current animal (i) and session
      idx_i <- idx[index2,]
      
      # find availability
      # availability: current animal (i) and session
      avail_i <- avail_idx[which(avail_idx$mouse_id == i  &  avail_idx$session == s),] # session rows for that animal when available
      available_index <- sapply(1:nrow(idx_i), function(x)  any(idx_i$t2[x] >= avail_i$t & idx_i$t1[x] <= avail_i$t ) ) # iterate over time "trials" in idx_i and see which contain availability timepoints 
      idx_i$available[which(available_index)] <- 1 # all trials that are in 
      idx$available[index2] <- idx_i$available # resave original idx
      
      # find stimulation
      # treatment: current animal (i) and session
      A_i <- A_idx[which(A_idx$mouse_id == i & A_idx$session == s),] # session rows for that animal
      stim_index <- sapply(1:nrow(idx_i), function(x)  any(idx_i$t2[x] >= A_i$t & idx_i$t1[x] <= A_i$t ) ) # iterate over time "trials" in idx_i and see which availability timepoints that fill in between trial start/stop
      idx_i$A[which(stim_index)] <- 1 # all trials that are in 
      idx$A[index2] <- idx_i$A # resave original idx
      
      rm(A_i, idx_i, index2, avail_i, available_index, stim_index)
    }
  }
  
  ##################################################################
  # remove trials with less than 3 consecutive timepoints
  ##################################################################
  if(post_threshold)   idx <- idx[pose_length >= min_trial_frames,]
  idx[, keyby = .(mouse_id, session), trial := seq(1, .N)]
  
  return(idx)
}



data_bin_conts <- function(idx, bin_size = 30, Hz = 1/30, max_len = 1800){
  
  # sequence only goes up until 1800!! So removes later portions of sessions that are longer
  seq1 <- seq(1, max_len/Hz + 1, by = bin_size/Hz)
  seq2 <- seq1[-1]
  seq1 <- seq1[-length(seq1)]
  
  res_lst <- vector("list", length = 336)
  ids <- unique(idx$mouse_id)
  cnt <- 0
  
  for(i in ids){
    print(i)
    # sessions for this animal
    index <- which(idx$mouse_id == i) # rows for this animal
    sessions <- unique(idx$session[index]) # unique sessions
    
    for(s in sessions){
      # iterate over sessions
      # indices for different datasets
      index2 <- which(idx$mouse_id == i        &      idx$session == s) # working dataset: current animal (i) and session
      idx_i <- idx[index2,]
      tg = unique(idx_i$target)
      sn = unique(idx_i$session_number)
      sr = unique(idx_i$session_repeat)
      
      # subset based on availability ****This depends on whether available is defined with or without feedback_status****
      idx_i <- idx_i[idx_i$available == 1,] # only use rows where the animal is available (did target pose)
      
      if(nrow(idx_i) > 0){
        # make sure there is some availability
        # find availability
        # availability: current animal (i) and session
        count_index <- sapply(1:length(seq1), function(x)  sum(seq1[x] <= idx_i$t1 & seq2[x] >= idx_i$t1 | 
                                                                 seq1[x] <= idx_i$t2 & seq2[x] >= idx_i$t2) ) # iterate over time "trials" in idx_i and see which contain start OR end timepoints (implicitly because idx_i only contains available rows)
        
        cnt <- cnt + 1
        # data-frame
        res_lst[[cnt]] <- data.frame(
          mouse_id = i, 
          session = s, 
          target = unique(idx_i$target), 
          bin_counts = sum(count_index), # this sums over counts so takes away bin structure (trial)
          total_counts = nrow(idx_i), # idx_i only has available rows so this is number of counts
          #bin = 1:length(count_index),
          session_number = unique(idx_i$session_number),
          session_repeat = unique(idx_i$session_repeat))
      }else{
        # no availability
        res_lst[[cnt]] <- data.frame(
          mouse_id = i, 
          session = s, 
          target = tg, 
          bin_counts = 0, # this sums over counts so takes away bin structure (trial)
          total_counts = 0, # idx_i only has available rows so this is number of counts
          #bin = 1:length(count_index),
          session_number = sn,
          session_repeat = sr)
      }
      
      rm(count_index, idx_i, index2)
    }
  }
  
  df <- do.call(rbind, res_lst)
  df <- df[order(df$mouse_id, df$session, df$target, df$bin),] # order rows
  
  
  return(df)
}

# makes into trial structure
data_bin_trials <- function(idx, bin_size = 30, Hz = 1/30, max_len = 1800){
  
  # sequence only goes up until 1800!! So removes later portions of sessions that are longer
  seq1 <- seq(1, max_len/Hz + 1, by = bin_size/Hz)
  bin_start <- seq(0, 1770, by = bin_size)
  seq2 <- seq1[-1]
  seq1 <- seq1[-length(seq1)]
  
  res_lst <- vector("list", length = 336)
  ids <- unique(idx$mouse_id)
  cnt <- 0
  
  for(i in ids){
    print(i)
    # sessions for this animal
    index <- which(idx$mouse_id == i) # rows for this animal
    sessions <- unique(idx$session[index]) # unique sessions
    
    for(s in sessions){
      # iterate over sessions
      # indices for different datasets
      index2 <- which(idx$mouse_id == i        &      idx$session == s) # working dataset: current animal (i) and session
      idx_i <- idx[index2,]
      tg = unique(idx_i$target)
      sn = unique(idx_i$session_number)
      sr = unique(idx_i$session_repeat)
      
      # subset based on availability ****This depends on whether available is defined with or without feedback_status****
      idx_i <- idx_i[idx_i$available == 1,] # only use rows where the animal is available (did target pose)
      
      if(nrow(idx_i) > 0){
        # make sure there is some availability
        # find availability
        # availability: current animal (i) and session
        count_index <- sapply(1:length(seq1), function(x)  sum(seq1[x] <= idx_i$t1 & seq2[x] >= idx_i$t1 | 
                                                                 seq1[x] <= idx_i$t2 & seq2[x] >= idx_i$t2) ) # iterate over time "trials" in idx_i and see which contain start OR end timepoints (implicitly because idx_i only contains available rows)
        
        cnt <- cnt + 1
        # data-frame
        res_lst[[cnt]] <- data.frame(
          mouse_id = i, 
          session = s, 
          target = unique(idx_i$target), 
          bin_counts = (count_index), # this sums over counts so takes away bin structure (trial)
          total_counts = nrow(idx_i), # idx_i only has available rows so this is number of counts
          bin = 1:length(count_index),
          bin_start = bin_start[1:length(count_index)],
          session_number = unique(idx_i$session_number),
          session_repeat = unique(idx_i$session_repeat))
      }else{
        # no availability
        res_lst[[cnt]] <- data.frame(
          mouse_id = i, 
          session = s, 
          target = tg, 
          bin_counts = 0, # this sums over counts so takes away bin structure (trial)
          total_counts = 0, # idx_i only has available rows so this is number of counts
          bin = 1:length(bin_start),
          bin_start = bin_start,
          session_number = sn,
          session_repeat = sr)
      }
      
      rm(count_index, idx_i, index2)
    }
  }
  
  df <- do.call(rbind, res_lst)
  df <- df[order(df$mouse_id, df$session, df$target, df$bin),] # order rows
  
  
  return(df)
}