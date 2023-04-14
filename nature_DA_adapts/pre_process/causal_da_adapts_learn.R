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
                      reward_DA = reward_DA - base_DA)

# join data together
dat <- as.data.frame( cbind(ids, ITIs, dat, DA_mat) )
# dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(iti_vec, ids, ITIs, DA_mat, photo, trial_num)


head(dat)


dat %>% as_tibble() %>%
  dplyr::filter(group == "slPlus") %>%
  group_by(stimState) %>%
  # summarise(m_cue = mean(cue_DA),
  #           m_rew = mean(reward_DA))
  ggplot(aes(x = as.factor(stimState), y = reward_DA)) +
  geom_boxplot()

# # make reward period DA into median to make IV indicator outcome below
# dat <- dat %>% as_tibble() %>%
#           dplyr::group_by(ids) %>%
#           dplyr::mutate(reward_median = median(reward_DA) )

################################################
# parametric, ATE, no-extrapolotation, Group 1
################################################
# E[X_{j+1} | X_{j} = 1, L_j, j, G = 1 ]
# cluster bootstrap for uncertainty
# GLMMM model wiith random intercepts
library(lme4)

g1 <- unique(dat$ids[dat$group == "slPlus"])
dat <- dat[order(dat$ids, dat$trial_num),] # make sure order is correct
id_list <- id_list2 <- vector(length = length(g1), "list")
da_list <- vector(length = length(g1)) # for previous trial's dopamine
dat$d_trt <- rep(NA, nrow(dat))
dat$d_rew_previous <- rep(NA, nrow(dat))

# find trials where x_{j-1} = 1
for(i in 1:length(g1)){
  id_num <- g1[i]
  idx <- which(dat$ids == id_num)
  reward_median <- median(dat$reward_DA[idx], na.rm = TRUE)
  id_list[[i]] <- idx[c(FALSE, sapply(2:length(idx), function(x) dat$lickState[ idx[x - 1]] == 1 ) )] # whether previous trial was a lick trial
  id_list2[[i]] <- 1 * sapply(2:length(id_list[[i]]), function(x) dat$stimState[ idx[x - 1]] == 1 ) # whether previous lick trial was a stim trial
  dat$d_trt[ id_list[[i]][-1] ] <- 1 * sapply(2:length(id_list[[i]]), function(x) dat$reward_DA[ idx[x - 1]] >= reward_median ) # go through each trial for which previous trial was a lick and save previous trial's reward period DA AUC
  dat$d_rew_previous[ id_list[[i]][-1] ] <- sapply(2:length(id_list[[i]]), function(x) dat$reward_DA[ idx[x - 1]] ) # go through each trial for which previous trial was a lick and save previous trial's reward period DA AUC
  id_list[[i]] <- id_list[[i]][-1] # remove first trial because even though it was a lick, it could not have had a laser/no-laser before
  }

id_list <- do.call(c, id_list)
id_list2 <- do.call(c, id_list2)

dat2 <- dat[id_list,]
dat2$stimPrev <- id_list2 # whether they licked on the previous trial

#### sample mean 
m <- dat2 %>% as_tibble %>%
      group_by(trial_num, stimPrev) %>%
      summarise(probs = mean(lickState))

cont1 <- data.frame(trial_num = m$trial_num[pred_mat$stimPrev == 1],  m$probs[pred_mat$stimPrev == 1])
cont0 <- data.frame(trial_num = m$trial_num[pred_mat$stimPrev == 0], m$probs[pred_mat$stimPrev == 0])
cont <- inner_join(cont1, cont0, by = "trial_num")
data.frame(cont, ate = cont$m.probs.pred_mat.stimPrev....1. - cont$m.probs.pred_mat.stimPrev....0.) %>%
  ggplot2::ggplot( aes(x = trial_num, y = ate)) +
  geom_line()

################


# dat2 <- dat2[dat2$trial_num <= 400,]

dat2$trial_num <- scale(dat2$trial_num)
# 
# mod1 <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + # + I(trial_num^2) +
#                       #stimPrev * trial_num + #stimPrev:I(trial_num^2) +
#                       (1|ids),
#                       family = "binomial",
#                       data = dat2)
# 
# 
# 
# # variance
# a = summary(mod1)
# b_var = a$vcov
# t_len <- length(trials)
# f_mat <- as.matrix(splines::bs(trials))  #  rep(c(TRUE, FALSE), t_len )
# contrast_mat <- cbind(t(replicate(t_len, c(0, 1, rep(0, ncol(f_mat))))), f_mat)
# var_est <- diag(contrast_mat %*% b_var %*% t(contrast_mat))
# preds <- plogis(contrast_mat %*% coef(git_gee) )
# ci_preds <- cbind(preds + 1.96 * sqrt(var_est),
#                   preds - 1.96 * sqrt(var_est))
# 
# plot(preds, type = "l", ylim = c(-0.6,1.2))
# lines(ci_preds[,1], col = "red")
# lines(ci_preds[,2], col = "red")
# abline(h = 0)
# 
# 
# # --------------------------------------
# 
# AIC(mod1)
# trials <- unique(dat2$trial_num)
# stim <- c(0, 1)
# pred_mat <- expand.grid(trials, stim, g1)
# colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
# preds <- predict(mod1, pred_mat, type = "response")
# pred_mat <- cbind(pred_mat, preds) %>% as_tibble() %>%
#                 group_by(trial_num, stimPrev) %>%
#                 summarise(probs = mean(preds)) %>%
#                 unique()
# 
# pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]
# 
# cont1 <- pred_mat$probs[pred_mat$stimPrev == 1]
# cont0 <- pred_mat$probs[pred_mat$stimPrev == 0]
# data.frame(trial = unique(pred_mat$trial_num), ate = cont1 - cont0) %>%
#             ggplot2::ggplot( aes(x = trial, y = ate)) +
#               geom_line()

dat2 %>% as_tibble() %>%
  group_by(ids) %>%
  unique() %>%
  summarise(m = n())

################################################################################
#  bootstrap LME
################################################################################
max_trial <- 600 # weird behavior after this
dat3 <- dat2[dat2$trial_num <= max_trial,]
trials <- unique(dat3$trial_num)
samp <- unique(dat3$ids)
stim <- c(0, 1)
pred_mat <- expand.grid(trials, stim, unique(samp))
boot_2nd_level <- FALSE # sample without replcement at 2nd level

# bootstrap

set.seed(1)
boots <- 200
res <- res_noIV <- matrix(nrow = length(trials), ncol = boots)
ids <- unique(dat3$ids)
n <- length(ids)

for(i in 1:boots){
  print(i)
  samp <- sample(ids, n, replace = TRUE)
  s1 <- sapply(samp, function(x) which(dat3$ids == x) )
  subj_id <- do.call(c, sapply(s1, function(x) sample(x, length(x), replace = boot_2nd_level) ))    
  
  # git_gee <- geeglm(lickState ~ stimPrev * splines::bs(trial_num),
  #                   id = ids,
  #                   family = "binomial",
  #                   scale.fix=TRUE,
  #                   corstr= "independence", #"exchangeable",
  #                   data = dat3[subj_id,] )
  mod1 <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + # + I(trial_num^2) +
                           #stimPrev * trial_num + #stimPrev:I(trial_num^2) +
                           (1|ids),
                         family = "binomial",
                         data = dat3[subj_id,])
  
  mod2 <- lme4::glmer(d_trt ~ stimPrev * splines::bs(trial_num) + # + I(trial_num^2) +
                        #stimPrev * trial_num + #stimPrev:I(trial_num^2) +
                        (1|ids),
                      family = "binomial",
                      data = dat3[subj_id,])
  
  # predict
  trials <- unique(dat3$trial_num)
  stim <- c(0, 1)
  # --------------------------------------------------
  pred_mat <- expand.grid(trials, stim, unique(samp))
  colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
  preds <- predict(mod1, pred_mat, type = "response")
  
  # ate numerator
  pred_mat <- cbind(pred_mat, preds) %>% as_tibble() %>%
    group_by(trial_num, stimPrev) %>%
    summarise(probs = mean(preds)) %>%
    unique()
  
  pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]
  
  cont1 <- pred_mat$probs[pred_mat$stimPrev == 1]
  cont0 <- pred_mat$probs[pred_mat$stimPrev == 0]
  # --------------------------------------------------
  
  # --------------------------------------------------
  # ate denominator
  trials <- unique(dat3$trial_num)
  stim <- c(0, 1)
  pred_mat <- expand.grid(trials, stim, unique(samp))
  colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
  
  preds2 <- predict(mod2, pred_mat, type = "response")
  
  pred_mat <- cbind(pred_mat, preds2) %>% as_tibble() %>%
    group_by(trial_num, stimPrev) %>%
    summarise(probs = mean(preds2)) %>%
    unique()
  
  pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]
  
  cont1_denom <- pred_mat$probs[pred_mat$stimPrev == 1]
  cont0_denom <- pred_mat$probs[pred_mat$stimPrev == 0]
  # --------------------------------------------------
  
  # IV
  res[,i] <- (cont1 - cont0) / (cont1_denom - cont0_denom)
  
  # no IV
  res_noIV[,i] <- (cont1 - cont0)
}

# standard errors
sd_vec <- apply(res, 1, sd)
sd_vec_noIV <- apply(res_noIV, 1, sd)

# means
mod1 <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + 
                      (1|ids),
                    family = "binomial",
                    data = dat3)

# means
mod2 <- lme4::glmer(d_trt ~ stimPrev * splines::bs(trial_num) + 
                      (1|ids),
                    family = "binomial",
                    data = dat3)

# predict
trials <- unique(dat3$trial_num)
stim <- c(0, 1)
pred_mat <- expand.grid(trials, stim, unique(samp))
colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
preds <- predict(mod1, pred_mat, type = "response")

# ate
pred_mat <- cbind(pred_mat, preds) %>% as_tibble() %>%
  group_by(trial_num, stimPrev) %>%
  summarise(probs = mean(preds)) %>%
  unique()

pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]

cont1 <- pred_mat$probs[pred_mat$stimPrev == 1]
cont0 <- pred_mat$probs[pred_mat$stimPrev == 0]

# no effect :()
ate <- cont1 - cont0

# --------------------------------------------------
# ate denominator
trials <- unique(dat3$trial_num)
stim <- c(0, 1)
pred_mat <- expand.grid(trials, stim, unique(samp))
colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
preds2 <- predict(mod2, pred_mat, type = "response")

pred_mat <- cbind(pred_mat, preds2) %>% as_tibble() %>%
  group_by(trial_num, stimPrev) %>%
  summarise(probs = mean(preds2)) %>%
  unique()
# pred_mat <- cbind(pred_mat, a=preds2/preds) %>% as_tibble() %>%
#   group_by(trial_num, stimPrev) %>%
#   summarise(probs = mean(a)) %>%
#   unique()

pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]

cont1_denom <- pred_mat$probs[pred_mat$stimPrev == 1]
cont0_denom <- pred_mat$probs[pred_mat$stimPrev == 0]
# --------------------------------------------------
jpeg("~/Desktop/NIMH Research/Causal/DA_adapts_learn_rate/IV_figs/iv_binary.jpg")
par(mfrow = c(1,2))
# IV
ate_iv <- ate / (cont1_denom - cont0_denom)
#ate_iv <- (cont1_denom - cont0_denom)
plot(ate_iv, type = "l", col = "red", ylim = c(-0.3,.35), main = "IV Binary")
lines(ate_iv + 1.96*sd_vec)
lines(ate_iv - 1.96*sd_vec)
abline(h = 0)

# no IV
plot(ate, type = "l", col = "red", ylim = c(-0.3,.15), main = "ITT")
lines(ate + 1.96*sd_vec_noIV)
lines(ate - 1.96*sd_vec_noIV)
abline(h = 0)
dev.off()


################################################################################
#  bootstrap LME - IV structural model (weight observations and use continous A)
################################################################################
max_trial <- 600 # weird behavior after this
dat3 <- dat2[dat2$trial_num <= max_trial,]
trials <- unique(dat3$trial_num)
samp <- unique(dat2$ids)
stim <- c(0, 1)
pred_mat <- expand.grid(trials, stim, unique(samp))

# bootstrap

set.seed(1)
boots <- 285
res <- res_noIV <- matrix(nrow = length(trials), ncol = boots)
ids <- unique(dat3$ids)
n <- length(ids)
boot_2nd_level <- FALSE # sample without replcement at 2nd level

for(i in 1:boots){
  print(i)
  samp <- sample(ids, n, replace = TRUE)
  s1 <- sapply(samp, function(x) which(dat3$ids == x) )
  subj_id <- do.call(c, sapply(s1, function(x) sample(x, length(x), replace = boot_2nd_level) ))    
  
  # git_gee <- geeglm(lickState ~ stimPrev * splines::bs(trial_num),
  #                   id = ids,
  #                   family = "binomial",
  #                   scale.fix=TRUE,
  #                   corstr= "independence", #"exchangeable",
  #                   data = dat3[subj_id,] )
  mod1 <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + # + I(trial_num^2) +
                        #stimPrev * trial_num + #stimPrev:I(trial_num^2) +
                        (1|ids),
                      family = "binomial",
                      data = dat3[subj_id,])
  
  # means
  mod2 <- lme4::lmer(d_rew_previous ~ stimPrev * splines::bs(trial_num) + 
                       (1|ids),
                     data = dat3[subj_id,])
  
  # predict
  trials <- unique(dat3$trial_num)
  stim <- c(0, 1)
  # --------------------------------------------------
  pred_mat <- expand.grid(trials, stim, unique(samp))
  colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
  preds <- predict(mod1, pred_mat, type = "response")
  preds2 <- predict(mod2, pred_mat)
  # --------------------------------------------------

  # ate
  pred_mat <- cbind(pred_mat, preds, preds2) %>% as_tibble() %>%
    group_by(trial_num, stimPrev) %>%
    # summarise(probs = mean(preds)) %>%
    unique()
  
  pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]
  
  cont1 <- pred_mat$preds[pred_mat$stimPrev == 1]
  cont0 <- pred_mat$preds[pred_mat$stimPrev == 0]
  
  cont1_denom <- pred_mat$preds2[pred_mat$stimPrev == 1]
  cont0_denom <- pred_mat$preds2[pred_mat$stimPrev == 0]
  
  # no effect :()
  ate <- cont1 - cont0
  ate_iv <- ate / (cont1_denom-cont0_denom)
  
  ate_mat <- cbind(ate_iv, ate, pred_mat$trial_num[pred_mat$stimPrev == 1]) 
  colnames(ate_mat) <- c("ate_iv", "ate", "trial_num")
  
  ate_mat <- ate_mat %>% 
    as_tibble() %>%
    group_by(trial_num) %>%
    summarise(ate_mean = mean(ate),
              ate_iv_mean = mean(ate_iv) ) %>%
    unique()
  # --------------------------------------------------
  # IV
  res[,i] <- ate_mat$ate_iv_mean
  
  # no IV
  res_noIV[,i] <- ate_mat$ate_mean
}

# standard errors
sd_vec <- apply(res[,1:boots], 1, sd)
sd_vec_noIV <- apply(res_noIV[,1:boots], 1, sd)

# means
mod1 <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + 
                      (1|ids),
                    family = "binomial",
                    data = dat3)

# means
mod2 <- lme4::lmer(d_rew_previous ~ stimPrev * splines::bs(trial_num) + 
                      (1|ids),
                    data = dat3)

# predict
trials <- unique(dat3$trial_num)
stim <- c(0, 1)
pred_mat <- expand.grid(trials, stim, unique(samp))
colnames(pred_mat) <- c("trial_num", "stimPrev", "ids")
preds <- predict(mod1, pred_mat, type = "response")
preds2 <- predict(mod2, pred_mat) # , type = "response"

# ate
pred_mat <- cbind(pred_mat, preds, preds2) %>% as_tibble() %>%
  group_by(trial_num, stimPrev) %>%
  # summarise(probs = mean(preds)) %>%
  unique()

pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]

cont1 <- pred_mat$preds[pred_mat$stimPrev == 1]
cont0 <- pred_mat$preds[pred_mat$stimPrev == 0]

cont1_denom <- pred_mat$preds2[pred_mat$stimPrev == 1]
cont0_denom <- pred_mat$preds2[pred_mat$stimPrev == 0]

# no effect :()
ate <- cont1 - cont0
ate_iv <- ate / (cont1_denom-cont0_denom)

ate_mat <- cbind(ate_iv, ate, pred_mat$trial_num[pred_mat$stimPrev == 1]) 
colnames(ate_mat) <- c("ate_iv", "ate", "trial_num")

ate_mat <- ate_mat %>% 
  as_tibble() %>%
  group_by(trial_num) %>%
  summarise(ate_mean = mean(ate),
            ate_iv_mean = mean(ate_iv) ) %>%
  unique()
# --------------------------------------------------

jpeg("~/Desktop/NIMH Research/Causal/DA_adapts_learn_rate/IV_figs/iv_continuous.jpg")
par(mfrow = c(1,2))
# IV
ate_iv <- ate_mat$ate_iv_mean
#ate_iv <- (cont1_denom - cont0_denom)
plot(ate_iv, type = "l", col = "red", ylim = c(-0.3,.15), main = "ATE IV Continuous")
lines(ate + 1.96*sd_vec)
lines(ate - 1.96*sd_vec)
abline(h = 0)

# no IV
ate_mean <- ate_mat$ate_mean
plot(ate_mean, type = "l", col = "red", ylim = c(-0.3,.15), main = "ITT")
lines(ate + 1.96*sd_vec_noIV)
lines(ate - 1.96*sd_vec_noIV)
abline(h = 0)
dev.off()


################################################################################
#  bootstrap GEE -- too slow!!
################################################################################
################
# GEE
################
library(geepack)
deg <- 3
knots <- c(seq(200, 800, by = 200))
git_gee <- geeglm(lickState ~ stimPrev * splines::bs(trial_num),
                  id = ids,
                  family = "binomial",
                  scale.fix=TRUE,
                  corstr= "independence", #"exchangeable",
                  data = dat2)

trials <- min(unique(dat2$trial_num)):max(unique(dat2$trial_num))
stim <- c(0, 1)
pred_mat <- rbind(cbind(trials, rep(0, nrow(x))), cbind(trials, rep(1, nrow(x))))
colnames(pred_mat) <- c("trial_num", "stimPrev")
pred_mat <- as.data.frame(pred_mat)
preds <- predict(git_gee, pred_mat, type = "response")

set.seed(1)
boots <- 500
res <- matrix(nrow = nrow(pred_mat), ncol = boots)
ids <- unique(dat2$ids)
n <- length(ids)

for(i in 1:boots){
  print(i)
  samp <- sample(ids, n, replace = TRUE)
  s1 <- sapply(samp, function(x) which(dat2$ids == x) )
  subj_id <- do.call(c, sapply(s1, function(x) sample(x, length(x), replace = TRUE) ))    
  
  # git_gee <- geeglm(lickState ~ stimPrev * splines::bs(trial_num),
  #                   id = ids,
  #                   family = "binomial",
  #                   scale.fix=TRUE,
  #                   corstr= "independence", #"exchangeable",
  #                   data = dat2[subj_id,] )
  git_gee <- lme4::glmer(lickState ~ stimPrev * splines::bs(trial_num) + # + I(trial_num^2) +
                        #stimPrev * trial_num + #stimPrev:I(trial_num^2) +
                        (1|ids),
                      family = "binomial",
                      data = dat2)
  
  # predict
  trials <- min(unique(dat2$trial_num)):max(unique(dat2$trial_num))
  stim <- c(0, 1)
  pred_mat <- rbind(cbind(trials, rep(0, nrow(x))), cbind(trials, rep(1, nrow(x))))
  colnames(pred_mat) <- c("trial_num", "stimPrev")
  pred_mat <- as.data.frame(pred_mat)
  res[,i] <- predict(git_gee, pred_mat, type = "response")
}

#plot(preds, type = "l")

# variance
a = summary(git_gee)
b_var = a$cov.unscaled
t_len <- length(trials)
f_mat <- as.matrix(splines::bs(trials))  #  rep(c(TRUE, FALSE), t_len )
contrast_mat <- cbind(t(replicate(t_len, c(0, 1, rep(0, ncol(f_mat))))), f_mat)
var_est <- diag(contrast_mat %*% b_var %*% t(contrast_mat))
preds <- plogis(contrast_mat %*% coef(git_gee) )
ci_preds <- cbind(preds + 1.96 * sqrt(var_est),
         preds - 1.96 * sqrt(var_est))

plot(preds, type = "l", ylim = c(-0.6,1.2))
lines(ci_preds[,1], col = "red")
lines(ci_preds[,2], col = "red")
abline(h = 0)




# trials <- min(unique(dat2$trial_num)):max(unique(dat2$trial_num))
# stim <- c(0, 1)
# x <- splines::bs(trials)
# pred_mat <- rbind(cbind(x, rep(0, nrow(x))), cbind(x, rep(1, nrow(x))))
# colnames(pred_mat) <- c(names(git_gee$coefficients)[3:5], "stimPrev")
# pred_mat <- as.data.frame(pred_mat)
# preds <- predict(git_gee, pred_mat, type = "response")

# ---------------------
pred_mat <- cbind(pred_mat, preds) %>% as_tibble() %>%
  group_by(trial_num, stimPrev) %>%
  summarise(probs = mean(preds)) %>%
  unique()

pred_mat <- pred_mat[order(pred_mat$trial_num, pred_mat$stimPrev),]

cont1 <- pred_mat$probs[pred_mat$stimPrev == 1]
cont0 <- pred_mat$probs[pred_mat$stimPrev == 0]
data.frame(trial = unique(pred_mat$trial_num), ate = cont1 - cont0) %>%
  ggplot2::ggplot( aes(x = trial, y = ate)) +
  geom_line()


m <- data.frame(trial_num = rep(1:1000, each = 2),
                stimPrev = rep(c(0,1), each = 2))
x = as.matrix(model.matrix( ~ stimPrev * splines::bs(trial_num), data = m))


ses <- sqrt(diag( x %*% tcrossprod(a, x) ))
# preds -- linear combo of betas




library(ggeffects)
x1 <- ggpredict(git_gee, pred_mat)
# plot(y = cont1 - cont0, x = unique(pred_mat$trial_num))

# pred_mat <- matrix(NA, ncol)

# dat2 %>% as_tibble() %>%
#   group_by(stimState, lickState) %>%
#   summarise(m = n())


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


# outcome for figure 5:
# As in the ACTR model, behaviour was biased in opposite directions for each contingency, with stimLick+ animals exhibiting lower and stimLick− animals exhibiting higher preparatory licking 
# groups: stim+, stim-, control (potentially could use uncalibrated DA in figure 6 in groups too)
# does it matter if instrument variesd in lickPlussStmPlus?