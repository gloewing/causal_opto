library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(refund)
library(dplyr)

setwd("~/Desktop/NIMH Research/Causal/DA_adapts_learn_rate")
matdat <- R.matlab::readMat("seshMerge.mat")
matdat <- matdat$seshMerge

photo_idx <- 7
opt_idx <- 4
session_idx <- 2
trial_idx <- 1
sess_idx <- c(1:6)

matdat[,1,]

nms <- rownames(matdat[,1,])
cue_nms <- nms[sess_idx] # names with same dimension about trial type
# length(matdat[1,,])

n <- 24

# session data
dat <- do.call(cbind, lapply(sess_idx, function(x) do.call(c, matdat[x,,])  ) ) # seq_along(cue_nms)
colnames(dat) <- cue_nms
ids <- do.call(c, lapply(seq_along(1:n), function(x) rep(x, length( matdat[1,,][[x]])  ) ))
dat <- as.data.frame( cbind(ids, dat) )
photo <- do.call(rbind, matdat[photo_idx,,] )
colnames(photo) <- paste0("photometry.", 1:ncol(photo) )
dat <- cbind(dat, photo)
