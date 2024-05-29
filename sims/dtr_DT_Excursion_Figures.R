# Figures and Tables for DTR -- Specific Excursion Effects

######################################################
# compare different K and betaVar -- no regularization
######################################################
setwd("/Users/Desktop/Research/msm_sims_normals/files_c") # msm_sims_dt2")
save_folder <- "/Causal/msm_hr/Final Sims/Figures"
wd <- "/Causal/msm_hr/Final Sims/Code/" 
source(paste0(wd,"saveFn.R"))
library(dplyr)
library(stringr)
library(latex2exp)
library(kableExtra)
library(tidyverse)

dd <- read.csv("/Desktop/Research/msm_sims_gauss/files_c/dtr_msm_n_4_tm_10_totSims_1000_boots_5000_c_agg2")
colnm = colnames(dd)
process <- TRUE
tVec <- c(10, 50, 500)
nVec <- c(4,6,8,10,30,100) 
nIncl <- c(6,10,30,100) # include in figures
nIncl2 <- c(10,30,100)
Tincl2 <- c(50,500) # timepoints to include
itrs <- numSims <- 1000 # number of simulation iterations
boots <- 5000
ls_ci <- ls_power <- ls_bias <- vector(length = length(nVec) * length(nVec), "list")

cnt <- 0
rmseMat <- matrix(nc = 4, nr = itrs)

# move all the "_c_" files to a different folder
# if(process)   move_files(list.files(), folder_name = "files_c")

for(t in tVec){
  for(n in nVec){
    # flNm <-  paste0("dtr_msm",
    #                 "_n_", n,
    #                 "_tm_", t,
    #                 "_totSims_", numSims,
    #                 "_boots_", boots, "_c_agg2")
    # if(file.exists(flNm)){
    #   
      flNm <-  paste0("dtr_msm",
                      "_n_", n,
                      "_tm_", t,
                      "_totSims_", numSims,
                      "_boots_", boots, "_c")
      
      # make sure outputted file does not exist
      fl_final <-  paste0("dtr_msm",
                      "_n_", n,
                      "_tm_", t,
                      "_totSims_", numSims,
                      "_boots_", boots, "_c.csv") # "_agg")#
      
      if(process & !file.exists(fl_final)){
        process_files(fileNm = flNm, itrs = numSims, colnm = colnm, rm = NULL) # "agg2"
      }   
      
      # check to see if file exists
      cnt <- cnt + 1
      flNm <-  paste0("dtr_msm",
                      "_n_", n,
                      "_tm_", t,
                      "_totSims_", numSims,
                      "_boots_", boots, "_c.csv") # "_agg")#
      
      if(!file.exists(flNm)){
        flNm <-  paste0("dtr_msm",
                        "_n_", n,
                        "_tm_", t,
                        "_totSims_", numSims,
                        "_boots_", boots, "_agg")#
      }   
      
      d <- read.csv(flNm) 
     
      
      # remove NAs
      ind <- apply(d, 1, function(x) all(is.na(x)))
      d <- d[!ind,]
      
      # CI coverage
      colIdx <- grep(paste0("^", "res_var"), colnm) # files with names that match
      ciMat <- d[,colIdx]
      
      # power
      colIdx <- grep(paste0("^", "power_var"), colnm) # files with names that match
      powerMat <- d[,colIdx]
      
      # bias
      colIdx <- grep(paste0("^", "bias_"), colnm) # files with names that match
      biasMat <- d[,colIdx]
      
      
      ls_ci[[cnt]] <- cbind(gather(ciMat), 
                              t,n)
       
      
      ls_power[[cnt]] <- cbind(gather(powerMat), 
                            t,n)
      
      ls_bias[[cnt]] <- cbind(gather(biasMat), 
                            t,n)
      
      d1 <- d
      rm(d)
    #}
  }
}

####################################################################################################
# set colors
myColors <- c("#ca0020", "#0868ac", "#525252", "#E69F00", "darkgreen", "darkgrey")

# excursion effect contrast matrix used in simulations
# c_mat <- rbind(c(0,0,1,0), # Blip
#                c(0,-1,1,0), # effect dissipation
#                c(1,0,0,0), # dose 0
#                c(1,0.5,0.5,0), # dose 1
#                rep(1,4) ) # dose 2

#########################
# 95% Coverage
#########################
dat <- do.call(rbind, ls_ci)
dat$t <- as.factor(dat$t)
dat$n <- as.factor(dat$n)

prefix <- "res_var_"
dat <- dat[complete.cases(dat),]

# sqrt(diag(var_array[1,,])) # HC0
# sqrt(diag(var_array[4,,])) # our estimator
# sqrt(diag(var_array[2,,])) # "HC2"
# sqrt(diag(var_array[3,,])) # "HC3"

# Blip and Effect Dissipation Excursion Effect ** Contrasts **
plt_rmse = 
  dat %>% tibble %>%  
  dplyr::mutate(var_idx = gsub( prefix, "", gsub("\\_coef.*","", key)), #substring(gsub(prefix, "", key) ,1,1), # extract variance
                coef_idx = gsub( ".*coef_", "", key)) %>% #str_sub(gsub(prefix, "", key) ,-1,-1) ) %>%
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl,
                var_idx %in% 1:4 # only HC0-HC3
                ) %>%
  dplyr::mutate(var_idx = paste0("HC", as.numeric(var_idx )),
                coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  dplyr::mutate(var_idx = plyr::revalue(var_idx, c("HC1" = "HC",
                                                   "HC4" = "Ours")),
                coef_idx = plyr::revalue(coef_idx, c("0" = "Blip",
                                                     "1" = "Dissipation",
                                                     "2" = "Dose 0",
                                                     "3" = "Dose 1",
                                                     "4" = "Dose 2") )) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n, coef_idx, var_idx) %>% 
  dplyr::summarize(my_mean = mean(as.numeric(value), na.rm = TRUE) ) %>% 
  ggplot(aes( y = my_mean, x = var_idx)) +
  facet_grid(t ~ n) +
  geom_point(aes(color = factor(coef_idx))) + 
  geom_hline(yintercept=0.95, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("95% CI Coverage" )+ 
  xlab("Robust Variance Estimator") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1))) +
  #guides(color= guide_legend(title=TeX('$\\beta$ Index')))  +
  guides(color= guide_legend(title="Excursion Effect"))

setwd(save_folder)
ggsave( "dtr_msm_DR_CI_excursion.pdf",
        plot = plt_rmse,
        width = 12,
        height = 8)

# short version
plt_rmse2 = 
  dat %>% tibble %>%  
  dplyr::mutate(var_idx = gsub( prefix, "", gsub("\\_coef.*","", key)), #substring(gsub(prefix, "", key) ,1,1), # extract variance
                coef_idx = gsub( ".*coef_", "", key)) %>% #str_sub(gsub(prefix, "", key) ,-1,-1) ) %>%
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl2,
                t %in% Tincl2,
                var_idx %in% 1:4 # only HC0-HC3
  ) %>%
  dplyr::mutate(var_idx = paste0("HC", as.numeric(var_idx )),
                coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  dplyr::mutate(var_idx = plyr::revalue(var_idx, c("HC1" = "HC",
                                                   "HC4" = "Ours")),
                coef_idx = plyr::revalue(coef_idx, c("0" = "Blip",
                                                     "1" = "Dissipation",
                                                     "2" = "Dose 0",
                                                     "3" = "Dose 1",
                                                     "4" = "Dose 2") )) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n, coef_idx, var_idx) %>% 
  dplyr::summarize(my_mean = mean(as.numeric(value), na.rm = TRUE) ) %>% 
  ggplot(aes( y = my_mean, x = var_idx)) +
  facet_grid(t ~ n) +
  geom_point(aes(color = factor(coef_idx))) + 
  geom_hline(yintercept=0.95, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("95% CI Coverage" )+ 
  xlab("Robust Variance Estimator") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1))) +
  #guides(color= guide_legend(title=TeX('$\\beta$ Index')))  +
  guides(color= guide_legend(title="Excursion Effect"))

#########################
# power
#########################
dat <- do.call(rbind, ls_power)
dat$t <- as.factor(dat$t)
dat$n <- as.factor(dat$n)

prefix <- "power_var_"
dat <- dat[complete.cases(dat),]

plt_power = 
  dat %>% tibble %>%  
  dplyr::mutate(var_idx = gsub( prefix, "", gsub("\\_coef.*","", key)), #substring(gsub(prefix, "", key) ,1,1), # extract variance
                coef_idx = gsub( ".*coef_", "", key)) %>% #str_sub(gsub(prefix, "", key) ,-1,-1) ) %>%
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl,
                var_idx %in% 1:4 # only HC0-HC3
  ) %>%
  dplyr::mutate(var_idx = paste0("HC", as.numeric(var_idx )),
                coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  dplyr::mutate(var_idx = plyr::revalue(var_idx, c("HC1" = "HC",
                                                   "HC4" = "Ours")),
                coef_idx = plyr::revalue(coef_idx, c("0" = "Blip",
                                                     "1" = "Dissipation",
                                                     "2" = "Dose 0",
                                                     "3" = "Dose 1",
                                                     "4" = "Dose 2") )) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n, coef_idx, var_idx) %>% 
  dplyr::summarize(my_mean = mean(as.numeric(value), na.rm = TRUE) ) %>% 
  ggplot(aes( y = my_mean, x = var_idx)) +
  facet_grid(t ~ n) +
  geom_point(aes(color = factor(coef_idx))) + 
  geom_hline(yintercept=0.95, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("Power" )+ 
  xlab("Robust Variance Estimator") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1))) +
  guides(color= guide_legend(title=TeX('$\\beta$ Index')))   

setwd(save_folder)
ggsave( "dtr_msm_DR_power_excursion.pdf",
        plot = plt_power,
        width = 12,
        height = 8)
########################################################

#########################
# bias
#########################
dat <- do.call(rbind, ls_bias)
dat$t <- as.factor(dat$t)
dat$n <- as.factor(dat$n)

prefix <- "bias_"
dat <- dat[complete.cases(dat),]

plt_bias = 
  dat %>% tibble %>%  
  dplyr::mutate(coef_idx = gsub( ".*bias_", "", key)) %>% 
  dplyr::mutate(coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  dplyr::mutate(coef_idx = plyr::revalue(coef_idx, c("0" = "Blip",
                                                     "1" = "Dissipation",
                                                     "2" = "Dose 0",
                                                     "3" = "Dose 1",
                                                     "4" = "Dose 2") )) %>%
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl,
                #coef_idx %in% 1:2 # blip and excursion effects
                ) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n, coef_idx) %>% 
  # dplyr::mutate(coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  ggplot(aes( y = value, x = coef_idx)) +
  facet_grid(t ~ n, scales = "free") +
  geom_boxplot(aes(color = factor(coef_idx)), lwd = 0.75, fatten = 0.5) + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(TeX('Relative Bias: $\\left (\\hat{\\beta} - \\beta \\right ) / \\beta$')) + 
  xlab("Coefficient Index") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1)),
         legend.position = "none") + 
  guides(color= guide_legend(title="Excursion Effect")) 
  

setwd(save_folder)
ggsave( "dtr_msm_DR_bias_excursion.pdf",
        plot = plt_bias,
        width = 12,
        height = 8)

# short version
plt_bias2 = 
  dat %>% tibble %>%  
  dplyr::mutate(coef_idx = gsub( ".*bias_", "", key)) %>% 
  dplyr::mutate(coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  dplyr::mutate(coef_idx = plyr::revalue(coef_idx, c("0" = "Blip",
                                                     "1" = "Dissipation",
                                                     "2" = "Dose 0",
                                                     "3" = "Dose 1",
                                                     "4" = "Dose 2") )) %>%
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl2,
                t %in% Tincl2,
  ) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n, coef_idx) %>% 
  # dplyr::mutate(coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  ggplot(aes( y = value, x = coef_idx)) +
  facet_grid(t ~ n, scales = "free") +
  geom_boxplot(aes(color = factor(coef_idx)), lwd = 0.75, fatten = 0.5) + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(TeX('Relative Bias: $\\left (\\hat{\\beta} - \\beta \\right ) / \\beta$')) + 
  xlab("Coefficient Index") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1)),
         legend.position = "none") + 
  guides(color= guide_legend(title="Excursion Effect")) 
########################################################
# combine plots
plt_comb <- ggpubr::ggarrange(plt_bias, plt_rmse, 
                              ncol=2, nrow=1, 
                              common.legend = TRUE, 
                              legend="bottom") 

setwd(save_folder)
ggsave( "dtr_msm_DR_combined.pdf",
        plot = plt_comb,
        width = 16,
        height = 6)

# short combine plots
plt_comb <- ggpubr::ggarrange(plt_bias2, plt_rmse2, 
                              ncol=2, nrow=1, 
                              common.legend = TRUE, 
                              legend="bottom") 

setwd(save_folder)
ggsave( "dtr_msm_DR_combined_short.pdf",
        plot = plt_comb,
        width = 12,
        height = 4)
#########################
# bias
#########################
dat <- do.call(rbind, ls_bias)
dat$t <- as.factor(dat$t)
dat$n <- as.factor(dat$n)

prefix <- "bias_"
dat <- dat[complete.cases(dat),]

plt_bias = 
  dat %>% tibble %>%  
  dplyr::mutate(coef_idx = gsub( ".*bias_", "", key)) %>% 
  dplyr::filter(value >= 0,
                value <= 1,
                n %in% nIncl) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(t, n) %>% 
  dplyr::mutate(coef_idx = as.factor(as.numeric(coef_idx)-1)) %>%
  ggplot(aes( y = value, x = t)) +
  # facet_grid(t ~ n, scales = "free") +
  facet_grid( ~ n, scales = "free") +
  geom_boxplot(aes(color = factor(coef_idx)), lwd = 0.75, fatten = 0.5) + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(TeX('Average Relative Bias: $1/p \\sum_j \\left (\\hat{\\beta_j} - \\beta_j \\right ) / \\beta_j$')) + 
  xlab("Number of Timepoints (T)") + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(face="bold", color="black", size = rel(1)),
         legend.position = "none")

setwd(save_folder)
ggsave( "dtr_msm_DR_bias_avg_excursion.pdf",
        plot = plt_bias,
        width = 12,
        height = 6)
