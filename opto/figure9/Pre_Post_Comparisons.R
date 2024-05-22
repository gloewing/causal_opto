# recreate Pre/Post analysis approach in Spontaneous DA paper

# Gabe Loewinger and Alex Levis 2-15-24
wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Data_processed/"
code_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/msm_hr/"
save_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/preProp_int/"
# learners
toml.list <- configr::read.config(file = "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoda_intermediate_results/closed_loop_learners.toml")
learner_ids <- toml.list[[1]]$learners

library(data.table)
library(dplyr)

control_ids <- read.csv(paste0(wd, "ctrl_ids.csv"))$mouse_id # control IDs
num_outcomes <- 6

# baseline responding
baseline <- data.table::fread(paste0(wd, "optoDA_pre.csv")) %>% 
  table.express::select(mouse_id, target, preCount) %>% 
  as_tibble() %>% unique()

# active session
dat <- data.table::fread(paste0(wd, "optoDA_active.csv")) %>% # does not use feedback_status criteria for Available (to be consistent with pre- data)
  table.express::mutate(opsin = 1*I(!mouse_id %in% control_ids)) %>% # treatment group only
  table.express::select(mouse_id, target, optoCount, opsin) %>% 
  as_tibble()

dat <- dplyr::left_join(dat, baseline, by = c("mouse_id", "target"))

targets <- unique(dat$target)[unique(dat$target) != -5] # target poses

# plot differences
myColors <- c("#ca0020", "#0868ac", "#525252", "#E69F00", "darkgreen", "darkgrey")

dat %>%
  dplyr::mutate(optoDiff = optoCount - preCount,
                Opto = ifelse(opsin == 1, "Opto", "Ctrl")) %>%
  dplyr::filter(target != -5) %>%
  ggplot(aes( y = optoDiff, x = factor(Opto))) +
  facet_grid( ~ factor(target)) +
  geom_boxplot(aes(color = factor(target))) + 
  # geom_point(aes(color = factor(coef_idx))) + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("Post - Pre Difference" )+ 
  xlab("Target Post") + 
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
  guides(color= guide_legend(title="Target"))   

# reformat to fit GLM
dat <- dat %>% 
  tidyr::pivot_longer(cols = c(optoCount, preCount), names_to = "prePost", values_to = "Y" ) %>%
  dplyr::mutate(Y = 2 * Y,
                prePost = ifelse(prePost == "preCount",0,1)) # turn average into total counts



# save data
coef_dim <- 4 # number of coefs for full factorial model
betaMat <- seMat <- matrix(NA, ncol = coef_dim, nrow = length(targets))
colnames(seMat) <- paste0("se_var_", 1:coef_dim)
colnames(betaMat) <- paste0("beta_var_", 1:coef_dim)

for(t in 1:length(targets)){
  # glm
  mod <- stats::glm(Y ~ opsin * prePost, 
                    data = dat[dat$target == targets[t],], 
                    family = poisson) 
  betaMat[t,] <- coef(mod)
  seMat[t,] <- as.numeric(sqrt(diag(sandwich::vcovCL(mod, type = "HC0", cluster = ~ mouse_id)))) # sandwich estimator
}


# process data for ggplot2
betas <- betaMat %>% 
  as_tibble() %>%
  dplyr::mutate(Outcome = 1:num_outcomes) %>%
  tidyr::pivot_longer(cols = starts_with("beta_var_"), 
                      names_to = "Coef", 
                      names_prefix = "beta_var_",
                      values_to = "beta") 

inference_table <- seMat %>% 
  as_tibble() %>%
  dplyr::mutate(Outcome = 1:num_outcomes) %>%
  tidyr::pivot_longer(cols = starts_with("se_var_"), 
                      names_to = "Coef", 
                      names_prefix = "se_var_",
                      values_to = "se") %>% 
  dplyr::right_join(betas, by = c("Outcome", "Coef")) %>%
  dplyr::filter(Coef == 4) # only need interaction

write.csv(inference_table, 
          "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/Pre-post/pre_post.csv",
          row.names = FALSE)

# treatment x opsin *Interaction*
myColors <- c("#ca0020", "darkgrey", "#0868ac", "#525252", "#E69F00", "darkgreen")

plt_int <-
  inference_table %>% 
  dplyr::rename(Pose = Outcome) %>% # for labeling
  ggplot(aes(x=Pose, y=beta, color=factor(Pose))) + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(latex2exp::TeX('Macro Longitudinal Treatment Effect:  $\\hat{\\beta}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:6) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) 
# 
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/Pre-post")
ggsave( "prePost.pdf",
        plot = plt_int,
        width = 6,
        height = 6)

rm(plt_int)


# ----------------------------------------------------------------
# -------------------------------
# GEE Table on their binned data (.parquet file)
# -------------------------------
# original data from authors
fl_nm <- "/Users/loewingergc/Desktop/Research/learning_timecourse_binsize-30.parquet"
dat_bin <- arrow::read_parquet(file = fl_nm) %>% 
  table.express::filter(session_number <= 2, # only pre-stim and stim days
                        rle == TRUE, # this is used to filter repeat rows
                        syllable == target_syllable) %>% # only consider counts for the target syllalbe
  dplyr::select(mouse_id, session_number, bin_start, target_syllable, count, area)

# control mice
ctrl_ids <- dat_bin %>% table.express::filter(area == "ctrl") %>% table.express::select(mouse_id) %>% unique()
ctrl_ids <- ctrl_ids$mouse_id

# learner mice
toml.list <- configr::read.config(file = "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoda_intermediate_results/closed_loop_learners.toml")
learner_ids <- toml.list[[1]]$learners

setDF(dat_bin)

baseline_avg <- dat_bin %>%
  table.express::filter(session_number %in% c(-1,0) ) %>%  # pre-stim sessions
  table.express::group_by(mouse_id, target_syllable, session_number) %>%
  dplyr::summarise(counts = sum(count)) %>% # calculate session-specific averages for each syllable, mouse
  table.express::group_by(mouse_id, target_syllable) %>%
  dplyr::summarise(counts_avg = mean(counts),
                   avg_baseline = mean(counts)) %>% # average across baseline sessions
  dplyr::mutate(time = 0)

opto_dat <- dat_bin %>%
  table.express::filter(session_number %in% 1:2) %>%  # pre-stim sessions
  table.express::group_by(mouse_id, target_syllable, session_number) %>%
  dplyr::summarise(counts = sum(count)) %>% # calculate session-specific averages for each syllable, mouse
  table.express::group_by(mouse_id, target_syllable) %>%
  dplyr::summarise(counts_avg = mean(counts),
                   avg_opto = mean(counts)) %>%
  dplyr::mutate(time = 1)

avg_table <- opto_dat %>%  # average across baseline sessions
  # join with baseline days and calculate difference
  dplyr::left_join(baseline_avg, by = c("mouse_id", "target_syllable")) %>%
  dplyr::mutate(optoDiff = avg_opto - avg_baseline,
                learner_orig = mouse_id %in% learner_ids,
                ctrl = I(mouse_id %in% ctrl_ids))

gee_dat <- rbind(baseline_avg, opto_dat) %>% 
  dplyr::mutate(ctrl = I(mouse_id %in% ctrl_ids),
                counts_avg = 2 * counts_avg) # turn avg counts into total counts

# iterate over targets
gee_dat$ctrl <- ifelse(gee_dat$ctrl, 0, 1)
gee_res <- matrix(NA, nrow = length(targets), ncol = 13)
colnames(gee_res) <- c("target", paste0("beta_", 0:3), paste0("se_", 0:3), paste0("signif_", 0:3))

for(t in 1:length(targets)){
  
  mod <- gee_dat %>%
    dplyr::filter(target_syllable == targets[t]) %>%
    #!mouse_id %in% paste0("dlight-chrimson-", 1:9)) %>% # ,
    glm(formula = counts_avg ~ ctrl*time, family = poisson) 
  gee_res[t,1] <- targets[t] # targets
  gee_res[t,2:5] <- coef(mod) # coefficients
  gee_res[t,6:9] <- as.numeric(sqrt(diag(sandwich::vcovCL(mod, type = "HC0", cluster = ~ mouse_id)))) # standard errors
  gee_res[t,10:13] <- I( (gee_res[t,2:5] + 1.96*gee_res[t,6:9]) *
                           (gee_res[t,2:5] - 1.96*gee_res[t,6:9]) > 0) # statistical significance
}


# plot their data's test

# process data for ggplot2
betas <- gee_res %>% 
  as_tibble() %>%
  dplyr::mutate(Outcome = target) %>%
  tidyr::pivot_longer(cols = starts_with("beta_"), 
                      names_to = "Coef", 
                      names_prefix = "beta_",
                      values_to = "beta") %>%
  dplyr::select(Outcome, Coef, beta)

inference_table <- gee_res %>% 
  as_tibble() %>%
  dplyr::mutate(Outcome = target) %>%
  tidyr::pivot_longer(cols = starts_with("se_"), 
                      names_to = "Coef", 
                      names_prefix = "se_",
                      values_to = "se") %>% 
  dplyr::select(Outcome, Coef, se) %>%
  dplyr::right_join(betas, by = c("Outcome", "Coef")) %>%
  dplyr::filter(Coef == 3) %>% # only need interaction
  dplyr::mutate(Outcome = 1:num_outcomes)

# treatment x opsin *Interaction*
myColors <- c("#ca0020", "darkgrey", "#0868ac", "#525252", "#E69F00", "darkgreen")

plt_int <-
  inference_table %>% 
  dplyr::rename(Pose = Outcome) %>% # for labeling
  ggplot(aes(x=Pose, y=beta, color=factor(Pose))) + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(latex2exp::TeX('Macro Longitudinal Treatment Effect (Binned):  $\\hat{\\beta}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:6) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1)),
         axis.title = element_text(face="bold", color="black", size=rel(1.2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) 
# 
# 
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/Pre-post")
ggsave( "prePost_bin.pdf",
        plot = plt_int,
        width = 6,
        height = 6)

