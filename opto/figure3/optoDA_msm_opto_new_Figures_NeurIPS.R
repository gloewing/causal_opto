# Dose x opto/no-opsin Interaction
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

delta <- 3
num_methods <- 3 # methods to compare: "original", "conditional", "ours
kVec <- 1
num_outcomes <- 6 # number of "poses"/outcomes
# import data for methods, calculate statistical significance for main effects at specific dose

# read in results for "our methods"
method_number <- 3
wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/Final/dose_treatment_dissipate_3/"
wd_ctrl <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/Final/dose_treatment_dissipate_3_ctrl/" # control group as baseline version of same model
learner <- NA # all animals
main_effect_dims <- 2:4 # coefficient indices of main effects
main_effect_group <- 5 # coefficient indices of main effect of group (null)
int_effect_dims <- 6:8 # coefficient indices of interactions effects
example_dose <- 3 # show significance of our method at dose = 3
resList <- resList2 <- resList3 <- resList4 <- vector(length = length(kVec), "list")
resList_ctrl <- vector(length = length(kVec), "list") # for control group as baseline

for(k in 1:length(kVec)){
  d <- read.csv( paste0(wd,"optoDA_Int_lrn_", learner, "_Y_Y", kVec[k]) ) # read in data
  
  # -----------------------------------------------------------------
  # main effect of laser 
  betaMat <- d[paste0("beta2_var_", main_effect_dims)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se2_var_", main_effect_dims)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  
  # main effect of laser 
  betas <- betaMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("beta2_var_"), 
                         names_to = "Dose", 
                         names_prefix = "beta2_var_",
                         values_to = "beta") 
  
  resList[[k]] <- seMat %>% 
                    as_tibble() %>%
                    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se2_var_"), 
                                 names_to = "Dose", 
                                 names_prefix = "se2_var_",
                                 values_to = "se") %>% 
                    dplyr::right_join(betas, by = c("Outcome", "Dose")) %>%
                    dplyr::mutate(Dose = delta - (as.numeric(Dose) - min(as.numeric(Dose))),
                                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  
  # -----------------------------------------------------------------
  # interactions
  betaMat <- d[paste0("beta2_var_", int_effect_dims)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se2_var_", int_effect_dims)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  # 
  # # interaction effect of laser 
  betas <- betaMat %>%
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("beta2_var_"),
                         names_to = "Dose",
                         names_prefix = "beta2_var_",
                         values_to = "beta") 
  
  resList2[[k]] <- seMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se2_var_"), 
                       names_to = "Dose", 
                       names_prefix = "se2_var_",
                       values_to = "se") %>% 
    dplyr::right_join(betas, by = c("Outcome", "Dose"))%>%
    dplyr::mutate(Dose = delta - (as.numeric(Dose) - min(as.numeric(Dose))),
                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  
  # -----------------------------------------------------------------
  # main effect of group 
  betaMat <- d[paste0("beta2_var_", main_effect_group)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se2_var_", main_effect_group)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  
  # main effect of laser 
  betas <- betaMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("beta2_var_"), 
                        names_to = "Dose", 
                        names_prefix = "beta2_var_",
                        values_to = "beta") 
  
  resList4[[k]] <- seMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se2_var_"), 
                        names_to = "Dose", 
                        names_prefix = "se2_var_",
                        values_to = "se") %>% 
    dplyr::right_join(betas, by = c("Outcome", "Dose")) %>%
    dplyr::mutate(Dose = delta - (as.numeric(Dose) - min(as.numeric(Dose))),
                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  
  # -----------------------------------------------------------------
  # conditional estimand
  # -----------------------------------------------------------------
  d <- read.csv( paste0(wd,"optoDA_lrn_", learner, "_Y_Y", kVec[k], "_conditional") ) # read in data
  
  # main effect of laser 
  betaMat <- d[paste0("beta_var_", 4)] # subset dimensions of interactions
  seMat <-   d[paste0("se_var_", 4)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  
  # main effect of laser 
  betas <- betaMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("beta_var_"), 
                        names_to = "Dose", 
                        names_prefix = "beta_var_",
                        values_to = "beta") 
  
  resList3[[k]] <- seMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se_var_"), 
                        names_to = "Dose", 
                        names_prefix = "se_var_",
                        values_to = "se") %>% 
    dplyr::right_join(betas, by = c("Outcome", "Dose")) %>%
    dplyr::mutate(Dose = as.numeric(Dose) - 1,
                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  
  # control group as baseline data of same model
  d <- read.csv( paste0(wd_ctrl,"optoDA_Int_lrn_", learner, "_Y_Y", kVec[k]) ) # read in data
  
  # -----------------------------------------------------------------
  # main effect of laser 
  betaMat <- d[paste0("beta2_var_", main_effect_dims)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se2_var_", main_effect_dims)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  
  # main effect of laser 
  betas <- betaMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("beta2_var_"), 
                        names_to = "Dose", 
                        names_prefix = "beta2_var_",
                        values_to = "beta") 
  
  resList_ctrl[[k]] <- seMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se2_var_"), 
                        names_to = "Dose", 
                        names_prefix = "se2_var_",
                        values_to = "se") %>% 
    dplyr::right_join(betas, by = c("Outcome", "Dose")) %>%
    dplyr::mutate(Dose = delta - (as.numeric(Dose) - min(as.numeric(Dose))),
                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  }

# ---------------------------------------------------------------------------
# Trials prior figure (lag)
myColors <- c("#ca0020", "darkgrey", "#0868ac", "#525252", "#E69F00", "darkgreen")

plt <-
  do.call(rbind, resList) %>%
  dplyr::filter(k == 1) %>%
  ggplot(aes(x=Dose, y=beta, color=factor(Outcome))) + 
  facet_grid(~ as.factor(Outcome), scales = "free_y") + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(latex2exp::TeX('Treatment Effect (Opto):  $\\{\\widehat{\\beta}_r\\}^3_{r=1}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:5) +  
  scale_y_continuous(breaks = c(-0.03, 0, 0.03)) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         # axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) +
  labs(x = "Trials Prior", tag = "A")   

setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
ggsave( "dissipate_trt_Int_maineffects.pdf",
        plot = plt,
        width = 6,
        height = 6)

# ---------------------------------------------------------------------------
# dose-response *Interaction*

plt_int <-
  do.call(rbind, resList2) %>%
  dplyr::filter(k == 1) %>%
  dplyr::mutate(beta = - beta) %>% # reverse sign so this corresponds to control as baseline
  ggplot(aes(x=Dose, y=beta, color=factor(Outcome))) + 
  facet_grid(~ as.factor(Outcome), scales = "free_y") + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  # ylab(latex2exp::TeX('Treatment Group Effect Modification:  $\\{\\widehat{\\beta}_r\\}^7_{r=5}$') ) + 
  ylab(latex2exp::TeX('Group Effect Modification:  $\\{\\widehat{\\beta}_r\\}^7_{r=5}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:5) +  
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         #axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) +
  labs(x = "Trials Prior", tag = "A")
# 
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
ggsave( "dissipate_trt_Int.pdf",
        plot = plt_int,
        width = 6,
        height = 6)


# ---------------------------------------------------------------------------
# Group Effect (zero vec treatment)
plt_grp <-
  do.call(rbind, resList4) %>%
  dplyr::filter(k == 1) %>%
  dplyr::mutate(beta = -beta) %>% # reverse sign so this corresponds to control as baseline
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
  ylab(latex2exp::TeX('Opto Group (No Stim) Effect:  $\\widehat{\\beta}_4$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:6) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         # axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) +
  labs(tag = "C")
# 
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
ggsave( "dissipate_trt_group.pdf",
        plot = plt_grp,
        width = 6,
        height = 6)
# ---------------------------------------------------------------------------
# main effect of lag in control group data (same model but different baseline for trt group)
plt_ctrl <-
  do.call(rbind, resList_ctrl) %>%
  dplyr::filter(k == 1) %>%
  ggplot(aes(x=Dose, y=beta, color=factor(Outcome))) + 
  facet_grid(~ as.factor(Outcome), scales = "free_y") + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(latex2exp::TeX("Treatment Effect (Control):  $\\{\\widehat{\\beta}'_r\\}_{r=1}^3$") ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:5) +  
  scale_y_continuous(breaks = c(-0.03, 0, 0.03)) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         # axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         # axis.title.x = element_blank(),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) +
  labs(x = "Trials Prior", tag = "B")   

setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
# ggsave( "dissipate_trt_Int_maineffects.pdf",
#         plot = plt,
#         width = 6,
#         height = 6)

# ---------------------------------------------------------------------------


# combined plot
plt_comb <- ggpubr::ggarrange(plt, plt_int, 
                              ncol=2, nrow=1, 
                              common.legend = TRUE, 
                              #widths = c(2,1),
                              legend="bottom") 

setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
ggsave( "dissipate_trt_Int_combined.pdf",
        plot = plt_comb,
        width = 10,
        height = 6)


# -----------------------------------------------------------
# standard analysis (not binned) -- GEE
# -----------------------------------------------------------
# see code for production of this table:
# /Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/Pre-post/Pre_Post_Comparisons.R

inference_table <- read.csv("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/Pre-post/pre_post.csv")

# treatment x opsin *Interaction*
plt_gee <-
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
  ylab(latex2exp::TeX('Macro Longitudinal Effect:  $\\widehat{\\gamma}_3$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:6) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         # axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose")) +
  labs(tag = "D")


# -----------------------------------------------------------
# conditional estimand
# -----------------------------------------------------------
plt_cond <-
  do.call(rbind, resList3) %>%
  dplyr::filter(k == 1) %>%
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
  ylab(latex2exp::TeX('Conditional Treatment Effect:  $\\widehat{\\alpha}_3$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = 1:6) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(1), face="bold"),
         axis.text.x=element_text(face="bold",color="black", size=rel(1.5)),
         # axis.title = element_text(face="bold", color="black", size=rel(1.25)),
         axis.title = element_text(face="bold", color="black", size=rel(1.5)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(1.5)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(1.5)),
         strip.text.y = element_text(face="bold", color="black", size = rel(1)),
         strip.background.x = element_blank(),
         legend.position = "bottom",
         strip.text.x = element_blank()) + 
  guides(color= guide_legend(title="Pose"))  +
  labs(tag = "B")

# combine plots
plt_comb <- ggpubr::ggarrange(plt, plt_cond, 
                              plt_ctrl, plt_grp,
                              plt_int, plt_gee,
                              ncol=2, nrow=3, 
                              common.legend = TRUE, 
                              #widths = c(2,1),
                              legend="bottom") 
# 
# # save figures
# setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
# ggsave( "dissipate_combined.pdf",
#         plot = plt_comb,
#         width = 10,
#         height = 15)



# Same as before but combine plots and make horizontal
plt_comb <- ggpubr::ggarrange(plt_int, 
                              plt_cond, plt_grp, plt_gee,
                              ncol=4, nrow=1, 
                              common.legend = TRUE, 
                              #widths = c(2,1),
                              legend="bottom") 

# save figures
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_dissipate_3/Figures")
ggsave( "dissipate_combined_horizontal_small.pdf",
        plot = plt_comb,
        width = 20,
        height = 6)

# main effect in control group
# main effect of group

