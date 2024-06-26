# dose model (no interactions)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

delta <- 5
num_methods <- 3 # methods to compare: "original", "conditional", "ours
kVec <- c(1,5,10)
num_outcomes <- 6 # number of "poses"/outcomes
resMat <- matrix(NA, nrow = num_outcomes * length(kVec), ncol = num_methods)
colnames(resMat) <- c("Classic    ", "Conditional  ", "Ours  ")
rownames(resMat) <- rep(1:num_outcomes, times = length(kVec)) # label which pose number (for now use 1:num_poses)

# import data for methods, calculate statistical significance for main effects at specific dose

# read in results for "standard methods"
method_number <- 1
wd <- "/Users/folder/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/preprocessed_data/"
orig_stats <- read.csv(paste0(wd, "allAnimals_stats.csv")) %>%
                dplyr::mutate(signif = p_val <= 0.05) %>% # calculate significance
                dplyr::select(target, signif)
resMat[,method_number] <- orig_stats$signif # repeats since does not have k

# read in results for "our methods"
method_number <- 3
wd <- "/Users/folder/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/Final/Dose/" 
learner <- NA # all animals
main_effect_dims <- 2:4 # coefficient indices of main effects
#int_effect_dims <- (delta + 3):(2*delta+2) # coefficient indices of interactions effects
example_dose <- 3 # show significance of our method at dose = 3
resList <- resList2 <- resList3 <- vector(length = length(kVec), "list")

for(k in 1:length(kVec)){
  d <- read.csv( paste0(wd,"optoDA_Int_lrn_", learner, "_Y_Y", kVec[k]) ) # read in data
  
  # -----------------------------------------------------------------
  # main effect of laser 
  betaMat <- d[paste0("beta2_var_", main_effect_dims)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se2_var_", main_effect_dims)] # subset dimensions of SEs
  
  # calculate Wald CIs
  high_CI <- betaMat + 1.96 * seMat
  low_CI <- betaMat - 1.96 * seMat
  
  # calculate significance for each pose
  signif_mat <- sapply(seq_along(main_effect_dims), 
                  function(x) I(high_CI[,x] > 0 & low_CI[,x] > 0 | 
                                high_CI[,x] < 0 & low_CI[,x] < 0) 
                      )
  
  row_range <- seq( (k-1)*num_outcomes +1, k*num_outcomes)
  resMat[row_range, method_number] <- signif_mat[,example_dose] # show significance at specified dose

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
                    dplyr::mutate(Dose = as.numeric(Dose) - 1,
                                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  
  # -----------------------------------------------------------------
  
  resList2[[k]] <- seMat %>% 
    as_tibble() %>%
    dplyr::mutate(Outcome = 1:num_outcomes) %>%
    tidyr::pivot_longer(cols = starts_with("se2_var_"), 
                       names_to = "Dose", 
                       names_prefix = "se2_var_",
                       values_to = "se") %>% 
    dplyr::right_join(betas, by = c("Outcome", "Dose")) %>%
    dplyr::mutate(Dose = as.numeric(Dose) - (2+delta),
                  k = kVec[k]) # artifact of column naming scheme
  # -----------------------------------------------------------------
  # conditional estimand
  # -----------------------------------------------------------------
  
  d <- read.csv( paste0(wd,"optoDA_lrn_", learner, "_Y_Y", kVec[k], "_conditional") ) # read in data
  
  # main effect of laser 
  betaMat <- d[paste0("beta_var_", 2)] # subset dimensions of main effects for delta--assuming glm specified with dose first
  seMat <-   d[paste0("se_var_", 2)] # subset dimensions of SEs
  
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
  
  }

# ---------------------------------------------------------------------------
# dose-response figure
myColors <- c("#ca0020", "darkgrey", "#0868ac", "#525252", "#E69F00", "darkgreen")

plt_full <-
  do.call(rbind, resList) %>%
  ggplot(aes(x=Dose, y=beta, color=factor(Outcome))) + 
  facet_grid(as.factor(k) ~ as.factor(Outcome), scales = "free_y") + 
  geom_line() +
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin=beta-(1.96*se), ymax=beta+(1.96*se)), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab(latex2exp::TeX('Dose-Response Treatment Effects:  $\\{\\widehat{\\beta}_r\\}^1_{r=3}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = c(1,2,3)) + 
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

# ALL k values: dose response curve, Delta = 5 for HR-MSM 
setwd("/Users/folder/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_MSM/Figures")
ggsave( "dose_response_full.pdf",
        plot = plt_full,
        width = 8,
        height = 12)

# ---------------------------------------------------------------------------
# k = 1 only
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
  ylab(latex2exp::TeX('Dose-Response Treatment Effects:  $\\{\\widehat{\\beta}_r\\}^1_{r=3}$') ) + 
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1),
                     breaks = c(1,2,3)) +  
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
  guides(color= guide_legend(title="Pose")) +
  labs(tag = "A")  

# conditional estimand
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
  ylab(latex2exp::TeX('Conditional Treatment Effect:  $\\widehat{\\alpha}_1$') ) + 
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
  guides(color= guide_legend(title="Pose")) +
  labs(tag = "B") 

# combined plot
plt_comb <- ggpubr::ggarrange(plt, plt_cond, 
                              ncol=2, nrow=1, 
                              common.legend = TRUE, 
                              widths = c(2,1),
                              legend="bottom") 

# k = 1 dose response curve, Delta = 5 for HR-MSM
setwd("/Users/folder/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_MSM/Figures")
ggsave( "dose_response_k1.pdf",
        plot = plt_comb,
        width = 8,
        height = 6)

