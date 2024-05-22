# Interaction with treatment and precision -- baseline responding
# HR-MSM for optoDA
# Gabe Loewinger and Alex Levis 1-25-24
# interaction with treatment/control

curium <- FALSE
if(curium){
  wd <- "/home/loewingergc/optoDA_preprocess/"
  code_wd <- "/home/loewingergc/msm_hr/"
  save_wd <- "/home/loewingergc/optoDA_results/"
  
  # learners
  toml.list <- configr::read.config(file = "/lscratch/SpontaneousBehaviour/optoda_intermediate_results/closed_loop_learners.toml")
  learner_ids <- toml.list[[1]]$learners
  
}else{
  wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Data_processed/"
  code_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/Figures/Final Figures/dose_treatment_precision/code/"
  save_wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/optoDA_results/Final/dose_treatment_precision/"
  # learners
  learner_ids <- read.csv(paste0(wd, "learners.csv"))$mouse_id 
}

delta <- 5
cov_type <- "HC0" # "HC3" is too computationally intensive
library(data.table)
library(dplyr)
control_ids <- read.csv(paste0(wd, "ctrl_ids.csv"))$mouse_id
myColors <- c("#ca0020", "#0868ac", "#525252", "#E69F00", "darkgreen", "darkgrey")

# read data
data.table::fread(paste0(wd, "optoDA_data.csv")) %>% 
  as_tibble() %>%
  dplyr::group_by(mouse_id, target) %>%
  dplyr::summarise(total_trials = n()) %>%
  dplyr::mutate(area = 1*I(!mouse_id %in% control_ids)) %>% # treatment/control
  dplyr::mutate(Opto = ifelse(area == 1, "Opto", "Ctrl")) %>%
  dplyr::filter(target != -5) %>%
  ggplot(aes( y = total_trials, x = factor(Opto))) +
  facet_grid( ~ factor(target)) +
  geom_boxplot(aes(color = factor(target))) + 
  # geom_point(aes(color = factor(coef_idx))) + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("Total Trials" )+ 
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
