# Jess Kunke, 2023
# This code generates various figures for the Facebook 100 simulations,
# including Figures 2-3 from the main paper and Figures 4-5 from the supplement

# If you haven't already generated the following data files:
#  - dratios.RData
#  - fb100_sims.RData
# please run the script fb100_simulations_run_once.R before you run this script.

library(knitr)
library(ggpubr)
library(scales)
library(R.matlab)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

# The absolute path to the facebook100 data folder
data.dir = "/your/absolute/path/here/ending/with/slash/"
# The absolute path to the folder containing dratios.RData and fb100_sims.RData
out.dir = "/your/absolute/path/here/ending/with/slash/"
# The absolute path to the folder where you wish to store figure output from
#  this script
fig.dir = "/your/absolute/path/here/ending/with/slash/"

# ---------------------------------------------------------------
# For all 1600 cases:
# - school name, id, and size (total population size)
# - hidden group name, id, and size
# - degree ratio for this choice of school and hidden group
load(paste0(out.dir, "dratios.RData"))

# summary(dratios$assort)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.09012  0.05237  0.18602  0.22921  0.34959  0.87839

school_inds = dratios$school_ind
hidden_inds = dratios$hidden_ind

load(paste0(out.dir, "fb100_sims.RData"))

# Two schools are called Texas and three are called UC; disambiguate them
# Texas: 80, 84 (appear as 63 and 64 in my list)
# UC: 33, 61, 64 (appear as 68-70 in my list)
# To verify:
# > i=63; unique(survey.data$School[(8000*(i-1)+1):(8000*i)])
# [1] "Texas"
# > i=64; unique(survey.data$School[(8000*(i-1)+1):(8000*i)])
# [1] "Texas"
# > i=68; unique(survey.data$School[(8000*(i-1)+1):(8000*i)])
# [1] "UC"
# > i=69; unique(survey.data$School[(8000*(i-1)+1):(8000*i)])
# [1] "UC"
# > i=70; unique(survey.data$School[(8000*(i-1)+1):(8000*i)])
dratios$school[(16*63+1):(16*64)] = "Texas2"
dratios$school[(16*68+1):(16*69)] = "UC2"
dratios$school[(16*69+1):(16*70)] = "UC3"

school.hidden.data$school[(16*63+1):(16*64)] = "Texas2"
school.hidden.data$school[(16*68+1):(16*69)] = "UC2"
school.hidden.data$school[(16*69+1):(16*70)] = "UC3"

survey.data$school[(8000*63+1):(8000*64)] = "Texas2"
survey.data$school[(8000*68+1):(8000*69)] = "UC2"
survey.data$school[(8000*69+1):(8000*70)] = "UC3"



# ---------------------------------------------------------
# Prep the data for making figures and tables

fb100.data = data.frame(
  survey.data %>%
  left_join(select(dratios, school, hidden, N, d.ratio, assort),
            by = c("school", "hidden")) %>%
  mutate(R.true = true.size/N,
         across(.cols = mle.size.est:mos.cens.est, .fns = ~.x / N),
         dratio.level = as.factor(ifelse(d.ratio < 0.8, "Low",
                                         ifelse(d.ratio > 1.2, "High",
                                                "Near 1"))),
         assort.level = as.factor(ifelse(assort>0, "Assortative", "Dissortative"))) %>%
  rename(MLE = mle.size.est, MLEknown = mle.size.est.known.d,
         PIMLE = pimle.size.est, PIMLEcens = pimle.cens.est,
         PIMLEknown = pimle.size.dtrue,
         MoS = mos.size.est, MoScens = mos.cens.est) %>%
  gather(key="Estimator", value="R.est", MLE:MoScens, factor_key=TRUE) %>%
  select(-c("true.size", "n.hidden")) %>%
  mutate(R.est = ifelse(is.infinite(R.est) | is.nan(R.est), NA, R.est),
         dratio.level = relevel(dratio.level, "Low")) %>%
  group_by(across(c(-R.est))) %>%
  # summarize(Nsim = n(), EffNsim = sum(!is.na(size)), # effective no. of sims, excluding NaNs and Infs which were set to NA
  #           Mean = mean(size, na.rm=TRUE), SD = sd(size, na.rm=TRUE)) %>%
  summarize(mean = mean(R.est, na.rm=TRUE), sd = sd(R.est, na.rm=TRUE)) %>%
  mutate(Bias = mean - R.true,
         AbsBias = abs(mean - R.true),
         RMSE = sqrt((mean - R.true)^2 + sd^2)) %>%
  mutate(BiasFrac = ifelse(abs(Bias)<1e-12, 0, Bias)/R.true,
         AbsBiasFrac = abs(Bias/R.true),
         SEFrac = sd/R.true,
         RMSEFrac = RMSE/R.true) %>%
  select(-c(mean)) %>%
  pivot_longer(cols=sd:RMSEFrac, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Estimator, values_from = Value) %>%
  mutate(smaller.PIMLEvMLE = ifelse(abs(abs(PIMLE) - abs(MLE)) < 10^-8,
                                    "Same", ifelse(abs(PIMLE) < abs(MLE),
                                                   "PIMLE", "MLE")),
         smaller.PIMLEvMoS = ifelse(abs(abs(PIMLE) - abs(MoS)) < 10^-8,
                                    "Same", ifelse(abs(PIMLE) < abs(MoS),
                                                   "PIMLE", "MoS")),
         smaller.MoSvMLE = ifelse(abs(abs(MoS) - abs(MLE)) < 10^-8,
                                    "Same", ifelse(abs(MoS) < abs(MLE),
                                                   "MoS", "MLE")),
         Metric = as.factor(Metric),
         PercDiff.PIMLE.MLE = abs(PIMLE-MLE),
         PercDiff.PIMLE.MoS = abs(PIMLE-MoS),
         PercDiff.MoS.MLE = abs(MoS-MLE)
  ) %>%
  rename(dRpR = MLE, dRpA = PIMLE, dApA = MoS)
)


# -------------------------------------------------
# Generate tables summarizing what proportion of the 
# time the PIMLE has lower bias, variance, and RMSE
# than the MLE and MoS
# -------------------------------------------------

summary(fb100.data$assort)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.09012  0.05237  0.18602  0.22921  0.34959  0.87839 
summary(filter(fb100.data, dratio.level=="Low")$assort)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03875 0.19495 0.40531 0.46113 0.76529 0.87839 
summary(filter(fb100.data, dratio.level=="Near 1")$assort)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.08150  0.04188  0.10132  0.17349  0.28948  0.87114 
summary(filter(fb100.data, dratio.level=="High")$assort)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.09012  0.06697  0.26466  0.25336  0.39086  0.63984 

summary(filter(fb100.data, Metric=="AbsBiasFrac", dratio.level=="Low")$dRpR)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2193  0.3061  0.3813  0.4049  0.4766  0.8136 
summary(filter(fb100.data, Metric=="AbsBiasFrac", dratio.level=="Low")$dRpA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00124 0.04175 0.08800 0.10237 0.13426 0.56095 


# counts for PIMLEvsMLE, for example, represent the number of cases
#   in which PIMLE < MLE for that metric (abs(bias), sd, or RMSE)
fb100.counts = fb100.data %>% 
  filter(Metric %in% c("AbsBias", "sd", "RMSE")) %>%
  group_by(dratio.level, Metric, assort.level) %>% 
  summarize(n = n(), PIMLEvsMLE = sum(dRpA<dRpR, na.rm = TRUE),
            PIMLEvsMoS = sum(dRpA<dApA, na.rm = TRUE),
            MoSvsMLE = sum(dApA<dRpR, na.rm = TRUE)
  )

fb100.summ = fb100.data %>% 
  filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac"), !is.na(smaller.PIMLEvMLE)) %>%
  group_by(dratio.level, Metric, smaller.PIMLEvMLE) %>% 
  summarize(n = n(),
            assort.mean = mean(assort),
            min.PIMLE.MLE = min(PercDiff.PIMLE.MLE, na.rm=TRUE),
            median.PIMLE.MLE = median(PercDiff.PIMLE.MLE, na.rm=TRUE),
            max.PIMLE.MLE = max(PercDiff.PIMLE.MLE, na.rm=TRUE)
  )
kable(fb100.summ, digits=3)




### Generate the figures
textbasesize = 18

gen.fb100.fig = function(data, prefix){
  plotdata = data %>%
    filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac")) %>%
    mutate(Metric = as.factor(ifelse(Metric=="AbsBiasFrac", "|Bias|/R",
                                     ifelse(Metric=="SEFrac", "SE/R", "RMSE/R")))) %>%
    mutate(Metric = relevel(relevel(Metric, "SE/R"), "|Bias|/R"),
           `Degree ratio` = d.ratio)
  
  if(prefix=="low_dratio_"){
    plotdata = filter(plotdata, d.ratio <= 0.8)
  }else if(prefix=="high_dratio_"){
    plotdata = filter(plotdata, d.ratio >= 1.2)
  }else if(prefix=="near1_dratio_"){
    plotdata = filter(plotdata, d.ratio > 0.8, d.ratio < 1.2)
  }
  
  pimle.mle = ggplot(plotdata, aes(x=dRpA, y=dRpR, color=`Degree ratio`)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    xlim(c(0, 1.2)) +
    ylim(c(0, 1.2)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1)
  print(pimle.mle)
  ggsave(paste0(fig.dir, prefix, "PIMLE_MLE.pdf"),
         width = 16, height = 5, units = "in")
  
  pimle.mos = ggplot(plotdata, aes(x=dRpA, y=dApA, color=`Degree ratio`)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    xlim(c(0, 1.2)) +
    ylim(c(0, 1.2)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1)
  print(pimle.mos)
  ggsave(paste0(fig.dir, prefix, "PIMLE_MoS.pdf"),
         width = 16, height = 5, units = "in")
  
  mos.mle = ggplot(plotdata, aes(x=dApA, y=dRpR, color=`Degree ratio`)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    xlim(c(0, 1.2)) +
    ylim(c(0, 1.2)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "bottom", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, legend.key.width = unit(1, 'in'))
  print(mos.mle)
  ggsave(paste0(fig.dir, prefix, "MoS_MLE.pdf"),
         width = 16, height = 7, units = "in")
  
  ggarrange(pimle.mle, pimle.mos, mos.mle,
            ncol=1, nrow=3, heights = c(5, 5, 6))
  ggsave(paste0(fig.dir, prefix, "combined.pdf"),
         width = 16, height = 17, units = "in")
}

# All cases (Figure 3 is the combined figure)
gen.fb100.fig(fb100.data, prefix="all_")
# Cases with low degree ratio (Figure 2 is the PIMLE_MLE figure)
gen.fb100.fig(fb100.data, prefix="low_dratio_")
# Cases with high degree ratio (Supplement Figure 5 is the PIMLE_MLE figure)
gen.fb100.fig(fb100.data, prefix="high_dratio_")
# Cases with degree ratio near 1 (Supplement Figure 4 is the PIMLE_MLE figure)
gen.fb100.fig(fb100.data, prefix="near1_dratio_")




