# This code generates Figure 1 in the analytical results section as well as
# Figures 1-3 in the online supplement, depending on how pHLval is set.

library(knitr)
library(ggpubr)
library(scales)
library(ggrastr)
library(formatR)
library(R.matlab)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

#### Specify inputs ####

# set pHLval to 0.01 for Figure 1 and Supplement Figures 1 and 3,
# set pHLval to 0.001 for Supplement Figure 2
pHLval = 0.01 #0.1, 0.05, 0.01
# textbasesize = 18
# width = 14 #inches
# height = 8 #inches
textbasesize = 8
width = 6 #inches
height = 3.4 #inches
# The absolute path to the folder in which you wish to save the figure output
fig.dir = "/your/absolute/path/here/ending/with/slash/"

# ----------- probably don't need to edit any code below this line ----------- #

#### Function to compute analytical approx. of bias and variance ####
# This function takes as input a list of parameter values and 
#  evaluates the approximate bias and variance expressions
#  presented in the main body of the paper
compute_nsum_probeinL = function(params, label){
  res = params %>%
    mutate(
      label = label,
      ER = (pHH == pHL & pLL == pHL),
      Assortative = (pHH > pHL & pLL > pHL),
      Dissortative = (pHH < pHL & pLL < pHL),
      Bias_RR = r*(r*a+1-r)/(r+(1-r)*a) - r,
      Bias_RA = r*(r*a+(1-r)/a) - r,
      Var_RR = r*(r*pHL + (1-r)*pLL)^2 *(r*pHH *(1-pHH) + (1-r)*pHL*(1-pHL))/ (n*N*(r*pHL+(1-r)*pLL)^4) + r^2*(r*pHH + (1-r)*pHL)^2 *(r*pHL*(1-pHL) + (1-r)*pLL*(1-pLL))/(n*NK*(r*pHL + (1-r)*pLL)^4),
      Var_RA = r*(r*pHH*(1-pHH) + (1-r)*pLL*(1-pHL))/(n*N*pHL^2) + r^2*(r*pHH^2*(1-pHL) + (1-r)*pHL^2*(1-pLL))/(n*NK*pHL^3)
    ) %>%
    mutate(
      Type = as.factor(ifelse(ER, "ER",
                              ifelse(Assortative, "Assort. non-ER",
                                     ifelse(Dissortative, "Dissortative",
                                            "Other")))),
      SE_RR = sqrt(Var_RR),
      SE_RA = sqrt(Var_RA),
      RMSE_RR = sqrt(Bias_RR^2 + Var_RR),
      RMSE_RA = sqrt(Bias_RA^2 + Var_RA)
    ) %>%
    mutate(
      Bias_RR = ifelse(abs(Bias_RR)<1e-12, 0, Bias_RR),
      Bias_RA = ifelse(abs(Bias_RA)<1e-12, 0, Bias_RA),
      Smaller_Bias = ifelse(abs(abs(Bias_RA) - abs(Bias_RR)) < 10^-8,
                            "Same", ifelse(abs(Bias_RA) < abs(Bias_RR),
                                           "RA", "RR")),
      Smaller_SE = ifelse(abs(SE_RA - SE_RR) < 10^-8,
                          "Same", ifelse(SE_RA < SE_RR, "RA", "RR")),
      Smaller_RMSE = ifelse(abs(RMSE_RA-RMSE_RR) < 10^-8,
                            "Same", ifelse(RMSE_RA < RMSE_RR,
                                           "RA", "RR")) #,
      # Sign_bias = ifelse(sign(Bias_RR)==sign(Bias_RA), "Same", "Different")
    ) %>%
    # # uncomment this code and comment out subsequent code if you wish to
    # # see these additional variables
    # mutate(
    #   Bias_Ratio = ifelse(Bias_RR==0 & Bias_RA==0, 1,
    #                       Bias_RA/Bias_RR),
    #   RMSE_Ratio = RMSE_RA/RMSE_RR
    # ) %>%
    # mutate(
    #   rH = r, # proportion of population in H (i.e. the prevalence r)
    #   rL = 1-r,
    #   EHdi = (r*N-1)*pHH + (1-r)*N*pHL, # expected degree in H
    #   ELdi = r*N*pHL + (N-r*N-1)*pLL # expected degree in L
    # ) %>%
    # mutate(Ed = rH*EHdi + rL*ELdi) %>% # expected degree
    # mutate(Edratio = EHdi/Ed) # expected degree ratio
    select(-c("ER":"Var_RA", "RMSE_RR":"RMSE_RA")) %>%
    pivot_longer(cols = c("Smaller_Bias":"Smaller_RMSE"), names_prefix = "Smaller_",
                 names_to = "Metric", values_to = "Smaller") %>%
    mutate(
      # Metric = ifelse(Metric=="Bias", "1Bias",
      #                 ifelse(Metric=="SE", "2SE", "3RMSE")),
      # # Metric = as.factor(Metric), Smaller = as.factor(Smaller),
      # Case = as.factor(paste0(label, Metric))
      Smaller = factor(Smaller, levels = c("RR", "RA", "Same")),
      Case = factor(paste0(Metric, " (", label, ")"),
                    levels = c("Bias (full)", "SE (full)", "RMSE (full)",
                              "Bias (restricted)",
                              "SE (restricted)", "RMSE (restricted)"))
    )
  return(res)
}


#### Compute full and restricted results ####
# "full" parameter range includes dissortative cases and r from 0.01 to 0.99
r = seq(0.01, 0.99, 0.01) #true prevalence of H
a = exp(seq(-4, 4, 0.1))
pHL = pHLval #between-group link probability
nN = 1000*c(5, 10, 50, 100, 500, 1000)
rK = c(0.01, 0.1, 0.2, 0.5, 0.8)

params.full = as.data.frame(expand_grid(r, a, pHL, nN, rK)) %>%
  mutate(pHH = a*pHL, pLL = a*pHL, n=50, N=nN/50) %>%
  mutate(NK = rK*N) %>%
  filter(pHH < 0.91, pLL < 0.91, NK < N)

full.data = compute_nsum_probeinL(params.full, "full")

# restricted range: just assortative and Erdos-Renyi cases with r 0.001-0.1
r = seq(0.001, 0.1, 0.001) #true prevalence of H
a = exp(seq(0, 4, 0.05)) # use -1 to 4? I keep only exp(seq(0, 4, 0.05)), but this ensures the factor levels are set up correctly for the variance
pHL = pHLval #between-group link probability
nN = 1000*c(5, 10, 50, 100, 500, 1000)
rK = c(0.01, 0.1, 0.2, 0.5, 0.8)

params.restr = as.data.frame(expand_grid(r, a, pHL, nN, rK)) %>%
  mutate(pHH = a*pHL, pLL = a*pHL, n=50, N=nN/50) %>%
  mutate(NK = rK*N) %>%
  filter(pHH < 0.91, pLL < 0.91, NK < N)

restr.data = compute_nsum_probeinL(params.restr, "restricted")

plot.data = rbind(full.data, restr.data) #%>%
  # filter(label == "full" | (label == "restricted" & a>=1))

#### Generate Figure 1 and supplemental plots #####

palette = c('#f0f0f0','#bdbdbd','#636363')

# Main paper, Figure 1 (make sure pHLval = 0.01)
figure1 = ggplot(filter(plot.data, nN==500000, rK==0.1),
                 aes(x = r, y = log(a), color = Smaller)) +
  facet_wrap(Case ~ ., scales="free") +
  scale_color_manual(values = palette) +
  rasterize(geom_point(cex=0.5), dpi=300, dev="ragg") +
  theme_minimal(base_size = textbasesize) +
  theme(legend.position="right", legend.spacing.y = unit(-0.2, 'cm')) +
  guides(color = guide_legend(title="Smaller:\n", byrow = TRUE))
print(figure1)
ggsave(paste0(fig.dir, "figure1_pHL", pHL, ".pdf"),
       width = width, height = height, units = "in")

# varying nN and pHL:
# Supplement Figure 1 (if pHLval = 0.01) or Supplement Figure 2 (if pHLval = 0.001)
comp.rmse.by.nN = ggplot(filter(plot.data, Metric=="RMSE",
                                label == "restricted"),
                         aes(x = r, y = log(a), color = Smaller)) +
  rasterize(geom_point(cex=0.5), dpi=300, dev="ragg") +
  facet_wrap(~ nN, labeller = label_both) +
  scale_color_manual(values = palette) +
  theme_minimal(base_size = textbasesize) +
  theme(panel.spacing.x = unit(1, "lines"), legend.spacing.y = unit(-0.2, 'cm')) +
  guides(color = guide_legend(title="Smaller RMSE:\n", byrow = TRUE))
print(comp.rmse.by.nN)
ggsave(paste0(fig.dir, "comp_RMSE_by_nN_pHL", pHL, ".pdf"),
       width = width, height = height, units = "in")

# Supplement Figure 3 (make sure pHLval = 0.01)
comp.rmse.by.rK = ggplot(filter(plot.data, Metric=="RMSE",
                                label == "restricted"),
                         aes(x = r, y = log(a), color = Smaller)) +
  rasterize(geom_point(cex=0.5), dpi=300, dev="ragg") +
  facet_wrap(~ rK, labeller = label_both) +
  scale_color_manual(values = palette) +
  theme_minimal(base_size = textbasesize) +
  theme(panel.spacing.x = unit(1, "lines"), legend.spacing.y = unit(-0.2, 'cm')) +
  guides(color = guide_legend(title="Smaller RMSE:\n", byrow = TRUE))
print(comp.rmse.by.rK)
ggsave(paste0(fig.dir, "comp_RMSE_by_rK_pHL", pHL, ".pdf"),
       width = width, height = height, units = "in")


