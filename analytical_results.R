# Jess Kunke, 2023
# This code generates Figure 1 in the analytical results section as well as
# Figures 1-3 in the online supplement.

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

textbasesize = 18
pHLval = 0.001 #0.1, 0.05, 0.01
# The absolute path to the folder in which you wish to save the figure output
fig.dir = "/your/absolute/path/here/ending/with/slash/"

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
      Bias_MLE = R*(R*a+1-R)/(R+(1-R)*a) - R,
      Bias_PIMLE = R*(R*a+(1-R)/a) - R,
      Var_MLE = R*(R*pHL + (1-R)*pLL)^2 *(R*pHH *(1-pHH) + (1-R)*pHL*(1-pHL))/ (n*N*(R*pHL+(1-R)*pLL)^4) + R^2*(R*pHH + (1-R)*pHL)^2 *(R*pHL*(1-pHL) + (1-R)*pLL*(1-pLL))/(n*NK*(R*pHL + (1-R)*pLL)^4),
      Var_PIMLE = R*(R*pHH*(1-pHH) + (1-R)*pLL*(1-pHL))/(n*N*pHL^2) + R^2*(R*pHH^2*(1-pHL) + (1-R)*pHL^2*(1-pLL))/(n*NK*pHL^3)
    ) %>%
    mutate(
      Type = as.factor(ifelse(ER, "ER",
                              ifelse(Assortative, "Assort. non-ER",
                                     ifelse(Dissortative, "Dissortative",
                                            "Other")))),
      SE_MLE = sqrt(Var_MLE),
      SE_PIMLE = sqrt(Var_PIMLE),
      RMSE_MLE = sqrt(Bias_MLE^2 + Var_MLE),
      RMSE_PIMLE = sqrt(Bias_PIMLE^2 + Var_PIMLE)
    ) %>%
    mutate(
      Bias_MLE = ifelse(abs(Bias_MLE)<1e-12, 0, Bias_MLE),
      Bias_PIMLE = ifelse(abs(Bias_PIMLE)<1e-12, 0, Bias_PIMLE),
      Smaller_Bias = ifelse(abs(abs(Bias_PIMLE) - abs(Bias_MLE)) < 10^-8,
                            "Same", ifelse(abs(Bias_PIMLE) < abs(Bias_MLE),
                                           "dRpA", "dRpR")),
      Smaller_SE = ifelse(abs(SE_PIMLE - SE_MLE) < 10^-8,
                          "Same", ifelse(SE_PIMLE < SE_MLE, "dRpA", "dRpR")),
      Smaller_RMSE = ifelse(abs(RMSE_PIMLE-RMSE_MLE) < 10^-8,
                            "Same", ifelse(RMSE_PIMLE < RMSE_MLE,
                                           "dRpA", "dRpR")) #,
      # Sign_bias = ifelse(sign(Bias_MLE)==sign(Bias_PIMLE), "Same", "Different")
    ) %>%
    # # uncomment this code and comment out subsequent code if you wish to
    # # see these additional variables
    # mutate(
    #   Bias_Ratio = ifelse(Bias_MLE==0 & Bias_PIMLE==0, 1,
    #                       Bias_PIMLE/Bias_MLE),
    #   RMSE_Ratio = RMSE_PIMLE/RMSE_MLE
    # ) %>%
    # mutate(
    #   rH = R, # proportion of population in H (i.e. the prevalence R)
    #   rL = 1-R,
    #   EHdi = (R*N-1)*pHH + (1-R)*N*pHL, # expected degree in H
    #   ELdi = R*N*pHL + (N-R*N-1)*pLL # expected degree in L
    # ) %>%
    # mutate(Ed = rH*EHdi + rL*ELdi) %>% # expected degree
    # mutate(Edratio = EHdi/Ed) # expected degree ratio
    select(-c("ER":"Var_PIMLE", "RMSE_MLE":"RMSE_PIMLE")) %>%
    pivot_longer(cols = c("Smaller_Bias":"Smaller_RMSE"), names_prefix = "Smaller_",
                 names_to = "Metric", values_to = "Smaller") %>%
    mutate(
      # Metric = ifelse(Metric=="Bias", "1Bias",
      #                 ifelse(Metric=="SE", "2SE", "3RMSE")),
      # # Metric = as.factor(Metric), Smaller = as.factor(Smaller),
      # Case = as.factor(paste0(label, Metric))
      Smaller = factor(Smaller, levels = c("dRpR", "dRpA", "Same")),
      Case = factor(paste0(Metric, " (", label, ")"),
                    levels = c("Bias (full)", "SE (full)", "RMSE (full)",
                              "Bias (restricted)",
                              "SE (restricted)", "RMSE (restricted)"))
    )
  return(res)
}


# --------------------------------------------------------------
# "full" parameter range includes dissortative cases and R from 0.01 to 0.99
R = seq(0.01, 0.99, 0.01) #true prevalence of H
a = exp(seq(-4, 4, 0.1))
pHL = pHLval #between-group link probability
nN = 1000*c(5, 10, 50, 100, 500, 1000)
rK = c(0.01, 0.1, 0.2, 0.5, 0.8)

params.full = as.data.frame(expand_grid(R, a, pHL, nN, rK)) %>%
  mutate(pHH = a*pHL, pLL = a*pHL, n=50, N=nN/50) %>%
  mutate(NK = rK*N) %>%
  filter(pHH < 0.91, pLL < 0.91, NK < N)

full.data = compute_nsum_probeinL(params.full, "full")

# restricted range: just assortative and Erdos-Renyi cases with R 0.001-0.1
R = seq(0.001, 0.1, 0.001) #true prevalence of H
a = exp(seq(0, 4, 0.05)) # use -1 to 4? I keep only exp(seq(0, 4, 0.05)), but this ensures the factor levels are set up correctly for the variance
pHL = pHLval #between-group link probability
nN = 1000*c(5, 10, 50, 100, 500, 1000)
rK = c(0.01, 0.1, 0.2, 0.5, 0.8)

params.restr = as.data.frame(expand_grid(R, a, pHL, nN, rK)) %>%
  mutate(pHH = a*pHL, pLL = a*pHL, n=50, N=nN/50) %>%
  mutate(NK = rK*N) %>%
  filter(pHH < 0.91, pLL < 0.91, NK < N)

restr.data = compute_nsum_probeinL(params.restr, "restricted")

plot.data = rbind(full.data, restr.data) #%>%
  # filter(label == "full" | (label == "restricted" & a>=1))

palette = "Paired"

figure1 = ggplot(filter(plot.data, nN==500000, rK==0.1),
                 aes(x = R, y = log(a), color = Smaller)) +
  # facet_grid(label ~ Metric, scales="free") +
  facet_wrap(Case ~ ., scales="free") +
  rasterize(geom_point(), dpi=300, dev="ragg") +
  scale_color_brewer(palette = palette) +
  theme_minimal(base_size = textbasesize) +
  theme(legend.position="right")
print(figure1)
ggsave(paste0(fig.dir, "figure1_pHL", pHL, ".pdf"),
       width = 14, height = 8, units = "in")

# Supplement Figure 1
comp.rmse.by.nN = ggplot(filter(plot.data, Metric=="RMSE",
                                label == "restricted"),
                         aes(x = R, y = log(a), color = Smaller)) +
  rasterize(geom_point(), dpi=300, dev="ragg") +
  facet_wrap(~ nN, labeller = label_both) +
  scale_color_brewer(palette = palette) +
  theme_minimal(base_size = textbasesize) + guides(color=guide_legend(title="Smaller RMSE:")) +
  theme(panel.spacing.x = unit(2, "lines"))
print(comp.rmse.by.nN)
ggsave(paste0(fig.dir, "comp_RMSE_by_nN_pHL", pHL, ".pdf"),
       width = 14, height = 8, units = "in")

# Supplement Figure 2
comp.rmse.by.rK = ggplot(filter(plot.data, Metric=="RMSE",
                                label == "restricted"),
                         aes(x = R, y = log(a), color = Smaller)) +
  rasterize(geom_point(), dpi=300, dev="ragg") +
  facet_wrap(~ rK, labeller = label_both) +
  scale_color_brewer(palette = palette) +
  theme_minimal(base_size = textbasesize) + guides(color=guide_legend(title="Smaller RMSE:")) +
  theme(panel.spacing.x = unit(2, "lines"))
print(comp.rmse.by.rK)
ggsave(paste0(fig.dir, "comp_RMSE_by_rK_pHL", pHL, ".pdf"),
       width = 14, height = 8, units = "in")


