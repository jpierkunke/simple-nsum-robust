# This code analyzes the results of the simulations and generates various 
# figures for the Facebook 100 simulations, including the figures from the main 
# paper and online supplement.

# If you haven't already generated the fb100_sims.RData file with simulation
# results, please run the script simulate_surveys_on_fb100_networks.R before
# running this script.

library(knitr)
library(ggpubr)
library(scales)
library(R.matlab)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

#### Specify inputs ####

set.seed(145)

# The absolute path to the facebook100 data folder
data.dir = "/your/absolute/path/here/ending/with/slash/"

# The absolute path to the folder where you wish to store fb100_sims.RData
# once this script creates it
out.dir = "/your/absolute/path/here/ending/with/slash/"

# The absolute path to the folder where you wish to store figures created by
#  this script
fig.dir = "/your/absolute/path/here/ending/with/slash/"


#### Load simulation results ####

# load probe and hidden group choices
load(paste0(data.dir, "all_probe_hidden_groups.RData"))

# load new sim results
load(paste0(out.dir, "school1.RData"))
case.data.new = case.data
individual.data.new = individual.data
survey.data.new = survey.data

for(i in 2:99){
  load(paste0(out.dir, "school", i, ".RData"))
  case.data.new = bind_rows(case.data.new, case.data)
  individual.data.new = bind_rows(individual.data.new, individual.data)
  survey.data.new = bind_rows(survey.data.new, survey.data)
}

load(paste0(out.dir, "school100.RData"))
case.data = bind_rows(case.data.new, case.data)
individual.data = bind_rows(individual.data.new, individual.data)
survey.data = bind_rows(survey.data.new, survey.data)
rm(case.data.new, individual.data.new, survey.data.new)

# Two schools are called Texas and three are called UC; disambiguate them
# unique(select(case.data, school, school.id))
# - Texas: 80, 84 (school.ids 63-64)
# - UC: 33, 61, 64 (school.ids 68-70)
case.data$school[case.data$school.id==64] = "Texas2"
case.data$school[case.data$school.id==69] = "UC2"
case.data$school[case.data$school.id==70] = "UC3"

names(case.data)
# [1] "case.id"   "school"    "school.id" "N"         "HH.edges"  "LL.edges"  "HL.edges"  "hidden"   
# [9] "hidden.id" "NH"        "r.true"    "d.ratio"   "a"         "assort"   

# to see results for all cases, not only the assortative ones,
# comment the following code between <----->
# <----->
assort.grps = hidden.grps %>%
  group_by(school.id) %>%
  slice_max(assort, n=10, with_ties = FALSE) %>%
  select(school.id, grp.name) %>%
  rename(hidden = grp.name)

case.data = assort.grps %>%
  left_join(case.data)
#now has 999 cases

individual.data = individual.data %>%
  filter(case.id %in% case.data$case.id)
# 12,082,391 rows

survey.data = survey.data %>%
  filter(case.id %in% case.data$case.id)
# 499,500 rows

# <----->

ncases = dim(case.data)[1]
ncases #1989 for all cases, 999 for assortative only
sum(case.data$r.true<0.1)/ncases
# 1


#### Compare and evaluate degree estimates ####

summary(individual.data$d.true)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   24.00   57.00   77.77  108.00 8246.00
summary(individual.data$d.R.est)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   22.72   59.18   85.93  119.65 8910.61 
summary(individual.data$d.A.est)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   22.64   59.02   85.70  119.57 9069.84  


# what this is telling us is that while both are roughly centered 
# around the true values, the A estimates have greater variance
sampids = sample(1:length(individual.data$case.id), 10000, replace=FALSE)
plot(individual.data$d.true[sampids], individual.data$d.R.est[sampids])
plot(individual.data$d.true[sampids], individual.data$d.A.est[sampids])
plot(individual.data$d.R.est[sampids], individual.data$d.A.est[sampids])

summary(abs(individual.data$d.R.est-individual.data$d.true)/individual.data$d.true)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.09103 0.19355 0.25924 0.33944 5.92845 
summary(abs(individual.data$d.A.est-individual.data$d.true)/individual.data$d.true)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.09297 0.19808 0.27306 0.35500 8.61467


# is it true that d.R.est = 0 if and only if d.R.est = 0? yes.
sum(individual.data$d.R.est == 0) == sum(individual.data$d.A.est == 0); sum(individual.data$d.R.est == 0) == sum(individual.data$d.R.est == 0 & individual.data$d.A.est == 0)
# how many estimates are 0? 274,140
sum(individual.data$d.R.est == 0)
# 274118
# what percent of estimates are 0s? about 2.3%
100*sum(individual.data$d.R.est == 0)/dim(individual.data)[1]
# 2.26874
# what are the true degrees like for people whose estimated degrees are zero?
summary(individual.data$d.true[individual.data$d.R.est == 0])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   2.681   3.000  63.000 

# what percent of y/d ratios are > 1 when the estimated degree is positive? 0.14%
# nonsense.yd = ((individual.data$y.iH/individual.data$d.R.est) > 1 & individual.data$d.R.est>0)
nonsense.yd = ((individual.data$y.iH > individual.data$d.R.est) & (individual.data$d.R.est>0))
100*sum(nonsense.yd)/dim(individual.data)[1]
# 0.1378287
hist(individual.data$y.iH[nonsense.yd]/individual.data$d.true[nonsense.yd])
summary(individual.data$y.iH[nonsense.yd]/individual.data$d.true[nonsense.yd])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02857 0.50000 0.66667 0.65711 0.81481 1.00000 

# indeed, people with y/d ratios > 1 tend to know many more people in the HTR group than average...
summary(individual.data$y.iH[nonsense.yd])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    6.00   15.00   22.16   29.00  625.00 
summary(individual.data$y.iH)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   1.247   1.000 625.000 

#... even though their degrees are smaller
summary(individual.data$d.true[nonsense.yd])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   11.00   23.00   35.66   47.00  664.00 
summary(individual.data$d.true)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   24.00   57.00   77.77  108.00 8246.00 


comp.deg.est = individual.data %>%
  # compute degree estimate percent error
  mutate(R.perc.err = (d.R.est-d.true)/d.true,
         A.perc.err = (d.A.est-d.true)/d.true) %>%
  # which estimator has the smaller percent error?
  # if difference is less than 0.5%, say they're roughly the same percent error
  mutate(smaller.perc.err = ifelse(abs(abs(R.perc.err)-abs(A.perc.err))<0.005, "Same",
                                   ifelse(abs(R.perc.err)<abs(A.perc.err), "R", "A")))

kable(table(comp.deg.est$smaller.perc.err)/nrow(comp.deg.est),
      caption = "Which degree estimator has smaller percent error?")
# |Var1 |      Freq|
# |:----|---------:|
# |A    | 0.4223141|
# |R    | 0.4984518|
# |Same | 0.0792342|
# R estimator tends to be better, but the R and A estimators are comparable

cor(individual.data$d.true, individual.data$d.R.est)
# [1] 0.9761899
cor(individual.data$d.true, individual.data$d.A.est)
# [1] 0.9742361
cor(individual.data$d.R.est, individual.data$d.A.est)
# [1] 0.9970101


fb100.data = data.frame(
  case.data %>%
    mutate(
      # Assign low, high, and near-1 labels to degree ratios
      dratio.level = as.factor(ifelse(d.ratio < 0.8, "Low",
                                      ifelse(d.ratio > 1.2, "High",
                                             "Near 1"))),
      # Assign labels to each case based on assortativity if a mix of cases are used
      assort.level = as.factor(ifelse(assort>0, "Assortative", "Dissortative"))
    ) %>%
    mutate(dratio.level = relevel(dratio.level, "Low")) %>%
    # join survey.data and case.data
    right_join(select(survey.data, -(survey.id:n.yd.gt1A)),
               by = c("case.id")) %>%
    # Pivot longer; after this step, nrow(fb100.data) should equal
    #   nrow(survey.data) times # estimators
    pivot_longer(cols = RR.r.est:AA.r.est4, names_to="Estimator", values_to="r.est") %>%
    group_by(across(c(-r.est))) %>%
    # compute mean and standard deviation of estimates across surveys for each estimator;
    # these are estimates of the mean and standard error of the estimators
    # after this step, nrow(fb100.data) should equal # cases times # estimators
    summarize(Mean = mean(r.est), SE = sd(r.est)) %>%
    # compute estimated Bias, RMSE
    mutate(Bias = Mean - r.true,
           AbsBias = abs(Mean - r.true),
           RMSE = sqrt((Mean - r.true)^2 + SE^2)) %>%
    # standardize bias, SE, RMSE by the true prevalence
    mutate(BiasFrac = ifelse(abs(Bias)<1e-12, 0, Bias)/r.true,
           AbsBiasFrac = abs(Bias/r.true),
           SEFrac = SE/r.true,
           RMSEFrac = RMSE/r.true) %>%
    select(-c(Mean)) %>%
    # at this point, # obs = # cases * # estimators
    # next step increases number of rows by factor of # performance metrics
    pivot_longer(cols=SE:RMSEFrac, names_to = "Metric", values_to = "Value") %>%
    # reduces number of rows by factor of 10 (# estimators)
    pivot_wider(names_from = Estimator, values_from = Value)
)

# look at the number/percent of zero-valued degrees by case
check.no0deg = (individual.data %>% group_by(case.id) %>% summarize(n = sum(d.R.est==0)))
no0deg = check.no0deg$case.id[check.no0deg$n == 0] # all cases have some zero-valued degree estimates
check.prop.0deg = individual.data %>% group_by(case.id) %>% summarize(r = mean(d.R.est==0))
summary(check.prop.0deg$r) #% of network with zero-valued est. degrees ranges from 0.4-5.3%
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004149 0.013536 0.020027 0.020666 0.026327 0.052506  

fb100.data = fb100.data %>%
  left_join(check.prop.0deg, by="case.id") %>%
  rename(r.0deg = r) %>%
  mutate(r.0deg.level = ifelse(r.0deg < summary(check.prop.0deg$r)[2], "1st Qu.",
                               ifelse(r.0deg < summary(check.prop.0deg$r)[3], "2nd Qu.",
                                      ifelse(r.0deg < summary(check.prop.0deg$r)[5], "3rd Qu.", "4th Qu."))))

summary(fb100.data$a)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1.579    7.253   17.148   28.861   36.264 1539.635 
hist(fb100.data$a)
hist(fb100.data$a[fb100.data$a<300])




#### Examine sensitivity of results to methods for handling NAs and Infs ####

# Check % diff in RMSE based on method for handling NAs and Infs
# What is the % difference in RMSE of the AA and RA estimators for each of
# the four ways to handle NAs and Infs?
# - for each case, compute the percent difference between the largest and
#   smallest RMSEs of the RA or AA estimator

# RA estimator:
fb100.data %>%
  select(-(AA.r.est1:AR.r.est), -RR.r.est, -(school:assort.level)) %>%
  pivot_longer(cols=starts_with("RA"), names_to = "Inf.NA.method", values_to = "Value") %>%
  filter(Metric =="RMSEFrac") %>%
  group_by(case.id, Metric) %>%
  summarize(
    max.diff = max(Value) - min(Value),
    sd = sd(Value),
    sd.frac = sd(Value)/max(Value)
  ) %>%
  pivot_wider(names_from = Metric, values_from = max.diff) %>%
  summary() 


# AA estimator:
fb100.data %>%
  select(-(AR.r.est:RR.r.est), -(school:assort.level)) %>%
  pivot_longer(cols=starts_with("AA"), names_to = "Inf.NA.method", values_to = "Value") %>%
  group_by(case.id, Metric) %>%
  summarize(
    max.diff = max(Value) - min(Value)
  ) %>%
  pivot_wider(names_from = Metric, values_from = max.diff) %>%
  summary()



# how do mean and max RMSEFrac compare between the four zero-degree-handling methods?

comp.methods = fb100.data %>%
  select(case.id, Metric, starts_with("RA")) %>%
  filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac")) %>%
  pivot_longer(
    cols=starts_with("RA"),
    names_to = c(".value", "Method"),
    names_pattern = "(.+)(.)"
  ) %>%
  group_by(Metric, Method) %>%
  summarize(
    min = min(RA.r.est),
    mean = mean(RA.r.est),
    max = max(RA.r.est)
  )

kable(comp.methods, digits=3)
# the methods don't make too huge of a difference
# |Metric      |Method |   min|  mean|   max|
# |:-----------|:------|-----:|-----:|-----:|
# |AbsBiasFrac |1      | 0.001| 0.263| 1.718|
# |AbsBiasFrac |2      | 0.000| 0.244| 1.651|
# |AbsBiasFrac |3      | 0.000| 0.241| 1.625|
# |AbsBiasFrac |4      | 0.000| 0.224| 1.559|
# |RMSEFrac    |1      | 0.044| 0.384| 1.860|
# |RMSEFrac    |2      | 0.043| 0.367| 1.792|
# |RMSEFrac    |3      | 0.044| 0.355| 1.767|
# |RMSEFrac    |4      | 0.043| 0.340| 1.699|
# |SEFrac      |1      | 0.042| 0.249| 1.405|
# |SEFrac      |2      | 0.041| 0.244| 1.343|
# |SEFrac      |3      | 0.039| 0.232| 1.279|
# |SEFrac      |4      | 0.039| 0.227| 1.262|

head(sort(case.data$NH), n=20)
# [1] 10 10 10 10 10 10 10 11 11 11 11 12 12 12 12 12 13 13 13 13


#### Choose one method for handling zero-valued degree estimates ####

fb100.data = fb100.data %>%
  # choose method 1 (set Infs to 1, exclude NAs)
  rename(RR = RR.r.est,
         AR = AR.r.est,
         RA = RA.r.est1,
         AA = AA.r.est1) %>%
  select(-ends_with("est2"), -ends_with("est3"), -ends_with("est4")) %>%
  # For each metric, compute whether the value is smaller for one estimator vs another
  mutate(
    smaller.RAvRR = ifelse(abs(abs(RA) - abs(RR)) < 10^-8,
                           "Same", ifelse(abs(RA) < abs(RR),
                                          "RA", "RR")),
    smaller.RAvAA = ifelse(abs(abs(RA) - abs(AA)) < 10^-8,
                           "Same", ifelse(abs(RA) < abs(AA),
                                          "RA", "AA")),
    smaller.RAvAR = ifelse(abs(abs(RA) - abs(AR)) < 10^-8,
                           "Same", ifelse(abs(RA) < abs(AR),
                                          "RA", "AR")),
    smaller.AAvRR = ifelse(abs(abs(AA) - abs(RR)) < 10^-8,
                           "Same", ifelse(abs(AA) < abs(RR),
                                          "AA", "RR")),
    smaller.AAvAR = ifelse(abs(abs(AA) - abs(AR)) < 10^-8,
                           "Same", ifelse(abs(AA) < abs(AR),
                                          "AA", "AR")),
    smaller.RRvAR = ifelse(abs(abs(RR) - abs(AR)) < 10^-8,
                           "Same", ifelse(abs(RR) < abs(AR),
                                          "RR", "AR")),
    Metric = as.factor(Metric),
    # these do not divide by r.true because RMSEFrac already does;
    # therefore look at ___Frac columns to assess percent difference
    PercDiff.RA.RR = RA-RR,
    PercDiff.RA.AA = RA-AA,
    PercDiff.RA.AR = RA-AR,
    PercDiff.AA.RR = AA-RR,
    PercDiff.AA.AR = AA-AR,
    PercDiff.RR.AR = RR-AR
  )


#### Examine assortativity, degree ratio, and true prevalence across the simulated cases ####

summary(fb100.data$assort)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01886 0.22914 0.28906 0.28160 0.33560 0.86687 

kable(fb100.data %>% group_by(dratio.level) %>% summarize(min = min(assort),
                                                          mean = mean(assort),
                                                          max = max(assort)))
# |dratio.level |       min|      mean|       max|
# |:------------|---------:|---------:|---------:|
# |Low          | 0.0559163| 0.2792333| 0.8668663|
# |High         | 0.0191208| 0.2932928| 0.5679451|
# |Near 1       | 0.0188558| 0.2715046| 0.5412684|


# true prevalences in these cases range from 0.10-9.9%, 2.2% on average
range(case.data$r.true)
# [1] 0.001010611 0.099428126
summary(case.data$r.true)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001011 0.006170 0.011650 0.021597 0.027194 0.099428


# number of cases whose degree ratios are low, high, or near 1
# (divide by 8 because there are 8 rows for each case, corresponding to the 8 metrics)
table(fb100.data$dratio.level)/8
# Low   High Near 1 
# 290    360    349 

fb100.data$nN = fb100.data$N*500
range(fb100.data$nN)
# [1]   384500 20777000


#### Evaluate and compare estimator performance ####
# Generate tables summarizing what proportion of the time the RA has lower
# bias, variance, and RMSE than the RR and AA


# Overall ranking of estimators by RMSE: what proportion of the time 
# does each estimator have the lowest RMSE?
fb100.rank = fb100.data %>%
  filter(Metric=="RMSEFrac") %>% # can change to other metrics
  select(AA:RR)

fb100.rank = fb100.rank %>%
  mutate(rank = colnames(fb100.rank)[apply(fb100.rank, 1, which.min)]) %>% # can change to which.max
  mutate(ranksmall = ifelse(abs(RR-AR) < 10^-8, "Same", ifelse(RR < AR, "RR", "AR")))

fb100.rank.table = fb100.rank %>%
  count(rank)
fb100.rank.table$prop = fb100.rank.table$n/dim(case.data)[1]

kable(fb100.rank.table)
# method 1: a slightly more modest result but RA still on top
# |rank |   n|      prop|
# |:----|---:|---------:|
# |AA   | 210| 0.2102102|
# |AR   | 182| 0.1821822|
# |RA   | 407| 0.4074074|
# |RR   | 200| 0.2002002|
# method 2:
# |rank |   n|      prop|
# |:----|---:|---------:|
# |AA   | 219| 0.2192192|
# |AR   | 174| 0.1741742|
# |RA   | 416| 0.4164164|
# |RR   | 190| 0.1901902|
# method 3:
# |rank |   n|      prop|
# |:----|---:|---------:|
# |AA   | 227| 0.2272272|
# |AR   | 166| 0.1661662|
# |RA   | 427| 0.4274274|
# |RR   | 179| 0.1791792|
# method 4: 
# |rank |   n|      prop|
# |:----|---:|---------:|
# |AA   | 239| 0.2392392|
# |AR   | 159| 0.1591592|
# |RA   | 432| 0.4324324|
# |RR   | 169| 0.1691692|


# Pairwise comparisons over all cases
fb100.counts = fb100.data %>% 
  filter(Metric %in% c("AbsBias", "SE", "RMSE")) %>%
  group_by(Metric) %>%
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))

kable(fb100.counts, digits=3)
# method 1
# |Metric  |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:-------|---:|------:|------:|------:|------:|------:|------:|
# |AbsBias | 999|  0.699|  0.636|  0.699|  0.667|  0.666|  0.504|
# |RMSE    | 999|  0.615|  0.671|  0.619|  0.584|  0.586|  0.523|
# |SE      | 999|  0.262|  0.730|  0.267|  0.253|  0.254|  0.646|
# method 2
# |Metric  |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:-------|---:|------:|------:|------:|------:|------:|------:|
# |AbsBias | 999|  0.725|  0.626|  0.727|  0.691|  0.692|  0.504|
# |RMSE    | 999|  0.637|  0.664|  0.636|  0.607|  0.608|  0.523|
# |SE      | 999|  0.267|  0.728|  0.270|  0.255|  0.259|  0.646|
# method 3
# |Metric  |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:-------|---:|------:|------:|------:|------:|------:|------:|
# |AbsBias | 999|  0.723|  0.619|  0.729|  0.695|  0.694|  0.504|
# |RMSE    | 999|  0.655|  0.659|  0.658|  0.627|  0.630|  0.523|
# |SE      | 999|  0.279|  0.731|  0.283|  0.269|  0.272|  0.646|
# method 4
# |Metric  |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:-------|---:|------:|------:|------:|------:|------:|------:|
# |AbsBias | 999|  0.751|  0.608|  0.750|  0.715|  0.719|  0.504|
# |RMSE    | 999|  0.672|  0.647|  0.672|  0.645|  0.646|  0.523|
# |SE      | 999|  0.286|  0.729|  0.290|  0.278|  0.283|  0.646|


# Pairwise comparisons by degree ratio
fb100.counts = fb100.data %>% 
  # filter(Metric %in% c("AbsBias", "SE", "RMSE")) %>%
  filter(Metric %in% c("RMSE")) %>%
  group_by(dratio.level, Metric) %>% 
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))


kable(fb100.counts, digits=3)
# |dratio.level |Metric |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:------------|:------|---:|------:|------:|------:|------:|------:|------:|
# |Low          |RMSE   | 290|  0.800|  0.641|  0.803|  0.762|  0.762|  0.421|
# |High         |RMSE   | 360|  0.783|  0.689|  0.792|  0.761|  0.767|  0.669|
# |Near 1       |RMSE   | 349|  0.287|  0.676|  0.287|  0.252|  0.252|  0.456|


# comparisons by assortativity; only informative if using all 1989 cases
# separates cases by whether they are among the 10 most assortative or 10 least assortative cases
fb100.data$assort.level2 = ifelse(fb100.data$assort >= 0.01, "Assortative", "Not Assortative")

fb100.counts = fb100.data %>%
  # filter(Metric %in% c("AbsBias", "SE", "RMSE")) %>%
  filter(Metric %in% c("RMSE")) %>%
  group_by(assort.level2, Metric) %>%
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))

kable(fb100.counts, digits=3)
# regardless of assortativity, RA tends to be better than AA
# otherwise, the other pairwise comparisons flip when assortativity changes
# |assort.level2   |Metric |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:---------------|:------|---:|------:|------:|------:|------:|------:|------:|
# |Assortative     |RMSE   | 999|  0.615|  0.671|  0.619|  0.584|  0.586|  0.523|
# |Not Assortative |RMSE   | 990|  0.187|  0.695|  0.192|  0.177|  0.176|  0.460|



# comparisons by true prevalence
# use thresholds 0.05, 0.1, 0.25, 0.5
fb100.r = fb100.data %>%
  mutate(r.level = ifelse(r.true < 0.01, "<0.01",
                          ifelse(r.true < 0.02, "<0.02",
                                 ifelse(r.true < 0.03, "<0.03",
                                        ifelse(r.true < 0.04, "<0.04", ">=0.04"))))) %>%
  filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac")) %>%
  group_by(Metric, r.level) %>% 
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))
kable(fb100.r, digits=3)
# method 1:
# |Metric      |r.level |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:-----------|:-------|---:|------:|------:|------:|------:|------:|------:|
# |AbsBiasFrac |<0.01   | 449|  0.733|  0.670|  0.733|  0.717|  0.719|  0.546|
# |AbsBiasFrac |<0.02   | 221|  0.638|  0.652|  0.633|  0.593|  0.597|  0.538|
# |AbsBiasFrac |<0.03   | 101|  0.693|  0.644|  0.703|  0.644|  0.644|  0.525|
# |AbsBiasFrac |<0.04   |  55|  0.509|  0.473|  0.509|  0.436|  0.436|  0.400|
# |AbsBiasFrac |>=0.04  | 173|  0.751|  0.572|  0.751|  0.717|  0.699|  0.370|

# |RMSEFrac    |<0.01   | 449|  0.617|  0.708|  0.626|  0.601|  0.601|  0.579|
# |RMSEFrac    |<0.02   | 221|  0.543|  0.674|  0.538|  0.507|  0.520|  0.552|
# |RMSEFrac    |<0.03   | 101|  0.644|  0.683|  0.634|  0.604|  0.604|  0.525|
# |RMSEFrac    |<0.04   |  55|  0.473|  0.491|  0.491|  0.400|  0.418|  0.400|
# |RMSEFrac    |>=0.04  | 173|  0.728|  0.618|  0.734|  0.682|  0.671|  0.376|

# |SEFrac      |<0.01   | 449|  0.439|  0.733|  0.445|  0.432|  0.432|  0.679|
# |SEFrac      |<0.02   | 221|  0.199|  0.738|  0.199|  0.190|  0.195|  0.588|
# |SEFrac      |<0.03   | 101|  0.050|  0.752|  0.059|  0.050|  0.050|  0.594|
# |SEFrac      |<0.04   |  55|  0.073|  0.545|  0.073|  0.036|  0.036|  0.600|
# |SEFrac      |>=0.04  | 173|  0.069|  0.757|  0.075|  0.058|  0.058|  0.676|

# comparisons by % network that has zero-valued degree estimates
fb100.counts = fb100.data %>% 
  # filter(Metric %in% c("AbsBias", "SE", "RMSE")) %>%
  filter(Metric %in% c("RMSE")) %>%
  group_by(r.0deg.level, Metric) %>% 
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))

kable(fb100.counts, digits=3)
# |r.0deg.level |Metric |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:------------|:------|---:|------:|------:|------:|------:|------:|------:|
# |1st Qu.      |RMSE   | 250|  0.664|  0.608|  0.664|  0.624|  0.632|  0.476|
# |2nd Qu.      |RMSE   | 240|  0.592|  0.683|  0.592|  0.583|  0.583|  0.621|
# |3rd Qu.      |RMSE   | 250|  0.644|  0.696|  0.656|  0.616|  0.620|  0.488|
# |4th Qu.      |RMSE   | 259|  0.560|  0.695|  0.564|  0.514|  0.510|  0.510|


# same but use thresholds skewed low
fb100.r = fb100.data %>%
  mutate(r.0deg.grp = ifelse(r.0deg < 0.01, "<0.01",
                          ifelse(r.0deg < 0.02, "<0.02",
                                 ifelse(r.0deg < 0.03, "<0.03", ">=0.03")))) %>%
  # filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac")) %>%
  filter(Metric %in% c("RMSEFrac")) %>%
  group_by(Metric, r.0deg.grp) %>% 
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))
kable(fb100.r, digits=3)
# |Metric   |r.0deg.grp |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:--------|:----------|---:|------:|------:|------:|------:|------:|------:|
# |RMSEFrac |<0.01      | 140|  0.693|  0.593|  0.686|  0.650|  0.657|  0.493|
# |RMSEFrac |<0.02      | 340|  0.591|  0.659|  0.594|  0.576|  0.576|  0.559|
# |RMSEFrac |<0.03      | 369|  0.623|  0.694|  0.634|  0.583|  0.588|  0.515|
# |RMSEFrac |>=0.03     | 150|  0.573|  0.713|  0.573|  0.540|  0.533|  0.487|


# comparisons by pHL level
fb100.data$pHL = fb100.data$HL.edges / (fb100.data$NH*(fb100.data$N - fb100.data$NH))
summary(fb100.data$pHL)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.404e-05 2.978e-03 6.074e-03 9.351e-03 1.167e-02 8.774e-02
fb100.data$pHL.level = ifelse(fb100.data$pHL < summary(fb100.data$pHL)[2], "1st Qu.",
                              ifelse(fb100.data$pHL < summary(fb100.data$pHL)[3], "2nd Qu.",
                                     ifelse(fb100.data$pHL < summary(fb100.data$pHL)[5],
                                            "3rd Qu.", "4th Qu.")))
# check:
table(fb100.data$pHL.level)/8
# 1st Qu. 2nd Qu. 3rd Qu. 4th Qu. 
#    249     250     250     250 

fb100.counts = fb100.data %>% 
  # filter(Metric %in% c("AbsBias", "SE", "RMSE")) %>%
  filter(Metric %in% c("RMSE")) %>%
  group_by(pHL.level, Metric) %>% 
  summarize(n = n(),
            RAvsRR = sum(RA<RR, na.rm = TRUE),
            RAvsAA = sum(RA<AA, na.rm = TRUE),
            RAvsAR = sum(RA<AR, na.rm = TRUE),
            AAvsRR = sum(AA<RR, na.rm = TRUE),
            AAvsAR = sum(AA<AR, na.rm = TRUE),
            RRvsAR = sum(RR<AR, na.rm = TRUE)
  ) %>%
  mutate(across(.cols = RAvsRR:RRvsAR, .fns = ~.x / n))

kable(fb100.counts, digits=3)
# just as we found with the theoretical results, RA does better the higher pHL is
# |pHL.level |Metric |   n| RAvsRR| RAvsAA| RAvsAR| AAvsRR| AAvsAR| RRvsAR|
# |:---------|:------|---:|------:|------:|------:|------:|------:|------:|
# |1st Qu.   |RMSE   | 249|  0.426|  0.711|  0.434|  0.406|  0.410|  0.482|
# |2nd Qu.   |RMSE   | 250|  0.656|  0.688|  0.656|  0.624|  0.628|  0.520|
# |3rd Qu.   |RMSE   | 250|  0.692|  0.648|  0.684|  0.652|  0.640|  0.540|
# |4th Qu.   |RMSE   | 250|  0.684|  0.636|  0.700|  0.652|  0.664|  0.548|




# examine degree estimator performance by degree ratio
deg.est.by.dratio = comp.deg.est %>%
  left_join(select(filter(fb100.data, Metric=="RMSEFrac"), case.id, dratio.level)) %>%
  group_by(dratio.level) %>%
  summarize(
    n = n(),
    nRbetter = sum(smaller.perc.err == "R"),
    nAbetter = sum(smaller.perc.err == "A")
  ) %>%
  mutate(across(.cols = nRbetter:nAbetter, .fns = ~.x / n))
kable(deg.est.by.dratio)
# it makes sense that the degree estimator performance doesn't depend much on degree ratio
# |dratio.level |       n|  nRbetter|  nAbetter|
# |:------------|-------:|---------:|---------:|
# |Low          | 2626489| 0.5002549| 0.4253085|
# |High         | 5559990| 0.5013795| 0.4188518|
# |Near 1       | 3895912| 0.4930579| 0.4252365|





#### Generate the figures ####

textbasesize = 10 #18

# estimated vs true degrees for supplemental figure
deg.max = ceiling(max(individual.data$d.A.est[sampids])/100)*100

deg.R.vs.true = ggplot(slice(individual.data, sampids), aes(x=d.true, y=d.R.est)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  xlim(c(0, deg.max)) + 
  ylim(c(0, deg.max)) +
  xlab("True degree") + ylab("R degree estimate") + 
  theme_minimal(base_size = textbasesize)
print(deg.R.vs.true)
ggsave(paste0(fig.dir, "deg_Rvstrue.pdf"),
       width = 3, height = 2, units = "in")

deg.A.vs.true = ggplot(slice(individual.data, sampids), aes(x=d.true, y=d.A.est)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  xlim(c(0, deg.max)) + 
  ylim(c(0, deg.max)) +
  xlab("True degree") + ylab("A degree estimate") + 
  theme_minimal(base_size = textbasesize)
print(deg.A.vs.true)
ggsave(paste0(fig.dir, "deg_Avstrue.pdf"),
       width = 3, height = 2, units = "in")



gen.fb100.fig = function(data, prefix){
  pt.size = 0.5
  xymax = 2.5
  palette = c('#bdbdbd','#000000','#636363')
  width = 6
  height = 2.2
  ht.with.lgnd = 3
  ht.overall = height*2 + ht.with.lgnd
  
  plotdata = data %>%
    filter(Metric %in% c("AbsBiasFrac", "SEFrac", "RMSEFrac")) %>%
    mutate(Metric = as.factor(ifelse(Metric=="AbsBiasFrac", "|Bias|/r",
                                     ifelse(Metric=="SEFrac", "SE/r", "RMSE/r")))) %>%
    mutate(Metric = relevel(relevel(Metric, "SE/r"), "|Bias|/r"),
           `Degree ratio` = dratio.level)
  
  if(prefix=="low_dratio_"){
    plotdata = filter(plotdata, d.ratio <= 0.8)
    # Figure 2
    Fig2 = ggplot(plotdata, aes(x=RA, y=RR)) +
      geom_point(cex=pt.size) +
      geom_abline(slope = 1, intercept = 0, color = "black") +
      xlim(c(0, xymax)) +
      ylim(c(0, xymax)) +
      facet_wrap(~ Metric) +
      theme_minimal(base_size = textbasesize) +
      theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
            aspect.ratio=1, plot.margin = margin(10, 20, 10, 10))
    print(Fig2)
    ggsave(paste0(fig.dir, prefix, "RA_RR_Fig2.pdf"),
           width = width, height = height, units = "in")
  }else if(prefix=="high_dratio_"){
    plotdata = filter(plotdata, d.ratio >= 1.2)
  }else if(prefix=="near1_dratio_"){
    plotdata = filter(plotdata, d.ratio > 0.8, d.ratio < 1.2)
  }
  
  RA.RR = ggplot(plotdata, aes(x=RA, y=RR, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10))
  print(RA.RR)
  ggsave(paste0(fig.dir, prefix, "RA_RR.pdf"),
         width = width, height = height, units = "in")

  RA.AA = ggplot(plotdata, aes(x=RA, y=AA, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10))
  print(RA.AA)
  ggsave(paste0(fig.dir, prefix, "RA_AA.pdf"),
         width = width, height = height, units = "in")

  AA.RR = ggplot(plotdata, aes(x=AA, y=RR, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "bottom", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10),
          legend.key.width = unit(0.6, 'in'))
  print(AA.RR)
  ggsave(paste0(fig.dir, prefix, "AA_RR.pdf"),
         width = width, height = ht.with.lgnd, units = "in")

  ggarrange(RA.RR, RA.AA, AA.RR,
            ncol=1, nrow=3, heights = c(2.2, 2.2, 2.9))
  ggsave(paste0(fig.dir, prefix, "combined.pdf"),
         width = width, height = ht.overall, units = "in")


  RA.AR = ggplot(plotdata, aes(x=RA, y=AR, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10))
  print(RA.AR)
  ggsave(paste0(fig.dir, prefix, "RA_AR.pdf"),
         width = width, height = height, units = "in")

  AA.AR = ggplot(plotdata, aes(x=AA, y=AR, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "none", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10))
  print(AA.AR)
  ggsave(paste0(fig.dir, prefix, "AA_AR.pdf"),
         width = width, height = height, units = "in")

  RR.AR = ggplot(plotdata, aes(x=RR, y=AR, color=`Degree ratio`)) +
    geom_point(cex=pt.size) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_color_manual(values = palette) +
    xlim(c(0, xymax)) +
    ylim(c(0, xymax)) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = textbasesize) +
    theme(legend.position = "bottom", panel.spacing = unit(0.5, "in"),
          aspect.ratio=1, plot.margin = margin(10, 20, 10, 10),
          legend.key.width = unit(0.6, 'in'))
  print(RR.AR)
  ggsave(paste0(fig.dir, prefix, "RR_AR.pdf"),
         width = width, height = ht.with.lgnd, units = "in")

  ggarrange(RA.AR, AA.AR, RR.AR,
            ncol=1, nrow=3, heights = c(2.2, 2.2, 2.9))
  ggsave(paste0(fig.dir, prefix, "combinedAR.pdf"),
         width = width, height = ht.overall, units = "in")
}

# All cases (Figure 3 is the combined figure)
gen.fb100.fig(fb100.data, prefix="all_")
# Cases with low degree ratio (Figure 2 is the RA_RR figure)
gen.fb100.fig(fb100.data, prefix="low_dratio_")
# Cases with high degree ratio (Supplement Figure 5 is the RA_RR figure)
gen.fb100.fig(fb100.data, prefix="high_dratio_")
# Cases with degree ratio near 1 (Supplement Figure 4 is the RA_RR figure)
gen.fb100.fig(fb100.data, prefix="near1_dratio_")








