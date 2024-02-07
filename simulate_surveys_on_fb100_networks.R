# This code simulates surveys from the Facebook 100 network data (main paper, 
# Section 5). Specifically,
#
#  - This script preprocesses the Traud et al (2011, 2012) data on 100
#    Facebook school networks from September 2005, identifies candidate probe
#    and hidden groups as any covariate value that constitutes 0.1-10% of the 
#    total school network and contains at least 10 people. For each school, the 
#    probe groups are taken to be the 20 largest candidate groups with 
#    assortativity coefficients whose magnitudes are < 0.1, and the hidden 
#    groups are chosen to be the ten most and ten least assortative groups with
#    prev between 0.1 and 10%, where less assortative means assortativity
#    coefficients close to or below 0.
#
#    ==> This output is stored in data.dir as all_probe_hidden_groups.RData.
#
#  - Next, this script uses this data to draw survey samples from each case
#    and compute the prevalence estimate for that case using each of the four 
#    estimators considered in the paper.
#
#    ==> This output is stored in out.dir as fb100_sims.RData.

library(knitr)
library(tictoc)
library(scales)
library(R.matlab)
library(magrittr)
library(assortnet)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

#### Specify inputs ####
n.surv = 500 # number of simulated surveys
n.frame = 500 # frame sample size for each survey

# The absolute path to the facebook100 data folder
data.dir = "/your/absolute/path/here/ending/with/slash/"

# The absolute path to the folder where you would like to store the data output
#  from this script
out.dir = "/your/absolute/path/here/ending/with/slash/"




#### Identify candidate groups, probe groups, hidden groups ####

filelist = list.files(path=data.dir, pattern="*.mat")

# Loop over the 100 school networks
for(school.id in seq_along(filelist)){
  # store case-level data (data about the school-hidden group combinations)
  # one row for each case (number of rows = n.cases.per.school * 100 schools)
  
  # extract school name from filename
  school.name = strsplit(filelist[school.id], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",
                         perl=TRUE)[[1]][1]
  print(paste("School", school.id, school.name))
  d <- readMat(paste0(data.dir,filelist[school.id]))
  linfo = as.data.frame(d$local.info)
  names(linfo) = c("status", "gender", "major", "minor", "dorm", "year", "hs")
  # set some years that don't make sense to 0 (and later to NA)
  linfo$year[linfo$year<0] = 0
  linfo$year[linfo$year>2010] = 0
  # If removing any rows/respondents, make sure to set N afterward
  N = dim(linfo)[1] # total population size of the network
  
  # Make the covariates into indicators for each possible level
  # e.g. names(covs) = status1, status2, ..., status5, gender1, gender2, major45, ...
  #  and the entries are all 0 or 1 for whether that respondent (row) is in that group (column)
  covs = linfo
  # make all columns factors
  covs = as.data.frame(lapply(covs[, 1:ncol(covs)], as.factor))
  covs = model.matrix(~ ., data=covs)[,-1] #make indicators from the levels
  
  # Set all 0s to NAs
  linfo[linfo==0] = NA
  # make all columns factors
  linfo = as.data.frame(lapply(linfo[, 1:ncol(linfo)], as.factor))
  
  ##### compute true degrees #####
  d.true = rowSums(d$A)
  
  ##### Identify candidate probe + hidden groups #####
  grp.sizes = colSums(covs)
  grp.prevs = grp.sizes/N
  n.edges = sum(d$A)/2
  # filter possible hidden/probe groups by prevalence and group size
  candidate.grps = data.frame(
    school.id = school.id,
    N = N,
    grp.id = 1:length(grp.sizes),
    grp.name = colnames(covs),
    grp.size = grp.sizes,
    grp.prev = grp.prevs) %>%
    filter(
      grp.prevs >= 0.001 & grp.prevs <= 0.1,
      grp.size >= 10) %>% # ensure there are at least 10 people in the group
    # compute assortativity and filter again
    mutate(
      HH.edges = unlist(map(grp.id,
                            function(x) sum(d$A[as.logical(covs[,x]), as.logical(covs[,x])])/2)),
      LL.edges = unlist(map(grp.id,
                            function(x) sum(d$A[!as.logical(covs[,x]), !as.logical(covs[,x])])/2)),
      HL.edges = unlist(map(grp.id,
                            function(x) sum(d$A[as.logical(covs[,x]), !as.logical(covs[,x])])))
    ) %>%
    mutate(
      f.same.grp = (HH.edges + LL.edges)/n.edges,
      C = ((HH.edges + HL.edges/2)/n.edges)^2 + ((LL.edges + HL.edges/2)/n.edges)^2
    ) %>%
    mutate(
      assort = (f.same.grp-C)/(1-C)
    ) %>%
    mutate(
      abs.assort = abs(assort)
    )
  save(candidate.grps, file = paste0(data.dir, "candidate_groups", school.id, ".RData"))
}

# read in results and decide on criteria for probe and hidden groups
load(paste0(data.dir, "candidate_groups1.RData"))
candidate.grps.all = candidate.grps

for(i in 2:99){
  load(paste0(data.dir, "candidate_groups", i, ".RData"))
  candidate.grps.all = bind_rows(candidate.grps.all, candidate.grps)
}

load(paste0(data.dir, "candidate_groups100.RData"))
candidate.grps = bind_rows(candidate.grps.all, candidate.grps)
rm(candidate.grps.all, i)

# choose 20 largest candidate groups with abs(assort) < 0.1
probe.grps = candidate.grps %>%
  filter(
    abs(assort) < 0.1
  ) %>%
  group_by(school.id) %>%
  slice_max(grp.prev, n=20, with_ties=FALSE)

# choose 10 groups with the highest assortativities
hidden.grps = setdiff(candidate.grps, probe.grps) %>%
  # first filter out probe groups
  # filter() %>%
  group_by(school.id) %>%
  slice_max(assort, n=10, with_ties = FALSE)

hidden.grps.lessassort = setdiff(candidate.grps, probe.grps) %>%
  # first filter out probe groups
  # filter() %>%
  group_by(school.id) %>%
  slice_min(assort, n=10, with_ties = FALSE)

range(hidden.grps.lessassort$assort)
# [1] -0.01289541  0.43038701

hidden.grps = union(hidden.grps, hidden.grps.lessassort)

# one school has only 9 cases, the other schools have 20 cases
table((hidden.grps %>% count(school.id))$n)
# indeed, the smallest number of candidate groups per school is 29, followed by 56
table((candidate.grps %>% count(school.id))$n)

save(probe.grps, hidden.grps, file = paste0(data.dir, "all_probe_hidden_groups.RData"))





#### Function to simulate surveys ####

simulate_surveys_on_fb100_networks = function(n.surv, n.frame,
                                              data.dir, out.dir){
  print("Inputs provided:")
  print(paste("n.surv =", n.surv))
  print(paste("n.frame =", n.frame))
  print(paste("data.dir =", data.dir))
  print(paste("out.dir =", out.dir))
  tic("Total run time for simulating surveys")
  set.seed(145)
  filelist = list.files(path=data.dir, pattern="*.mat")
  case.id=0
  
  # Load probe and hidden group choices
  load(paste0(data.dir, "all_probe_hidden_groups.RData"))
  
  # Loop over the 100 school networks
  for(this.school.id in seq_along(filelist)){
    # store case-level data (data about the school-hidden group combinations)
    # one row for each case (number of rows = n.cases.per.school * 100 schools)
    case.data = data.frame(
      case.id = integer(),
      school = character(),
      # there are 2 schools named Texas and 3 named UC,
      # so school.id disambiguates them
      school.id = integer(),
      N = integer(),
      HH.edges = integer(),
      LL.edges = integer(),
      HL.edges = integer(),
      hidden = character(),
      hidden.id = integer(),
      NH = numeric(),
      r.true = numeric(),
      d.ratio = numeric(),
      a = numeric(),
      assort = numeric()
    )
    
    # store individual people's data
    # number of rows = sum of 100 network sizes times n.cases.per.school; sum(case.data$N)
    individual.data = data.frame(
      case.id = integer(),
      is.hidden = logical(),
      y.iH = integer(),
      d.true = integer(),
      d.R.est = numeric(),
      d.A.est = numeric()
    )
    
    # store data from each survey
    # number of rows = n.surv surveys * n.cases.per.school * 100 schools
    survey.data = data.frame(
      case.id = integer(),
      survey.id = integer(),
      n.hidden = integer(),
      n.deg0 = integer(),
      n.yd.inf = integer(),
      n.yd.na = integer(),
      n.yd.gt1R = integer(),
      n.yd.gt1A = integer(),
      RR.r.est = numeric(),
      AR.r.est = numeric(),
      RA.r.est1 = numeric(),
      RA.r.est2 = numeric(),
      RA.r.est3 = numeric(),
      RA.r.est4 = numeric(),
      AA.r.est1 = numeric(),
      AA.r.est2 = numeric(),
      AA.r.est3 = numeric(),
      AA.r.est4 = numeric()
    )
    
    # Extract school name from filename
    school.name = strsplit(filelist[this.school.id], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",
                   perl=TRUE)[[1]][1]
    print(paste("School", this.school.id, school.name))
    tic(paste("Run time for School", this.school.id, school.name))
    
    ##### Read in and preprocess school network data #####
    tic("Run time to read network")
    d <- readMat(paste0(data.dir,filelist[this.school.id]))
    toc()
    linfo = as.data.frame(d$local.info)
    names(linfo) = c("status", "gender", "major", "minor", "dorm", "year", "hs")
    # set some years that don't make sense to 0 (and later to NA)
    linfo$year[linfo$year<0] = 0
    linfo$year[linfo$year>2010] = 0
    # If removing any rows/respondents, make sure to set N afterward
    N = dim(linfo)[1] # total population size of the network
    
    # Make the covariates into indicators for each possible level
    # e.g. names(covs) = status1, status2, ..., status5, gender1, gender2, major45, ...
    #  and the entries are all 0 or 1 for whether that respondent (row) is in that group (column)
    covs = linfo
    # make all columns factors
    covs = as.data.frame(lapply(covs[, 1:ncol(covs)], as.factor))
    covs = model.matrix(~ ., data=covs)[,-1] #make indicators from the levels
    
    # Set all 0s to NAs
    linfo[linfo==0] = NA
    # make all columns factors
    linfo = as.data.frame(lapply(linfo[, 1:ncol(linfo)], as.factor))
    
    #compute true degrees
    d.true = rowSums(d$A)
    
    ##### Identify the probe groups + hidden group #####
    curr.hidden.grps = filter(hidden.grps, school.id == this.school.id)
    n.cases.per.school = dim(curr.hidden.grps)[1]
    # probe groups (same for all cases for a given school)
    curr.probe.grps = filter(probe.grps, school.id == this.school.id)
    probe.groups = curr.probe.grps$grp.name
    probe.ids = curr.probe.grps$grp.id
    Njs = curr.probe.grps$grp.size
    Nj.sum = sum(Njs)
    n.probe = length(probe.groups)
    
    # Loop over each case (each choice of hidden group)
    for(j in seq_along(curr.hidden.grps$grp.id)){
      cat(paste0("..", j)) # progress indicator
      case.id = case.id + 1
      
      # hidden group
      hidden.group = curr.hidden.grps$grp.name[j]
      hidden.id = curr.hidden.grps$grp.id[j]
      NH = curr.hidden.grps$grp.size[j]
      r.true = curr.hidden.grps$grp.prev[j]
      
      # ID membership in H and compute degree ratio
      # indices for people in H; length(is.hidden) = # people in H
      is.hidden = as.logical(covs[,hidden.id])
      # degrees of hidden people
      d.hidden = d.true[which(is.hidden)]
      # degree ratio (hidden/frame, and here frame = total pop)
      d.ratio = mean(d.hidden) / mean(d.true)
      # Frame ratio is 1 (frame = general population)
      # True positive rate is 1 (no simulated transmission/response error)
      
      # Compute number of edges and estimate a #####
      HH.edges = curr.hidden.grps$HH.edges[j]
      LL.edges = curr.hidden.grps$LL.edges[j]
      HL.edges = curr.hidden.grps$HL.edges[j]
      poss.HH.edges = NH*(NH-1)/2
      poss.LL.edges = (N-NH)*(N-NH-1)/2
      poss.HL.edges = NH*(N-NH)
      a = mean(c(HH.edges/poss.HH.edges, LL.edges/poss.LL.edges)) / (HL.edges/poss.HL.edges)
      
      ##### Update case.data #####
      case.data = case.data %>%
        add_row(case.id = case.id,
                school = school.name,
                school.id = this.school.id,
                N = N,
                HH.edges = HH.edges,
                LL.edges = LL.edges,
                HL.edges = HL.edges,
                hidden = hidden.group,
                hidden.id = j,
                NH = NH,
                r.true = r.true,
                d.ratio = d.ratio,
                a = a,
                assort = curr.hidden.grps$assort[j])
    
      ##### Compute ARD and estimated degrees for everyone in the network #####
      # this uses all people in the network (d$A), not just respondents (d$A[resp_id,])
      # dim = network size x number of groups (hidden + probe)
      # each element is Y_ij, the number of people that respondent i knows in group j
      ard.all = as.matrix(d$A %*% as.matrix(covs[,c(probe.ids, hidden.id)]))
      
      # compute degree estimates for all people in the population
      # using ARD responses for the probe groups only
      # Ratio-of-average degree estimates: N * sum_j y_ij / sum_j N_j
      d.R.est = N * rowSums(ard.all[,-(n.probe+1)]) / Nj.sum
      # Average-of-ratio degree estimates: N * mean_j (y_ij / N_j)
      d.A.est = N * rowSums(ard.all[,-(n.probe+1)] %*% diag(1/Njs))/n.probe
      y.iH = ard.all[, (n.probe+1)]
      
      ##### Update individual.data #####
      individual.data = individual.data %>%
        add_row(case.id = case.id,
                is.hidden = is.hidden,
                y.iH = y.iH,
                d.true = d.true,
                d.R.est = d.R.est,
                d.A.est = d.A.est)
      
      ##### Draw n.surv many surveys from the school network under this choice of hidden group and probe groups #####
      for(i in 1:n.surv){
        ##### Draw a survey sample of n.frame many people from the frame pop #####
        frame.id = sample(1:N, n.frame, replace = FALSE)
        # the ARD and degrees for this sample are just a subset of everyone's
        d.R.sample = d.R.est[frame.id]
        d.A.sample = d.A.est[frame.id]
        d.true.sample = d.true[frame.id]
        y.iH.sample = y.iH[frame.id]
        # Compute n_H (n.hidden), the number of people in the sample that are hidden
        n.hidden = sum(is.hidden[frame.id])
        # record how many degrees were estimated to be zero
        n.deg0 = sum(d.R.sample == 0)
        
        ##### Compute prevalence estimates for H using this survey sample #####
        y.iH.sum = sum(y.iH.sample)
        d.R.sum = sum(d.R.sample)
        d.A.sum = sum(d.A.sample)
        
        # RR (MLE) estimate: sum_i y_iu / sum_i d_i, with R degree estimates
        RR.r.est = y.iH.sum / d.R.sum
        
        # AR estimate: sum_i y_iu / sum_i d_i, with A degree estimates
        AR.r.est = y.iH.sum / d.A.sum
        
        # for the RA and AA estimators, we will have to handle zero-valued degree estimates
        ratios.dR = y.iH.sample / d.R.sample
        ratios.dA = y.iH.sample / d.A.sample
        # record how many NAs and Infs occurred in this sample
        n.yd.inf = sum(is.infinite(ratios.dR))
        n.yd.na = sum(is.na(ratios.dR))
        n.yd.gt1R = sum(ratios.dR[!is.infinite(ratios.dR)]>1,  na.rm=TRUE)
        n.yd.gt1A = sum(ratios.dA[!is.infinite(ratios.dR)]>1,  na.rm=TRUE)
        # censor ratios > 1 to 1 even when est. degree > 0
        ratios.dR = ifelse(ratios.dR>1 & !is.infinite(ratios.dR), 1, ratios.dR)
        ratios.dA = ifelse(ratios.dA>1 & !is.infinite(ratios.dA), 1, ratios.dA)
        
        # Approach 1: set Inf to 1, exclude NA
        # censor individual ratios >1 to be 1
        # we found all cases of ratios > 1 were Infs that occurred because
        # the denominator (estimated degree) was 0
        ratios.dR.Infto1 = ifelse(is.infinite(ratios.dR), 1, ratios.dR)
        ratios.dA.Infto1 = ifelse(is.infinite(ratios.dA), 1, ratios.dA)
        RA.r.est1 = mean(ratios.dR.Infto1, na.rm = T)
        AA.r.est1 = mean(ratios.dA.Infto1, na.rm = T)
        
        # Approach 2: set Inf to 1 and NA to 0
        ratios.dR.both = ifelse(is.na(ratios.dR.Infto1), 0, ratios.dR.Infto1)
        ratios.dA.both = ifelse(is.na(ratios.dA.Infto1), 0, ratios.dA.Infto1)
        RA.r.est2 = mean(ratios.dR.both, na.rm = T)
        AA.r.est2 = mean(ratios.dA.both, na.rm = T)
        
        # Approach 3: exclude both Inf (#/0) and NA (0/0)
        ratios.dR.InftoNA = ifelse(is.infinite(ratios.dR), NA, ratios.dR)
        ratios.dA.InftoNA = ifelse(is.infinite(ratios.dA), NA, ratios.dA)
        RA.r.est3 = mean(ratios.dR.InftoNA, na.rm = T)
        AA.r.est3 = mean(ratios.dA.InftoNA, na.rm = T)
        
        # Approach 4: set NA to 0, exclude Inf
        ratios.dR.NAto0 = ifelse(is.na(ratios.dR.InftoNA), 0, ratios.dR.InftoNA)
        ratios.dA.NAto0 = ifelse(is.na(ratios.dA.InftoNA), 0, ratios.dA.InftoNA)
        RA.r.est4 = mean(ratios.dR.NAto0, na.rm = T)
        AA.r.est4 = mean(ratios.dA.NAto0, na.rm = T)
        
        ##### Add row to survey.data #####
        survey.data = survey.data %>%
          add_row(case.id = case.id,
                  survey.id = i,
                  n.hidden = n.hidden,
                  n.deg0 = n.deg0,
                  n.yd.inf = n.yd.inf,
                  n.yd.na = n.yd.na,
                  n.yd.gt1R = n.yd.gt1R,
                  n.yd.gt1A = n.yd.gt1A,
                  RR.r.est = RR.r.est,
                  AR.r.est = AR.r.est,
                  RA.r.est1 = RA.r.est1,
                  RA.r.est2 = RA.r.est2,
                  RA.r.est3 = RA.r.est3,
                  RA.r.est4 = RA.r.est4,
                  AA.r.est1 = AA.r.est1,
                  AA.r.est2 = AA.r.est2,
                  AA.r.est3 = AA.r.est3,
                  AA.r.est4 = AA.r.est4)
      } # end of loop to draw surveys
    } # end of loop over probe groups
    cat("\n")
    print(paste("Ended with case id", case.id))
    toc()
    save(case.data, individual.data, survey.data,
         n.surv, n.frame,
         file = paste0(out.dir, "school", this.school.id, ".RData"))
  } # end of loop over school networks
  toc()
} # end of survey simulating function


#### Run the simulation function once to generate the data for all cases ####
simulate_surveys_on_fb100_networks(n.surv, n.frame,
                                   data.dir, out.dir)
system("say The simulations have finished")




