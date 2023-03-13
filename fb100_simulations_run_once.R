# Jess Kunke, 2023
# This code generates the facebook simulation data. You only need to run this once.
# Specifically:
#  - First, this script preprocesses the Traud et al (2011, 2012) data on 100
#    Facebook school networks from September 2005, identifies candidate probe
#    and hidden groups as the 16 largest groups based on covariate values,
#    computes the sizes and degree ratios of these groups, and computes the 
#    assortativityof each combination (which we will call a case) of school, 
#    hidden group, and 15 probe groups.
#    ==> This output is stored as dratios.RData.
#  - Next, this script uses this data to draw survey samples from each of these
#    1600 cases and compute the prevalence estimate for each case using each of
#    the estimators considered in the paper.
#    ==> This output is stored as fb100_sims.RData.

library(knitr)
library(scales)
library(R.matlab)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

# The absolute path to the facebook100 data folder
data.dir = "/your/absolute/path/here/ending/with/slash/"
# The absolute path to the folder where you would like to store the data output
#  from this script
out.dir = "/your/absolute/path/here/ending/with/slash/"

# ---------------------------------------------------------------
# Part 1: generate the dratios.RData file
# read in all the school data and generate a data set with the
# following variables for each of the 1600 school-hidden group cases we
# consider:
  # - school name, id, and size (total population size)
  # - hidden group name, id, and size
  # - degree ratio for this choice of school and hidden group
set.seed(145)
n.grps = 16 # use the n.grps many largest groups as hidden + probe groups
filelist = list.files(path=data.dir, pattern="*.mat")

dratios = data.frame(school = character(),
                     school_ind = integer(),
                     N = integer(),
                     hidden = character(),
                     hidden_ind = integer(),
                     true.size = integer(),
                     d.ratio = numeric(),
                     assort = numeric())

for(school.ind in 1:length(filelist)){
  print(paste0("School: ", school.ind))
  # School name
  tmp = strsplit(filelist[school.ind], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",
                 perl=TRUE)[[1]]
  school.name = tmp[1]
  print(school.name)
  # networks$school[i] = tmp[1]
  # Read in file; this is the time intensive step in this loop (based on using
  # tictoc- the other steps were 0 sec, while this one tended to be 0.1-1 sec)
  print("Reading matrix data...")
  d <- readMat(paste0(data.dir,filelist[school.ind]))
  print("Done reading data.")
  
  linfo = as.data.frame(d$local.info)
  names(linfo) = c("status", "gender", "major", "minor", "dorm", "year", "hs")
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
  
  ## Use the 16 largest groups
  groups.ind = order(colSums(covs), decreasing = TRUE)[1:n.grps]
  groups = colnames(covs[,groups.ind])
  groups.true.sizes = colSums(covs[,groups.ind])
  # compute true degrees
  d.true = rowSums(d$A)
  
  for(j in 1:n.grps){ # for each choice of hidden group...
    # name and index of hidden group, chosen to be #j of the 16 largest groups
    hidden.group = groups[j]
    hidden.ind = groups.ind[j]
    # names and indices of probe groups (all the other 15 largest groups)
    probe.groups = groups[-j]
    probe.inds = groups.ind[-j]
    # Nu = true hidden pop size
    Nu = groups.true.sizes[j]
    #sum(covs[,hidden.ind])
    # sum of N_js (the number of people in each probe group)
    Njs = groups.true.sizes[-j]
    Nj.sum = sum(Njs) #sum(colSums(covs[,probe.inds]))
    
    # inds for people in the hidden group; length(is.hidden) = # of people in hidden group
    is.hidden = as.logical(covs[,hidden.ind])
    # degrees of hidden people
    d.hidden = d.true[which(is.hidden)]
    # degree ratio (hidden/frame, and here frame = total pop)
    d.ratio = mean(d.hidden) / mean(d.true)
    
    # compute assortativity coefficient (normalized modularity)
    n.edges = sum(d$A)/2 #length(d$A@x)/2
    # A.sum.same.grp = sum(d$A[is.hidden, is.hidden]) + sum(d$A[!is.hidden, !is.hidden])
    actual.edges.same.grp = 0.5*(sum(d$A[is.hidden, is.hidden]) + 
                                   sum(d$A[!is.hidden, !is.hidden]))
    deg.prod = 0
    for(g in 1:length(d.true[is.hidden])){
      deg.prod = deg.prod +
        sum(d.true[is.hidden][g]*d.true[is.hidden][-g]) #-d.true[is.hidden][g]^2
    }
    for(g in 1:length(d.true[!is.hidden])){
      deg.prod = deg.prod +
        sum(d.true[!is.hidden][g]*d.true[!is.hidden][-g]) #-d.true[!is.hidden][g]^2
    }
    exp.edges.same.grp = 0.25*deg.prod/n.edges
    assort = (actual.edges.same.grp-exp.edges.same.grp)/(n.edges-exp.edges.same.grp)
    
    dratios = dratios %>%
      add_row(school = school.name,
              school_ind = school.ind,
              N = N,
              hidden = hidden.group,
              hidden_ind = j,
              true.size = Nu,
              d.ratio = d.ratio,
              assort = assort)
  }
}

save(dratios, file = paste0(out.dir, "dratios.RData"))
system("say d ratios have been stored")




# ---------------------------------------------------------------
# Part 2: Simulate surveys from the Facebook 100 network data

# Uncomment this line if you want to reload dratios.RData
# load(paste0(out.dir, "dratios.RData"))

# Note:
# - It appears that the larger networks (more nodes and/or more edges) are 
#   actually more likely to have zeroes in the numerator and/or denominator.
# - The NaNs occur when both (1) the sum of how many people each person knows 
#   in the target group (part of the numerator) and (2) the sum (over 
#   respondents and over probe groups) of how many people each person knows in 
#   the probe groups (the denominator) are zero.
# - Inf occurs when only the denominator (quantity 2 above) is zero.

# This function expects a vector of school indices and a vector of hidden
#   indices of the same length, corresponding to the school-hidden group 
#   combinations of interest
include_PIMLE_with_true_degrees = function(school_inds, hidden_inds, data.dir,
                                           
                                           outdatafilename){
  set.seed(145)
  n.grps = 16 # use the n.grps many largest groups as hidden + probe groups
  n.surv = 500 # number of simulated surveys
  n.frame = 500 # frame sample size for each survey
  filelist = list.files(path=data.dir, pattern="*.mat")
  
  # no. rows = length(school_inds) = length(hidden_inds)
  school.hidden.data = data.frame(school = character(),
                                  N = integer(),
                                  hidden = character(),
                                  true.size = integer(),
                                  d.ratio = numeric(),
                                  n.frame = integer())
  
  # no. rows = sum of network sizes N for each school network (but if
  #  m hidden group choices for a given school, need to have m*N rows
  #  for that school because many of these variables also depend on
  #  hidden group)
  individual.data = data.frame(school = character(),
                               hidden = character(),
                               true.size = integer(),
                               is.hidden = logical(),
                               y.iu = integer(),
                               d.true = integer(),
                               d.killworth.est = numeric(),
                               d.mos.est = numeric())
  
  survey.data = data.frame(school = character(),
                           hidden = character(),
                           true.size = integer(),
                           n.hidden = integer(),
                           mle.size.est = numeric(),
                           mle.size.est.known.d = numeric(),
                           pimle.size.est = numeric(),
                           pimle.cens.est = numeric(),
                           pimle.size.dtrue = numeric(),
                           mos.size.est = numeric(),
                           mos.cens.est = numeric())
  
  for(ind in 1:length(school_inds)){
    school.ind = school_inds[ind]
    print(paste0(ind, " School: ", school.ind))
    # School name
    tmp = strsplit(filelist[school.ind], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])",
                   perl=TRUE)[[1]]
    school.name = tmp[1]
    print(school.name)
    # networks$school[i] = tmp[1]
    # Read in file; this is the time intensive step in this loop (based on using
    # tictoc- the other steps were 0 sec, while this one tended to be 0.1-1 sec)
    d <- readMat(paste0(data.dir,filelist[school.ind]))
    
    linfo = as.data.frame(d$local.info)
    names(linfo) = c("status", "gender", "major", "minor", "dorm", "year", "hs")
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
    
    ## Use the 16 largest groups
    groups.ind = order(colSums(covs), decreasing = TRUE)[1:n.grps]
    groups = colnames(covs[,groups.ind])
    groups.true.sizes = colSums(covs[,groups.ind])
    # compute true degrees
    d.true = rowSums(d$A)
    
    j = hidden_inds[ind]
    # name and index of hidden group, chosen to be #j of the 16 largest groups
    hidden.group = groups[j]
    hidden.ind = groups.ind[j]
    print(paste0("Hidden ind ", j, ": ", hidden.group, " ", hidden.ind))
    # names and indices of probe groups (all the other 15 largest groups)
    probe.groups = groups[-j]
    probe.inds = groups.ind[-j]
    # Nu = true hidden pop size
    Nu = groups.true.sizes[j] #sum(covs[,hidden.ind])
    # sum of N_js (the number of people in each probe group)
    Njs = groups.true.sizes[-j]
    Nj.sum = sum(Njs) #sum(colSums(covs[,probe.inds]))
    
    # inds for people in the hidden group; length(is.hidden) = # of people in hidden group
    is.hidden = as.logical(covs[,hidden.ind])
    # degrees of hidden people
    d.hidden = d.true[which(is.hidden)]
    # degree ratio (hidden/frame, and here frame = total pop)
    d.ratio = mean(d.hidden) / mean(d.true)
    ## Frame ratio is 1 (frame = total pop)
    ## True positive rate is 1
    ## so Feehan correction factor is just d.ratio
    
    # this uses all people in the network (d$A), not just respondents (d$A[resp_id,])
    # dim = nrespondents x ngroups (hidden + probe)
    # each element is y_ik, the number of people that respondent i knows in group k
    ard.all = as.matrix(d$A %*% as.matrix(covs[,c(probe.inds, hidden.ind)]))
    
    # compute degree estimates for all people in the population
    # Killworth degree estimate; originally called degree.est
    #   N * sum_j y_ij / sum_j N_j
    d.killworth.est = N * rowSums(ard.all[,-16]) / Nj.sum
    # MoS degree estimate
    #   N * mean_j (y_ij / N_j)
    d.mos.est = N * rowSums(ard.all[,-16] %*% diag(1/Njs))/15
    
    # save individual data for this choice of hidden group
    # for.this.hidden.grp = cbind(rep(hidden.group, length(Nu)), is.hidden,
    #                             ard.all[,16], d.true,
    #                             d.killworth.est, d.mos.est)
    # individual.data = rbind(individual.data, for.this.hidden.grp)
    individual.data = individual.data %>%
      add_row(school = school.name,
              hidden = rep(hidden.group, N),
              true.size = Nu,
              is.hidden = is.hidden,
              y.iu = ard.all[,16],
              d.true = d.true,
              d.killworth.est = d.killworth.est,
              d.mos.est = d.mos.est)
    
    # For each simulated survey...
    for(i in 1:n.surv){
      # draw a sample of sample-size many people from the frame pop
      frame.id = sample(1:N, n.frame, replace = FALSE)
      # the ARD for this sample are just a subset of everyone's ARD (ard.all)
      ard.sample = ard.all[frame.id,]
      d.killworth.sample = d.killworth.est[frame.id]
      d.mos.sample = d.mos.est[frame.id]
      d.true.sample = d.true[frame.id]
      # number of people in the sample that are hidden
      n.hidden = sum(is.hidden[frame.id])
      
      # MLE estimate: N * sum_i y_iu / sum_i d_i
      mle.size.est = N * sum(ard.sample[,16]) / sum(d.killworth.sample)
      
      # MLE estimate if we knew the degrees exactly
      mle.size.est.known.d = N * sum(ard.sample[,16]) / sum(d.true.sample)
      
      # PIMLE estimate: N * mean(y_iu / d_i)
      pimle.vec = ard.sample[,16] / d.killworth.sample # individual ratios y_iu/d_i
      pimle.size.est = N * mean(pimle.vec, na.rm = T) # average indiv ratio * N
      
      # PIMLE with true degrees (assume we know them)
      pimle.size.dtrue = N * mean(ard.sample[,16]/d.true.sample, na.rm = T)
      
      # PIMLE censored estimate (set Infs and any ratios > 1 to 1)
      ## We can censor the ratio at 1 (helps with variance too)
      # i.e. any individual ratios y_iu/d_i above 1 are set to 1
      # when degree is estimated to be 0, then pimle.vec = Inf
      # this changes that value to 1.
      pimle.vec[pimle.vec > 1] = 1
      pimle.cens.est = N * mean(pimle.vec, na.rm = T)
      
      # MoS estimate: N * mean_i (y_iu / d_i)
      mos.vec = ard.sample[,16]/d.mos.sample
      mos.size.est = N * mean(mos.vec, na.rm = T)
      
      # MoS censored estimate (set Infs and any ratios > 1 to 1)
      mos.vec[mos.vec > 1] = 1
      mos.cens.est = N * mean(mos.vec, na.rm = T)
      
      # store results in survey.data
      # for.this.survey = cbind(n.hidden, mle.size.est,
      #                         mle.size.est.known.d, pimle.size.est,
      #                         pimle.cens.est, mos.size.est)
      # survey.data = rbind(survey.data, for.this.survey)
      survey.data = survey.data %>%
        add_row(school = school.name,
                hidden = hidden.group,
                true.size = Nu,
                n.hidden = n.hidden,
                mle.size.est = mle.size.est,
                mle.size.est.known.d = mle.size.est.known.d,
                pimle.size.est = pimle.size.est,
                pimle.cens.est = pimle.cens.est,
                pimle.size.dtrue = pimle.size.dtrue,
                mos.size.est = mos.size.est,
                mos.cens.est = mos.cens.est)
    }
    school.hidden.data = school.hidden.data %>%
      add_row(school = school.name,
              N = N,
              hidden = hidden.group,
              true.size = Nu,
              d.ratio = d.ratio,
              n.frame = n.frame)
  }
  
  save(school.hidden.data, individual.data, survey.data,
       file = paste0(out.dir, outdatafilename))
}

dratios = dratios %>%
  mutate(`%N` = num(100*true.size/N, digits=1), d.ratio = num(d.ratio, digits=2)) %>%
  relocate(`%N`, .after = true.size) %>%
  rename(School = school, `School index` = school_ind,
         `Hidden group` = hidden, `Hidden group index` = hidden_ind,
         `True hidden group size` = true.size, `Degree ratio` = d.ratio)

school_inds = dratios$`School index`
hidden_inds = dratios$`Hidden group index`

# # run this once to generate the data for all 1600 cases
# (100 schools x 16 hidden/probe group sets for each school)
include_PIMLE_with_true_degrees(school_inds, hidden_inds, data.dir,
                                "fb100_sims.RData")
system("say The simulations have finished")


