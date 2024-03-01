#This script is to calculate effect value and its CI (not adjusting expsoure).

#Method: RMediation
#exp, med, out
#N_M_b, N_M_se disorders to mediators
#M_I_b, M_I_se mediators to outcome
#b_total, se_total disorders to outcome

setwd(".")

library(data.table)
library(dplyr)
library(RMediation)

#Import candidate mediators
candidate <- fread(input = "Candidate_neuro_mediators_ibs.csv") %>% data.frame()

#Calculate indirect effects
medci_table = function(dat, name.b1, name.se1, name.b2, name.se2, type="asymp") { 
  
  estimate <- c() #estimates
  se.indirect <- c()
  lci <- c()
  uci <- c()
  
  #dat refers to the data.frame that contains beta and se for multiple pairs (exposure-outcome)
  for (i in 1:nrow(dat)) {
    
    mu.x <- dat[i,name.b1]
    
    mu.y <- dat[i,name.b2]
    
    se.x <- dat[i,name.se1]
    
    se.y <- dat[i,name.se2]
      
    medci(mu.x = mu.x, mu.y = mu.y, se.x = se.x, se.y = se.y, type = type) -> result
    
    estimate <- c(estimate, result$Estimate)
    se.indirect <- c(se.indirect, result$SE)
    lci <- c(lci, result$`95% CI`[1])
    uci <- c(uci, result$`95% CI`[2])
  }
  
  #output
  dat %>% mutate(estimate.delta=estimate) %>% mutate(se.delta=se.indirect) %>% mutate(lci.delta=lci) %>% mutate(uci.delta=uci) -> dat
  
  return(dat)
}

#Calculate proportion of indirect effects
prop_ci_table = function(dat, name.b1, name.se1, name.b2, name.se2, name.b.total, name.se.total, type="asymp") { 
  
  estimate <- c() #estimates
  se.prop <- c()
  lci <- c()
  uci <- c()
  
  #dat refers to the data.frame that contains beta and se for multiple pairs (exposure-outcome)
  for (i in 1:nrow(dat)) {
    
    mu.x <- dat[i,name.b1]
    
    mu.y <- dat[i,name.b2]
    
    mu.z <- dat[i,name.b.total]
    
    se.x <- dat[i,name.se1]
    
    se.y <- dat[i,name.se2]
    
    se.z <- dat[i,name.se.total]
    
    ci(mu=c(b1=mu.x,b2=mu.y,b3=mu.z), Sigma = c(se.x,0,0,se.y,0,se.z), quant = ~b1*b1/b3, type = type) -> result
    
    estimate <- c(estimate, result$Estimate)
    se.prop <- c(se.prop, result$SE[1,1])
    lci <- c(lci, result$`97.5% CI`[1])
    uci <- c(uci, result$`97.5% CI`[2])
  }
  
  #output
  dat %>% mutate(prop.estimate.delta=estimate) %>% mutate(prop.se.delta=se.prop) %>% mutate(prop.lci.delta=lci) %>% mutate(prop.uci.delta=uci) -> dat
  
  return(dat)
}

candidate <- medci_table(dat = candidate, name.b1="N_M_b",name.se1 = "N_M_se",name.b2 = "M_I_b",name.se2 = "M_I_se")

candidate <- prop_ci_table(dat = candidate, name.b1="N_M_b",name.se1 = "N_M_se",name.b2 = "M_I_b",name.se2 = "M_I_se", name.b.total="b_total", name.se.total="se_total")

write.csv(candidate, file = "Not_Adjusted_Effect_CI_candidate.csv")
