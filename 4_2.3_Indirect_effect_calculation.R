#This R script is to provide the confidence interval for indirect effect via mediators.

#To improve efficiency, first calculate CI for each chain and save it locally and subsequently combine

#Method 1: mediation::mediate
#Load R packages
library(data.table)
library(dplyr)
library(parallel)
library(RMediation)

#File of significant after exposure_adjusted
dat.input <- fread(input = "./Adjusted_Significant_Neuro_Mediators_IBS.csv") %>% data.frame()

#read chain
Sig_chain <- dat.input$chain

#Initialize cores
cl <- makeCluster(40)

#Indirect effect via mediators calculation
indirect_effect_mediation = function(chain) {

  #load packages
  library(data.table)
  library(dplyr)
  library(bruceR)
  
  #read mv.dat
  mv.hardat <- fread(input = paste("./exposure_adjusted_harmonise_mvdat/mvdat_",chain,".csv",sep = "")) %>% data.frame()
  
  #set a data.frame for output
  sig_adj_ci <- data.frame(matrix(nrow = 1)) %>% select(-matrix.nrow...1.) %>% mutate(chain=chain) %>%
    mutate(BETA_indirect_mediation=NA) %>% mutate(LCI_indirect_mediation=NA) %>% mutate(UCI_indirect_mediation=NA) %>% mutate(Pval_indirect_mediation=NA) %>%
    mutate(Prop_BETA_indirect_mediation=NA) %>% mutate(Prop_LCI_indirect_mediation=NA) %>% mutate(Prop_UCI_indirect_mediation=NA) %>% mutate(Prop_Pval_indirect_mediation=NA)
  
  #effects calculation
  model.m <- lm(beta.mediator ~ -1 + beta.exposure, weights = 1/se.outcome^2, data = mv.hardat)
  model.y <- lm(beta.outcome ~ -1 + beta.exposure + beta.mediator, weights = 1/se.outcome^2, data = mv.hardat)
  mediation::mediate(model.m = model.m, model.y = model.y, sims = 3000, treat = "beta.exposure", mediator = "beta.metabolite", boot = T) -> Process.results

  #Write
  sig_adj_ci[1,2:5] <- c(Process.results$d1,Process.results$d1.ci,Process.results$d1.p)
  sig_adj_ci[1,6:9] <- c(Process.results$n1,Process.results$n1.ci[1:2],Process.results$n1.p)

  
  write.csv(sig_adj_ci, file = paste("./indirect_effect_single/mediation_",chain,"_indirect_effect_ci.csv",sep = ""),row.names = F)
  print("Saved")
}

parLapply(cl, Sig_chain, indirect_effect_PROCESS)

stopCluster(cl)

sig_adj_ci_total <- data.frame(c())

for (chain in Sig_chain) {
  
  dat <- fread(input = paste("./indirect_effect_single/mediation_",chain,"_indirect_effect_ci.csv",sep = "")) %>% data.frame()
  
  sig_adj_ci_total <- rbind(sig_adj_ci_total, dat)
  
  write.csv(sig_adj_ci_total, file = "./Adjusted_mediation_indirect_effect_ci.csv",row.names = F)
}

#read chain
dat.input %>% left_join(sig_adj_ci_total, by=c("chain"="chain")) -> dat.output

write.csv(dat.output, file = "./Indirect_effect_CI.csv") 

#Method 2: Product or Difference by bootstrap
indirect_effect_2 = function(chain) {

  #load packages
  library(data.table)
  library(dplyr)
  library("AER")
  library ("boot")
  
  #read mv.dat
  mv.hardat <- fread(input = paste("./exposure_adjusted_harmonise_mvdat/mvdat_",chain,".csv",sep = "")) %>% data.frame()
  
  #set a data.frame for output
  sig_adj_ci <- data.frame(matrix(nrow = 1)) %>% select(-matrix.nrow...1.) %>% mutate(chain=chain) %>%
    mutate(BETA_indirect_product=NA) %>% mutate(LCI_indirect_product=NA) %>% mutate(UCI_indirect_product=NA) %>%
    mutate(BETA_indirect_product_proportion=NA) %>% mutate(LCI_indirect_product_proportion=NA) %>% mutate(UCI_indirect_product_proportio=NA) %>%
    mutate(BETA_indirect_difference=NA) %>% mutate(LCI_indirect_difference=NA) %>% mutate(UCI_indirect_difference=NA) %>%
    mutate(BETA_indirect_difference_proportion=NA) %>% mutate(LCI_indirect_difference_proportion=NA) %>% mutate(UCI_indirect_difference_proportion=NA)
  
  ## 1 - Code for two-step MR - Product of coefficients method
  # Total effect of Exposure 1 (Neuro) on binary Y (IBS)
  total_out <- ivreg(beta.outcome ~ beta.exposure, weights = 1/(se.outcome^2), data = mv.hardat)
  
  # Effect of X on M (metabolite or microbe, Exposure 2)
  exp_mediator <- ivreg(beta.mediator ~  beta.exposure, weights = 1/(se.mediator^2), data = mv.hardat)
  
  # Effect of M (Exposure 2) on Y
  mediator_out <- ivreg(beta.outcome ~ beta.mediator + beta.exposure, weights = 1/(se.outcome^2), data = mv.hardat)
  
  # Indirect effect (including bootstrapped confidence intervals)
  indirect_effect_product <- exp_mediator$coef[2]*mediator_out$coef[2]
  
  # Write into results data.frame
  sig_adj_ci[1,2] <- indirect_effect_product
  
  set.seed(1234)
  indirect_product <- function(data, indices) { 
    sample <- data[indices,]
    exp_mediator <- ivreg(beta.mediator ~  beta.exposure, weights = 1/(se.mediator^2), data = sample)
    mediator_out <- ivreg(beta.outcome ~ beta.mediator + beta.exposure, weights = 1/(se.outcome^2), data = sample)
    return(exp_mediator$coef[2]*mediator_out$coef[2])
  }
  
  print("Running product boot")
  boot_indirect_product <- boot(mv.hardat, indirect_product, R=1000)
  print("Boot done!")
  
  boot_ci=boot.ci(boot.out= boot_indirect_product, type = c("norm"))
  sig_adj_ci[1,3:4] <- boot_ci$normal[2:3]
  
  # Proportion mediated (including bootstrapped confidence intervals)
  proportion_effect_product <- indirect_effect_product/total_out$coef[2]
  sig_adj_ci[1,5] <-  proportion_effect_product
  
  set.seed(1234)
  proportion_product <- function(data, indices) { 
    sample <- data[indices,]
    total_out <- ivreg(beta.outcome ~ beta.exposure, weights = 1/(se.outcome^2), data = sample)
    exp_mediator <- ivreg(beta.mediator ~  beta.exposure, weights = 1/(se.mediator^2), data = sample)
    mediator_out <- ivreg(beta.outcome ~ beta.mediator + beta.exposure, weights = 1/(se.outcome^2), data = sample)
    return(exp_mediator$coef[2]*mediator_out$coef[2]/total_out$coef[2])
  }
  
  print("Running product_proportion boot")
  boot_proportion_product <- boot(mv.hardat, proportion_product, R=1000)
  print("Boot done!")
  
  boot_ci=boot.ci(boot.out= boot_proportion_product, type = c("norm"))
  sig_adj_ci[1,6:7] <- boot_ci$normal[2:3]
  
  ## 2 - Code for two-step MR - Difference in coefficients method
  total_out <- ivreg(beta.outcome ~ beta.exposure, weights = 1/(se.outcome^2), data = mv.hardat)
  
  # Direct effect of X on Y controlling for M
  direct <- ivreg(beta.outcome ~ beta.exposure + beta.mediator, weights = 1/(se.outcome^2), data = mv.hardat)
  
  # Indirect effect (including bootstrapped confidence intervals)
  indirect_effect_difference <- total_out$coef[2]-direct$coef[2]
  sig_adj_ci[1,8] <-  indirect_effect_difference
  
  set.seed(1234)
  indirect_difference <- function(data, indices) { 
    sample <- data[indices,]
    total_out <- ivreg(beta.outcome ~ beta.exposure, weights = 1/(se.outcome^2), data = sample)
    direct <- ivreg(beta.outcome ~ beta.exposure + beta.mediator, weights = 1/(se.outcome^2), data = sample)
    return(total_out$coef[2]-direct$coef[2])
  }
  print("Running difference boot")
  boot_indirect_difference <- boot(mv.hardat,indirect_difference,R=1000)
  print("Boot done!")
  
  boot_ci=boot.ci(boot.out= boot_indirect_product, type = c("norm"))
  sig_adj_ci[1,9:10] <- boot_ci$normal[2:3]
  
  # Proportion mediated (including bootstrapped confidence intervals)
  proportion_effect_difference <- indirect_effect_difference/total_out$coef[2]
  sig_adj_ci[1,11] <-  proportion_effect_difference
  
  set.seed(1234)
  proportion_difference <- function(data, indices) { 
    sample <- data[indices,]
    total_out <- ivreg(beta.outcome ~ beta.exposure, weights = 1/(se.outcome^2), data = sample)
    direct <- ivreg(beta.outcome ~ beta.exposure + beta.mediator, weights = 1/(se.outcome^2), data = sample)
    return((total_out$coef[2]-direct$coef[2])/total_out$coef[2])
  }
  
  print("Running difference_proportion boot")
  boot_proportion_difference <- boot(mv.hardat, proportion_difference, R=1000)
  print("Boot done!")
  
  boot_ci=boot.ci(boot.out= boot_proportion_difference, type = c("norm"))
  sig_adj_ci[1,12:13] <- boot_ci$normal[2:3]
  
  
  write.csv(sig_adj_ci, file = paste("./indirect_effect_single/2_",chain,"_indirect_effect_ci.csv",sep = ""),row.names = F)
  print("Saved")
}

parLapply(cl, Sig_chain, indirect_effect_2)

stopCluster(cl)

sig_adj_ci_total <- data.frame(c())

for (chain in Sig_chain) {
  
  dat <- fread(input = paste("./indirect_effect_single/2_",chain,"_indirect_effect_ci.csv",sep = "")) %>% data.frame()
  
  sig_adj_ci_total <- rbind(sig_adj_ci_total, dat)
  
  write.csv(sig_adj_ci_total, file = "./Adjusted_2_significant_indirect_effect_ci.csv",row.names = F)
}

#read chain
dat.output %>% left_join(sig_adj_ci_total, by=c("chain"="chain")) -> dat.output

write.csv(dat.output, file = "./Indirect_effect_CI.csv") 

##Method 3: delta method
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

dat.output <- medci_table(dat.output, name.b1="N_M_b",name.se1 = "N_M_se", name.b2 = "BETA_Met_Out_adjusted", name.se2 = "SE_Met_Out_adjusted")

dat.output <- prop_ci_table(dat.output, 
                               name.b1="N_M_b",name.se1 = "N_M_se", 
                               name.b2 = "BETA_Met_Out_adjusted", name.se2 = "SE_Met_Out_adjusted",
                               name.b.total="b_total",name.se.total="se_total")

write.csv(dat.output, file = "./Indirect_effect_CI.csv") 

