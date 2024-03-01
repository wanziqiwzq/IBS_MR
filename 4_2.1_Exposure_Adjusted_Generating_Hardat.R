#This R script is to generate mvmr dat based on re-clumped SNP.

#Having saved re_clumped SNP for all chains(exposure_mediators_outcome)

#Multithread
library(parallel)
library(data.table)
library(dplyr)

#initialise cores
cl <- makeCluster(40)

#Import the file of bidirectional two-sample MR
bi_ts_mr <- fread(input = "./Candidate_neuro_mediators_ibs.csv") %>% data.frame()

#To improve the efficiency, run per exp_out dat
pair_exp_out <- c(paste(bi_ts_mr[,"Exposure_Neuro"],bi_ts_mr[,"Outcome_IBS"],sep=" "))

pair_exp_out <- unique(pair_exp_out)

print(pair_exp_out)

#Set a data.frame for beta_adjusted
beta_adjust = function(exp_out) {

  #Load R packages
  library(data.table)
  library(plyr)
  library(dplyr)
  library(TwoSampleMR)
  
  exp_out <- strsplit(exp_out, " ")
  
  #Exposure = exp, Outcome = out, Mediator = med
  exp <- exp_out[[1]][1]
  out <- exp_out[[1]][2]
  
  #Step 1: Get all mediators for this exp_out chain
  #Import the file of bidirectional two-sample MR
  bi_ts_mr <- fread(input = "./Candidate_neuro_mediators_ibs.csv") %>% data.frame()
  
  bi_ts_mr[which(bi_ts_mr[,"Exposure_Neuro"]==exp & bi_ts_mr[,"Outcome_IBS"]==out),3] -> mediators
  
  #Skip tasks done
  exist <- list.files(path = "./exposure_adjusted_harmonise_mvdat",pattern = paste(paste("mvdat",exp,sep = "_"),"*",out,sep = ""))
  if(length(exist)==length(mediators)) {
    print("mvdat data exists")
    next}
  
  #Set the environment where the clump data is
  setwd(".")
  
  #Load fread and format function
  source("./0_fread_adj_format.R")
  
  #dat.exp      #dat.out      #dat.med
  extract_snp = function(dataset, asformat, snps=NULL, type, clump = clump_com){
    
    if(grepl("GCST",asformat)) {asformat <- "SERO"} 
    
    if(grepl("family.",asformat) | grepl("order.",asformat) | grepl("genus.",asformat) | grepl("class.",asformat) | grepl("phylum.",asformat)) {
      asformat <- "micro"}
    
    f=getFunction(paste("format",asformat,sep = "_"))
    dataset %>% f(snps = snps, type = type) -> dataset
    return(dataset)

  }
  
  #Import outcome ss
  dat.out <- fread_IBS(out)
  
  #Import exposure ss
  dat.exp <- fread_adj(exp) 
  
  for (med in mediators) {   #For each mediator
    
    #Skip task done
    exist <- list.files(path = "./exposure_adjusted_harmonise_mvdat",pattern = paste("mvdat",exp,med,out,sep = "_"))
    if(length(exist)==1) {
      print("mvdat data exists")
      next}
    
    #Import mediator ss
    path_to_med <- list.files(path = c("./metabolites","./MiBiGen"), pattern = med, full.names = T)
    
    dat.med <- fread(input = path_to_med) %>% data.frame() %>% mutate(Pheno=med)
    
    #Import clump_com
    clump_com <- fread(input = paste("./Clump_5e_6_new/exposure_adjusted/",paste(exp,med,out,sep = "_"),".csv",sep = "")) %>% data.frame()
    
    #Extract SNP from outcome data (save this to harmonise)
    extract_snp(dat.out, asformat = out, snps = clump_com$SNP, type = "outcome") -> snp_out
    print("SNP in outcome")
    
    #Prepare mv.dat
    snp_out %>% dplyr::select(SNP,beta.outcome,se.outcome) -> mv.dat
    
    #Extract SNP from exposure dat and harmonise with outcome and output beta_se into mv.dat
    extract_snp(dat.exp, asformat = exp, snps = clump_com$SNP, type = "exposure") %>% harmonise_data(outcome_dat = snp_out) %>%
      subset(mr_keep==T) %>% dplyr::select(SNP,beta.exposure,se.exposure) %>% inner_join(mv.dat, by=c("SNP"="SNP")) -> mv.dat
    print("SNP in exposure")
    
    #Extract SNP from mediator dat and harmonise with outcome and output beta_se into mv.dat
    extract_snp(dat.med, asformat = med, snps = clump_com$SNP, type = "exposure") %>% harmonise_data(outcome_dat = snp_out) %>%
      subset(mr_keep==T) %>% mutate(beta.mediator=beta.exposure) %>% mutate(se.mediator=se.exposure) %>%
      dplyr::select(SNP,beta.mediator,se.mediator) %>% inner_join(mv.dat, by=c("SNP"="SNP")) -> mv.dat
    print("SNP in mediators")
    
    #Discard any SNP with NA
    for (p in 1:nrow(mv.dat)) {
      for (q in 1:ncol(mv.dat)) {
        if (is.na(mv.dat[p,q])) {
          print(c(p,q))
          mv.dat <- mv.dat[-p,]
        }
      }
    }
    
    #Save the mv.dat locally
    write.csv(mv.dat, file = paste(paste("./exposure_adjusted_harmonise_mvdat/mvdat_",exp,sep = ""),med,paste(out,".csv",sep = ""),sep = "_"))
    
    print(head(mv.dat))
  } #Screen all mediators for pair_exp_out

}

parLapply(cl, pair_exp_out, beta_adjust)

stopCluster(cl)

