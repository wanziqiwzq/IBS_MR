#This R script is to generate harmonized data for of-interest disorders (exposure) - metabolites (outcome) and metabolites (exposure) - IBS (outcome).

#Note that the clump data was saved in Clump_5e_6_new

#Note that this harmonized data was not adjusted (not involved MVMR).

#Load parallel packages
library(parallel)

#Multithread
cl <- makeCluster(40) #set the number of clusters per service

#Set metabolites group
metabolites <- paste("GCST90",199621:201020,sep = "")

Hardat_metabolites = function(j) {
  
  #Load packages
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  
  #Loading format functions and fread functions
  source("./0_fread_adj_format.R")
  print("Loaded fread_adj_format.R")
  
  #Step 1: of-interest disorders (exposure) - metabolites (outcome)
  
  #Set diseases
  forward.exp <- c("ADHD","BIP","MDD","HOA","INSOM") #Identified significant
  
  #Set format function for metabolites
  format_SERO=function(dataset, snps = NULL, type = "exposure") {
    dataset <- format_data(dataset, snps = snps, type = type,
                           phenotype_col = "Pheno",
                           snp_col = "variant_id",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           eaf_col = "effect_allele_frequency",
                           beta_col = "beta", 
                           se_col = "standard_error",
                           pval_col = "p_value",
                           chr_col = "chromosome", pos_col = "base_pair_location",
                           samplesize_col = "samplesize")
    return(dataset)
  } #From GCST90200422
  
  #Generating harmonised data per metabolite outcome
  #Read outcome data
  dat <- fread(input = paste("./metabolites/",j,"_buildGRCh38.tsv.gz",sep = "")) %>% data.frame()
  print(j)
  
  for (i in forward.exp) {
    
    print(i)
    
    #If harmonised data exists, skip this script
    exist <- list.files(path = "./twosampleMR_metabolites_harmonise/",pattern = paste(i,j,sep = "_"))
    if(length(exist)>0) {
      print("Harmonised data exists")
      next}
    
    #If not existing, run the jobs
    print(paste("Harmonising",i,j,"..."))
    
    #Import clumped data
    clump <- fread(input = paste("./Clump_5e_6_new/",i,".csv",sep = ""))
    clump <- data.frame(clump)
    
    #Harmonise
    dat %>% format_SERO(snps=clump$SNP, type="outcome") %>% 
      harmonise_data(exposure_dat=clump) %>%
      mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
      mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat

    #Save the harmonised data
    write.csv(har.dat,file = paste("./twosampleMR_metabolites_harmonise/",paste(i,j,sep = "_"),".csv",sep = ""))
  }
  
  #Step 2: metabolites (exposure) - IBS (outcome)
  
  #Set IBS SS
  forward.out <- c("IBS_3")
  
  #Read outcome data
  dat.out <- fread_IBS(forward.out)
  
  #If harmonised data exists, skip this script
  exist <- list.files(path = "./twosampleMR_metabolites_harmonise/",pattern = paste(j,forward.out,sep = "_"))
  if(length(exist)>0) {
    print("Harmonised data exists")
    next}
  
  #If not existing, run the jobs
  print(paste("Harmonising",j,forward.out,"..."))
  
  #Import clumped data
  clump <- fread(input = paste("./Clump_5e_6_new/metabolites/",j,".csv",sep = "")) %>% data.frame()
  print("Successfully read clump data for metabolites")
  
  #Harmonise
  
  dat.out %>% format(snps=clump$SNP, type="outcome", asformat=forward.out) %>% 
    harmonise_data(exposure_dat=clump) %>%
    mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
    mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
  
  #Save the harmonised data
  write.csv(har.dat,file = paste("./twosampleMR_metabolites_harmonise/",paste(j,forward.out,sep = "_"),".csv",sep = ""))

}

parLapply(cl, metabolites, Hardat_metabolites)

stopCluster(cl)

