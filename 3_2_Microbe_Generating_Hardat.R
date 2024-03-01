#This R script is to generate harmonized data for of-interest disorders (exposure) - microbes (outcome) and microbes (exposure) - IBS (outcome).

#Note that the clump data was saved in Clump_5e_6_new

#Note that this harmonized data was not adjusted (not involved MVMR).

#Load parallel packages
library(parallel)

#Multithread
cl <- makeCluster(40)

#Set microbiome group
micros <- list.files(path = "./MiBiGen", pattern = "*.txt.gz")
micros <- gsub(pattern = ".summary.txt.gz", replacement = "", micros)

#Set function
Hardat_microbiome = function(j) {
  
  #Load packages
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  
  #Loading format functions and fread functions
  source("./0_fread_adj_format.R")
  print("Loaded fread_adj_format.R")
  
  #Step 1: of-interest disorders (exposure) - microbes (outcome)
  
  #Set diseases
  forward.exp <- c("ADHD","BIP","MDD","HOA","INSOM") #Identified significant
  
  #Set format function for metabolites
  format_micro=function(dataset, snps = NULL, type = "exposure") {
    dataset <- format_data(dataset, snps = snps, type = type,
                           phenotype_col = "bac",
                           snp_col = "rsID",
                           effect_allele_col = "eff.allele",
                           other_allele_col = "ref.allele",
                           #eaf_col = "effect_allele_frequency",
                           beta_col = "beta", 
                           se_col = "SE",
                           pval_col = "P.weightedSumZ",
                           chr_col = "chr", pos_col = "bp",
                           samplesize_col = "N")
    return(dataset)
  } #From class.Actinobacteria.id.419
  
  #Generating harmonised data per microbiome outcome
  #Read outcome data
  dat <- fread(input = paste("./MiBiGen/",j,".summary.txt.gz",sep = "")) %>% data.frame()
  print(j)
  
  for (i in forward.exp) {
    
    print(i)
    
    #If harmonised data exists, skip this script
    exist <- list.files(path = "./twosampleMR_microbiome_harmonise/",pattern = paste(i,j,sep = "_"))
    if(length(exist)>0) {
      print("Harmonised data exists")
      next}
    
    #If not existing, run the jobs
    print(paste("Harmonising",i,j,"..."))
    
    #Import clumped data
    clump <- fread(input = paste("./Clump_5e_6_new/",i,".csv",sep = ""))
    clump <- data.frame(clump)
    
    #Harmonise
    #Ensure there is SNP for harmonising
    dat %>% dplyr::filter(rsID %in% clump$SNP) -> sub.snp #Note here rsID is from microbiome ss in an identical format.
    if(nrow(sub.snp)==0) {
      
      #No SNP for harmonising
      #Copy har.dat from same exposure
      format.file <- list.files(path = "./twosampleMR_microbiome_harmonise", pattern = paste(i,"_",sep = ""), full.names = T)
      har.dat <- fread(input = format.file[1]) %>% data.frame()
      har.dat[1,] <- NA
      har.dat <- har.dat[1,]
        
    } else {
      
      #There is snp.
      dat %>% format_micro(snps=clump$SNP, type="outcome") %>% 
        harmonise_data(exposure_dat=clump) %>%
        mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
        mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
    }

    #Save the harmonised data
    write.csv(har.dat,file = paste("./twosampleMR_microbiome_harmonise/",paste(i,j,sep = "_"),".csv",sep = ""))
  }
  
  #Step 2: microbes (exposure) - IBS (outcome)
  
  #Set IBS SS
  forward.out <- c("IBS_3")
  
  #Read outcome data
  dat.out <- fread_IBS(forward.out)
  
  #If harmonised data exists, skip this script
  exist <- list.files(path = "./twosampleMR_microbiome_harmonise/",pattern = paste(j,forward.out,sep = "_"))
  if(length(exist)>0) {
    print("Harmonised data exists")
    next}
  
  #If not existing, run the jobs
  print(paste("Harmonising",j,forward.out,"..."))
  
  #Import clumped data
  clump <- fread(input = paste("./Clump_5e_6_new/microbiome/",j,".summary.txt.gz.csv",sep = "")) %>% data.frame()
  print("Successfully read clump data for microbiome")
  
  #Ensure having at least 1 SNP
  har.dat <-  tryCatch(
    {
      dat.out %>% format(snps=clump$SNP, type="outcome", asformat=forward.out) %>% 
        harmonise_data(exposure_dat=clump) %>%
        mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
        mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
    },
    error=function(e) {
      format.file <- list.files(path = "./twosampleMR_microbiome_harmonise", pattern = paste(j,"_",sep = ""), full.names = T)
      har.dat <- fread(input = format.file[1]) %>% data.frame()
      har.dat[1,] <- NA
      har.dat <- har.dat[1,]
    }
  )
  
  #Save the harmonised data
  write.csv(har.dat,file = paste("./twosampleMR_microbiome_harmonise/",paste(j,forward.out,sep = "_"),".csv",sep = ""))
}

parLapply(cl, micros, Hardat_microbiome)

stopCluster(cl)






