#This R script is to re-clump SNP for exposure_adjusted mediation effect.

#Copy this file 5 times and re-name them to 4_2_Reclump_1 to 4_2_Reclump_5 (to improve efficiency).

#Step 1 of whole mediation analysis: Combine SNP and Re-clump

library(data.table)
library(plyr)
library(dplyr)
library(ieugwasr)
library(scriptName)

#Set the environment where the clump data is
setwd("./Clump_5e_6_new")

#Import the file of bidirectional two-sample MR
bi_ts_mr <- fread(input = "./Candidate_neuro_mediators_ibs.csv") %>% data.frame()

#mutate chain
Chain <- c(paste(bi_ts_mr[,"Exposure_Neuro"],bi_ts_mr[,"Mediator"],bi_ts_mr[,"Outcome_IBS"],sep=" "))

print(Chain)

filename <- current_filename()
print(filename)

filename %>% gsub(pattern = "4_2_Reclump_", replacement = "") %>% gsub(pattern = ".R", replacement = "") %>% as.numeric() -> n

for (chain in Chain[(5*n-4):(5*n)]) {
  
  chain <- strsplit(chain, " ")
  
  #Exposure = j, Mediator = k, Outcome = i
  j <- chain[[1]][1]
  k <- chain[[1]][2]
  i <- chain[[1]][3]
  print(paste(j,k,i))
  
  #If harmonised data exists, skip this script
  exist <- list.files(path = "./Clump_5e_6_new/exposure_adjusted",pattern = paste(j,k,i,sep = "_"))
  if(length(exist)>0) {
    print("Re-clump data exists")
    next}
  
  #Clump fread function
  clump_fread = function(exposure) {
    
    path_to_clump <- list.files(path = c(".","./metabolites","./microbiome"), pattern = exposure, full.names = T)
    
    clump <- fread(input = path_to_clump) %>% data.frame()
    
    return(clump)
  }
  
  #import clump_dat
  clump_exp <- clump_fread(j)
  clump_mediator <- clump_fread(k)
  
  #Combine two sets of clumped SNP
  clump_com <- rbind.fill(clump_exp,clump_mediator)
  
  #adjust the exposure and further clump (to ensure independence)
  clump_com %>% mutate(exposure=j) %>% mutate(id.exposure=clump_exp[1,"id.exposure"]) -> clump_com
  
  #Re-clump function
  clump_new = function(dataset) {
    ieugwasr::ld_clump_local(dat = dplyr::tibble(rsid=dataset$SNP, pval=dataset$pval.exposure),
                             clump_kb = 5000, clump_r2 = 0.001, clump_p = 1,
                             bfile = "./1kg.v3/EUR",
                             plink_bin = "./plink_linux_x86_64_20230116/plink") -> snpafterclump
    clumpset <- subset(dataset, SNP %in% snpafterclump$rsid)
    return(clumpset)
  }
  
  #re-clump
  clump_com <- clump_com %>% clump_new() %>% mutate(chr_pos=paste(chr.exposure,pos.exposure,sep = ":"))
  
  #Save it locally
  write.csv(clump_com, file = paste("./Clump_5e_6_new/exposure_adjusted/",paste(j,k,i,sep = "_"),".csv",sep = ""), row.names = F)
  
  print(paste(j,k,i,"Re-clump Done"))
}
