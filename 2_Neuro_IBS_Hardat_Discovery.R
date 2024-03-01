#This R script is for generating harmonised data for bidirectional Two-Sample MR between IBS and neurological/psychological diseases.

#The interpretation is NOT on this script. (See 2_Neuro_IBS_MR_analysis.R and 2_Neuro_IBS_MR_PRESSO.R)

#IBS SS is from Finngen for discovery phase.

#There were 23 of-interest diseases.

##Some notions:
#ALZ and PD ss do not contain rsids therefore they need Munge to find rsids by chr_pos information
#IBS_Finn is from hg38, thus cannot simply find information by matching rs information by their position when analyzing IBS_Fin and ALZ/PD.

#Install R packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(MungeSumstats)
library("SNPlocs.Hsapiens.dbSNP155.GRCh37") #Loading this package locally if using laptop
library("BSgenome.Hsapiens.1000genomes.hs37d5") #Loading this package locally if using laptop

#Set environment
setwd(".")

#Loading format functions and fread functions
source("0_fread_adj_format.R")
print("Loaded fread_adj_format.R")

#Set 23 disorders
forward.exp <- c("ALZ","PD","ADHD","ANERV","ANXIETY","ASD","BIP","MDD","PANIC","OCD","HOA",
                 "TS","SAPNEA","SCZ","PTSD",
                 "CHRON","INSOM","DSLEEP","SDURA","NAP",
                 "STROKE","IANE","LS")

#Set IBS SS
forward.out <- c("IBS_3")

#Generating harmonised data per IBS outcome
for (j in forward.out) { #IBS as outcome
  #Read outcome data
  dat <- fread_IBS(j)
  
  for (neuro in forward.exp) {
   
    #If harmonised data exists, skip this script
    exist <- list.files(path = "./Neuro_IBS/harmonise",pattern = paste(neuro,j,sep = "_"))
    if(length(exist)>0) {next}
    
    #If not existing, run the jobs
    print(paste("Harmonising",neuro,j,"..."))
    
    #Import clumped data
    clump <- fread(input = paste("./Clump_5e_6_new/",neuro,".csv",sep = "")) %>% as.data.frame()
    
    #Harmonise
    dat %>% format(snps=clump$SNP, type="outcome", asformat=j) %>% 
      harmonise_data(exposure_dat=clump) %>%
      mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
      mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
    
    #Save the harmonised data
    write.csv(har.dat,file = paste("./Neuro_IBS/harmonise/",paste(neuro,j,sep = "_"),".csv",sep = ""))
  }
}

##Generating harmonised data per neuro outcome
for (neuro in forward.exp) {
  #Read outcome data
  dat <- fread_adj(neuro)
  
  for (j in forward.out) { #IBS as exposure
    #If harmonised data exists, skip this script
    exist <- list.files(path = "./Neuro_IBS/harmonise",pattern = paste(j,neuro,sep = "_"))
    if(length(exist)>0) {next}
    
    #If not existing, run the jobs
    print(paste("Harmonising",j,neuro,"..."))
    
    #Import clumped data
    clump <- fread(input = paste("./Clump_5e_6_new/",j,".csv",sep = "")) %>% as.data.frame()
    
    #Harmonise
    if(neuro %in% c("PD","ALZ")) {
      if(j == "IBS_3") {
        #Converting hg38 to hg37
        clump %>% dplyr::select(SNP,other_allele.exposure,effect_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure) %>% 
          MungeSumstats::format_sumstats(ref_genome = "GRCh37", dbSNP = 155, return_data = T, return_format = "data.table") -> clump_hg37
        
        #combine it to previous clump data
        clump_hg37 %>% mutate(chr_pos=paste(CHR,BP,sep = ":")) -> clump_hg37
        
      } 
      
      #Extracting clump$SNP in neuro_ss by identical chr_pos
      print("Extracting clump$SNP in neuro_ss by identical chr_pos")
      dat %>% subset(chr_pos %in% clump_hg37$chr_pos) %>% 
        left_join(clump_hg37[,c("SNP","chr_pos")], by=c("chr_pos"="chr_pos")) %>% 
        format(snps=clump$SNP, type="outcome", asformat=neuro) %>% 
        harmonise_data(exposure_dat=clump) %>%
        mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
        mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
      
    } else {
      
      dat %>% format(snps=clump$SNP, type="outcome", asformat=neuro) %>% 
        harmonise_data(exposure_dat=clump) %>%
        mutate(outcome.exposure.p=pval.outcome/pval.exposure>1e4) %>%
        mutate(mr_keep=mr_keep & outcome.exposure.p) -> har.dat
    }
    
    #Save the harmonised data
    write.csv(har.dat,file = paste("./Neuro_IBS/harmonise/",paste(j,neuro,sep = "_"),".csv",sep = ""))
  }
}

