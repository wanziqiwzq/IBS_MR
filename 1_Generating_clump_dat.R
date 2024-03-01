#This R script is for generating clumped SNP for 23 of-interest disorders and 3 IBS ss.

##Some notions:
#ALZ and PD ss do not contain rsids therefore they need Munge to find rsids by chr_pos information
#IBS_Finn(IBS_3) is from hg38, thus CHR:BP is different from other ss.
#IBS_2 is from Bellygenes (by request to Dr.Mauro D’Amato)

#Set environment
setwd("/home/wan/IBS") #delete this during submission

#Install R packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ieugwasr)
library(MungeSumstats)
library("SNPlocs.Hsapiens.dbSNP155.GRCh37") #Loading this package locally if using laptop
library("BSgenome.Hsapiens.1000genomes.hs37d5") #Loading this package locally if using laptop

#Set clump function
clump_new = function(dataset) {
  ieugwasr::ld_clump_local(dat = dplyr::tibble(rsid=dataset$SNP, pval=dataset$pval.exposure),
                           clump_kb = 5000, clump_r2 = 0.001, clump_p = 1,
                           bfile = "./1kg.v3/EUR", #download this from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
                           plink_bin = "./plink_linux_x86_64_20230116/plink") -> snpafterclump
  clumpset <- subset(dataset, SNP %in% snpafterclump$rsid)
  return(clumpset)
}

#Import format and fread functions
source("0_fread_adj_format.R")

#Set 23 of-interest disorders
forward.exp <- c("ALZ","PD","ADHD","ANERV","ANXIETY","ASD","BIP","MDD","PANIC","OCD","HOA",
                 "TS","SAPNEA","SCZ","PTSD",
                 "CHRON","INSOM","DSLEEP","SDURA","NAP",
                 "STROKE","IANE","LS")

#Set 3 IBS Summary statistics
forward.out <- c("IBS_1","IBS_2","IBS_3")

#Do the following jobs to generate clumped data for each Summary statistics
for (i in c(forward.out,forward.exp)) {
  
  print(i)
  
  #For re-running R script, clumped data do not run.
  a=list.files(path = "./Clump_5e_6_new", pattern = i)
  if(length(a)==1) {next
    print("Clump data exists, Skip")}
  
  #read the corresponding ss
  if(i %in% forward.out) { 
    dat <- fread_IBS(i)
    } else{
    dat <- fread_adj(i) 
    }
  print("Freading the whole SS done")
  
  #give the corresponding format function
  f=getFunction(paste("format",i,sep = "_"))
  
  #Find rs id and format, clump
  if(i %in% c("ALZ","PD")) { #There is no ids for SNPs in these two SS. Find SNP ids according to their positions.
    if(i=="ALZ") {
      dat %>% subset(p_value<5e-06) -> dat
      
      print("Find rsids")
      dat %>% dplyr::select(chromosome,base_pair_location,other_allele,effect_allele,effect_allele_frequency,beta,standard_error,p_value) %>% 
        MungeSumstats::format_sumstats(ref_genome = "GRCh37", dbSNP = 155, return_data = T, return_format = "data.table") -> dat_rs
      }
    if(i=="PD") {
      dat %>% subset(p<5e-06) -> dat
      
      print("Find rsids")
      dat %>% dplyr::select(CHR,POS,A2,A1,freq,b,se,p) %>% 
        MungeSumstats::format_sumstats(ref_genome = "GRCh37", dbSNP = 155, return_data = T, return_format = "data.table") -> dat_rs
    }
    #Add the position into the converted dat_rs
    dat_rs %>% mutate(chr_pos=paste(CHR,BP,sep = ":")) -> dat_rs
    
    #Combine the information in the original dataset with current one by identical chr_pos
    print("Combine the information by identical chr_pos")
    dat %>% subset(chr_pos %in% dat_rs$chr_pos) %>% 
      left_join(dat_rs[,c("chr_pos","SNP")], by=c("chr_pos"="chr_pos")) %>% 
      f() %>% clump_new() -> clump
  
  } else{ #format and clump
    f(dat) %>% subset(pval.exposure<5e-6) %>% clump_new() -> clump}
  
  #Save clump data locally
  write.csv(clump, file = paste("./Clump_5e_6_new/",i,".csv",sep = ""))
}
