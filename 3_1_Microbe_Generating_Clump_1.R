#This R script is to generate clumped data for microbiome

#This is an example script. Copy this to make 21 scripts and rename them by changing the last number from 1 to 21.

#Save into /home/wan/IBS/Clump_5e_6_new/microbiome

#To improve efficiency, define 10 microbiome as a group. (Clump_data can not be run in parallel function).

#Set main environment
setwd("/home/wan/IBS")

#Install R packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ieugwasr)
library(scriptName)

#Set clump function
clump_new = function(dataset) {
  ieugwasr::ld_clump_local(dat = dplyr::tibble(rsid=dataset$SNP, pval=dataset$pval.exposure),
                           clump_kb = 5000, clump_r2 = 0.001, clump_p = 1,
                           bfile = "./1kg.v3/EUR",
                           plink_bin = "./plink_linux_x86_64_20230116/plink") -> snpafterclump
  clumpset <- subset(dataset, SNP %in% snpafterclump$rsid)
  return(clumpset)
}

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

#This function is for ss of multiple metabolitess in a same format.
format=function(dataset, asformat="micro", snps=NULL, type="exposure") {
  tryCatch(
    {
      f=getFunction(paste("format",asformat,sep = "_"))
      return(f(dataset, snps = snps, type = type))
    },
    error=function(e) {
      message("No format function provided. Please set corresponding format function")
    }
  )
}

#Get all microbiome ss
micros <- list.files(path = "/home/wan/IBS/MiBiGen", pattern = "*.txt.gz")

#Set this group number of 10 microbiome from filename
filename <- current_filename()
print(filename)
filename %>% gsub(pattern = "3_1_Microbe_1_Generating_Clump_", replacement = "") %>% gsub(pattern = ".R", replacement = "") %>% as.numeric() -> n

#Set the group file name
if(n==21) {
  micros <- micros[(10*n-9):(10*n+1)]
} else {
  micros <- micros[(10*n-9):(10*n)]
}

#Do the following jobs to generate clumped data for each ss
for (i in micros) {
  
  print(i)
  
  #For re-running R script, clumped data do not run.
  a=list.files(path = "/home/wan/IBS/Clump_5e_6_new/microbiome", pattern = i)
  if(length(a)==1) {next
    print("Clump data exists, Skip")}
  
  #read the corresponding ss
  dat <- fread(input = paste("/home/wan/IBS/MiBiGen/",i,sep = ""),header = T)
  print(paste(i,"Done fread"))
  
  #format and clump
  dat %>% subset(P.weightedSumZ<5e-6) %>% format() %>% clump_new() -> clump
  print(paste(i,"has",nrow(clump),"SNP after clumping")) 
  
  #Save clump data locally
  write.csv(clump, file = paste("./Clump_5e_6_new/microbiome/",i,".csv",sep = ""))
}


