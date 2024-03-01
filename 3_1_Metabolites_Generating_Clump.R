#This R script is to generate clumped data for metabolites

#Save into /home/wan/IBS/Clump_5e_6_new/metabolites

#Set main environment
setwd("/home/wan/IBS")

#Install R packages
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ieugwasr)

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

#This function is for ss of multiple metabolitess in a same format.
format=function(dataset, asformat="SERO", snps=NULL, type="exposure") {
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

#Set metabolites
metabolites <- paste("GCST90",199621:201020,sep = "")

#Do the following jobs to generate clumped data for each ss
for (i in metabolites) {
  
  print(i)
  
  #For re-running R script, clumped data do not run.
  a=list.files(path = "./Clump_5e_6_new/metabolites", pattern = i)
  if(length(a)==1) {next
    print("Clump data exists, Skip")}
  
  #read the corresponding ss
  dat <- fread(input = paste("./metabolites/",i,"_buildGRCh38.tsv.gz",sep = ""),header = T)
  dat %>% mutate(Pheno=i) -> dat
  print(paste(i,"Done fread"))
  
  #format and clump
  dat %>% subset(p_value<5e-6) %>% format() %>% clump_new() -> clump
  print(paste(i,"has",nrow(clump),"SNP after clumping")) 
  
  #Save clump data locally
  write.csv(clump, file = paste("./Clump_5e_6_new/metabolites/",i,".csv",sep = ""))
}

