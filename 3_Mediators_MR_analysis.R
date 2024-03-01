#This R script is to analyse twosample har.dat between metabolites and neuro/ibs.

#To incorporate the sample size from original studies to ss that do not included sample size.
#To calculate R2 and F value
#To adjust P value 

library(TwoSampleMR)
library(dplyr)
library(data.table)

#Set a directory where har.dat is
setwd(".")

#Set diseases (exposure)
forward.exp <- c("ADHD","BIP","MDD","HOA","INSOM")

#Set 3 IBS SS
forward.out <- c("IBS_3")

#Set metabolites
metabolites <- paste("GCST90",199621:201020,sep = "")

#Set microbes
microbiome <- list.files(path = "./MiBiGen", pattern = "*.txt.gz")
microbiome <- gsub(pattern = ".summary.txt.gz", replacement = "", microbiome)

#Add information to previous file OR Create a data frame to record results
mr_results <- tryCatch(
  {mr_results <- fread(input = "./Mediators_NeuroIBS_mr_not_presso.csv")
  mr_results <- data.frame(mr_results)
  },
  error=function(e) {
    mr_results <- data.frame(matrix(ncol = 17))
    names(mr_results) <- c("Exposure","Outcome","N_SNP","R2","F","Pleiotropy","Heterogeneity","b","se",
                           "mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode",
                           "Independent","P_ivw_unadjust","P_adjust")
    return(mr_results)
  }
)

#Give an initial row
k=1

#Import SS_info
SS_info <- read.csv(file = "./0_Database_list.csv")

#Set analyse function
mr_interpretation = function(i,j) {
  
  #Set metabolites
  metabolites <- paste("GCST90",199621:201020,sep = "")
  
  #Set microbes
  microbiome <- list.files(path = "./MiBiGen", pattern = "*.txt.gz")
  microbiome <- gsub(pattern = ".summary.txt.gz", replacement = "", microbiome)
  
  #Test if file exists
  exist <- list.files(pattern = paste("^",i,"_",j,".csv",sep = ""), path = c("./twosampleMR_metabolites_harmonise","./twosampleMR_microbiome_harmonise"), full.names = T)
  if(length(exist)==0) {
    mr_results[k,c("Exposure","Outcome","N_SNP")] <- c(i,j,"Har.dat not found")
    print("Harmonised data do not exist")
    return(mr_results) }
  
  #Import har.dat IBS as outcome
  print("Read harmonised data")
  har.dat <- fread(input = exist)
  har.dat <- data.frame(har.dat)
  mr_results[k,c("Exposure","Outcome","N_SNP")] <- c(i,j,nrow(subset(har.dat,mr_keep==T)))
  
  #Some har.dat have 0 SNP therefore do not need mr interpretation
  if(mr_results[k,"N_SNP"]==0) {
    print("Harmonised data has 0 SNP")
    return(mr_results)
  }
  
  #Provide sample size if not existing
  if(i %in% metabolites) {exposure_info_row=which(SS_info$Abbreviations=="Metabolites")
  } else {
    if (i %in% microbiome) {exposure_info_row=which(SS_info$Abbreviations=="Microbiome")
    } else {
      exposure_info_row=which(SS_info$Abbreviations==i) }
  }
  
  if(!("samplesize.exposure" %in% names(har.dat))) {
    print(SS_info$Sample.size[exposure_info_row])
    har.dat %>% mutate(samplesize.exposure=as.integer(SS_info$Sample.size[exposure_info_row])) -> har.dat
  }
  
  #Recalculate R2 if eaf provided
  if(NA %in% har.dat$eaf.exposure) {
    print("No EAF, R2 from get_r_from_bsen")
    har.dat %>% mutate(R2=get_r_from_bsen(b=abs(beta.exposure),se=se.exposure, n=samplesize.exposure)^2) -> har.dat
  } else { #if eaf provided
    print("Have EAF, R2 from eaf and beta")
    har.dat %>% dplyr::mutate(R2=2*eaf.exposure*(1-eaf.exposure)*beta.exposure^2) -> har.dat
  }
  
  #Recalculate F value
  R2=sum(subset(har.dat,mr_keep==TRUE)$R2)
  nsnp=nrow(subset(har.dat,mr_keep==TRUE))
  samplesize=max(subset(har.dat,mr_keep==TRUE)$samplesize.exposure)
  F=(samplesize-nsnp-1)/nsnp*R2/(1-R2)
  
  #Save new R2 and F into mr_results
  print(paste(R2,F))
  mr_results[k,c("R2","F")] <- c(R2, F)
  
  #Save pleiotropy and heterogeneity
  #har.dat has only 1 SNP could not be tested 
  if(mr_results[k,"N_SNP"]==1) {
    mr_results[k,c("Pleiotropy","Heterogeneity")] <- NA
  } else {
    pleio <- mr_pleiotropy_test(har.dat)
    ivw_hetero <- mr_heterogeneity(har.dat) %>% dplyr::filter(method=="Inverse variance weighted")
    mr_results[k,c("Pleiotropy","Heterogeneity")] <- c(pleio$pval,ivw_hetero$Q_pval)
    print(paste(pleio$pval,ivw_hetero$Q_pval))
  }
  
  #Save beta and se
  b_se <- mr(dat = har.dat, method_list = "mr_ivw")[1,c("b","se")]
  mr_results[k,c("b","se")] <- b_se[1,c("b","se")]
  
  #Save OR_CI_P by different methods
  mr <- generate_odds_ratios(mr_res = mr(har.dat, method_list = c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")))
  mr %>% mutate(OR_CI_P=paste("OR",or,"CI",or_lci95,or_uci95,"P",pval,sep = "_")) -> mr
  
  #Some har.dat do not have enough number of SNP for all the methods, need supplement NA
  p=nrow(mr)
  for (method in c("Inverse variance weighted","Inverse variance weighted (multiplicative random effects)","Weighted median","MR Egger","Weighted mode")) {
    if(method %in% mr$method) {next}
    p=p+1
    mr[p,1:5] <- c(mr[1,1:4],method)
  }
  mr_results[k,c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")] <- mr$OR_CI_P
  mr_results[k,"P_ivw_unadjust"] <- mr$pval[1]
  
  #Identify overlap
  if(j %in% metabolites) {outcome_info_row=which(SS_info$Abbreviations=="Metabolites")} else {
    if(j %in% microbiome) {outcome_info_row=which(SS_info$Abbreviations=="Microbiome")} else {
      outcome_info_row=which(SS_info$Abbreviations==j) }
  }
  
  if(FALSE %in% grepl("Large",SS_info$UKB[c(exposure_info_row,outcome_info_row)])) {
    mr_results[k,"Independent"] <- "OK"
  } else {mr_results[k,"Independent"] <- "Largely_ovelapped"}
  print(mr_results[k,])
  
  #Save mr_results on time
  write.csv(mr_results, file="./Mediators_NeuroIBS_mr_not_presso.csv", fileEncoding = "UTF-8",row.names = FALSE)
  print(paste(i,j,"interpretation done"))
  return(mr_results)
}

#Do the following for every pair of exposure and outcome
for (j in c(metabolites, microbiome)) {
  
  #metabolites as outcome
  for (i in forward.exp) {
  
    print(paste(i,j))
    
    #Find if i_j pair exists
    i_row=which(mr_results$Exposure==i & mr_results$Outcome==j)
    if(length(i_row)==0) { #Do not find this pair row, run the below
      runscript=1} else {
        m=1
        runscript=0
        while (m<ncol(mr_results)) {
          if(is.na(mr_results[i_row,m])) {
            runscript=1 #Find NA in previous row, run the below
            break
          } else {
            m=m+1
          }
        }
      }
    
    #Run or not
    if(runscript==0) {
      print("Skipping") } else {
      mr_results <- mr_interpretation(i,j)
      print(paste(i,j,"MR analysis Done"))}
    
    k=k+1 }
  
  #metabolites as exposure
  print(paste(j,forward.out))
  
  #Find if j_i pair exists
  j_row=which(mr_results$Exposure==j & mr_results$Outcome==forward.out)
  if(length(j_row)==0) { #Do not find this pair row, run the below
    runscript=1} else {
      m=1
      runscript=0
      while (m<ncol(mr_results)) {
        if(is.na(mr_results[j_row,m])) {
          runscript=1 #Find NA in previous row, run the below
          break
        } else {
          m=m+1
        }
      }
    }
  
  #Run or not
  if(runscript==0) {
    print("Skipping") } else {
      mr_results <- mr_interpretation(j,forward.out)
      print(paste(j,forward.out,"MR analysis Done"))
    }
  
  k=k+1
  
  print(paste(j,"Done!"))
}
