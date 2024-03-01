#This R script is to give subsequent data modification for TwoSample mr results between metabolites and neuro/ibs

#Including

#adjusted p value by multiple testing.

#Select mediators with unadjusted significance (neuro-mediators-ibs).

#Load packages
library(data.table)
library(dplyr)
library(TwoSampleMR)

#Set where the csv is
setwd("/home/wan/IBS")

#Import Mediators_NeuroIBS_MR_after_presso.csv
M_NorI_MR_after_presso <- fread(input = "Mediators_NeuroIBS_MR_after_presso.csv") %>% data.frame()

#Set mediators group
metabolites <- paste("GCST90",199621:201020,sep = "")

#Set microbes
microbiome <- list.files(path = "/home/wan/IBS/MiBiGen", pattern = "*.txt.gz")
microbiome <- gsub(pattern = ".summary.txt.gz", replacement = "", microbiome)

#Set diseases (exposure)
forward.exp <- c("ADHD","BIP","MDD","HOA","INSOM")

#Set 3 IBS SS
forward.out <- c("IBS_3")

#Task 0: only 1 SNP, supply b and se
mr_one_snp <- function(exposure, outcome, path) {
  
  #find har.dat
  exist <- list.files(path = path, pattern = paste("^",exposure,"_",outcome,".csv",sep = ""), full.names = T)
  
  #import har.dat
  har.dat <- fread(input = exist) %>% data.frame()
  
  #Ensure only 1 SNP and has 1 SNP for use
  if(nrow(subset(har.dat,mr_keep==T))!=1) {
    result <- "0 or at least 2 SNP, not suitable for this method"
    return(result)
  }
  
  result <- generate_odds_ratios(mr(dat = har.dat, method_list = "mr_wald_ratio"))
  result %>% mutate(OR_CI_P=paste("OR",or,"CI",or_lci95,or_uci95,"P",pval,sep = "_")) -> result
  return(result)
}

for (i in 1:nrow(M_NorI_MR_after_presso)) {
  
  if(M_NorI_MR_after_presso[i,"N_SNP"]==1) {
    
    #Run MR analysis Wald ratio
    har.path <- c("./twosampleMR_metabolites_harmonise","./twosampleMR_microbiome_harmonise")
    exposure <- M_NorI_MR_after_presso[i,"Exposure"]
    outcome <- M_NorI_MR_after_presso[i,"Outcome"]
    har.run <- mr_one_snp(exposure = exposure, outcome = outcome, path = har.path)

    #Save results into the summary result table
    M_NorI_MR_after_presso[i,c("b","se","mr_ivw","P_ivw_unadjust")] <- har.run[1,c("b","se","OR_CI_P","pval")]
  }
}

write.csv(M_NorI_MR_after_presso, file = "Mediators_NeuroIBS_MR_after_presso.csv", row.names = F)

#Task 1: Adjusted p value by multiple testing
print("Task 1: Adjusted p value by multiple testing")

for (i in list(forward.exp,"IBS_3")) {
  
  #chang list to character
  print(i)
  
  range_exp_out <- c(metabolites,microbiome,i)
  
  #adjusted p value by multiple testing
  M_NorI_MR_after_presso %>% dplyr::filter((Exposure %in% range_exp_out) & (Outcome %in% range_exp_out)) %>%
    mutate(P_adjust=p.adjust(P_ivw_unadjust, method = "fdr")) -> bTS_MR
  
  if("ADHD" %in% i) { i <- "neuro"}
  
  assign(paste("bTS_MR","mediators",i,sep = "_"),bTS_MR)
}

#Task 2: select significance mediators

#Rule: mediating mediators should at least have significant relationship with both neuro and ibs

#Simplify: according to the bTS_MR_Mediators_neuro, set neuro as exposure, test ibs dat

#set exposure_outcome with significant p-adjusted
#bTS_MR_mediators_neuro %>% dplyr::filter(P_adjust<0.05) -> dat

#Scenario: No significant pair!
print("select significance mediators")
bTS_MR_mediators_neuro %>% dplyr::filter(P_ivw_unadjust<0.05) -> dat

#Set the data.frame to record significant pair
Neuro_Mediators_IBS <- data.frame(matrix(ncol = 37))
names(Neuro_Mediators_IBS) <- c("Exposure_Neuro","Mediator",paste("N_M",c("N_SNP","R2","F","Pleiotropy","Heterogeneity","b","se",
                                                                          "mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode",
                                                                          "Independent","P_ivw_unadjust","P_adjust","MRPRESSO_result","Pair"), sep = "_"),
                                  "Outcome_IBS",paste("M_I",c("N_SNP","R2","F","Pleiotropy","Heterogeneity","b","se",
                                                                "mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode",
                                                                "Independent","P_ivw_unadjust","P_adjust","MRPRESSO_result","Pair"), sep = "_"))
nmi_row=1

#Test the other side
for(j in 1:nrow(dat)) {
  
  #Neuro-Mediators-IBS
  if(dat[j,"Outcome"] %in% c(metabolites,microbiome)) {
    
    #Test IBS summary statistics
    for (ibs in forward.out) {
      
      #Get ibs dat
      ibs_dat <- get(paste("bTS_MR","mediators",ibs,sep = "_"))
      
      #find the corresponding row
      ibs_dat %>% dplyr::filter(Exposure==dat[j,"Outcome"]) -> test_info
      
      #Test the p-value
      if(test_info$P_ivw_unadjust<0.05) {
        
        #Write into results data frame
        Neuro_Mediators_IBS[nmi_row,] <- c(dat[j,],test_info[1,2:19])
        
        #output
        print(Neuro_Mediators_IBS[nmi_row,])
        
        #save it locally
        write.csv(Neuro_Mediators_IBS, file = "Candidate_neuro_mediators_ibs.csv")
        
        nmi_row=nmi_row+1
        
      } #End of test 
      
      print("Test Done")
      
    } #End of test
    
  } #End of scenerio
  
} #End of selection

#Task 3: Combine total effects into results
neuro_ibs <- fread(input = "/home/wan/IBS/Significant_Neuro_IBS.csv") %>% data.frame()

Neuro_Mediators_IBS %>% left_join(neuro_ibs, by=c("Exposure_Neuro"="Exposure","Outcome_IBS"="Outcome"), keep=T) -> Neuro_Mediators_IBS

#save it locally
write.csv(Neuro_Mediators_IBS, file = "Candidate_neuro_mediators_ibs.csv")

