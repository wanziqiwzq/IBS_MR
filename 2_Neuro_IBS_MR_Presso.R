#This R script is to do sensitivity analysis for har.dat with significant pleiotropy.

#Set the environment where mr_results save
setwd(".")

#Loading packages
library(data.table)
library(dplyr)
library(parallel)

#Import mr results of bidirectional two-sample MR
N_I_unpresso <- fread(input = "Neuro_IBS_mr_not_presso.csv", header = T) %>% data.frame()
print(head(N_I_unpresso))

#Add MRPRESSO results in a new column
N_I_unpresso %>% mutate(MRPRESSO_result=NA) %>% mutate(Pair=paste(Exposure,Outcome,sep = "_")) -> N_I_unpresso
print(head(N_I_unpresso))

#Select har.dat having significant pleiotropy
N_I_unpresso %>% dplyr::filter(Pleiotropy<0.05) -> N_I_need_presso
print(N_I_need_presso$Pair)

#Multithread
cl <- makeCluster(20)

Pair_Presso <- function(pair) {
  
  #Load packages
  library(data.table)
  library(dplyr)
  library(MRPRESSO)
  library(TwoSampleMR)
  
  #Import mr results of bidirectional two-sample MR
  N_I_unpresso <- fread(input = "Neuro_IBS_mr_not_presso.csv", header = T) %>% data.frame()
  
  #Add MRPRESSO results in a new column
  N_I_unpresso %>% mutate(MRPRESSO_result=NA) %>% mutate(Pair=paste(Exposure,Outcome,sep = "_")) -> N_I_unpresso
  print(head(N_I_unpresso))
  
  #Read har.dat by pair name
  har.dat <- fread(input = paste("./Neuro_IBS/harmonise/",pair,".csv",sep = "")) %>% data.frame()
  print("Fread har.dat Done")
  
  #Run MRPRESSO
  print("Running MR-presso")
  presso <- run_mr_presso(har.dat, NbDistribution = 3000)
  
  #Test MRPRESSO P value
  presso_global <- presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
  print(presso_global)
  
  #Find the right row for the pair
  pair_row=which(N_I_unpresso$Pair==pair)
  
  #If Significant pleiotropy with mr-egger and mr-presso
  if(presso_global<0.05) {
    
    if(is.na(presso[[1]]$`Main MR results`[2,"P-value"])) {
      
      #Senerio 1: Significant pleiotropy with mr-egger and mr-presso but could not be solved by testing outliners
      print("hardat had significant pleiotropy that could not reduced by mr-presso")
      presso_results <- paste(presso_global,"significant pleiotropy, no outliners detected") 
      
    } else{
      
      #Senerio 2: Significant pleiotropy with mr-egger and mr-presso and outliners found
      outliners <- subset(presso[[1]]$`MR-PRESSO results`$`Outlier Test`,Pvalue<0.1) #Identify outliners
      har.dat <- har.dat[-as.integer(rownames(outliners)),] #Delete outliners
      print("New har.dat")
      
      #Set presso results
      presso_results <-  paste(presso_global,"significant pleiotropy, outliners deleted")
      
      #Recalculate F value
      R2=sum(subset(har.dat,mr_keep==TRUE)$R2)
      nsnp=nrow(subset(har.dat,mr_keep==TRUE))
      samplesize=max(subset(har.dat,mr_keep==TRUE)$samplesize.exposure)
      F=(samplesize-nsnp-1)/nsnp*R2/(1-R2)
      print(paste(R2, F))
      
      #Re-test Pleiotropy Heterogeneity
      pleio <- mr_pleiotropy_test(har.dat)
      ivw_hetero <- mr_heterogeneity(har.dat) %>% dplyr::filter(method=="Inverse variance weighted")
      print(paste(pleio$pval, ivw_hetero$Q_pval))
      
      #Save beta and se
      b_se <- mr(dat = har.dat, method_list = "mr_ivw")[1,c("b","se")]
      
      #Save OR_CI_P by different methods
      mr <- generate_odds_ratios(mr_res = mr(har.dat, method_list = c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")))
      mr %>% mutate(OR_CI_P=paste("OR",or,"CI",or_lci95,or_uci95,"P",pval,sep = "_")) -> mr
      
      #Write new MR results
      N_I_unpresso[pair_row,c("N_SNP","R2","F","Pleiotropy","Heterogeneity","b_total","se_total",
                              "mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode",
                              "P_ivw_unadjust")] <- c(nsnp, R2, F, pleio$pval, ivw_hetero$Q_pval,b_se[1,c("b","se")],mr$OR_CI_P, mr$pval[1])
      
      #Note that this parallel function could not synchronize to write the same file
      #Parallel generate a new environment without previous files
  
      #Save the new har.dat
      write.csv(har.dat, file = paste("./Neuro_IBS/harmonise/",pair,"_PRESSO.csv",sep = ""), row.names = FALSE)
      print("New har.dat Saved")
        
    }
  } else { 
    
    ##Senerio 3: Not having significant pleiotropy with MRPRESSO
    print("hardat did not have significant pleiotropy with global test pval")
    presso_results <-  paste(presso_global,"No significant pleiotropy by MRPRESSO")
    
  }
  
  #Add PRESSO into the column
  N_I_unpresso[pair_row,"MRPRESSO_result"] <- presso_results
  
  #Save
  write.csv(N_I_unpresso[pair_row,], file = paste("./Neuro_IBS/PRESSO_summary_",pair,".csv",sep = ""), row.names = FALSE)
  
  #Output
  return(c(pair, presso_results))
 
}

parLapply(cl, N_I_need_presso$Pair, Pair_Presso)

stopCluster(cl)

#Read local files
for (pair in N_I_need_presso$Pair) {
  
  #Find the right row for the pair
  pair_row=which(N_I_unpresso$Pair==pair)
  
  #Read after presso results
  afterpresso <- fread(input = paste("./Neuro_IBS/PRESSO_summary_",pair,".csv",sep = "")) %>% data.frame()
  
  #Import result into original results
  N_I_unpresso[pair_row,] <- afterpresso[1,]
  
  #Save
  write.csv(N_I_unpresso, file = "Neuro_IBS_MR_after_presso.csv", row.names = FALSE)
  print("Combine Done")
  
}
