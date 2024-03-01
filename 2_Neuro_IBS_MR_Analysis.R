##This R script is to analyze Neuro_IBS bidirectional TwoSampleMR harmonise data

#To incorporate the sample size from orginal studies to ss that do not included sample size.
#To calculate R2 and F value
#To run MRPresso for data with significant pleiotropy
#To adjust P value 

library(TwoSampleMR)
library(dplyr)
library(data.table)

#Set a directory where har.dat is
setwd("/home/wan/IBS/Neuro_IBS/harmonise")

#Set 23 neuro diseases
forward.exp <- c("ALZ","PD","ADHD","ANERV","ANXIETY","ASD","BIP","MDD","PANIC","OCD","HOA",
                 "TS","SAPNEA","SCZ","PTSD",
                 "CHRON","INSOM","DSLEEP","SDURA","NAP",
                 "STROKE","IANE","LS")

#Set 3 IBS SS
forward.out <- c("IBS_3")

#Add information to previous file OR Create a data frame to record results
mr_results <- tryCatch(
  {mr_results <- fread(input = "/home/wan/IBS/Neuro_IBS_mr_not_presso.csv")
   mr_results <- data.frame(mr_results)
  },
  error=function(e) {
    mr_results <- data.frame(matrix(ncol = 17))
    names(mr_results) <- c("Exposure","Outcome","N_SNP","R2","F","Pleiotropy","Heterogeneity","b_total","se_total",
                           "mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode",
                           "Independent","P_ivw_unadjust","P_adjust")
    return(mr_results)
  }
)

#Give an initial row
k=1

#Import SS_info
SS_info <- read.csv(file = "/home/wan/IBS/0_Database_list.csv")

#Do the following for every pair of exposure and outcome
for (neuro in forward.exp) {
  
  for (j in forward.out) {
    
    print(paste(neuro,j))
    #Find if neuro_j pair exists
    neuro_row=which(mr_results$Exposure==neuro & mr_results$Outcome==j)
    if(length(neuro_row)==0) { 
      
      #Do not find this pair row, run the below
      runscript=1
      
      } else { 
        
      #Find the start from existing results
      m=1
      runscript=0
      
      while (m<ncol(mr_results)) {
        
        if(is.na(mr_results[neuro_row,m])) {
          
          runscript=1 #Find NA in previous row, run the below
          break
          
        } else {
          m=m+1
        }
      }
    }
    
    #Run or not
    if(runscript==0) {
      k=k+1
      print("Skipping")
      next}

    #Import har.dat IBS as outcome
    print("Read harmonised data")
    har.dat <- fread(input = paste(paste(neuro,j,sep = "_"),".csv",sep = ""))
    har.dat <- data.frame(har.dat)
    mr_results[k,c("Exposure","Outcome","N_SNP")] <- c(neuro,j,nrow(subset(har.dat,mr_keep==T)))
    
    #Provide sample size if not existing
    exposure_info_row=which(SS_info$Abbreviations==neuro)
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
    pleio <- mr_pleiotropy_test(har.dat)
    ivw_hetero <- mr_heterogeneity(har.dat) %>% dplyr::filter(method=="Inverse variance weighted")
    mr_results[k,c("Pleiotropy","Heterogeneity")] <- c(pleio$pval,ivw_hetero$Q_pval)
    print(paste(pleio$pval,ivw_hetero$Q_pval))
    
    #Save beta and se
    b_se <- mr(dat = har.dat, method_list = "mr_ivw")[1,c("b","se")]
    mr_results[k,c("b_total","se_total")] <- b_se[1,c("b","se")]
    
    #Save OR_CI_P by different methods
    mr <- generate_odds_ratios(mr_res = mr(har.dat, method_list = c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")))
    mr %>% mutate(OR_CI_P=paste("OR",or,"CI",or_lci95,or_uci95,"P",pval,sep = "_")) -> mr
    mr_results[k,c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")] <- mr$OR_CI_P
    mr_results[k,"P_ivw_unadjust"] <- mr$pval[1]
    
    #Identify overlap
    outcome_info_row=which(SS_info$Abbreviations==j)
    if(FALSE %in% grepl("Large",SS_info$UKB[c(exposure_info_row,outcome_info_row)])) {
      mr_results[k,"Independent"] <- "OK"
    } else {mr_results[k,"Independent"] <- "Largely_ovelapped"}
    print(mr_results[k,])
    
    k=k+1
    
    #Save mr_results on time
    write.csv(mr_results, file="/home/wan/IBS/Neuro_IBS_mr_not_presso.csv", fileEncoding = "UTF-8",row.names = FALSE) 
  }
}

#Do the following for every pair of exposure and outcome
for (ibs in forward.out) {
  for (neuro in forward.exp) {
    
    print(paste(ibs,neuro))
    #Find if ibs_neuro pair exists
    ibs_row=which(mr_results$Exposure==ibs & mr_results$Outcome==neuro)
    if(length(ibs_row)==0) { #Do not find this pair row, run the below
      runscript=1} else {
        m=1
        runscript=0
        while (m<ncol(mr_results)) {
          if(is.na(mr_results[ibs_row,m])) {
            runscript=1 #Find NA in previous row, run the below
            break
          } else {
            m=m+1
          }
        }
      }
    
    #Run or not
    if(runscript==0) {
      k=k+1
      print("Skipping")
      next}
    
    #Import har.dat IBS as outcome
    print("Read harmonised data")
    har.dat <- fread(input = paste(paste(ibs,neuro,sep = "_"),".csv",sep = ""))
    har.dat <- data.frame(har.dat)
    mr_results[k,c("Exposure","Outcome","N_SNP")] <- c(ibs,neuro,nrow(subset(har.dat,mr_keep==T)))
    
    #Provide sample size if not existing
    ##Not adding sample size before R2 is to enable calculation of R2 by get_r_from_bsen function
    exposure_info_row=which(SS_info$Abbreviations==ibs)
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
    pleio <- mr_pleiotropy_test(har.dat)
    ivw_hetero <- mr_heterogeneity(har.dat) %>% dplyr::filter(method=="Inverse variance weighted")
    mr_results[k,c("Pleiotropy","Heterogeneity")] <- c(pleio$pval,ivw_hetero$Q_pval)
    print(paste(pleio$pval,ivw_hetero$Q_pval))
    
    #Save beta and se
    b_se <- mr(dat = har.dat, method_list = "mr_ivw")[1,c("b","se")]
    mr_results[k,c("b_total","se_total")] <- b_se[1,c("b","se")]
    
    #Save OR_CI_P by different methods
    mr <- generate_odds_ratios(mr_res = mr(har.dat, method_list = c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")))
    mr %>% mutate(OR_CI_P=paste("OR",or,"CI",or_lci95,or_uci95,"P",pval,sep = "_")) -> mr
    mr_results[k,c("mr_ivw","mr_ivw_mre","mr_weighted_median","mr_egger_regression","mr_weighted_mode")] <- mr$OR_CI_P
    mr_results[k,"P_ivw_unadjust"] <- mr$pval[1]
    
    #Identify overlap
    outcome_info_row=which(SS_info$Abbreviations==neuro)
    if(FALSE %in% grepl("Large",SS_info$UKB[c(exposure_info_row,outcome_info_row)])) {
      mr_results[k,"Independent"] <- "OK"
    } else {mr_results[k,"Independent"] <- "Largely_ovelapped"}
    print(mr_results[k,])
    
    k=k+1
    
    #Save mr_results on time
    write.csv(mr_results, file="/home/wan/IBS/Neuro_IBS_mr_not_presso.csv", fileEncoding = "UTF-8",row.names = FALSE) 
  }
}

mr_results$P_adjust <- p.adjust(mr_results$P_ivw_unadjust, method = "fdr") 
write.csv(mr_results, file="/home/wan/IBS/Neuro_IBS_mr_not_presso.csv", fileEncoding = "UTF-8",row.names = FALSE) #Table S1

mr_results %>% dplyr::filter(P_ivw_unadjust < 0.05)  #Identify significant causal relationships
write.csv(subset(mr_results, P_ivw_unadjust < 0.05), file="/home/wan/IBS/Significant_Neuro_IBS.csv", fileEncoding = "UTF-8",row.names = FALSE) #Table 1
