#This R script contains functions for freading and format

library(data.table)
library(dplyr)

#Note the directory is for Linux, not for /Volume
fread_adj=function(neuro) {
  
  if(neuro=="ADHD") {
    dat <- fread(input = "./Neuro/ADHD2022_iPSYCH_deCODE_PGC.meta.gz", header = T)
    dat %>% mutate(Beta=log(OR)) %>% mutate(Pheno="ADHD") -> dat
    }
 
  if(neuro=="ANERV"){
    dat <- fread(input = "./Neuro/pgcAN2.2019-07.vcf.tsv.gz", header = T)
    dat %>% mutate(Pheno="Anorexia nervosa") -> dat
  }
  
  if(neuro=="ANXIETY"){
    dat <- fread(input = "./Neuro/anxiety.meta.full.cc.tbl.gz", header = T)
    dat %>% mutate(Pheno="Anxiety") -> dat
  }
  
  if(neuro=="ASD"){
    dat <- fread(input = "./Neuro/iPSYCH-PGC_ASD_Nov2017.gz", header = T)
    dat %>% mutate(BETA=log(OR)) %>% mutate(Pheno="Autism spectrum disorder") -> dat
  }
  
  if(neuro=="BIP"){
    dat <- fread(input = "./Neuro/daner_bip_pgc3_nm_noukbiobank.gz", header = T, fill = T)
    dat %>% mutate(Beta=log(OR)) %>% mutate(Pheno="BIP") -> dat
  }
  
  if(neuro=="MDD"){
    dat <- fread(input = "./Neuro/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz", header = T)
    dat %>% mutate(Pheno="Major depreesion") %>% mutate(BETA=log(OR)) -> dat
  }
  
  if(neuro=="PTSD"){
    dat <- fread(input = "./Neuro/pts_eur_freeze2_overall.results.gz", header = T)
    dat %>% mutate(Beta=log(OR)) %>% mutate(Pheno="PTSD")-> dat
  }
  
  if(neuro=="SCZ"){
    dat <- fread(input = "./Neuro/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", header = T)
    dat %>% mutate(Pheno="SCZ") -> dat
  }
  
  if(neuro=="HOA"){
    dat <- fread(input = "./Neuro/hoarding2022.vcf.tsv.gz", header = T)
    dat %>% mutate(Pheno="Hoarding syndrome") -> dat
  }
  
  if(neuro=="SAPNEA"){
    dat <- fread(input = "./Neuro/chen_and_cade_et_al_2018_ahi_3_percent_european_americans.txt", header = T)
    dat %>% mutate(Pheno="Sleep apnea") -> dat
  }
  
  if(neuro=="TS"){
    dat <- fread(input = "./Neuro/TS_Oct2018.gz", header = T)
    dat %>% mutate(BETA=log(OR)) %>% mutate(Pheno="Tourette syndrome")-> dat
  }
  
  if(neuro=="OCD"){
    dat <- fread(input = "./Neuro/ocd_aug2017.gz", header = T)
    dat %>% mutate(Pheno="OCD") %>% mutate(BETA=log(OR))  -> dat
  }
  
  if(neuro=="PANIC"){
    dat <- fread(input = "./Neuro/pgc-panic2019.vcf.tsv.gz", header = T)
    dat %>% mutate(Pheno="Panic") -> dat
  }
  
  if(neuro=="SDURA"){
    dat  <- fread(input = "./Neuro/sleepdurationsumstats.txt", header = T)
    dat  %>% mutate(Pheno="Sleep duration") -> dat
  }

    if(neuro=="CHRON"){
    dat <- fread(input = "./Neuro/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz", header = T)
    dat %>% mutate(Pheno="Chronotype") -> dat
  }
  
  if(neuro=="NAP"){
    dat <- fread(input = "./Neuro/bolt_453K_geneticEurCluster.Nap_noBMI_VM.bgen.stats_HRC_1KG_LDSC_cleaned.txt", header = T)
    dat %>% mutate(Pheno="Napping") -> dat
  }
  
  if(neuro=="DSLEEP"){
    dat <- fread(input = "./Neuro/Saxena.fullUKBB.DaytimeSleepiness_adjBMI.sumstats.txt", header = T)
    dat %>% mutate(Pheno="Daytime Sleepness") -> dat
  }
  
  if(neuro=="INSOM"){
    dat <- fread(input = "./Neuro/Saxena_fullUKBB_Insomnia_summary_stats.txt", header = T)
    dat %>% mutate(Pheno="Insomnia") -> dat
  }
  
  if(neuro=="STROKE"){
    dat <- fread(input = "./Neuro/MEGASTROKE.1.AS.EUR.out", header = T)
    dat %>% mutate(Pheno="Stroke") -> dat
  }
  
  if(neuro=="LS"){
    dat <- fread(input = "./Neuro/LacunarStroke-GWAS-TraylorM-et-al-European-04122020.txt.gz", header = T)
    dat %>% mutate(Pheno="Lacunar stroke") ->dat
  }
  
  if(neuro=="IANE"){
    dat <- fread(input = "./Neuro/IA.GWAS.BakkerMK.2020.sumstats.Stage_1_excludingUKBB.txt.gz", header = T)
    dat %>% mutate(Pheno="Intracranial aneurysm") -> dat
  }
  
  if(neuro=="ALZ"){
    dat <- fread(input = "./Neuro/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz", header = T)
    dat %>% mutate(chr_pos=paste(chromosome,base_pair_location, sep = ":")) %>% mutate(Pheno="ALZ") -> dat
  }
  
  if(neuro=="PD") {
    dat <- fread(input = "./Neuro/nallsEtAl2019_excluding23andMe_allVariants.tab", header = T)
    dat %>% mutate(chr_pos=gsub("chr","",SNP)) %>% dplyr::select(-SNP) %>% mutate(Pheno="Parkinson Disease") -> dat
    as.data.table(do.call(rbind,strsplit(dat$chr_pos, split = ":")))-> dat_chr_pos
    dat %>% mutate(CHR=dat_chr_pos$V1) %>% mutate(POS=dat_chr_pos$V2) -> dat
  }
  return(dat)
}

#Reading IBS summary statistics
fread_IBS=function(ibs){
  if(ibs=="IBS_1") {
    print("Freading IBS_1 ...")
    dat <- fread(input = "BonfiglioF_29626450.txt.gz", header = T) #From UKB
    dat %>% mutate(Beta=log(OR)) %>% mutate(Beta_SE=abs(log(OR)/qnorm(P/2))) %>% mutate(Pheno="IBS_1")-> dat
  }
  
  if(ibs=="IBS_2"){
    print("Freading IBS_2 ...")
    dat <- fread(input = "METAANALYSIS1.TBL.pos", header = T) #From Bellygenes
    dat %>% mutate(Pheno="IBS_2") -> dat
  }
  
  if(ibs=="IBS_3"){
    print("Freading IBS_3 ...")
    dat <- fread(input = "finngen_R9_K11_IBS.gz", header = T) #Finnland
    dat %>% mutate(Pheno="IBS_3") -> dat
  }
  return(dat)
}

#format_functions for each of-interest disorder
format_ADHD=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset, snps = snps, type = type,
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "FRQ_A_38691",
                         beta_col = "Beta", 
                         se_col = "SE",
                         pval_col = "P", 
                         ncase_col = "Nca",ncontrol_col = "Nco",
                         phenotype_col = "Pheno",
                         chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_ANERV=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset, snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "ID",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "PVAL",
                         chr_col = "CHR", pos_col = "BP",
                         ncase_col = "NCAS", ncontrol_col = "NCON")
  return(dataset)
}
format_ANXIETY=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset, snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNPID",
                         effect_allele_col = "Allele1",
                         other_allele_col = "Allele2",
                         eaf_col = "Freq1",
                         beta_col = "Effect", 
                         se_col = "StdErr",
                         pval_col = "P.value",
                         chr_col = "CHR", pos_col = "BP",
                         samplesize_col = "TotalN")
  return(dataset)
}
format_ASD=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset, snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_BIP=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset, snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "FRQ_A_40463",
                         beta_col = "Beta", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP",
                         ncase_col = "Nca", ncontrol_col = "Nco")
  return(dataset)
}
format_MDD=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset,snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "FRQ_A_45396",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP",
                         ncase_col = "Nca", ncontrol_col = "Nco")
  return(dataset)
}
format_OCD=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset,snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_PANIC=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset,snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "ID",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "FCAS",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "PVAL",
                         ncase_col = "NCAS", ncontrol_col = "NCON",
                         chr_col = "#CHROM", pos_col = "POS")
  return(dataset)
}
format_PTSD=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset,snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "FRQ_A_23212",
                         beta_col = "Beta", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP",
                         ncase_col = "Nca", ncontrol_col = "Nco")
  return(dataset)
}
format_SCZ=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "ID",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       eaf_col = "FCAS",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "PVAL",
                       ncase_col = "NCAS", ncontrol_col = "NCON",
                       chr_col = "CHROM", pos_col = "POS")
  return(dataset)
}
format_SAPNEA=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       eaf_col = "A1FREQ",
                       beta_col = "BETA", 
                       se_col = "BETA_SE",
                       pval_col = "P",
                       chr_col = "CHR", pos_col = "BP",
                       samplesize_col = "N")
  return(dataset)
}
format_TS=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "P",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_HOA=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       snp_col = "ID",
                       phenotype_col = "Pheno",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       eaf_col = "FCAS",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "PVAL",
                       ncase_col = "NCAS", ncontrol_col = "NCON",
                       chr_col = "CHROM", pos_col = "POS")
  return(dataset)
}

format_SDURA=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       eaf_col = "A1FREQ",
                       beta_col = "BETA_SLEEPDURATION", 
                       se_col = "SE_SLEEPDURATION",
                       pval_col = "P_SLEEPDURATION",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_CHRON=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       eaf_col = "A1FREQ",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "P_BOLT_LMM",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_INSOM=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       eaf_col = "A1FREQ",
                       beta_col = "BETA_INSOMNIA", 
                       se_col = "SE_INSOMNIA",
                       pval_col = "P_INSOMNIA",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_DSLEEP=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "ALLELE1",
                       other_allele_col = "ALLELE0",
                       eaf_col = "A1FREQ",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "P",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_NAP=function(dataset, snps = NULL, type = "exposure") {
  dataset <- format_data(dataset,snps = snps, type = type,
                         phenotype_col = "Pheno",
                         snp_col = "SNP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "EAF",
                         beta_col = "BETA", 
                         se_col = "SE",
                         pval_col = "P",
                         chr_col = "CHR", pos_col = "BP",
                         samplesize_col = "N")
  return(dataset)
}

format_STROKE=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "MarkerName",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       eaf_col = "Freq1",
                       beta_col = "Effect", 
                       se_col = "StdErr",
                       pval_col = "P-value")
  return(dataset)
}
format_IANE=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       snp_col = "SNP",
                       effect_allele_col = "A_EFF",
                       other_allele_col = "A_NONEFF",
                       eaf_col = "Freq_EFF",
                       beta_col = "BETA", 
                       se_col = "SE",
                       pval_col = "P",
                       chr_col = "CHR", pos_col = "BP",
                       samplesize_col = "N",
                       phenotype_col = "Pheno")
  return(dataset)
}
format_LS=function(dataset, snps = NULL, type = "exposure") {
  dataset<-format_data(dataset, snps = snps, type = type,
                       snp_col = "SNP",
                       effect_allele_col = "EA",
                       other_allele_col = "OA",
                       beta_col = "Beta", 
                       se_col = "SE",
                       pval_col = "pval",
                       chr_col = "CHR", pos_col = "BP",
                       phenotype_col = "Pheno")
  return(dataset)
}

format_PD=function(dataset, snps = NULL, type = "exposure"){
  as.data.table(do.call(rbind,strsplit(dataset$chr_pos, split = ":")))-> dat_chr_pos
  dataset %>% mutate(CHR=dat_chr_pos$V1) %>% mutate(POS=dat_chr_pos$V2) ->dataset
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       eaf_col = "freq",
                       beta_col = "b", 
                       se_col = "se",
                       pval_col = "p",
                       ncase_col = "N_cases",ncontrol_col = "N_controls",
                       chr_col = "CHR", pos_col = "POS")
}

format_ALZ=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       snp_col = "SNP",
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       eaf_col = "effect_allele_frequency",
                       beta_col = "beta", 
                       se_col = "standard_error",
                       pval_col = "p_value", samplesize_col = "N",
                       phenotype_col = "Pheno",
                       chr_col = "chromosome", pos_col = "base_pair_location")
}
format_IBS_1=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       beta_col = "Beta",
                       se_col = "Beta_SE",
                       pval_col = "P",
                       eaf_col = "EAF",
                       effect_allele_col = "A1",
                       other_allele_col = "A2", samplesize_col = "N",
                       chr_col = "CHR",
                       pos_col = "POS")
  return(dataset)
}
format_IBS_2=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "SNP",
                       beta_col = "Effect",
                       se_col = "StdErr",
                       pval_col = "P-value",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       ncase_col = "N_CASE", ncontrol_col = "N_CONTROL",
                       samplesize_col = "N_TOTAL",
                       eaf_col = "Freq1",
                       chr_col = "CHR", pos_col = "BP")
  return(dataset)
}
format_IBS_3=function(dataset, snps = NULL, type = "exposure"){
  dataset<-format_data(dataset, snps = snps, type = type,
                       phenotype_col = "Pheno",
                       snp_col = "rsids",
                       beta_col = "beta",
                       se_col = "sebeta",
                       pval_col = "pval",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       eaf_col = "af_alt",
                       chr_col = "#chrom", pos_col = "pos")
  return(dataset)
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
} #From GCST90200422 Serotonin levels

#This function is for ss of multiple metabolites in a same format.
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

#Set format function for microbiome
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