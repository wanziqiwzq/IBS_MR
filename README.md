# IBS_MR
Customized codes for IBS MR

We posted our customized R codes for discovery for an unpublished paper "Metabolome-wide Mendelian randomization identifies ursodeoxycholate as a negative mediator of the effect of insomnia on irritable bowel disease" for replication.

Part 1: Generate independent instrument variables
Run script 1_Generating_clump_dat.R to generate independent instrument variables that will be used later.

Part 2: Bidirectional two-sample MR between IBS and of-interest disorders
Run scripts started with 2_ in a row to get MR results and their sensitivity analyses.
Read 0_Database_list.csv to check disorders of interest.

Part 3: Mediation analysis for two-step MR 
Run scripts started with 3_ in a row to get MR results and their sensitivity analyses.
Read 0_Database_list.csv to check metabolome and microbiome as the mediator pool.

Part 4: Calculate indirect effects via mediators
Script 4_1 selects candidate mediators that receive causal relations of exposure which have causal relations on outcomes.
Script 4_2 calculates exposure-adjusted effects of mediators on the causal relations of exposures on outcomes (by three methods).
Script 4_3 calculates not-adjusted effects of mediators on the causal relations of exposures on outcomes (by asymmetrical-delta method).
