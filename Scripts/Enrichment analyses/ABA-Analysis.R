################################################################################
#Allen Brain Atlas - Analysis
#Zillah Daysh
#2025
################################################################################
#R version 4.1.3 

#load packages
library(ABAEnrichment)

################################################################################
#IBD/SLE PPIN gene set analysis
################################################################################

#KAT8/INO80E gene set from Multimorbid3D
KAT8_INO80E_PPIN = c("KAT8","MSL1","MSL2","MSL3","KANSL1","KANSL2","KANSL3",
                     "RBBP5","PHF20","WDR5","RUVBL2","KAT5","H4C6","MCRS1",
                     "EP400","MORF4L1","H3C12","TRRAP","YEATS4","UCHL5",
                     "RUVBL1","ACTR5","INO80B","INO80C","TFPT","NFRKB",
                     "INO80","INO80E","PSMD3","ADRM1")

IBD_SLE_input = data.frame(KAT8_INO80E_PPIN, is_candidate = 1)

res_IBD_SLE_devel = aba_enrich(IBD_SLE_input, 
                               dataset='5_stages', 
                               cutoff_quantiles = c(0.5,0.6,0.7,0.8,0.9),
                               n_randsets = 10000)

sig_res_IBD_SLE_1 = res_IBD_SLE_devel_1[c(1:27),c(1,3,4)]

IBD_SLE_results = res_IBD_SLE_devel$results
IBD_SLE_results_genes = res_IBD_SLE_devel$genes
IBD_SLE_thresholds = res_IBD_SLE_devel$cutoffs

write.csv(IBD_SLE_results, "results_IBD_SLE_genes_devel.csv")

expression_chrom_remod_genes_devel = get_expression(structure_ids = IBD_SLE_results$structure_id,
                                                    gene_ids = IBD_SLE_results_genes$KAT8_INO80E_PPIN,
                                                    dataset = '5_stages')

save(expression_chrom_remod_genes_devel, file = "expression_chrom_remod_genes_devel.RData")

################################################################################
#T2DM PPIN gene set analysis
################################################################################

KAT8_PPIN = c("KAT8","MSL1","MSL2","MSL3","KANSL1","KANSL2","KANSL3",
              "RBBP5","PHF20","WDR5","RUVBL2","KAT5","H4C6","MCRS1",
              "EP400","MORF4L1","H3C12","H3C13","CBX3", "EHMT2", "CHD3",
              "TRRAP","YEATS4","INO80B","SETD1A","INO80","INO80E", "RBBP4")

BIN1_PPIN = c("BIN1", "MYC", "FBXW7", "SYNJ1", "ITSN1", "CDC42", "TNK2",
              "DNM2", "HCLS1", "CTTN", "SRC", "BCAR1", "RIN3", "RAB5B",
              "PIK3C3", "BECN1", "PRKN", "SNCA", "SHGL2", "TRRAP", "BRCA1",
              "HIF1A", "EP300", "CBL", "RAPGEF1", "CD2AP", "SH3GL2", "ARHGAP1",
              "UBE2L3")

ACE_PPIN = c("ACE", "KNG1", "KLKB1", "F12", "SERPING1", "MASP2")

APOC2_PPIN = c("APOC2", "APOE", "LPL", "LRP1", "HSP90B1", "CALR", "TAPBP")

PLEKHM1_PPIN = c("PLEKHM1", "MAP1LC3C", "GABARAP", "CALCOCO2", "RB1CC1",
                 "ATG16L1", "ATG7", "ATG3", "ATG4B", "OPTN", "BNIP3",
                 "BNIP3L", "RHEB", "GABARAPL2","TBK1", "BCL2L1", "BCL2",
                 "BECN1", "NBR1", "ULK2", "RPTOR", "IKBKE")


KAT8_T2DM_input = data.frame(KAT8_PPIN, is_candidate = 1)
BIN1_T2DM_input = data.frame(BIN1_PPIN , is_candidate = 1)
ACE_T2DM_input = data.frame(ACE_PPIN , is_candidate = 1)
APOC2_T2DM_input = data.frame(APOC2_PPIN , is_candidate = 1)
PLEKHM1_T2DM_input = data.frame(PLEKHM1_PPIN , is_candidate = 1)


res_KAT8_T2DM_devel = aba_enrich(KAT8_T2DM_input, dataset='5_stages',
                                 cutoff_quantiles = c(0.5,
                                                      0.6,
                                                      0.7,
                                                      0.8,
                                                      0.9),
                                 n_randsets = 10000)

res_BIN1_T2DM_devel = aba_enrich(BIN1_T2DM_input, dataset='5_stages',
                                 cutoff_quantiles = c(0.5,
                                                      0.7,
                                                      0.9),
                                 n_randsets = 10000)

res_ACE_T2DM_devel = aba_enrich(ACE_T2DM_input, dataset='5_stages',
                                cutoff_quantiles = c(0.5,
                                                     0.6,
                                                     0.7,
                                                     0.8,
                                                     0.9), 
                                n_randsets = 10000)

res_APOC2_T2DM_devel = aba_enrich(APOC2_T2DM_input, dataset='5_stages',
                                  cutoff_quantiles = c(0.5,
                                                       0.7,
                                                       0.9), 
                                  n_randsets = 10000)

res_PLEKHM1_T2DM_devel = aba_enrich(PLEKHM1_T2DM_input, dataset='5_stages',
                                    cutoff_quantiles = c(0.5,
                                                         0.6,
                                                         0.7,
                                                         0.8,
                                                         0.9),
                                    n_randsets = 10000)

#enrichment results

results_KAT8_T2DM = res_KAT8_T2DM_devel$results
results_BIN1_T2DM = res_BIN1_T2DM_devel$results
results_ACE_T2DM = res_ACE_T2DM_devel$results
results_APOC2_T2DM = res_APOC2_T2DM_devel$results
results_PLEKHM1_T2DM = res_PLEKHM1_T2DM_devel$results

write.csv(results_KAT8_T2DM, "ABA_enrichment_results_KAT8_T2DM.csv")
write.csv(results_ACE_T2DM, "ABA_enrichment_results_ACE_T2DM.csv")
write.csv(results_APOC2_T2DM, "ABA_enrichment_results_APOC2_T2DM.csv")
write.csv(results_BIN1_T2DM_11, "ABA_enrichment_results_BIN1_T2DM.csv")
write.csv(results_PLEKHM1_T2DM, "ABA_enrichment_results_PLEKHM1_T2DM.csv")


#expression data for enriched gene sets

expression_KAT8_T2DM = get_expression(structure_ids = results_KAT8_T2DM_10000$structure_id,
                                      gene_ids = KAT8_PPIN,
                                      dataset = '5_stages')
save(expression_KAT8_T2DM, file = "expression_KAT8_T2DM.RData")


expression_BIN1_T2DM = get_expression(structure_ids = results_BIN1_T2DM$structure_id,
                                      gene_ids = BIN1_PPIN,
                                      dataset = '5_stages')
save(expression_BIN1_T2DM, file = "expression_BIN1_T2DM.RData")

expression_PLEKHM1_T2DM = get_expression(structure_ids = results_PLEKHM1_T2DM$structure_id,
                                         gene_ids = PLEKHM1_PPIN,
                                         dataset = '5_stages')
save(expression_PLEKHM1_T2DM, file = "expression_PLEKHM1_T2DM.RData")

################################################################################
#Depression gene set analysis
################################################################################

depression_genes = c("LINC02210", "ZSCAN26", "ZNF603P", "ZSCAN31", "KANSL1", 
                     "KANSL1-AS1","LRRC37A4P", "RP11-259G18.3", "MAPK8IP1P2", 
                     "RP11-259G18.1", "MAPT-AS1", "CELF1", "LRRC37A17P", "WNT3",
                     "ZNF165", "RP11-707O23.1", "DND1P1", "MAPK8IP1P1",
                     "MAPT-IT1", "RP11-669E14.4", "MAPT", "CTB-39G8.3", "NUP160", 
                     "RAPSN","ARHGAP27", "ZKSCAN3", "MADD", "UBXN2A", "PLEKHM1",
                     "CHRNE", "KAT8", "LACTB")

depression_input = data.frame(depression_genes, is_candidate = 1)

res_depression_devel = aba_enrich(depression_input,
                                  dataset='5_stages', 
                                  cutoff_quantiles = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                  n_randsets = 10000)

results_depression = res_depression_devel$results

write.csv(results_depression, "ABA_enrichment_results_depression.csv")
