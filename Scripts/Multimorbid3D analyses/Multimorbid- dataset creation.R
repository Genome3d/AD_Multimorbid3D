################################################################################
##Multimorbid3D Analysis
#Zillah Daysh
#2025
################################################################################
#R version 4.4.1

#load packages
library(tidyverse)
library(forcats)
library(scales)
library(viridis)
library(pheatmap)
library(gprofiler2)
library(wesanderson)

#data from Multimorbid3D pipeline
#https://github.com/Genome3d/multimorbid3D

################################################################################
#Adult Cortex Multimorbid
################################################################################

#read in significant traits 
adult_cortex_multimorbid_sig_traits_bt = read.delim("brain_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

adult_cortex_multimorbid_sig_traits_bt$Level = fct_recode(adult_cortex_multimorbid_sig_traits_bt$Level, 
                                                          "L0" = "level0",
                                                          "L1" = "level1",
                                                          "L2" = "level2",
                                                          "L3" = "level3",
                                                          "L4" = "level4")

#get extended info on gene-snps involved in significant traits - level by level
adult_cortex_level_0_gene_traits = read.delim("brain_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in%  (adult_cortex_multimorbid_sig_traits_bt %>%
                        filter(Level == "L0") %>%
                        pull(Trait)))

adult_cortex_level_1_gene_traits = read.delim("brain_morbid/level1_sig_interactions.txt") %>%
  filter(trait %in% (adult_cortex_multimorbid_sig_traits_bt %>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

adult_cortex_level_2_gene_traits = read.delim("brain_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in%  (adult_cortex_multimorbid_sig_traits_bt %>%
                        filter(Level == "L2") %>%
                        pull(Trait)))

adult_cortex_level_3_gene_traits = read.delim("brain_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in%  (adult_cortex_multimorbid_sig_traits_bt %>%
                        filter(Level == "L3") %>%
                        pull(Trait)))

adult_cortex_level_4_gene_traits = read.delim("brain_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in%  (adult_cortex_multimorbid_sig_traits_bt %>%
                        filter(Level == "L4") %>%
                        pull(Trait)))

adult_cortex_level_0_gene_traits$level = "L0"
adult_cortex_level_1_gene_traits$level = "L1"
adult_cortex_level_2_gene_traits$level = "L2"
adult_cortex_level_3_gene_traits$level = "L3"
adult_cortex_level_4_gene_traits$level = "L4"

#join into one dataset
adult_cortex_all_gene_traits = list(adult_cortex_level_0_gene_traits,
                                    adult_cortex_level_1_gene_traits,
                                    adult_cortex_level_2_gene_traits,
                                    adult_cortex_level_3_gene_traits,
                                    adult_cortex_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Aorta Multimorbid
################################################################################

aorta_multimorbid_sig_traits_bt = read.delim("aorta_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

aorta_multimorbid_sig_traits_bt$Level = fct_recode(aorta_multimorbid_sig_traits_bt$Level, 
                                                   "L0" = "level0",
                                                   "L1" = "level1",
                                                   "L2" = "level2",
                                                   "L3" = "level3",
                                                   "L4" = "level4")

aorta_level_0_gene_traits = read.delim("aorta_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (aorta_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))

aorta_level_1_gene_traits = read.delim("aorta_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (aorta_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

aorta_level_2_gene_traits = read.delim("aorta_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (aorta_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

aorta_level_3_gene_traits = read.delim("aorta_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (aorta_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

aorta_level_4_gene_traits = read.delim("aorta_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (aorta_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))

#info regarding # of genes, traits, and eQTLs at each level and total
aorta_level_0_gene_traits$level = "L0"
aorta_level_1_gene_traits$level = "L1"
aorta_level_2_gene_traits$level = "L2"
aorta_level_3_gene_traits$level = "L3"
aorta_level_4_gene_traits$level = "L4"

aorta_all_gene_traits = list(aorta_level_0_gene_traits,
                             aorta_level_1_gene_traits,
                             aorta_level_2_gene_traits,
                             aorta_level_3_gene_traits,
                             aorta_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Coronary artery Multimorbid
################################################################################

coronary_multimorbid_sig_traits_bt = read.delim("coronary_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

coronary_multimorbid_sig_traits_bt$Level = fct_recode(coronary_multimorbid_sig_traits_bt$Level, 
                                                      "L0" = "level0",
                                                      "L1" = "level1",
                                                      "L2" = "level2",
                                                      "L3" = "level3",
                                                      "L4" = "level4")

coronary_level_0_gene_traits = read.delim("coronary_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (coronary_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))

coronary_level_1_gene_traits = read.delim("coronary_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (coronary_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

coronary_level_2_gene_traits = read.delim("coronary_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (coronary_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

coronary_level_3_gene_traits = read.delim("coronary_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (coronary_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

coronary_level_4_gene_traits = read.delim("coronary_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (coronary_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))

#info regarding # of genes, traits, and eQTLs at each level and total
coronary_level_0_gene_traits$level = "L0"
coronary_level_1_gene_traits$level = "L1"
coronary_level_2_gene_traits$level = "L2"
coronary_level_3_gene_traits$level = "L3"
coronary_level_4_gene_traits$level = "L4"

coronary_all_gene_traits = list(coronary_level_0_gene_traits,
                                coronary_level_1_gene_traits,
                                coronary_level_2_gene_traits,
                                coronary_level_3_gene_traits,
                                coronary_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Fetal Cortex Multimorbid
################################################################################
fetal_cortex_multimorbid_sig_traits_bt = read.delim("fetal_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

fetal_cortex_multimorbid_sig_traits_bt$Level = fct_recode(fetal_cortex_multimorbid_sig_traits_bt$Level, 
                                                          "L0" = "level0",
                                                          "L1" = "level1",
                                                          "L2" = "level2",
                                                          "L3" = "level3",
                                                          "L4" = "level4")

fetal_cortex_level_0_gene_traits = read.delim("fetal_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (fetal_cortex_multimorbid_sig_traits_bt %>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


fetal_cortex_level_1_gene_traits = read.delim("fetal_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (fetal_cortex_multimorbid_sig_traits_bt %>%
                       filter(Level == "L1") %>%
                       pull(Trait)))


fetal_cortex_level_2_gene_traits = read.delim("fetal_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in%  (fetal_cortex_multimorbid_sig_traits_bt %>%
                        filter(Level == "L2") %>%
                        pull(Trait)))


fetal_cortex_level_3_gene_traits = read.delim("fetal_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in%  (fetal_cortex_multimorbid_sig_traits_bt%>%
                        filter(Level == "L3") %>%
                        pull(Trait)))

fetal_cortex_level_4_gene_traits = read.delim("fetal_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in%  (fetal_cortex_multimorbid_sig_traits_bt%>%
                        filter(Level == "L4") %>%
                        pull(Trait)))

#info regarding # of genes, traits, and eQTLs at each level and total
fetal_cortex_level_0_gene_traits$level = "L0"
fetal_cortex_level_1_gene_traits$level = "L1"
fetal_cortex_level_2_gene_traits$level = "L2"
fetal_cortex_level_3_gene_traits$level = "L3"
fetal_cortex_level_4_gene_traits$level = "L4"

fetal_cortex_all_gene_traits = list(fetal_cortex_level_0_gene_traits,
                                    fetal_cortex_level_1_gene_traits,
                                    fetal_cortex_level_2_gene_traits,
                                    fetal_cortex_level_3_gene_traits,
                                    fetal_cortex_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Liver Multimorbid
################################################################################
liver_multimorbid_sig_traits_bt = read.delim("liver_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

liver_multimorbid_sig_traits_bt$Level = fct_recode(liver_multimorbid_sig_traits_bt$Level, 
                                                   "L0" = "level0",
                                                   "L1" = "level1",
                                                   "L2" = "level2",
                                                   "L3" = "level3",
                                                   "L4" = "level4")

liver_level_0_gene_traits = read.delim("liver_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (liver_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


liver_level_1_gene_traits = read.delim("liver_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (liver_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

liver_level_2_gene_traits = read.delim("liver_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (liver_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

liver_level_3_gene_traits = read.delim("liver_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (liver_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

liver_level_4_gene_traits = read.delim("liver_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (liver_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))

#info regarding # of genes, traits, and eQTLs at each level and total
liver_level_0_gene_traits$level = "L0"
liver_level_1_gene_traits$level = "L1"
liver_level_2_gene_traits$level = "L2"
liver_level_3_gene_traits$level = "L3"
liver_level_4_gene_traits$level = "L4"

liver_all_gene_traits = list(liver_level_0_gene_traits,
                             liver_level_1_gene_traits,
                             liver_level_2_gene_traits,
                             liver_level_3_gene_traits,
                             liver_level_4_gene_traits) %>%
  reduce(full_join)



################################################################################
#Lung Multimorbid
################################################################################

lung_multimorbid_sig_traits_bt = read.delim("lung_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

lung_multimorbid_sig_traits_bt$Level = fct_recode(lung_multimorbid_sig_traits_bt$Level, 
                                                  "L0" = "level0",
                                                  "L1" = "level1",
                                                  "L2" = "level2",
                                                  "L3" = "level3",
                                                  "L4" = "level4")

lung_level_0_gene_traits = read.delim("lung_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (lung_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


lung_level_1_gene_traits = read.delim("lung_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (lung_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

lung_level_2_gene_traits = read.delim("lung_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (lung_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

lung_level_3_gene_traits = read.delim("lung_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (lung_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

lung_level_4_gene_traits = read.delim("lung_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (lung_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))


#info regarding # of genes, traits, and eQTLs at each level and total
lung_level_0_gene_traits$level = "L0"
lung_level_1_gene_traits$level = "L1"
lung_level_2_gene_traits$level = "L2"
lung_level_3_gene_traits$level = "L3"
lung_level_4_gene_traits$level = "L4"

lung_all_gene_traits = list(lung_level_0_gene_traits,
                            lung_level_1_gene_traits,
                            lung_level_2_gene_traits,
                            lung_level_3_gene_traits,
                            lung_level_4_gene_traits) %>%
  reduce(full_join)


################################################################################
#Skin unexposed Multimorbid
################################################################################

skin_unexposed_multimorbid_sig_traits_bt = read.delim("skin_unexposed_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

skin_unexposed_multimorbid_sig_traits_bt$Level = fct_recode(skin_unexposed_multimorbid_sig_traits_bt$Level, 
                                                            "L0" = "level0",
                                                            "L1" = "level1",
                                                            "L2" = "level2",
                                                            "L3" = "level3",
                                                            "L4" = "level4")

skin_unexposed_level_0_gene_traits = read.delim("skin_unexposed_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (skin_unexposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


skin_unexposed_level_1_gene_traits = read.delim("skin_unexposed_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (skin_unexposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

skin_unexposed_level_2_gene_traits = read.delim("skin_unexposed_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (skin_unexposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

skin_unexposed_level_3_gene_traits = read.delim("skin_unexposed_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (skin_unexposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

skin_unexposed_level_4_gene_traits = read.delim("skin_unexposed_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (skin_unexposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))


#info regarding # of genes, traits, and eQTLs at each level and total
skin_unexposed_level_0_gene_traits$level = "L0"
skin_unexposed_level_1_gene_traits$level = "L1"
skin_unexposed_level_2_gene_traits$level = "L2"
skin_unexposed_level_3_gene_traits$level = "L3"
skin_unexposed_level_4_gene_traits$level = "L4"

skin_unexposed_all_gene_traits = list(skin_unexposed_level_0_gene_traits,
                                      skin_unexposed_level_1_gene_traits,
                                      skin_unexposed_level_2_gene_traits,
                                      skin_unexposed_level_3_gene_traits,
                                      skin_unexposed_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Skin exposed Multimorbid
################################################################################

skin_exposed_multimorbid_sig_traits_bt = read.delim("skin_exposed_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

skin_exposed_multimorbid_sig_traits_bt$Level = fct_recode(skin_exposed_multimorbid_sig_traits_bt$Level, 
                                                          "L0" = "level0",
                                                          "L1" = "level1",
                                                          "L2" = "level2",
                                                          "L3" = "level3",
                                                          "L4" = "level4")


skin_exposed_level_0_gene_traits = read.delim("skin_exposed_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (skin_exposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


skin_exposed_level_1_gene_traits = read.delim("skin_exposed_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (skin_exposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

skin_exposed_level_2_gene_traits = read.delim("skin_exposed_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (skin_exposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

skin_exposed_level_3_gene_traits = read.delim("skin_exposed_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (skin_exposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

skin_exposed_level_4_gene_traits = read.delim("skin_exposed_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (skin_exposed_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))


#info regarding # of genes, traits, and eQTLs at each level and total
skin_exposed_level_0_gene_traits$level = "L0"
skin_exposed_level_1_gene_traits$level = "L1"
skin_exposed_level_2_gene_traits$level = "L2"
skin_exposed_level_3_gene_traits$level = "L3"
skin_exposed_level_4_gene_traits$level = "L4"

skin_exposed_all_gene_traits = list(skin_exposed_level_0_gene_traits,
                                    skin_exposed_level_1_gene_traits,
                                    skin_exposed_level_2_gene_traits,
                                    skin_exposed_level_3_gene_traits,
                                    skin_exposed_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Whole Blood Multimorbid
################################################################################

blood_multimorbid_sig_traits_bt = read.delim("blood_morbid/significant_enrichment_bootstrap.txt") %>%
  filter(sim_pval <= 0.05)%>%
  rename(Level = "level",
         '# sceQTL' = "trait_eqtls",
         Trait = "trait") %>%
  arrange(Trait)

blood_multimorbid_sig_traits_bt$Level = fct_recode(blood_multimorbid_sig_traits_bt$Level, 
                                                   "L0" = "level0",
                                                   "L1" = "level1",
                                                   "L2" = "level2",
                                                   "L3" = "level3",
                                                   "L4" = "level4")


blood_level_0_gene_traits = read.delim("blood_morbid/level0_sig_interactions.txt") %>%
  filter(trait %in% (blood_multimorbid_sig_traits_bt%>%
                       filter(Level == "L0") %>%
                       pull(Trait)))


blood_level_1_gene_traits = read.delim("blood_morbid/level1_sig_interactions.txt")%>%
  filter(trait %in% (blood_multimorbid_sig_traits_bt%>%
                       filter(Level == "L1") %>%
                       pull(Trait)))

blood_level_2_gene_traits = read.delim("blood_morbid/level2_sig_interactions.txt")%>%
  filter(trait %in% (blood_multimorbid_sig_traits_bt%>%
                       filter(Level == "L2") %>%
                       pull(Trait)))

blood_level_3_gene_traits = read.delim("blood_morbid/level3_sig_interactions.txt")%>%
  filter(trait %in% (blood_multimorbid_sig_traits_bt%>%
                       filter(Level == "L3") %>%
                       pull(Trait)))

blood_level_4_gene_traits = read.delim("blood_morbid/level4_sig_interactions.txt")%>%
  filter(trait %in% (blood_multimorbid_sig_traits_bt%>%
                       filter(Level == "L4") %>%
                       pull(Trait)))


#info regarding # of genes, traits, and eQTLs at each level and total
blood_level_0_gene_traits$level = "L0"
blood_level_1_gene_traits$level = "L1"
blood_level_2_gene_traits$level = "L2"
blood_level_3_gene_traits$level = "L3"
blood_level_4_gene_traits$level = "L4"

blood_all_gene_traits = list(blood_level_0_gene_traits,
                             blood_level_1_gene_traits,
                             blood_level_2_gene_traits,
                             blood_level_3_gene_traits,
                             blood_level_4_gene_traits) %>%
  reduce(full_join)

################################################################################
#Total tissue summaries
#General summaries
################################################################################
adult_cortex_all_gene_traits$tissue = "Adult_cortex"
fetal_cortex_all_gene_traits$tissue = "Fetal_cortex"
aorta_all_gene_traits$tissue = "Aorta"
coronary_all_gene_traits$tissue = "Coronary"
liver_all_gene_traits$tissue = "Liver"
lung_all_gene_traits$tissue = "Lung"
skin_unexposed_all_gene_traits$tissue = "Skin_unexposed"
skin_exposed_all_gene_traits$tissue = "Skin_exposed"
blood_all_gene_traits$tissue = "Whole_blood"

multimorbid_all_tissues_gene_traits = list(adult_cortex_all_gene_traits,
                                           fetal_cortex_all_gene_traits,
                                           aorta_all_gene_traits,
                                           coronary_all_gene_traits,
                                           liver_all_gene_traits,
                                           lung_all_gene_traits,
                                           skin_unexposed_all_gene_traits,
                                           skin_exposed_all_gene_traits,
                                           blood_all_gene_traits) %>%
  reduce(full_join)



write.csv(multimorbid_all_tissues_gene_traits, "multimorbid_all_tissues_gene_traits.csv")

