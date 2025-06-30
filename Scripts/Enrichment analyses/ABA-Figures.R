################################################################################
#Allen Brain Atlas - Gene set enrichment
#Zillah Daysh
#2025
################################################################################
#R version 4.4.1 

#load packages
library(tidyverse)
library(pheatmap)
library(ggh4x)

################################################################################
#IBD/SLE Analysis
################################################################################
#Enrichment figure

results_IBD_SLE_genes_devel = read.csv("results_IBD_genes_devel.csv")

#filter data to only keep significant brain regions
results_IBD_SLE_filtered = results_IBD_SLE_genes_devel %>%
  group_by(structure) %>%
  filter(!all(n_significant == "0")) %>%
  mutate(age_category = as.factor(age_category)) %>%
  ungroup()

#give full names to developmental stages
results_IBD_SLE_filtered$age_definition = (results_IBD_SLE_filtered$age_category) %>%
  fct_recode("Prenatal" = "1",
             "Infant" = "2",
             "Child" = "3",
             "Adolescent" = "4",
             "Adult" = "5")

palette = c(
  "#B32357",  '#35A7FF', "#D55E00",  "#6ACC65", "#911eb4", "#F781BF", "#009E73", "#FF6F61",
  "#0072B2",  "#b66dff", "#FFC000", "#117733",    "#00CED1", "#808000",  
  "#F0E442",  "#FF8C00",   "#D33682", 
  "#4682B4",  "#000000")

#Figure 3-3A : enrichment plot
ggplot(results_IBD_SLE_filtered, aes(y = structure, x = min_FWER, colour = age_definition)) + 
  geom_point(size = 3)+
  geom_vline(xintercept = 0.05, lty = 2,lwd = 0.75)+
  theme_bw(base_size = 16)+
  labs(colour = "Age Category", y = "Brain Structure", x = "FWER")+
  scale_color_manual(values  = palette) + 
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1.0, 0.05),
                     limits = c(0,1))

################################################################################
#Expression figure
load("expression_chrom_remod_genes_devel.RData")

#get brain regions (structures) enriched for gene set
IBD_SLE_ABA_structures = results_IBD_SLE_genes_devel %>%
  filter(!n_significant == 0) %>%
  select(structure_id, structure)

#get expression data for significant age group (prenatal)  
IBD_SLE_expression = expression_chrom_remod_genes_devel$age_category_1

#filter dataset
IBD_SLE_expression = IBD_SLE_expression %>%
  as.data.frame() %>%
  filter(rownames(IBD_SLE_expression) %in% IBD_SLE_ABA_structures$structure_id) %>% #filter to only keep "significant" brain regions
  rownames_to_column() %>%
  rename("structure_id" = "rowname")

#get full structure names
IBD_SLE_expression = left_join(IBD_SLE_expression, IBD_SLE_ABA_structures, by = "structure_id")

rownames(IBD_SLE_expression) = IBD_SLE_expression$structure
IBD_SLE_expression = IBD_SLE_expression[,-c(1,26)]

#Figure 3-3B : expression heatmap
pheatmap(IBD_SLE_expression,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = TRUE,
         cluter_columns = TRUE)

################################################################################
#T2DM analysis
################################################################################
#Enrichment figure

results_BIN1_T2DM = read.csv("ABA_enrichment_results_BIN1_T2DM.csv")
results_KAT8_T2DM = read.csv("ABA_enrichment_results_KAT8_T2DM.csv")
results_PLEKHM1_T2DM = read.csv("ABA_enrichment_results_PLEKHM1_T2DM.csv")

results_BIN1_T2DM$PPIN = "BIN1"
results_KAT8_T2DM$PPIN = "KAT8"
results_PLEKHM1_T2DM$PPIN = "PLEKHM1"

#join datasets
T2DM_full_results = list(results_BIN1_T2DM,
                         results_KAT8_T2DM,
                         results_PLEKHM1_T2DM) %>%
  reduce(full_join)

#filter to keep structures that have significant enrichment
T2DM_full_results = T2DM_full_results %>%
  group_by(PPIN, structure) %>%
  filter(!all(n_significant == "0")) %>%
  mutate(age_category = as.factor(age_category)) %>%
  ungroup()

#get full name
T2DM_full_results$age_definition = (T2DM_full_results$age_category) %>%
  fct_recode("Prenatal" = "1",
             "Infant" = "2",
             "Child" = "3",
             "Adolescent" = "4",
             "Adult" = "5")

palette = c(
  "#B32357",  '#35A7FF', "#D55E00",  "#6ACC65", "#911eb4", "#F781BF", "#009E73", "#FF6F61",
  "#0072B2",  "#b66dff", "#FFC000", "#117733",    "#00CED1", "#808000",  
  "#F0E442",  "#FF8C00",   "#D33682", 
  "#4682B4",  "#000000")


#setting design for plot
design = c(1,2,3)

#Figure 3-9
ggplot(T2DM_full_results, aes(y = structure, x = min_FWER, colour = age_definition)) + 
  geom_point(size = 3)+
  geom_vline(xintercept = 0.05, lty = 2,lwd = 0.5)+
  theme_bw(base_size = 16)+
  labs(colour = "Age Category", y = "Brain Structure", x = "FWER")+
  scale_color_manual(values  = palette) + 
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1.0, 0.05),
                     limits = c(0, 1))+
  facet_manual(~PPIN,
               scales = "free",
               design = design,
               heights = c(9,6,4))

################################################################################
#Expression figure

#BIN1

T2D_BIN1_ABA_structures = results_BIN1_T2DM %>%
  filter(!n_significant == 0) %>%
  select(structure_id, structure, age_category)


load("expression_BIN1_T2DM.RData")
BIN1_brain_expression_1 = expression_BIN1_T2DM$age_category_1
BIN1_brain_expression_2 = expression_BIN1_T2DM$age_category_2
BIN1_brain_expression_4 = expression_BIN1_T2DM$age_category_4
BIN1_brain_expression_5 = expression_BIN1_T2DM$age_category_5


BIN1_T2D_genes_fetal = BIN1_brain_expression_1  %>%
  as.data.frame() %>%
  filter(rownames(BIN1_brain_expression_1) %in% (T2D_BIN1_ABA_structures %>%
                                                   filter(age_category == "1") %>%
                                                   pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


BIN1_T2D_genes_fetal = left_join(
  BIN1_T2D_genes_fetal, (T2D_BIN1_ABA_structures %>%
                           filter(age_category == "1")), by = "structure_id") %>%
  mutate(structure = paste("Prenatal - ", structure))

##########
BIN1_T2D_genes_infant = BIN1_brain_expression_2  %>%
  as.data.frame() %>%
  filter(rownames(BIN1_brain_expression_2) %in% (T2D_BIN1_ABA_structures %>%
                                                   filter(age_category == "2") %>%
                                                   pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


BIN1_T2D_genes_infant = left_join(
  BIN1_T2D_genes_infant, (T2D_BIN1_ABA_structures %>%
                            filter(age_category == "2")), by = "structure_id")%>%
  mutate(structure = paste("Infant - ", structure))


##########
BIN1_T2D_genes_adolescent = BIN1_brain_expression_4  %>%
  as.data.frame() %>%
  filter(rownames(BIN1_brain_expression_4) %in% (T2D_BIN1_ABA_structures %>%
                                                   filter(age_category == "4") %>%
                                                   pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


BIN1_T2D_genes_adolescent = left_join(
  BIN1_T2D_genes_adolescent, (T2D_BIN1_ABA_structures %>%
                                filter(age_category == "4")), by = "structure_id") %>%
  mutate(structure = paste("Adolescent - ", structure))


##########
BIN1_T2D_genes_adult = BIN1_brain_expression_5  %>%
  as.data.frame() %>%
  filter(rownames(BIN1_brain_expression_5) %in% (T2D_BIN1_ABA_structures %>%
                                                   filter(age_category == "5") %>%
                                                   pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


BIN1_T2D_genes_adult = left_join(
  BIN1_T2D_genes_adult, (T2D_BIN1_ABA_structures %>%
                           filter(age_category == "5")), by = "structure_id") %>%
  mutate(structure = paste("Adult - ", structure))

BIN1_full_expression = list(BIN1_T2D_genes_fetal, 
                            BIN1_T2D_genes_infant,
                            BIN1_T2D_genes_adolescent,
                            BIN1_T2D_genes_adult) %>%
  purrr::reduce(full_join)%>%
  pivot_longer(cols = c(2:26), names_to = "gene", values_to = "ave_expression") %>%
  mutate(ave_expression = ifelse(ave_expression > 100, 
                                 100, 
                                 ave_expression)) %>%
  pivot_wider(names_from = "gene", values_from = "ave_expression") %>%
  as.data.frame()


rownames(BIN1_full_expression) = BIN1_full_expression$structure
BIN1_full_expression = BIN1_full_expression[,-c(1:3)]

#Figure 3-10A
pheatmap(BIN1_full_expression,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE)

################################################################################
#KAT8 

T2D_KAT8_ABA_structures = results_KAT8_T2DM %>%
  filter(!n_significant == 0) %>%
  select(structure_id, structure)

load("expression_KAT8_T2DM.RData")
KAT8_brain_expression = expression_KAT8_T2DM$age_category_1

KAT8_T2D_genes = KAT8_brain_expression  %>%
  as.data.frame() %>%
  filter(rownames(KAT8_brain_expression ) %in% T2D_KAT8_ABA_structures $structure_id) %>%
  rownames_to_column()%>%
  rename("structure_id" = "rowname") %>%
  pivot_longer(cols = c(2:22), names_to = "gene", values_to = "ave_expression") %>%
  mutate(ave_expression = ifelse(ave_expression > 100, 
                                 100, 
                                 ave_expression)) %>%
  pivot_wider(names_from = "gene", values_from = "ave_expression")

KAT8_T2D_genes_1 = left_join(
  KAT8_T2D_genes, T2D_KAT8_ABA_structures, by = "structure_id") %>%
  as.data.frame()

rownames(KAT8_T2D_genes_1) = paste("Prenatal -",KAT8_T2D_genes_1$structure)
KAT8_T2D_genes_1  = KAT8_T2D_genes_1[,-c(1,23)]

#Figure 3-10B
pheatmap(KAT8_T2D_genes_1,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 12)

################################################################################
#PLEKHM1

T2D_PLEKHM1_ABA_structures = results_PLEKHM1_T2DM %>%
  filter(!n_significant == 0) %>%
  select(structure_id, structure, age_category)


load("expression_PLEKHM1_T2DM.RData")
PLEKHM1_brain_expression_1 = expression_PLEKHM1_T2DM$age_category_1
PLEKHM1_brain_expression_2 = expression_PLEKHM1_T2DM$age_category_2
PLEKHM1_brain_expression_3 = expression_PLEKHM1_T2DM$age_category_3
PLEKHM1_brain_expression_4 = expression_PLEKHM1_T2DM$age_category_4



PLEKHM1_T2D_genes_fetal = PLEKHM1_brain_expression_1  %>%
  as.data.frame() %>%
  filter(rownames(PLEKHM1_brain_expression_1) %in% (T2D_PLEKHM1_ABA_structures %>%
                                                      filter(age_category == "1") %>%
                                                      pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


PLEKHM1_T2D_genes_fetal = left_join(
  PLEKHM1_T2D_genes_fetal, (T2D_PLEKHM1_ABA_structures %>%
                              filter(age_category == "1")), by = "structure_id") %>%
  mutate(structure = paste("Prenatal - ", structure))


##########
PLEKHM1_T2D_genes_infant = PLEKHM1_brain_expression_2  %>%
  as.data.frame() %>%
  filter(rownames(PLEKHM1_brain_expression_2) %in% (T2D_PLEKHM1_ABA_structures %>%
                                                      filter(age_category == "2") %>%
                                                      pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


PLEKHM1_T2D_genes_infant = left_join(
  PLEKHM1_T2D_genes_infant, (T2D_PLEKHM1_ABA_structures %>%
                               filter(age_category == "2")), by = "structure_id") %>%
  mutate(structure = paste("Infant - ", structure))


##########
PLEKHM1_T2D_genes_adolescent = PLEKHM1_brain_expression_4  %>%
  as.data.frame() %>%
  filter(rownames(PLEKHM1_brain_expression_4) %in% (T2D_PLEKHM1_ABA_structures %>%
                                                      filter(age_category == "4") %>%
                                                      pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


PLEKHM1_T2D_genes_adolescent = left_join(
  PLEKHM1_T2D_genes_adolescent, (T2D_PLEKHM1_ABA_structures %>%
                                   filter(age_category == "4")), by = "structure_id") %>%
  mutate(structure = paste("Adolescent - ", structure))


##########
PLEKHM1_T2D_genes_child = PLEKHM1_brain_expression_3  %>%
  as.data.frame() %>%
  filter(rownames(PLEKHM1_brain_expression_3) %in% (T2D_PLEKHM1_ABA_structures %>%
                                                      filter(age_category == "3") %>%
                                                      pull(structure_id)))%>%
  rownames_to_column()%>%
  rename("structure_id" ="rowname")


PLEKHM1_T2D_genes_child = left_join(
  PLEKHM1_T2D_genes_child, (T2D_PLEKHM1_ABA_structures %>%
                              filter(age_category == "3")), by = "structure_id") %>%
  mutate(structure = paste("Child - ", structure))


PLEKHM1_full_expression = list(PLEKHM1_T2D_genes_fetal, 
                               PLEKHM1_T2D_genes_infant,
                               PLEKHM1_T2D_genes_child, 
                               PLEKHM1_T2D_genes_adolescent) %>%
  purrr::reduce(full_join) %>%
  pivot_longer(cols = c(2:21), names_to = "gene", values_to = "ave_expression") %>%
  mutate(ave_expression = ifelse(ave_expression > 100, 
                                 100, 
                                 ave_expression)) %>%
  pivot_wider(names_from = "gene", values_from = "ave_expression") %>%
  as.data.frame()


rownames(PLEKHM1_full_expression) = PLEKHM1_full_expression$structure
PLEKHM1_full_expression = PLEKHM1_full_expression[,-c(1:3)]

#Figure 3-10C
pheatmap(PLEKHM1_full_expression,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 12)

