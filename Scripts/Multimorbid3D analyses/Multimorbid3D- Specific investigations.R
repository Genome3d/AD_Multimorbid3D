################################################################################
#Multimorbid investigation into specific PPINs
#Zillah Daysh
#2025
################################################################################
#R version 4.4.1

#load packages
library(tidyverse)
library(biomaRt)

################################################################################

#Full multimorbid results - from Multimorbid- dataset creation R file
multimorbid_all_tissues_gene_traits = read.csv("multimorbid_all_tissues_gene_traits.csv")

#datasets that show links between L0 genes to L1 to L4 genes
adult_cortex_ppin = read.delim("brain_morbid/graph_filtered.txt")
fetal_cortex_ppin = read.delim("fetal_morbid/graph_filtered.txt")
aorta_ppin = read.delim("aorta_morbid/graph_filtered.txt")
coronary_ppin = read.delim("coronary_morbid/graph_filtered.txt")
liver_ppin = read.delim("liver_morbid/graph_filtered.txt")
lung_ppin = read.delim("lung_morbid/graph_filtered.txt")
skin_exposed_ppin = read.delim("skin_exposed_morbid/graph_filtered.txt")
skin_unexposed_ppin = read.delim("skin_unexposed_morbid/graph_filtered.txt")
blood_ppin = read.delim("blood_morbid/graph_filtered.txt")


# List of tissue PPINs
tissue_ppins = list(
  aorta = aorta_ppin,
  liver = liver_ppin,
  adult_cortex = adult_cortex_ppin,
  fetal_cortex = fetal_cortex_ppin,
  lung = lung_ppin, 
  coronary = coronary_ppin, 
  skin_exposed = skin_exposed_ppin, 
  skin_unexposed = skin_unexposed_ppin, 
  whole_blood = blood_ppin
)

# List of genes
genes = c("LACTB", 
          "INO80E",
          "KAT8", 
          "ACE", 
          "APOC2",
          "BIN1",
          "EARS2",
          "PLEKHM1")

# Placeholder for all tissue pathways and gene lists
all_tissue_pathways = list()
all_tissue_gene_lists = list()

##############extract out L0 genes of interest and their PPINs##################

# Iterate over each tissue
for(tissue in names(tissue_ppins)) {
  
  ppin = tissue_ppins[[tissue]]
  tissue_pathways = list()
  tissue_gene_lists = list()
  
  # Iterate over each gene
  for(gene in genes) {
    
    if(!gene %in% ppin$geneA) next  # Skip if gene is not in the tissue PPIN
    
    df_level0 = ppin %>% 
      filter(geneA == gene)
    
    L0_interaction = unique(df_level0$geneB)
    
    gene_pathways = data.frame() 
    
    # Iterate over interactions
    for (j in L0_interaction) {
      
      df_level1 = ppin %>% 
        filter(geneA == j)
      
      L1_interaction = unique(df_level1$geneB)
      
      for(k in L1_interaction) {
        
        df_level2 = ppin %>%
          filter(geneA == k)
        
        L2_interaction = unique(df_level2$geneB)
        
        for(l in L2_interaction) {
          
          df_level3 = ppin %>%
            filter(geneA == l)
          
          L3_interaction = unique(df_level3$geneB)
          
          for(m in L3_interaction) {
            
            df_level4 = ppin %>%
              filter(geneA == m)
            
            L4_interaction = unique(df_level4$geneB)
            
            gene_pathways = rbind(gene_pathways,
                                  data.frame(L0 = gene,
                                             L1 = j,
                                             L2 = k,
                                             L3 = l,
                                             L4 = m))  
          }
        }
      }
    }
    
    tissue_pathways[[gene]] = rbind(tissue_pathways[[gene]], gene_pathways)
    
    # Extract unique genes with their level info (L0 to L4) for each pathway
    gene_level_info = data.frame(
      gene = c(gene_pathways$L0, gene_pathways$L1, gene_pathways$L2, gene_pathways$L3, gene_pathways$L4),
      level = c(rep("L0", length(gene_pathways$L0)),
                rep("L1", length(gene_pathways$L1)),
                rep("L2", length(gene_pathways$L2)),
                rep("L3", length(gene_pathways$L3)),
                rep("L4", length(gene_pathways$L4)))
    )
    tissue_gene_lists[[gene]] = rbind(tissue_gene_lists[[gene]], unique(gene_level_info))
  }
  
  # Store the pathways and gene lists for this tissue into the master list
  all_tissue_pathways[[tissue]] = tissue_pathways
  all_tissue_gene_lists[[tissue]] = tissue_gene_lists

}

################################################################################
#INO80E
################################################################################

#extract significant traits 
#PPINs should be the same regardless of tissue 
#so filtering by one of the tissues where the gene has been observed in 
Sig_PPIN_traits_INO80E = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
         paste(all_tissue_gene_lists$skin_exposed$INO80E$gene, all_tissue_gene_lists$skin_exposed$INO80E$level)) %>%
  dplyr::select(gene,level,trait,snp,tissue) %>%
  distinct() %>%
  arrange(gene)

#get pathways linking genes  
INO80E_ppin_pathways = all_tissue_pathways$skin_exposed$INO80E %>%
  filter(L1 %in% (Sig_PPIN_traits_INO80E %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_INO80E %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_INO80E %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_INO80E %>%
                      filter(level == "L4") %>%
                      pull(gene)))

INO80E_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_INO80E$gene) %>%
  mutate(PPIN = "INO80E")


################################################################################
#LACTB
################################################################################

Sig_PPIN_traits_LACTB = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$skin_exposed$LACTB$gene, all_tissue_gene_lists$skin_exposed$LACTB$level)) %>%
  dplyr::select(gene,level,trait,snp, tissue) %>%
  distinct() %>%
  arrange(gene)

LACTB_ppin_pathways = all_tissue_pathways$skin_exposed$LACTB %>%
  filter(L1 %in% (Sig_PPIN_traits_LACTB %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_LACTB %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_LACTB %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_LACTB %>%
                      filter(level == "L4") %>%
                      pull(gene)))

LACTB_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_LACTB$gene) %>%
  mutate(PPIN = "LACTB")

################################################################################
#KAT8
################################################################################

Sig_PPIN_traits_KAT8 = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$skin_exposed$KAT8$gene, all_tissue_gene_lists$skin_exposed$KAT8$level)) %>%
  dplyr::select(gene,level,trait,snp, tissue) %>%
  distinct() %>%
  arrange(gene)

KAT8_ppin_pathways = all_tissue_pathways$skin_exposed$KAT8 %>%
  filter(L1 %in% (Sig_PPIN_traits_KAT8 %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_KAT8 %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_KAT8 %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_KAT8 %>%
                      filter(level == "L4") %>%
                      pull(gene)))

KAT8_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_KAT8$gene) %>%
  mutate(PPIN = "KAT8")


################################################################################
#ACE
################################################################################

Sig_PPIN_traits_ACE = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$skin_exposed$ACE$gene, all_tissue_gene_lists$skin_exposed$ACE$level)) %>%
  dplyr::select(gene,level,trait,snp,tissue) %>%
  distinct() %>%
  arrange(gene)

ACE_ppin_pathways = all_tissue_pathways$skin_exposed$ACE %>%
  filter(L1 %in% (Sig_PPIN_traits_ACE %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_ACE %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_ACE %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_ACE %>%
                      filter(level == "L4") %>%
                      pull(gene)))

ACE_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_ACE$gene) %>%
  mutate(PPIN = "ACE")
################################################################################
#APOC2
################################################################################

Sig_PPIN_traits_APOC2 = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$skin_unexposed$APOC2$gene, all_tissue_gene_lists$skin_unexposed$APOC2$level)) %>%
  dplyr::select(gene,level,trait,snp,tissue) %>%
  distinct() %>%
  arrange(gene)

APOC2_ppin_pathways = all_tissue_pathways$skin_exposed$APOC2 %>%
  filter(L1 %in% (Sig_PPIN_traits_APOC2 %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_APOC2 %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_APOC2 %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_APOC2 %>%
                      filter(level == "L4") %>%
                      pull(gene)))

APOC2_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_APOC2$gene) %>%
  mutate(PPIN = "APOC2")
################################################################################
#BIN1
################################################################################

Sig_PPIN_traits_BIN1 = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$liver$BIN1$gene, all_tissue_gene_lists$liver$BIN1$level)) %>%
  dplyr::select(gene,level,trait,snp, tissue) %>%
  distinct() %>%
  arrange(gene)

BIN1_ppin_pathways = all_tissue_pathways$liver$BIN1 %>%
  filter(L1 %in% (Sig_PPIN_traits_BIN1 %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_BIN1 %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_BIN1 %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_BIN1 %>%
                      filter(level == "L4") %>%
                      pull(gene)))

BIN1_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_BIN1$gene) %>%
  mutate(PPIN = "BIN1")
################################################################################
#EARS2
################################################################################

Sig_PPIN_traits_EARS2 = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$skin_exposed$EARS2$gene, all_tissue_gene_lists$skin_exposed$EARS2$level)) %>%
  dplyr::select(gene,level,trait,snp, tissue) %>%
  distinct() %>%
  arrange(gene)

EARS2_ppin_pathways = all_tissue_pathways$skin_exposed$EARS2 %>%
  filter(L1 %in% (Sig_PPIN_traits_EARS2 %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_EARS2 %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_EARS2 %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_EARS2 %>%
                      filter(level == "L4") %>%
                      pull(gene)))

EARS2_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_EARS2$gene) %>%
  mutate(PPIN = "EARS2")
################################################################################
#PLEKHM1
################################################################################

Sig_PPIN_traits_PLEKHM1 = multimorbid_all_tissues_gene_traits %>%
  filter(paste(gene,level) %in%
           paste(all_tissue_gene_lists$lung$PLEKHM1$gene, all_tissue_gene_lists$lung$PLEKHM1$level)) %>%
  dplyr::select(gene,level,trait,snp, tissue) %>%
  distinct() %>%
  arrange(gene)

PLEKHM1_ppin_pathways = all_tissue_pathways$lung$PLEKHM1 %>%
  filter(L1 %in% (Sig_PPIN_traits_PLEKHM1 %>%
                    filter(level == "L1") %>%
                    pull(gene)) |
           L2 %in% (Sig_PPIN_traits_PLEKHM1 %>%
                      filter(level == "L2") %>%
                      pull(gene)) |
           L3 %in% (Sig_PPIN_traits_PLEKHM1 %>%
                      filter(level == "L3") %>%
                      pull(gene))|
           L4 %in% (Sig_PPIN_traits_PLEKHM1 %>%
                      filter(level == "L4") %>%
                      pull(gene)))

PLEKHM1_multimorbid_traits = multimorbid_all_tissues_gene_traits %>%
  filter(gene %in% Sig_PPIN_traits_PLEKHM1$gene) %>%
  mutate(PPIN = "PLEKHM1")

################################################################################
#Combine for summary figures / tables 

joint_gene_multimorbid = list(KAT8_multimorbid_traits,
                              INO80E_multimorbid_traits,
                              LACTB_multimorbid_traits,
                              ACE_multimorbid_traits,
                              APOC2_multimorbid_traits,
                              BIN1_multimorbid_traits,
                              EARS2_multimorbid_traits,
                              PLEKHM1_multimorbid_traits) %>%
  purrr::reduce(full_join)

joint_gene_multimorbid = joint_gene_multimorbid[-1]
write.csv(joint_gene_multimorbid, "joint_gene_multimorbid.csv")

################################################################################
#FOR PPIN FIGURES 
################################################################################
#Look at files for particular PPIN: 
#Sig_PPIN_traits_*insert PPIN name* 
#*insert PPIN name*_ppin_pathways

##########For full PPIN
#use all genes within PPIN pathways that link to L1 to L4 genes
#that are regulated by trait-associated sceQTLs

##########For trait-specific PPINs##
#Find trait-specific genes in Sig_PPIN_traits_*insert PPIN name* 
#Then look in *insert PPIN name*_ppin_pathways to find specific pathways linking
#L0 gene to trait-sceQTL regulated gene

#input that list into STRING web tool
#https://string-db.org/
#change settings to only above 0.9 and removing text mining 
#change to level of confidence lines
#save 
#figures created in adobe illustrator from output from STRING

################################################################################
#Depression investigation - no PPINs as all genes L0 or L1
################################################################################

depression = multimorbid_all_tissues_gene_traits %>%
  filter(grepl("depress", trait, ignore.case = TRUE))

depression_matrix = as.data.frame.matrix(table(depression$gene, 
                                               depression$trait)) %>% 
  mutate_if(is.numeric, ~ ifelse(. > 50, 50, .))

#Figure 3-12
pheatmap(t(depression_matrix),
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = TRUE,
         cluter_columns = TRUE,
         cellwidth = 15,
         cellheight = 15,
         fontsize = 13)

depression_genes = depression %>%
  dplyr::select(gene) %>%
  distinct() 

write.csv(depression_genes, "depression_genes.csv")

#gene mart for annotating risk genes
gene_mart=useMart(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl")

#use ensembl to get location of risk genes
depression_genes_location = getBM(attributes=c('hgnc_symbol',
                                               'chromosome_name',
                                               'start_position',
                                               'end_position',
                                               'band'),
                                  filters = ("hgnc_symbol"),
                                  values=list(depression_genes$gene),
                                  mart=gene_mart) 

#Find L0 genes for L1 depression genes
adult_cortex_ppin %>%
  filter(geneB == "KANSL1")

skin_exposed_ppin %>%
  filter(geneB == "RAPSN")

skin_exposed_ppin %>%
  filter(geneB == "UBXN2A")
