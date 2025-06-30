################################################################################
#Script for enrichment analysis - Gut, Immune, Pancreatic cells
#Zillah Daysh
#2025
################################################################################
#R version 4.4.1

#load libraries
library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(celldex)
library(scRNAseq)

#######################################################################s#########
#Gut data analysis
################################################################################
#import data
#data from : https://www.gutcellatlas.org/ 
#Space-Time Gut Cell Atlas
# Full single cell RNA-seq dataset of 428K intestinal cells from fetal,
#pediatric, adult donors, and up to 11 intestinal regions. 
full_gut_data = readH5AD("Full_obj_log_counts_soupx_v2.h5ad")
#rename assay to norm_counts 
names(assays(full_gut_data))[1] = "norm_counts"

#access metadata
gut_rowdata = rowData(full_gut_data)
gut_coldata = colData(full_gut_data)

#gut gene set
chromatin_remodelling_genes = c("KAT8","MSL1","MSL2","MSL3","KANSL1","KANSL2","KANSL3",
                                "RBBP5","PHF20","WDR5","RUVBL2","KAT5","H4C6","MCRS1",
                                "EP400","MORF4L1","H3C12","TRRAP","YEATS4","UCHL5",
                                "RUVBL1","ACTR5","INO80B","INO80C","TFPT","NFRKB",
                                "INO80","INO80E","PSMD3","ADRM1")

#set gene set and subsetting variables
gene_set_gut = chromatin_remodelling_genes
cell_types_gut = unique(colData(full_gut_data)$category)
tissue = unique(colData(full_gut_data)$Region)
#select only healthy individuals (dataset also has Crohn's paediatric data)
age = unique(colData(full_gut_data)$Diagnosis)[c(1,3,4)] %>%
  as.vector()

#set expression cutoffs for incremental analyses
expression_cutoffs = c(0.5, 0.6, 0.7, 0.8, 0.9)

#list for results to aggregate into
results_full_gut = list()

#number of permutations for FWER calculation
n_permutations = 1000



for (a in age) {
  #subset by age
  age_data = full_gut_data[,colData(full_gut_data)$Diagnosis == a] 
  
for (t in tissue) {
  #subset by tissue
  tissue_data = age_data[, colData(age_data)$Region== t] 
  
  if (length(tissue_data) == 0) next
  
  #repeat for each  quantile
  for (q in expression_cutoffs) { 
    
    for (ct in cell_types_gut) {
      #keep columns corresponding certain cell type
      specific_cell_type = colnames(tissue_data)[colData(tissue_data)$category== ct] 
      
      if (length(specific_cell_type) == 0) next
      
      #get RNA seq matrix data for cell type
      cell_expression_matrix = assay(tissue_data, "norm_counts")[, specific_cell_type]
      
      #get background gene set for that age:tissue:cell combo
      background_cell_genes = rownames(cell_expression_matrix)
      
      #get average expression for genes within cell type
      global_cell_expression = rowMeans(cell_expression_matrix)
      
      #establish quantile cutoff thresholds based on mean gene expression
      threshold = quantile(global_cell_expression, probs = q)
      
      #get names of expressed genes which are above given threshold
      expressed_genes = names(global_cell_expression[global_cell_expression > threshold])
      
      #get length of genes within gene set also expressed above threshold
      overlap = length(intersect(gene_set_gut, expressed_genes))
      
      if (overlap == 0) next
      
      #calculate p value from hypergeometric test
      pval = phyper(
        q = overlap - 1,
        m = length(expressed_genes),
        n = length(background_cell_genes) - length(expressed_genes), #non expressed genes
        k = length(gene_set_gut),
        lower.tail = FALSE
      )
      
      #calculate 1000 permuted p values for validation
      permuted_pvals = replicate(n_permutations, {
        permuted_set = sample(background_cell_genes, length(gene_set_gut))
        perm_overlap = length(intersect(permuted_set, expressed_genes))
        phyper(
          q = perm_overlap - 1,
          m = length(expressed_genes),
          n = length(background_cell_genes) - length(expressed_genes),
          k = length(permuted_set),
          lower.tail = FALSE
        )
      })
      
      #get FWER
      fwer = mean(permuted_pvals <= pval)
      
      #store results
      results_full_gut[[paste0(a, "_", t, "_", ct, "_", q)]] <- data.frame(
        age = a,
        tissue = t,
        cell_type = ct,
        cutoff = q,
        threshold = threshold,
        pval = pval,
        fwer = fwer,
        overlap = length(intersect(gene_set_gut, expressed_genes)))
      
    }
  }
}
}

# Combine all results
results_full_gut = do.call(rbind, results_full_gut)

# Multiple testing correction
results_full_gut$adj_fwer = p.adjust(results_full_gut$fwer, method = "BH")

#arrange by cutoff
results_full_gut = results_full_gut %>%
  arrange(-cutoff)

#save
write.csv(results_full_gut, "enrichment_results_gut_full_final.csv")

#results_full_gut = read.csv("enrichment_results_gut_full_final.csv")

#filter to highest cutoff, only keep age:tissue:cell type combos that are significant
results_full_gut_filtered = results_full_gut %>%
  filter(cutoff == "0.9") %>%
  dplyr::select(age, tissue, cell_type, adj_fwer) %>%
  group_by(cell_type) %>%
  filter(!all(adj_fwer > 0.05)) %>%
  ungroup()

#rename ages
results_full_gut_filtered$Age = factor(results_full_gut_filtered$age) %>%
  fct_recode("Prenatal" = "fetal",
             "Adult" = "Healthy adult", 
             "Paediatric" = "Pediatric healthy") %>%
 fct_relevel(c("Prenatal", "Paediatric", "Adult"))

#rename tissues
results_full_gut_filtered$tissue = factor(results_full_gut_filtered$tissue) %>%
  fct_recode("Appendix" = "APD",
             "Large Intestine" = "LargeInt",
             "Lymph Node" = "lymph node", 
             "Rectum" = "REC", 
             "Small Intestine" = "SmallInt")

#set colour palette
palette = c(
  "#B32357",  '#35A7FF', "#D55E00",  "#6ACC65", "#911eb4", "#F781BF", "#009E73", "#FF6F61",
  "#0072B2",  "#b66dff", "#FFC000", "#117733",    "#00CED1", "#808000",  
  "#F0E442",  "#FF8C00",   "#D33682", 
  "#4682B4",  "#000000")

#Figure 3-4A
ggplot(results_full_gut_filtered , aes(y = cell_type, x = adj_fwer, colour = tissue)) + 
  geom_point(size = 3)+
  geom_vline(xintercept = 0.05, lty = 2)+
  theme_bw(base_size = 16)+
  labs(colour = "Gut Region", y = "Cell type", x = "FWER")+
  scale_color_manual(values  = palette)+
  facet_wrap(~Age) + 
  scale_x_continuous(breaks = c(0.05, 0.25, 0.45, 0.65), 
                     limits = c(0,0.7))

################################################################################
#Gut gene expression figure 3-4B
################################################################################

#Obtain gene expression data for each age-tissue-cell type combination for genes in PPIN
overlapping_gut_genes = intersect(chromatin_remodelling_genes, rownames(full_gut_data))
gut_full_gene_subset = full_gut_data[overlapping_gut_genes,]

#list for results
gut_expression_results = list()


for (a in age) {
  #subset by age
  age_data = gut_full_gene_subset[,colData(gut_full_gene_subset)$Diagnosis == a]
  
  for (t in tissue) {
    #subset by tissue
    tissue_data = age_data[, colData(age_data)$Region == t]
    
      for (ct in cell_types_gut) {
        #select cell type within tissue
        specific_cell_type = colnames(tissue_data)[colData(tissue_data)$category == ct]
        
        if (length(specific_cell_type) == 0) next
        
        #get RNA seq matrix  for cell type
        cell_expression_matrix = assay(tissue_data, "norm_counts")[, specific_cell_type]
        
        #get mean expression for genes within cell type
        global_cell_expression = rowMeans(cell_expression_matrix)
        
        #save results
        gut_expression_results[[paste0(a, "_", t, "_", ct)]] <- data.frame(
          age = a,
          tissue = t,
          cell_type = ct,
          gene = names(global_cell_expression),
          expression = global_cell_expression)
      }
  }
}

gut_expression_results = do.call(rbind, gut_expression_results)

#write.csv(gut_expression_results,"full_gut_expression.csv")

#rename age
gut_expression_results$age = factor(gut_expression_results$age) %>%
  fct_recode("Prenatal" = "fetal",
             "Adult" = "Healthy adult", 
             "Paediatric" = "Pediatric healthy") %>%
  fct_relevel(c("Prenatal", "Paediatric", "Adult"))

#rename tissue
gut_expression_results$tissue = factor(gut_expression_results$tissue) %>%
  fct_recode("Appendix" = "APD",
             "Large Intestine" = "LargeInt",
             "Lymph Node" = "lymph node", 
             "Rectum" = "REC", 
             "Small Intestine" = "SmallInt")

#age-tissue-cell type combo
gut_expression_results$rownames = paste(gut_expression_results$age,
                                        gut_expression_results$tissue,
                                        gut_expression_results$cell_type, 
                                             sep = " - ")

#age-tissue-cell type combo significantly enriched from filtered results
results_full_gut_filtered$rownames = paste(results_full_gut_filtered$Age,
                                           results_full_gut_filtered$tissue,
                                           results_full_gut_filtered$cell_type, 
                                           sep = " - ")

#filter full expression for age-tissue-cell type combo significantly enriched
expression_matrix_1_summary = gut_expression_results %>%
  filter(rownames %in% (results_full_gut_filtered %>%
                          filter(adj_fwer < 0.05) %>%
                          pull(rownames))) %>%
  mutate(expression = ifelse(expression > 1, 1, expression))

#re-arrange data for plotting
expression_matrix_final_plot = expression_matrix_1_summary[4:6] %>%
  pivot_wider(names_from = "gene", values_from = "expression") %>%
  as.data.frame()

rownames(expression_matrix_final_plot) = expression_matrix_final_plot$rownames

expression_matrix_final_plot_1 = expression_matrix_final_plot[2:29] 

#Figure 3-4B
pheatmap(expression_matrix_final_plot_1,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 14)

################################################################################
#DICE - Immune cell analysis
################################################################################
#access data using celldex
DICE_data = DatabaseImmuneCellExpressionData(cell.ont = "all")

#convert to SingleCellExperiment format
DICE_sce = SingleCellExperiment(
  assays = assays(DICE_data),
  rowData = rowData(DICE_data),
  colData = colData(DICE_data),
  metadata = metadata(DICE_data)
)

chromatin_remodelling_genes = c("KAT8","MSL1","MSL2","MSL3","KANSL1","KANSL2","KANSL3",
                                "RBBP5","PHF20","WDR5","RUVBL2","KAT5","H4C6","MCRS1",
                                "EP400","MORF4L1","H3C12","TRRAP","YEATS4","UCHL5",
                                "RUVBL1","ACTR5","INO80B","INO80C","TFPT","NFRKB",
                                "INO80","INO80E","PSMD3","ADRM1")

#set gene set, background set and subsetting variable
gene_set_immune = chromatin_remodelling_genes
background_genes_DICE = rownames(DICE_sce)  
cell_types_DICE = unique(colData(DICE_sce)$label.fine)

#expression quantile cutoffs
expression_cutoffs = c(0.5,0.6,0.7,0.8,0.9)
n_permutations = 1000
results_DICE = list()

for (q in expression_cutoffs) {
  #for each cutoff quantile
  
  for (ct in cell_types_DICE) {
    #subset by cell type
    specific_cell_type = colnames(DICE_sce)[colData(DICE_sce)$label.fine == ct]
    if (length(specific_cell_type) == 0) next
    
    cell_expression_matrix = assay(DICE_sce, "logcounts")[, specific_cell_type]
    
    #get average gene expression over cell type
    global_cell_expression = rowMeans(cell_expression_matrix)
    
    #set expression threshold based on average expression distribution
    threshold = quantile(global_cell_expression, probs = q)
    
    #keep genes above this threshold
    expressed_genes = names(global_cell_expression[global_cell_expression > threshold])
    
    #get number of genes within gene set that overlaps at each threshold
    overlap = length(intersect(gene_set_immune, expressed_genes))
    
    if (overlap == 0) next
    
    #calculate hpyergeometric distribution p value
    pval = phyper(
      q = overlap - 1,
      m = length(expressed_genes),
      n = length(background_genes_DICE) - length(expressed_genes),
      k = length(gene_set_immune),
      lower.tail = FALSE
    )
    
    #permute p value for FWER
    permuted_pvals = replicate(n_permutations, {
      permuted_set = sample(background_genes_DICE, length(gene_set_immune))
      perm_overlap = length(intersect(permuted_set, expressed_genes))
      phyper(
        q = perm_overlap - 1,
        m = length(expressed_genes),
        n = length(background_genes_DICE) - length(expressed_genes),
        k = length(gene_set_immune),
        lower.tail = FALSE
      )
    })
    
    fwer = mean(permuted_pvals <= pval)
    
    #store results
    results_DICE[[paste0(ct, "_", q)]] = data.frame(
      cell_type = ct,
      cutoff = q,
      threshold = threshold,
      p_value = pval,
      fwer = fwer,
      overlap = overlap
    )
  }
}


# Combine all results
results_DICE = do.call(rbind, results_DICE)
# Multiple testing correction
results_DICE$adj_fwer = p.adjust(results_DICE$fwer, method = "BH")

#rearrange for highest cutoff first
results_DICE  = results_DICE  %>%
  arrange(-cutoff)

#write.csv(results_DICE, "enrichment_results_DICE.csv")
#results_DICE = read.csv("enrichment_results_DICE.csv")

#filter to keep only highest cutoff
full_results_DICE_filtered = results_DICE %>%
  filter(cutoff == "0.9")

palette = c(
  "#B32357",  '#35A7FF', "#D55E00",  "#6ACC65", "#911eb4", "#F781BF", "#009E73", "#FF6F61",
  "#0072B2",  "#b66dff", "#FFC000", "#117733",    "#00CED1", "#808000",  
  "#F0E442",  "#FF8C00",   "#D33682", 
  "#4682B4",  "#000000")

#Figure 3-6A
ggplot(full_results_DICE_filtered, aes(y = cell_type, x = adj_fwer)) + 
  geom_point(size = 3,colour = "#B32357")+
  geom_vline(xintercept = 0.05, lty = 2)+
  theme_bw(base_size = 16)+
  labs(y = "Cell type", x = "FWER")


################################################################################
#DICE gene expression Figure 3-6B
################################################################################

#genes in gene set present in DICE dataset
DICE_PPIN_overlap= intersect(chromatin_remodelling_genes, rownames(DICE_sce))
#filter dataset to only include these genes
DICE_overlap_sce = DICE_sce[DICE_PPIN_overlap,]

cell_types_DICE = unique(colData(DICE_sce)$label.fine)

expression_matrix_DICE = list()

for (ct in cell_types_DICE) {
      #select cell type within tissue
      specific_cell_type = colnames(DICE_overlap_sce)[colData(DICE_overlap_sce)$label.fine == ct]
      
      #get RNA seq matrix data for cell type
      cell_expression_matrix = assay(DICE_overlap_sce, "logcounts")[, specific_cell_type]
      
      #get mean expression for genes within cell type
      global_cell_expression = rowMeans(cell_expression_matrix)
      
      #store results
      expression_matrix_DICE[[paste0(ct)]] = data.frame(
        cell_type = ct,
        gene = names(global_cell_expression),
        expression = global_cell_expression)
    }


expression_matrix_DICE = do.call(rbind, expression_matrix_DICE)

#rearrange data
expression_matrix_DICE_summary = expression_matrix_DICE %>%
  group_by(cell_type, gene) %>%
  summarise(av_expression = mean(expression)) %>%
  mutate(raw_expression = 2^av_expression,
         raw_expression = ifelse(raw_expression> 100, 100, raw_expression)) %>%
  ungroup()
  
expression_matrix_DICE_final_plot = expression_matrix_DICE_summary %>%
  dplyr::select(!av_expression) %>%
  pivot_wider(names_from = "gene", values_from = "raw_expression") %>%
  as.data.frame()

rownames(expression_matrix_DICE_final_plot) = expression_matrix_DICE_final_plot$cell_type

expression_matrix_DICE_final_plot = expression_matrix_DICE_final_plot[2:29] 

#Figure 3-6B
pheatmap(expression_matrix_DICE_final_plot,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE,
         fontsize = 14)


################################################################################
#T2D analysis
################################################################################

#access data using scRNA-seq
pancreas_rpkm_data = fetchDataset("xin-pancreas-2016", "2023-12-19")

#get metadata
pancreas_coldata = colData(pancreas_rpkm_data)
pancreas_rowdata = rowData(pancreas_rpkm_data)

#change rownames to gene names 
rownames(pancreas_rpkm_data) = pancreas_rowdata$symbol

#filter data to only healthy individuals for consistency with previous analysis
pancreas_rpkm_data_healthy = pancreas_rpkm_data[, colData(pancreas_rpkm_data)$condition == "Healthy"]

#gene sets 
BIN1_PPIN = c("BIN1", "MYC", "FBXW7", "SYNJ1", "ITSN1", "CDC42", "TNK2",
              "DNM2", "HCLS1", "CTTN", "SRC", "BCAR1", "RIN3", "RAB5B",
              "PIK3C3", "BECN1", "PRKN", "SNCA", "SHGL2", "TRRAP", "BRCA1",
              "HIF1A", "EP300", "CBL", "RAPGEF1", "CD2AP", "SH3GL2", "ARHGAP1",
              "UBE2L3")

KAT8_PPIN = c("KAT8","MSL1","MSL2","MSL3","KANSL1","KANSL2","KANSL3",
              "RBBP5","PHF20","WDR5","RUVBL2","KAT5","H4C6","MCRS1",
              "EP400","MORF4L1","H3C12","H3C13","CBX3", "EHMT2", "CHD3",
              "TRRAP","YEATS4","INO80B","SETD1A","INO80","INO80E", "RBBP4")

PLEKHM1_PPIN = c("PLEKHM1", "MAP1LC3C", "GABARAP", "CALCOCO2", 
                 "RB1CC1", "ATG16L1",
                 "ATG7", "ATG3", "ATG4B", "OPTN", "BNIP3", "BNIP3L",
                 "RHEB", "GABARAPL2",
                 "TBK1", "BCL2L1", "BCL2", "BECN1",
                 "NBR1", "ULK2", "RPTOR", "IKBKE")

T2D_gene_sets = list(BIN1 = BIN1_PPIN,
                     KAT8 = KAT8_PPIN,
                     PLEKHM1 = PLEKHM1_PPIN)

 
cell_types_pancreas = unique(colData(pancreas_rpkm_data_healthy)$'cell.type')

expression_cutoffs = c(0.5, 0.6, 0.7,0.8, 0.9)
results_pancreas = list()
n_permutations = 1000

for (g in names(T2D_gene_sets)){
  #subset by PPIN gene set
  gene_set = T2D_gene_sets[[g]]
  
  for (q in expression_cutoffs) {
    #for each cutoff quantile:
    
    for (ct in cell_types_pancreas) {
      #subset by cell type
      specific_cell_type = colnames(pancreas_rpkm_data_healthy)[colData(pancreas_rpkm_data_healthy)$cell.type == ct] 
      
      if (length(specific_cell_type) == 0) next
      
      cell_expression_matrix = assay(pancreas_rpkm_data_healthy, "rpkm")[, specific_cell_type]
      
      background_genes_pancreas = rownames(cell_expression_matrix) 
      
      #get average expression over cell
      global_cell_expression = rowMeans(cell_expression_matrix)
      
      #get expression threshold based on average expression distribution
      threshold = quantile(global_cell_expression, probs = q) 
      
      #find genes that are above this threshold
      expressed_genes = names(global_cell_expression[global_cell_expression > threshold])
      
      #number of genes in gene set above threshold
      overlap = length(intersect(gene_set, expressed_genes))
      
      if (overlap == 0) next
      
      #calculate hypergeometric distribution p value
      pval = phyper(
        q = overlap - 1,
        m = length(expressed_genes),
        n = length(background_genes_pancreas) - length(expressed_genes),
        k = length(gene_set),
        lower.tail = FALSE
      )
      
      #permute p value to calculate FWER
      permuted_pvals = replicate(n_permutations, {
        permuted_set = sample(background_genes_pancreas, length(gene_set))
        perm_overlap = length(intersect(permuted_set, expressed_genes))
        phyper(
          q = perm_overlap - 1,
          m = length(expressed_genes),
          n = length(background_genes_pancreas) - length(expressed_genes),
          k = length(gene_set),
          lower.tail = FALSE
        )
      })
      
      fwer = mean(permuted_pvals <= pval)
      
      #store results
      results_pancreas[[paste0(g, "_", ct, "_", q)]] = data.frame(
        PPIN = g,
        cell_type = ct,
        cutoff = q,
        threshold = threshold,
        pval = pval,
        fwer = fwer,
        overlap = length(intersect(gene_set, expressed_genes)))
    }
  }
}


# Combine all results
results_pancreas = do.call(rbind, results_pancreas)

# Multiple testing correction
results_pancreas$adj_fwer = p.adjust(results_pancreas$fwer, method = "BH")

#remove contaminated cells 
results_pancreas_final = results_pancreas %>%
  filter(!grepl("contaminated", cell_type)) 

results_pancreas_final = results_pancreas_final %>%
  arrange(-cutoff)


#write.csv(results_pancreas_final, "enrichment_results_pancreas_final.csv")
#results_pancreas = read.csv("enrichment_results_pancreas_final.csv")

#subset to keep highest cutoff
results_pancreas_subset = results_pancreas_final %>%
  filter(cutoff == "0.9")

#rename cell type
results_pancreas_subset$cell_name = factor(results_pancreas_subset$cell_type) %>%
  fct_recode("Pancreatic Polypeptide cells" = "PP",
             "Alpha cells" = "alpha", 
             "Beta cells" = "beta", 
             "Delta cells" = "delta")

palette = c(
  "#B32357",  '#35A7FF', "#D55E00",  "#6ACC65", "#911eb4", "#F781BF", "#009E73", "#FF6F61",
  "#0072B2",  "#b66dff", "#FFC000", "#117733",    "#00CED1", "#808000",  
  "#F0E442",  "#FF8C00",   "#D33682", 
  "#4682B4",  "#000000")

#Figure 3-11A
ggplot(results_pancreas_subset, aes(y = cell_name, x = adj_fwer, colour = PPIN)) + 
  geom_point(size = 3)+
  geom_vline(xintercept = 0.05, lty = 2)+
  theme_bw(base_size = 16)+
  labs(colour = "PPIN gene set", y = "Cell type", x = "FWER")+
  scale_color_manual(values  = palette) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1.0, 0.05),
                     limits = c(0,0.25))

################################################################################
#Pancreas gene expression figure 3-11B to 3-11D
################################################################################

#get overlapping gene set
pancreas_BIN1_overlap = intersect(BIN1_PPIN, rownames(pancreas_rpkm_data_healthy))

#subset data by overlapping gene set
pancreas_BIN1 = pancreas_rpkm_data_healthy[pancreas_BIN1_overlap,]

BIN1_pancreas_expression_results = list()

#loop to get expression dta from each cell type 
    for (ct in cell_types_pancreas) {
      #subset by cell type
      specific_cell_type = colnames(pancreas_BIN1)[colData(pancreas_BIN1)$cell.type == ct] 
      
      if (length(specific_cell_type) == 0) next
      
      cell_expression_matrix = assay(pancreas_BIN1, "rpkm")[, specific_cell_type]
        
        #get mean expression for genes within cell type
        global_cell_expression = rowMeans(cell_expression_matrix) 
         
        #save results
        BIN1_pancreas_expression_results[[paste0(ct)]] <- data.frame(
          cell_type = ct,
          gene = names(global_cell_expression),
          ave_expression = global_cell_expression) %>%
          mutate(ave_expression = ifelse(ave_expression > 100, 100, ave_expression))%>%
          pivot_wider(names_from = "gene", values_from = "ave_expression")
      }


# Combine all results
BIN1_pancreas_expression_results = do.call(rbind, BIN1_pancreas_expression_results)

#remove contaminated cells 
BIN1_pancreas_expression_results = BIN1_pancreas_expression_results %>%
  filter(!grepl("contaminated", cell_type)) %>%
  as.data.frame() %>%
  arrange(cell_type)

#rename cell type
BIN1_pancreas_expression_results$cell_type = factor(BIN1_pancreas_expression_results$cell_type) %>%
  fct_recode("Pancreatic Polypeptide cells" = "PP",
             "Alpha cells" = "alpha", 
             "Beta cells" = "beta", 
             "Delta cells" = "delta")

rownames(BIN1_pancreas_expression_results) = BIN1_pancreas_expression_results$cell_type

BIN1_pancreas_expression_results = BIN1_pancreas_expression_results %>%
  arrange(cell_type) %>%
  select(!cell_type)

#Figure 3-11B
pheatmap(BIN1_pancreas_expression_results,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 14)




################################################################################
pancreas_KAT8_overlap = intersect(KAT8_PPIN, rownames(pancreas_rpkm_data_healthy))

pancreas_KAT8 = pancreas_rpkm_data_healthy[pancreas_KAT8_overlap,]

KAT8_pancreas_expression_results = list()

for (ct in cell_types_pancreas) {
  #subset by cell type
  specific_cell_type = colnames(pancreas_KAT8)[colData(pancreas_KAT8)$cell.type == ct] 
  
  if (length(specific_cell_type) == 0) next
  
  cell_expression_matrix = assay(pancreas_KAT8, "rpkm")[, specific_cell_type]
  
  #get mean expression for genes within cell type
  global_cell_expression = rowMeans(cell_expression_matrix) 
  
  #save results
  KAT8_pancreas_expression_results[[paste0(ct)]] <- data.frame(
    cell_type = ct,
    gene = names(global_cell_expression),
    ave_expression = global_cell_expression) %>%
    mutate(ave_expression = ifelse(ave_expression > 100, 100, ave_expression))%>%
    pivot_wider(names_from = "gene", values_from = "ave_expression")
}


# Combine all results
KAT8_pancreas_expression_results = do.call(rbind, KAT8_pancreas_expression_results)

#remove contaminated cells 
KAT8_pancreas_expression_results = KAT8_pancreas_expression_results %>%
  filter(!grepl("contaminated", cell_type)) %>%
  as.data.frame() %>%
  arrange(cell_type)

#rename cell type
KAT8_pancreas_expression_results$cell_type = factor(KAT8_pancreas_expression_results$cell_type) %>%
  fct_recode("Pancreatic Polypeptide cells" = "PP",
             "Alpha cells" = "alpha", 
             "Beta cells" = "beta", 
             "Delta cells" = "delta")

rownames(KAT8_pancreas_expression_results) = KAT8_pancreas_expression_results$cell_type

KAT8_pancreas_expression_results = KAT8_pancreas_expression_results %>%
  arrange(cell_type) %>%
  select(!cell_type)

#Figure 3-11C
pheatmap(KAT8_pancreas_expression_results,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 14)

################################################################################
pancreas_PLEKHM1_overlap = intersect(PLEKHM1_PPIN, rownames(pancreas_rpkm_data_healthy))

pancreas_PLEKHM1 = pancreas_rpkm_data_healthy[pancreas_PLEKHM1_overlap,]

PLEKHM1_pancreas_expression_results = list()

for (ct in cell_types_pancreas) {
  #subset by cell type
  specific_cell_type = colnames(pancreas_PLEKHM1)[colData(pancreas_PLEKHM1)$cell.type == ct] 
  
  if (length(specific_cell_type) == 0) next
  
  cell_expression_matrix = assay(pancreas_PLEKHM1, "rpkm")[, specific_cell_type]
  
  #get mean expression for genes within cell type
  global_cell_expression = rowMeans(cell_expression_matrix) 
  
  #save results
  PLEKHM1_pancreas_expression_results[[paste0(ct)]] <- data.frame(
    cell_type = ct,
    gene = names(global_cell_expression),
    ave_expression = global_cell_expression) %>%
    mutate(ave_expression = ifelse(ave_expression > 100, 100, ave_expression))%>%
    pivot_wider(names_from = "gene", values_from = "ave_expression")
}


# Combine all results
PLEKHM1_pancreas_expression_results = do.call(rbind, PLEKHM1_pancreas_expression_results)

#remove contaminated cells 
PLEKHM1_pancreas_expression_results = PLEKHM1_pancreas_expression_results %>%
  filter(!grepl("contaminated", cell_type)) %>%
  as.data.frame() %>%
  arrange(cell_type)

#rename cell type
PLEKHM1_pancreas_expression_results$cell_type = factor(PLEKHM1_pancreas_expression_results$cell_type) %>%
  fct_recode("Pancreatic Polypeptide cells" = "PP",
             "Alpha cells" = "alpha", 
             "Beta cells" = "beta", 
             "Delta cells" = "delta")

rownames(PLEKHM1_pancreas_expression_results) = PLEKHM1_pancreas_expression_results$cell_type

PLEKHM1_pancreas_expression_results = PLEKHM1_pancreas_expression_results %>%
  arrange(cell_type) %>%
  select(!cell_type)

#Figure 3-11D
pheatmap(PLEKHM1_pancreas_expression_results,
         color=rev(hcl.colors(80,"Reds 3")),
         cluster_rows = FALSE,
         cluter_columns = TRUE, 
         fontsize = 14)
