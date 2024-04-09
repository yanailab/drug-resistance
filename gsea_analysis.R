library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)
library(pals)
library(fgsea)

setwd("/path/kura_analysis")

load("kura_all_combined.rda")
load("gene_sets.rda")


diff_expression_gsea <- function(seurat_obj,
                                 cluster_column,
                                 ctrl = F,
                                 ctrl_id = c("Sample_id"), rest_no_ctrl = T) {
  # calculate differential expression for cells that define modules 
  # and don't filter for pvalue or anything to be able to rank the genes for
  # GSEA analysis. Returns a list with dataframes of diff. expression
  
  # get clusters with signature names
  clusters = sort(as.vector(unique(seurat_obj@meta.data[, cluster_column])))
  # numbers of modules present
  modules = unique(clusters)
  # list to keep the output
  list_diff_exp = c()
  for(module in modules) {
    # get names of the signatures for the module
    # cells in module or cluster
    cells_module = rownames(subset(seurat_obj@meta.data, RNA_snn_res.0.3 == module))
    #genes_module = genes_module[genes_module > 0.1]
    # keep log_ratios for all signatures of a module
    
    # this compares all module cells against ctrl (of course exluding ctrl from module)  
    if((ctrl == T) & (rest_no_ctrl == F)) {
      # take ctrl cells
      rest_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident %in% ctrl_id, ])
      # subset cells module that are not ctrl
      cells_module = subset(cells_module, !(cells_module %in% rest_cells))
      # this compares module cells against all others, excluding ctrls  
    } else if((ctrl == T) & (rest_no_ctrl == T)) {
      # take ctrl cells
      ctrl_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident %in% ctrl_id, ])
      # subset cells module that are not ctrl
      all_cells = colnames(seurat_obj@data)
      all_cells = subset(all_cells, !(all_cells %in% ctrl_cells)) # exclude Ctrl
      cells_module = subset(cells_module, !(cells_module %in% ctrl_cells)) # exclude Ctrl 
      # rest cells have no ctrl and no cells in module
      rest_cells = subset(all_cells, !(all_cells %in% cells_module))
      # this compares all module cells against all other cells
    } else {
      all_cells = rownames(seurat_obj@meta.data)
      #print(head(all_cells))
      #print("-----")
      # rest cells have no cells in module
      rest_cells = subset(all_cells, !(all_cells %in% cells_module))
      #print(head(cells_module))
      #print(head(rest_cells))
    }
    # test diff expression in cells_module vs rest: The default comparison
    markers = FindMarkers(seurat_obj, ident.1 = cells_module, ident.2 = rest_cells,
                          min.cells.feature = 10,
                          test.use = "MAST", only.pos = F, min.pct = -Inf,
                          logfc.threshold = -Inf, min.diff.pct = -Inf)
    list_diff_exp[[module]] <- markers
  }
  return(list_diff_exp)
}


list_diff_exp_all = diff_expression_gsea(kura_all_combined)
save(list_diff_exp_all, 
     file = "list_diff_exp_all_for_gsea.rda")

load("list_diff_exp_all_for_gsea.rda")


#########################################
# Run fgsea on each of the module genes #
#########################################

library(fgsea)

rank_genes <- function(markers) {
  # remove crap
  markers = markers[!is.na(markers$p_val), ]
  # get log10 and check for inf values
  markers_log_pval = -log10(markers$p_val)
  n_inf = length(markers_log_pval[is.infinite(markers_log_pval)])
  # -log10 pval goes only until 320, otherwise it becomes infinite, so add 1
  # to each of them in the order to make them all different. GSEA does not like
  # equal rankings. If there are infs, fix them.
  if(n_inf > 0) {
    new_inf = c()
    for(i in n_inf:1) { 
      new_inf = c(new_inf, 320 + i) 
    }
    markers_log_pval[1:n_inf] <- new_inf
  }
  markers_rank = sign(markers$avg_log2FC) * (markers_log_pval)
  names(markers_rank) <- rownames(markers)
  markers_rank = sort(markers_rank, decreasing = T)
  return(markers_rank)
}


run_fgsea <- function(list_diff_exp, gene_sets) {
  # run fgsea for each of the dataframes contained in list_diff_exp and
  # store in a list
  
  list_fgsea = list()
  # use the counter to catch names for list output
  for(i in 1:length(list_diff_exp)) {
    ranked_genes = rank_genes(list_diff_exp[[i]])
    fgseaRes <- fgsea(pathways = gene_sets, 
                      stats    = ranked_genes,
                      minSize  = 10,
                      nperm = 10000,
                      maxSize  = 500)
    list_fgsea[[i]] = fgseaRes
  }
  names(list_fgsea) <- names(list_diff_exp)
  return(list_fgsea)
}


fgsea_modules_pdx_per = run_fgsea(list_diff_exp_all_pdx_per, gene_sets_hall_kegg)


###########################
# Creating tables to plot #
###########################

pathways_to_select = c("HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                       "REACTOME_KERATINIZATION", "KEGG_ADHERENS_JUNCTION",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                       "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_TRANSLATION", 
                       "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                       "KEGG_FERROPTOSIS", "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
                       "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                       "REACTOME_METABOLISM_OF_CARBOHYDRATES", 
                       "REACTOME_METABOLISM_OF_NUCLEOTIDES"
)

# store NES scores
nes_data = data.frame()
for(p in pathways_to_select) {
  nes_scores = (sapply(fgsea_modules_no_ctrl, function(x) x[pathway == p, ][, c(NES)]))
  nes_data = rbind(nes_data, nes_scores)
}
colnames(nes_data) <- names(fgsea_modules_no_ctrl)
rownames(nes_data) <- pathways_to_select
# reorder modules
nes_data = nes_data[, c(3,1,2,4,5,6)]
nes_data$pathway <- rownames(nes_data)
# long format for ggplot
nes_data_long = melt(nes_data, id.vars = c("pathway"), value.name = "NES",
                     variable.name = "Module")

# store pvals scores
pval_data = data.frame()
for(p in pathways_to_select) {
  pvals = (sapply(fgsea_modules_no_ctrl, function(x) x[pathway == p, ][, c(pval)]))
  pval_data = rbind(pval_data, pvals)
}
colnames(pval_data) <- names(fgsea_modules_no_ctrl)
rownames(pval_data) <- pathways_to_select
# reorder modules
pval_data = pval_data[, c(3,1,2,4,5,6)]
pval_data$pathway <- rownames(pval_data)
# long format for ggplot
pval_data_long = melt(pval_data, id.vars = c("pathway"), value.name = "pval",
                      variable.name = "Module")
# merging both data (adding the column)
nes_pval_data_long = nes_data_long
nes_pval_data_long$pval <- pval_data_long$pval
# reorder y axis
nes_pval_data_long$pathway <- factor(nes_pval_data_long$pathway, 
                                     levels = rev(pathways_to_select))

# Plot pathways as bubble plot
ggplot(nes_pval_data_long, aes(x = Module, y = pathway, size = -log10(pval), fill = NES)) + 
  geom_point(shape = 21) +
  scale_size_continuous(range = c(1, 8), breaks = c(0, 1, 2, 3)) +
  scale_fill_gradientn(colours = pals::coolwarm(100)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10))