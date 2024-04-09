library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(infercnv)

merge_exp <- function(exp_matrix_list) {
  # Merge multiple expression matrices by row names (by=0)
  # exp_matrix_list: List of expression matrices dataframes.
  # Returns an expression matrix with all cells (columns) merged. It is
  # important that they have different column names, like _"sampleID"
  # function calling example: merge_exp(list(a, b, c))
  
  merged_data = Reduce(function(x, y) transform(merge(x, y, by=0, all=TRUE), 
                                                row.names=Row.names, Row.names=NULL),
                       exp_matrix_list)
  # Remove eventual NAs and replace by 0
  merged_data[is.na(merged_data)] <- 0
  return(merged_data)
}


setwd("/path/inferCNV_kura")

# annotation state for all
cell_ann_all = read.table("cell_ann_all.tsv")
cell_ann_t1 = read.table("cell_ann_t1.tsv")
cell_ann_t2.5 = read.table("cell_ann_t2.5.tsv")
cell_ann_t5 = read.table("cell_ann_t5.tsv")
cell_ann_t10 = read.table("cell_ann_t10.tsv")
cell_ann_t20 = read.table("cell_ann_t20.tsv")
cell_ann_t40 = read.table("cell_ann_t40.tsv")
cell_ann_t80 = read.table("cell_ann_t80.tsv")
cell_ann_t160 = read.table("cell_ann_t160.tsv")
cell_ann_t320 = read.table("cell_ann_t320.tsv")

cell_ann_list = list(cell_ann_t1, cell_ann_t2.5, cell_ann_t5, cell_ann_t10,
                     cell_ann_t20, cell_ann_t40, cell_ann_t80, cell_ann_t160,
                     cell_ann_t320)
# remove C from annotation
cell_ann_list = lapply(cell_ann_list, function(x) x = x[x$V2 != "C", ])
cell_ann_all = cell_ann_all[cell_ann_all$V2 != "C", ]

# get gene order
gene_ordering = read.table("gene_ordering_file.txt")
# read Exome shared genes to plot the exact same genes from exome analysis
#genes_cnv_shared_exome = as.character(read.table("genes_cnv_exome_shared.txt")$V1)

# read matrices
infercnv_t1 = readRDS("t1/run.final.infercnv_obj")
infercnv_t2.5 = readRDS("t2_5/run.final.infercnv_obj")
infercnv_t5 = readRDS("t5/run.final.infercnv_obj")
infercnv_t10 = readRDS("t10/run.final.infercnv_obj")
infercnv_t20 = readRDS("t20/run.final.infercnv_obj")
infercnv_t40 = readRDS("t40/run.final.infercnv_obj")
infercnv_t80 = readRDS("t80/run.final.infercnv_obj")
infercnv_t160 = readRDS("t160/run.final.infercnv_obj")
infercnv_t320 = readRDS("t320/run.final.infercnv_obj")


infercnv_obj_list = list(infercnv_t1@expr.data, infercnv_t2.5@expr.data, infercnv_t5@expr.data,
                         infercnv_t10@expr.data, infercnv_t20@expr.data, infercnv_t40@expr.data,
                         infercnv_t80@expr.data, infercnv_t160@expr.data, infercnv_t320@expr.data)
# shared genes in all
all_genes = lapply(infercnv_obj_list, rownames)
genes_shared = Reduce(intersect, all_genes)

# subset infercnv with only treated cells and common genes in all
for(i in 1:length(infercnv_obj_list)){
  cells = as.character(cell_ann_list[[i]]$V1)
  infercnv_obj_list[[i]] <- infercnv_obj_list[[i]][genes_shared, cells]
}


##################################################
# Trying sampling cells and plot as single-cells #
##################################################

# sample X % of cells in each group
sample_cells <- function(exp_cnv_matrix, fraction = 0.3) {
  # sample cells to plot cnvs
  n_cells = round(length(colnames(exp_cnv_matrix)) * fraction)
  sample_cells = sample(colnames(exp_cnv_matrix), n_cells)
  return(exp_cnv_matrix[, sample_cells])
}


prepare_data_plot <- function(infercnv_obj_list_sample_ind, cell_ann_list_ind,
                              genes_shared, color_cluster_ind, color_sample_ind) {
  
  # prepare data for individual samples to plot, the thing ready to plot?
  infercnv_obj_list_sample = lapply(list(infercnv_obj_list_sample_ind),
                                    sample_cells, fraction = 0.3)
  
  # merge list of dataframes into one single dataframe
  infercnv_obj_sample_merged = merge_exp(infercnv_obj_list_sample)
  infercnv_obj_sample_merged = infercnv_obj_sample_merged[genes_shared, ]
  
  cell_ann_all = cell_ann_list_ind
  rownames(cell_ann_all) <- as.character(cell_ann_all$V1)
  print(length(cell_ann_all))
  cell_ann_states = cell_ann_all[colnames(infercnv_obj_sample_merged), ]
  print(dim(infercnv_obj_sample_merged))
  sample_state_ann = stringr::str_split(as.character(cell_ann_states$V2), pattern = "_")
  # extracting samples
  sample_ann = factor(sapply(sample_state_ann, "[[", 1))
  # extracting clusters
  state_ann = factor(sapply(sample_state_ann, "[[", 2))
  
  # chr annotation
  rownames(gene_ordering) <- gene_ordering$V1
  gene_order_shared = gene_ordering[genes_shared, ]
  chr_ann = factor(gene_order_shared$V2, levels = c(paste("chr",1:22, sep = "")))
  # add a number in the end to make unique names
  rand_number = c(1:length(colnames(infercnv_obj_sample_merged)))
  # create new cell names that will be sorted after dendrogram! this helps to keep the
  # sample and state order according to the time points.
  new_col_names_for_sort = paste(as.character(sample_ann), as.character(state_ann), 
                                 1:length(colnames(infercnv_obj_sample_merged)), sep = "_")
  # add the new cell names
  colnames(infercnv_obj_sample_merged) <- new_col_names_for_sort
  # sample state annotation
  sample_state_ann_heat = rowAnnotation(df = data.frame(sample = sample_ann, state = state_ann), 
                                        simple_anno_size = unit(2, "mm"),
                                        col = list(state = color_cluster_ind,
                                                   sample = color_sample_ind),
                                        show_annotation_name = F,
                                        annotation_legend_param = list(title_gp = gpar(fontsize = 6),
                                                                       title = "State",
                                                                       labels_gp = gpar(fontsize = 6)))
  # plot heatmap
  col_fun = colorRamp2(c(0.8, 1, 1.2), c("#000A87", "#ffffff", "#900512")) # this works well
  
  library(dendextend)
  dend = as.dendrogram(hclust(dist(t(infercnv_obj_sample_merged))), type = "average") # can be used as cluster_rows and columns arg instead of T
  # this does the nice trick to sort cells based on their names that I redefined by sample and state!
  # this orders strings based on numbers! and keeps the clustering constraints
  dend <- rotate(dend, stringr::str_sort(labels(dend), numeric = T))
  
  heatmap_cnv = Heatmap(t(infercnv_obj_sample_merged), cluster_rows = dend, row_dend_reorder = F,
                        cluster_columns = F,
                        show_row_names = F, show_column_names = F, show_row_dend = F,
                        use_raster = T, left_annotation = sample_state_ann_heat,
                        col = col_fun, column_split = chr_ann, column_gap = unit(1, "mm"),
                        border = TRUE, column_title_side = "bottom", 
                        column_title_gp = gpar(fontsize = 6),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                    title = "Expression",
                                                    labels_gp = gpar(fontsize = 6),
                                                    legend_height = unit(2, "cm")))
  return(heatmap_cnv)  
}


cluster_colors = list(c("0" = "#caafd4", "1" = "#af88bf", "2" = "#cb71ef", "3" = "#d261ff"),
                      c("0" = "#7c5c8a", "1" = "#824D99", "2" = "#8a39ac", "3" = "#9717cf"),
                      c("0" = "#4E78C4", "1" = "#2c6ee8"),
                      c("0" = "#57A2AC", "1" = "#84bbc2", "2" = "#37bccd", "3" = "#12d8f3"),
                      c("0" = "#87ab82", "1" = "#70cb62", "2" = "#61e04d"),
                      c("0" = "#b9a75b", "1" = "#d0b643", "2" = "#e8c52c"),
                      c("0" = "#c68353", "1" = "#dd7f33", "2" = "#ff791a"),
                      c("0" = "#b43e3c", "1" = "#e3514f", "2" = "#e4100c"),
                      c("0" = "#421510", "1" = "#74241b", "2" = "#8e1c10"))


sample_colors = list(c("T1" = "#B997C7"), 
                     c("T2.5" = "#824D99"),
                     c("T5" = "#4E78C4"),
                     c("T10" = "#57A2AC"),
                     c("T20" = "#7EB875"),
                     c("T40" = "#D0B541"),
                     c("T80" = "#E67F33"),
                     c("T160" = "#CE2220"), 
                     c("T320" = "#80271E"))

#infercnv_obj_list_sample = lapply(infercnv_obj_list, sample_cells, fraction = 0.3)


list_heatmaps = list()

for(i in 1:length(infercnv_obj_list_sample)) {
  heat_cnv = prepare_data_plot(infercnv_obj_list[[i]], 
                               cell_ann_list[[i]],
                               genes_shared, 
                               cluster_colors[[i]],
                               sample_colors[[i]]
  )
  list_heatmaps[[i]] <- heat_cnv
}


draw(list_heatmaps[[4]])

#########################################
# plot as infercnv subcluster structure #
#########################################

cluster_colors = list(c("0" = "#caafd4", "1" = "#af88bf", "2" = "#cb71ef", "3" = "#d261ff"),
                      c("0" = "#7c5c8a", "1" = "#824D99", "2" = "#8a39ac", "3" = "#9717cf"),
                      c("0" = "#4E78C4", "1" = "#2c6ee8"),
                      c("0" = "#57A2AC", "1" = "#84bbc2", "2" = "#37bccd", "3" = "#12d8f3"),
                      c("0" = "#87ab82", "1" = "#70cb62", "2" = "#61e04d"),
                      c("0" = "#b9a75b", "1" = "#d0b643", "2" = "#e8c52c"),
                      c("0" = "#c68353", "1" = "#dd7f33", "2" = "#ff791a"),
                      c("0" = "#b43e3c", "1" = "#e3514f", "2" = "#e4100c"),
                      c("0" = "#421510", "1" = "#74241b", "2" = "#8e1c10"))


sample_colors = list(c("T1" = "#B997C7"), 
                     c("T2.5" = "#824D99"),
                     c("T5" = "#4E78C4"),
                     c("T10" = "#57A2AC"),
                     c("T20" = "#7EB875"),
                     c("T40" = "#D0B541"),
                     c("T80" = "#E67F33"),
                     c("T160" = "#CE2220"), 
                     c("T320" = "#80271E"))

infercnv_list = list(infercnv_t1, infercnv_t2.5, infercnv_t5,
                     infercnv_t10, infercnv_t20, infercnv_t40,
                     infercnv_t80, infercnv_t160, infercnv_t320)

list_heatmaps = list()
for(i in 1:length(infercnv_list)) {
  # excludes ctrl cells
  cells_sample_order = c()
  state_ann = c()
  for(j in 1:(length(infercnv_list[[i]]@tumor_subclusters$hc) - 1)) {
    cluster_id = names(infercnv_list[[i]]@tumor_subclusters$hc[j])
    cluster_id = stringr::str_split(cluster_id, pattern = "_")
    # get cluster id
    cluster_id = cluster_id[[1]][2]
    #sample_id = cluster_id[[1]][1]
    #print(cluster_id)
    # get cell order
    cells_ind = infercnv_list[[i]]@tumor_subclusters$hc[[j]]$labels
    # subset 30% of cells
    indices_rand = sample(1:(abs(length(cells_ind)*0.3)))
    print(length(indices_rand))
    cells_ind = cells_ind[indices_rand]
    
    state_ann_ind = rep(cluster_id, length(cells_ind))
    cells_sample_order = c(cells_sample_order, cells_ind)
    state_ann = c(state_ann, state_ann_ind)
  }
  
  infer_cnv_sub = infercnv_list[[i]]@expr.data[, cells_sample_order]
  rownames(gene_ordering) <- gene_ordering$V1
  gene_order_shared = gene_ordering[rownames(infer_cnv_sub), ]
  chr_ann = factor(gene_order_shared$V2, levels = c(paste("chr",1:22, sep = "")))
  
  sample_state_ann_heat = rowAnnotation(df = data.frame(state = state_ann), 
                                        simple_anno_size = unit(2, "mm"),
                                        col = list(state = cluster_colors[[i]]),
                                        show_annotation_name = F,
                                        annotation_legend_param = list(title_gp = gpar(fontsize = 6),
                                                                       title = "Cluster",
                                                                       labels_gp = gpar(fontsize = 6)))
  
  heatmap_cnv = Heatmap(t(infer_cnv_sub), cluster_rows = F, row_dend_reorder = F,
                        cluster_columns = F,
                        show_row_names = F, show_column_names = F, show_row_dend = F,
                        use_raster = T, left_annotation = sample_state_ann_heat,
                        col = col_fun, column_split = chr_ann, column_gap = unit(1, "mm"),
                        border = TRUE, column_title_side = "bottom", 
                        column_title_gp = gpar(fontsize = 6),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                    title = "Expression",
                                                    labels_gp = gpar(fontsize = 6),
                                                    legend_height = unit(2, "cm")))
  list_heatmaps[[i]] <- heatmap_cnv
  
}
