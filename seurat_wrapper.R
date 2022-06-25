###################################
# Seurat v2.3.1 wrapper functions #
###################################

library(Seurat)
library(magrittr)
library(dplyr)
library(RColorBrewer)

source("singleCell_QC.R")

# Step 1: Create seurat object
create_seurat_object <- function(exp_matrix, names.field = 2, names.delim = "_",
                                 proj_name = "Sample") {
  # creates a new seurat object
  # With percentage of mitochondrial and ribosomal genes assigned
  
  obj = CreateSeuratObject(counts = exp_matrix, 
                           project = proj_name, names.field = names.field,
                           names.delim = names.delim)
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(obj@assays$RNA@counts), value = TRUE)
  #ribo.genes <- grep(pattern = "^RPL", x = rownames(x = kura@data), value = TRUE)
  percent.mito = colSums(as.matrix(obj@assays$RNA@counts)[mito.genes, ]) / colSums(as.matrix(obj@assays$RNA@counts))
  #exp_ = exp_matrix[!(rownames(exp_matrix) %in% ribo_genes), ]
  percent.ribo = colSums(as.matrix(obj@assays$RNA@counts)[rownames(obj@assays$RNA@counts) %in% ribo_genes, ]) / colSums(as.matrix(obj@assays$RNA@counts))
  obj <- AddMetaData(obj, metadata = percent.mito, 
                     col.name = "percent.mito")
  obj <- AddMetaData(obj, metadata = percent.ribo, 
                     col.name = "percent.ribo")
  
  return(obj)
}

# Step 2: Add custom normalized data
add_custom_norm_data <- function(obj, exp_norm) {
  # assign normalized data to Seurat (when custom normalized data is performed)
  
  obj@assays$RNA@data <- as(as.matrix(exp_norm), "dgCMatrix") 
  return(obj)
}


# Step 1: Create seurat object, V 2.0
create_seurat_object <- function(exp_matrix, names.field = 2, names.delim = "_",
                                     proj_name = "Sample") {
  # creates a new seurat object
  # With percentage of mitochondrial and ribosomal genes assigned
  
  obj = CreateSeuratObject(exp_matrix, 
                           project = proj_name, names.field = names.field,
                           names.delim = names.delim)
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(obj@raw.data), value = TRUE)
  #ribo.genes <- grep(pattern = "^RPL", x = rownames(x = kura@data), value = TRUE)
  percent.mito = colSums(as.matrix(obj@raw.data[mito.genes, ])) / colSums(as.matrix(obj@raw.data))
  #exp_ = exp_matrix[!(rownames(exp_matrix) %in% ribo_genes), ]
  percent.ribo = colSums(as.matrix(obj@raw.data[rownames(obj@raw.data) %in% ribo_genes, ])) / colSums(as.matrix(obj@raw.data))
  obj <- AddMetaData(obj, metadata = percent.mito, 
                     col.name = "percent.mito")
  obj <- AddMetaData(obj, metadata = percent.ribo, 
                     col.name = "percent.ribo")
  
  return(obj)
}

# Step 2: Add custom normalized data
add_custom_norm_data <- function(obj, exp_norm) {
  # assign normalized data to Seurat (when custom normalized data is performed)
  
  obj@data <- as(as.matrix(exp_norm), "dgCMatrix") 
  return(obj)
}


# Run cell cycle analysis
cell_cycle <- function(obj) {
  # runs Seurat cell cycle classification
  
  # Read in a list of cell cycle markers, from Tirosh et al, 2015
  cc.genes <- readLines(con = "/home/gu/Dropbox/projects/kuramochi_olaparib_evolution/datasets_external/regev_lab_cell_cycle_genes.txt")
  
  # We can segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]
  
  # Assign cell cycle scores
  obj <- CellCycleScoring(obj, s.genes = s.genes, g2m.genes = g2m.genes, 
                          set.ident = TRUE)
  return(obj)
}


get_g1_cells <- function(obj) {
  # get cells in G1 (G2M + S score < -0.1)
  # Use SubsetData function from Seurat to subset those cells from the original
  # object. Ex: SubsetData(obj, cells.use = g1_cells)
  g1_cells = rownames(obj@meta.data[(obj@meta.data$Phase == "G1"), ])
  return(g1_cells)
}


get_g2m_cells <- function(obj) {
  # get cells in G1 (G2M + S score < -0.1)
  # Use SubsetData function from Seurat to subset those cells from the original
  # object. Ex: SubsetData(obj, cells.use = g1_cells)
  g2m_cells = rownames(obj@meta.data[(obj@meta.data$Phase == "G2M"), ])
  return(g2m_cells)
}


get_s_cells <- function(obj) {
  # get cells in G1 (G2M + S score < -0.1)
  # Use SubsetData function from Seurat to subset those cells from the original
  # object. Ex: SubsetData(obj, cells.use = g1_cells)
  s_cells = rownames(obj@meta.data[(obj@meta.data$Phase == "S"), ])
  return(s_cells)
}


#get_top_pc_genes <- function(seurat_object, pc = 1, n = 50) {
  # get PC X genes from seurat object, sorted in Decreasing order, so both top N
  # positive and top N negative ones are more meaningful.
  
#  top_up = names(head(sort(seurat_object@dr$pca@gene.loadings[, pc], decreasing = T), n))
#  top_down = names(head(sort(seurat_object@dr$pca@gene.loadings[, pc], decreasing = F), n))
  
#  pc_genes = unique(c(top_up, top_down))
#  return(pc_genes)
#}


get_top_pc_genes_old <- function(seurat_object, pc = 1, n = 40) {
  # get PC X genes from seurat object, sorted in Decreasing order, so both top N
  # positive and top N negative ones are more meaningful. Works for seurat v > 3.0
  top_up = names(head(sort(seurat_object@dr$pca@gene.loadings[, pc], decreasing = T), n))
  top_down = names(head(sort(seurat_object@dr$pca@gene.loadings[, pc], decreasing = F), n))
  
  pc_genes = unique(c(top_up, top_down))
  return(pc_genes)
}


get_top_diff_genes <- function(seurat_diff_exp, n = 100) {
  # get top N differentially expressed genes returned by Seurat differential
  # expression table. The matrix is already sorted by p-value, so it is just
  # getting the head() of N elements for each cluster/group.
  
  diff_genes = c()
  clusters = levels(unique(seurat_diff_exp$cluster))
  for(cluster in clusters) {
    # get genes for that cluster. If number of rows/diff expressed genes
    # in that cluster is lower, get all of them
    if(nrow(seurat_diff_exp) < n) {
      genes = seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene
    } else {
      genes = head(seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene, n)
    }
    
    diff_genes = c(diff_genes, genes)
  
  }
  return(diff_genes)
}


summary_exp_by_cluster <- function(norm_exp_matrix, metadata, cluster_col,
                                   threshold = 0.01, suffix = "S1") {
  # for an expression matrix, return the average normalized expression 
  # across all cells, the sum of raw counts, the frequency of cells in which
  # the gene is expressed above threshold for raw or normalized expression.
  # also store the standard deviation.
  # this is better used for averaging already normalized data
  
  # iterate over each Cluster ID
  # create a list of vectors e and create a dataframe from it
  
  # List of of dataframes
  all_cluster_data = list()
  for(cluster in levels(as.factor(metadata[, cluster_col]))) {
    # get cell ids for a particular cluster
    cell_ids = rownames(metadata[which(metadata[, cluster_col] == cluster), ])
    # get avg norm exp
    avg_norm = round(rowMeans(norm_exp_matrix[, cell_ids]), 5)
    # get standard deviation of avg_norm
    std_norm = round(apply(norm_exp_matrix[, cell_ids], 1, sd), 5)
    
    # frequency of expression
    calc_freq <- function(vector) {
      n = sum(vector > threshold)
      freq = round(n / length(vector), 3)
      return(freq)
    }
    
    freq = apply(norm_exp_matrix[, cell_ids], 1, calc_freq)
    # creating a dataframe for the cluster
    cluster_data = data.frame(avg_norm = avg_norm, freq_norm = freq, 
                              std_norm = std_norm)
    # fixing colnames
    colnames(cluster_data) <- paste(paste(colnames(cluster_data), suffix,
                                          sep = "_"), cluster, sep = "_c")
    # add rownames for further merging dataframes from all clusters
    rownames(cluster_data) <- rownames(norm_exp_matrix)
    # appending 
    all_cluster_data = rlist::list.append(all_cluster_data, cluster_data)
  }
  # merging data from all clusters from the list to use with merge_exp()
  all_data = merge_exp(all_cluster_data)
  return(all_data)
}  


summary_exp_by_cluster_counts <- function(raw_exp_matrix, metadata, cluster_col,
                                          suffix = "S1") {
  # for an expression matrix, return the average normalized expression 
  # across all cells, the sum of raw counts, the frequency of cells in which
  # the gene is expressed above threshold for raw or normalized expression.
  # also store the standard deviation.
  # this is better used for averaging already normalized data
  
  # iterate over each Cluster ID
  # create a list of vectors e and create a dataframe from it
  
  # List of of dataframes
  all_cluster_data = list()
  for(cluster in levels(as.factor(metadata[, cluster_col]))) {
    # get cell ids for a particular cluster
    cell_ids = rownames(metadata[which(metadata[, cluster_col] == cluster), ])
    # get avg norm exp
    sum_exp = rowSums(raw_exp_matrix[, cell_ids])
    # creating a dataframe for the cluster
    cluster_data = data.frame(sum_exp = sum_exp)
    # fixing colnames
    colnames(cluster_data) <- paste(paste(colnames(cluster_data), suffix,
                                          sep = "_"), cluster, sep = "_c")
    # add rownames for further merging dataframes from all clusters
    rownames(cluster_data) <- rownames(raw_exp_matrix)
    # appending 
    all_cluster_data = rlist::list.append(all_cluster_data, cluster_data)
  }
  # merging data from all clusters from the list to use with merge_exp()
  all_data = merge_exp(all_cluster_data)
  return(all_data)
}  


avg_freq_exp_by_cluster <- function(norm_exp_matrix, metadata, cluster_col,
                                    threshold = 0.01, sample_id = "C") {
  # for an expression matrix, return the average normalized expression 
  # across all cells, the frequency of cells in which
  # the gene is expressed above threshold for raw or normalized expression.
  # this is better used for averaging already normalized data.
  # Returns a dataframe with Gene, Sample, cluster, avg_clust expression, freq_clust
  # expression above threshold.
  
  # Final dataframe
  result_all = data.frame(gene = c(), sample = c(), clusters = c(), 
                          avg_clust = c(), freq_clust = c())
  
  for(cluster in levels(as.factor(metadata[, cluster_col]))) {
    # get cell ids for a particular cluster
    cell_ids = rownames(metadata[which(metadata[, cluster_col] == cluster), ])
    # get avg norm exp
    avg_clust = round(rowMeans(norm_exp_matrix[, cell_ids]), 5)
    # get standard deviation of avg_norm
    #std_norm = round(apply(norm_exp_matrix[, cell_ids], 1, sd), 5)
    
    # frequency of expression
    calc_freq <- function(vector) {
      n = sum(vector > threshold)
      freq = round(n / length(vector), 3)
      return(freq)
    }
    
    freq_clust = apply(norm_exp_matrix[, cell_ids], 1, calc_freq)
    
    result = data.frame(gene = rownames(norm_exp_matrix), 
                        sample = rep(sample_id, nrow(norm_exp_matrix)),
                        clusters = rep(cluster, nrow(norm_exp_matrix)),
                        avg_clust = as.numeric(avg_clust),
                        freq_clust = as.numeric(freq_clust))
    
    result_all = rbind(result_all, result)
   
  }
  #result_all$clusters <- paste("c", result_all$clusters, sep = "")
  return(result_all)
}  


get_top_diff_genes_to_list <- function(seurat_diff_exp, 
                                       sample_name = "sample", all = T, n = 100) {
  # get top N differentially expressed genes returned by Seurat differential
  # expression table. The matrix is already sorted by p-value, so it is just
  # getting the head() of N elements for each cluster/group. You can run for each
  # matrix (sample) then merge all into a single list using rlist::list.merge
  
  diff_genes = list()
  clusters = levels(unique(seurat_diff_exp$cluster))
  list_counter = 0
  for(cluster in clusters) {
    # get genes for that cluster. If number of rows/diff expressed genes
    # in that cluster is lower, get all of them
    
    if(all == T) {
      genes = seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene
    } else {
      if(nrow(seurat_diff_exp) < n) {
        genes = seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene
      } else {
        genes = head(seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene, n)
      }
    }
    
    list_counter = list_counter + 1
    diff_genes = rlist::list.append(diff_genes, genes)
    names(diff_genes)[list_counter] <- paste(sample_name, 
                                             paste("c", cluster, sep = ""),
                                             sep = "_") 
    #diff_genes = c(diff_genes, genes)
  }
  return(diff_genes)
}


get_cluster_genes <- function(gene_markers_tables, clusters_to_select,
                              column = "cluster", unique = T) {
  # This is useful to select genes from particular clusters to further subset
  # expression matrices or use the genes for pseudotime or clustering.
  # arg1: gene_markers_tables: A list of seurat diff. expression table for each
  #       sample. Ex: list(kura_ctrl_g1_markers, kura_T1_g1_markers)
  # arg2: clusters_to_select: A list of clusters ids for each of the diff.
  #       expression tables to select. Ex: (list(c(0, 1), c(1, 2), etc.)). 
  #       It must be in the same clusters in the same order you pass your tables.
  # arg3 (optional): Number of genes to select (top p-value) or avg_fold?
  #      do this for later
  # !TODO: Make it more flexible to select genes from particular avg_fold or pval filter
  # like adding another filtering steps if true.
  # output: A vector of genes representing each selected clusters.
  
  genes = c()
  counter = 1
  
  for(table in gene_markers_tables) { 
    filtered_genes = dplyr::filter(table, UQ(as.name(column)) %in% clusters_to_select[[counter]])$gene
    genes = c(genes, filtered_genes)
    counter = counter + 1
  }
  if(unique == T) {
    return(unique(genes))  
  } else {
    return(genes)
  }
}


get_cluster_cells <- function(cells_metadata_tables, clusters_to_select, 
                              column = "res.0.3") {
  # This is useful to select cells from particular clusters to further subset
  # expression matrices.
  
  # arg1: gene_markers_tables: A list of seurat diff. expression table for each
  #       sample. Ex: list(kura_ctrl_g1_markers, kura_T1_g1_markers)
  # arg2: clusters_to_select: A list of clusters ids for each of the diff.
  #       expression tables to select. Ex: (list(c(0, 1), c(1, 2), etc.)). 
  #       It must be in the same clusters in the same order you pass your tables.
  # arg3: which column to evaluate. Default is cluster from res.0.3
  #       The trick here to transform a string to an argument is done by 
  #       UQ(as.name(column)) in dplyr.
  # output: A vector of genes representing each selected clusters.
  
  cells = c()
  counter = 1
  
  for(table in cells_metadata_tables) { 
    # we need to make rownames as a column before filtering because dplyr removes it
    table$cells <- rownames(table)
    filtered_cells = dplyr::filter(table, UQ(as.name(column)) %in% clusters_to_select[[counter]])$cells
    cells = c(cells, filtered_cells)
    counter = counter + 1
  }
  return(cells)  
}


create_pheno_data <- function(cells_metadata_tables, selected_cells) {
  # create a pheno_data object for monocle input. The pheno_data contains columns
  # from metadata object from seurat, but with the selected cells subset.
  
  # arg1: Same as get_cluster_cells. A list of metadata seurat objects.
  # arg2: vector of cells to select from. Must contain the union of all cells
  #       wanted to merge the metadata.
  # output: A pheno_data dataframe with orig_ident, cluster, ident_cluster cols.
  
  ident = c()
  cluster = c()
  ident_cluster = c()
  
  for(table in cells_metadata_tables) { 
    # we need to make rownames as a column before filtering because dplyr removes it
    filtered_table = table[rownames(table) %in% selected_cells, ]
    orig_ident = as.character(filtered_table$orig.ident)
    clusters = as.character(filtered_table$res.0.3)
    ident_clusters = paste(orig_ident, clusters, sep = "_c")
    # filling the vectors
    ident = c(ident, orig_ident) 
    cluster = c(cluster, clusters)
    ident_cluster = c(ident_cluster, ident_clusters)
  }
  pheno_data = data.frame(orig_ident = ident,
                          cluster = cluster, 
                          ident_cluster = ident_cluster)
  return(pheno_data)  
}


genes_to_list <- function(gene_markers_tables, sample_names, 
                          column = "cluster") {
  # This is useful to create a list of vectors for several diff. exp tables.
  # like Seurat Diff. expression output. Uniting a list of tables and clusters
  # arg1: gene_markers_tables: A list of seurat diff. expression table for each
  #       sample. Ex: list(kura_ctrl_g1_markers, kura_T1_g1_markers)
  
  genes_list = list()
  counter_table = 1
  element = 1
  
  for(table in gene_markers_tables) {
    # get the cluster labels
    clusters = levels(unique(table$cluster))
    
    for(cl in clusters) {
      #print(cluster)
      filtered_genes = dplyr::filter(table, UQ(as.name(column)) == cl)$gene
      #filtered_genes = table[table$cluster == cluster, ]$gene
      genes_list = rlist::list.append(genes_list, filtered_genes)
      # fix how to give names to list elements and figure out the counter
      names(genes_list)[element] <- paste(sample_names[counter_table], cl, 
                                          sep = "_c")
      element = element + 1
    }
    counter_table = counter_table + 1
  }
  return(genes_list)
}


get_top_diff_genes_to_list <- function(seurat_diff_exp, 
                                       sample_name = "sample", all = T, n = 100) {
  # get top N differentially expressed genes returned by Seurat differential
  # expression table. The matrix is already sorted by p-value, so it is just
  # getting the head() of N elements for each cluster/group. You can run for each
  # matrix (sample) then merge all into a single list using rlist::list.merge
  
  diff_genes = list()
  clusters = levels(unique(seurat_diff_exp$cluster))
  list_counter = 0
  for(cluster in clusters) {
    # get genes for that cluster. If number of rows/diff expressed genes
    # in that cluster is lower, get all of them
    
    if(all == T) {
      genes = seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene
    } else {
      if(nrow(seurat_diff_exp) < n) {
        genes = seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene
      } else {
        genes = head(seurat_diff_exp[seurat_diff_exp$cluster == cluster, ]$gene, n)
      }
    }
    
    list_counter = list_counter + 1
    diff_genes = rlist::list.append(diff_genes, genes)
    names(diff_genes)[list_counter] <- paste(sample_name, 
                                             paste("c", cluster, sep = ""),
                                             sep = "_") 
    #diff_genes = c(diff_genes, genes)
  }
  return(diff_genes)
}


genes_in_clusters <- function(genes, list_genes) {
  # genes: should contain the universe of all genes that should be looked for.
  # for a vector of genes, creates a matrix of its presence (1) absence (0) in
  # the selected clusters. This will help to find genes exclusive or shared
  # between clusters and states.
  
  data_output = data.frame()
  
  for(gene in genes) {
    gene_presence = c()
    
    for(cluster in list_genes) {
      gene_status = 0
      if(gene %in% cluster) {
        gene_status = 1
      }
      gene_presence = c(gene_presence, gene_status)
    }
    data_output = rbind(data_output, gene_presence)
  }
  rownames(data_output) <- genes
  colnames(data_output) <- names(list_genes)
  return(data_output)
}


