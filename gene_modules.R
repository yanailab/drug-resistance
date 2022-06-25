library(Seurat)
library(bayesbio)
library(rlist)
library(parallel)

#######################
# Loading Seurat data #
#######################

# These Are Seurat objects initialized and with Scaling by umis, mitochondrial,
# and ribosomal genes regressed out; Variable genes and PCA. The script that
# originated are the ones from each sample.


load("kura_all_g1_together_seurat.rda")
load("genes_to_remove.rda")


find_markers_parallel <- function(seurat_obj, cluster_list, ncores=4) {
  # from a seurat object and a list of cell clusters, test the markers for
  # each cell cluster.
  # Input: seurat_obj; clust_list (output of create_cell_clusters())
  # Returns a list of dataframes of markers for each cluster that should be
  # filtered out later.
  
  # x is the variable for each cell cluster (an element of the list)
  #register(SnowParam(ncores))
  find_marker_wrapper <- function(x, seurat_obj){
    Seurat::FindMarkers(seurat_obj, ident.1 = x, 
                        min.pct = 0.2, min.cells.gene  = 10,
                        test.use = "MAST", only.pos = T, print.bar = F)
  }
  
  #register(param)
  Markers <- mclapply(cluster_list, find_marker_wrapper,
                      seurat_obj = seurat_obj, mc.cores = ncores)
  #BPPARAM=MulticoreParam(ncores))
  return(Markers)
}


cluster_cells <- function(seurat_obj, genes_to_remove, n = 1000) {
  # Use hierarchical clustering on all cells of seurat obj and variable genes
  # using 1 - pearson correlation as distance
  # input Seurat 2 object
  # Returns hclust object.
  
  # select genes and center the expression, remove mito, ribo genes
  vargene_data = subset(seurat_obj@hvg.info, !(rownames(seurat_obj@hvg.info) %in% genes_to_remove))
  var_genes = head(rownames(vargene_data), n)
  expression_center = t(scale(t(seurat_obj@data[rownames(seurat_obj@data) %in% var_genes, ]), 
                              center = T, scale = F))
  # calculate distance
  expression_dist = as.dist(1 - cor(expression_center))
  # hierarchical clustering
  cluster = hclust(expression_dist, method="average")
  return(cluster)
}


get_clusters <- function(cluster, start, end) {
  # get all possible clusters from a hclust obj
  # Input: cluster: hclust object
  # start: number of K for initial and K for end, getting a range of clusters
  # typically start is 2 and End is the number of units (cells) clustered.
  
  clusters = cutree(cluster, start:end)
  return(clusters)
}


create_cell_clusters <- function(cluster_table, min = 10, max = 0.85) {
  # cluster your data before, than use cuttree(cluster, k = x:y), with this
  # output is the input for this function. We go over the possible clusters
  # filtering according to min and maximum members per cluster and detect the
  # useful clusters, returning as a list of clusters and cells (members)
  # belonging to it.
  # cluster_table: cutree output
  
  cluster_table = t(cluster_table)
  total_cells = ncol(cluster_table) # total number of cells
  cluster_list = list()
  #i = 1
  for(row in 1:nrow(cluster_table)) {
    # create a frequency table from clusters for each row
    freq_cluster = as.data.frame(table(cluster_table[row, ]))
    # filter clusters based on parameters
    filtered_clusters = freq_cluster[(freq_cluster$Freq >= min) & (freq_cluster$Freq <= (max*total_cells)), ]
    valid_clusters = as.vector(filtered_clusters$Var1)
    # for each valid cluster, get the corresponding cells and append to a list
    for(c in valid_clusters) {
      cells = names(which(cluster_table[row, ] == c))
      cluster_list <- list.append(cluster_list, cells)
    }
  }
  # filter to remove duplicated clusters
  cluster_list = unique(cluster_list)
  cluster_names = paste("c", 1:length(cluster_list), sep = "")
  names(cluster_list) <- cluster_names
  return(cluster_list)
}


create_signatures <- function(seurat_obj, cluster_list, logfold = 0.58,
                              padj = 0.001, n_sig = 50) {
  # takes expression matrix with all cells, test differential expression for
  # all clusters in cluster_list and create a signature according to those
  # fold, padj and size criteria.
  # returns a list of signatures, each containing the genes associated with
  # each valid cluster
  
  signatures_list = list()
  # for here for each cluster
  for(i in 1:length(cluster_list)) {
    cluster_name = names(cluster_list[i])
    cluster_cells = cluster_list[[i]]
    signature = FindMarkers(seurat_obj, ident.1 = cluster_cells,
                            min.pct = 0.2, min.cells.gene  = 10,
                            test.use = "MAST", only.pos = T, print.bar = F)
    # test if element is not null
    if(!is.null(signature)) {
      # if not null, test if the number of diff. expressed genes > n_sig
      if(nrow(signature) >= n_sig) {
        # filtering genes by fold change and padj
        signature = signature[(signature$avg_logFC >= logfold) & (signature$p_val_adj <= padj), ]
        genes_signature = rownames(signature)
        # if signature has more than n_sig members, keep it
        if(length(genes_signature) > n_sig) {
          signatures_list[[cluster_name]] <- genes_signature
          #signatures_list = rlist::list.append(signatures_list, genes_signature)
        }
      }
    }
  }
  # give it a name using the cluster name for this list
  return(signatures_list)
}


filter_cluster_signatures <- function(cluster_list, signatures_list,
                                       jaccard_th = 0.75) {
  # take clusters of cells and signatures of those clusters and compare the
  # pairwise jaccard index across cell clusters to select the ones with bigger
  # size of cells cluster for those which overlap > 0.75
  # cluster_list: output of create_cell_clusters()
  # signatures_list: output of create_signatures()
  # Returns: vector with ids of clusters to keep. Useful to filter out 
  # cluster_list and signatures_list for further analysis: for example 
  # filtering signatures for jaccard clustering.
  
  # subsetting the clusters with valid signatures
  cluster_list = cluster_list[names(signatures_list)]
  # pairwise cluster jaccard overlap matrix
  cluster_overlap = sapply(cluster_list, 
                           function(x) sapply(cluster_list,
                                              function(y) jaccardSets(x,y)))
  # get the comparisons which have jaccard > 0.75
  clusters_to_remove = c()
  for(i in 1:nrow(cluster_overlap)) {
    # name of the cluster in row i
    cluster_test = rownames(cluster_overlap)[i]
    # cluster pairs that have high overlap (> 0.75) excluding itself ( < 1)
    row_jacc = names(which((cluster_overlap[i, ] > jaccard_th) & (cluster_overlap[i, ] < 1)))
    # if any pair has overlap > than jaccard_th, test to remove it.
    if(length(row_jacc) >= 1) {
      # size of the signature for cluster_test
      cluster_sig_size = length(signatures_list[[cluster_test]])
      # size of the signatures of the pairs
      sig_sizes = sapply(signatures_list[row_jacc], length)
      # if cluster not in merge to remove then test
      if(!(cluster_test %in% clusters_to_remove)) {
        # size of the cluster (n of cells) for cluster_test
        cluster_cell_size = length(cluster_list[[cluster_test]])
        # size of the signatures of the pairs
        cell_sizes = sapply(cluster_list[row_jacc], length)
        c_names = c(cluster_test, row_jacc) # join names
        c_sizes = c(cluster_cell_size, cell_sizes) # join sizes
        c_to_remove = c_names[-which.max(c_sizes)] # remove the smaller ones
        clusters_to_remove = append(clusters_to_remove, c_to_remove) # add to remove
      }
    }
  }
  
  all_clusters = names(cluster_list)
  clusters_to_remove = unique(clusters_to_remove)
  clusters_to_keep = subset(all_clusters, !(all_clusters %in% clusters_to_remove))
  return(clusters_to_keep)
}


################################################################################
# Running functions and get outputs                                            #
################################################################################


# 1. Cluster cells
cat("Step1: Clustering cells\n")
cell_cluster = cluster_cells(kura_all, genes_to_remove) 

#plot(cell_cluster, label = F)

# 2. Cutree
cat("Step2: Getting clusters\n")
cluster_size = ncol(kura_all@data)
clusters = get_clusters(cell_cluster, start = 2, end = cluster_size)

# 3. Create cell clusters
cat("Step3: Creating cell clusters\n")
cluster_list = create_cell_clusters(clusters, min = 30, max = 0.80)

# 4. Find all markers for clusters
cat("Step4: Finding all markers for clusters: parallel\n")
signatures = find_markers_parallel(kura_all, cluster_list[1:4], ncores = 4)


# 5. create signatures
cat("Step5: Creating signatures\n")
signature_list = create_signatures(kura_all, cluster_list)


# 6. Filter clusters/signatures
cat("Step6: Filtering clusters and signatures\n")
filtered_clusters = filter_cluster_signatures(cluster_list, signature_list,
                                              jaccard_th = 0.75)
cluster_list = cluster_list[filtered_clusters]
signature_list = signature_list[filtered_clusters]

# 7. Saving the objects
cat("Step 7: Saving the objects\n")
save(cell_cluster, clusters, signatures,
     cluster_list, signature_list, file = "neftel_signatures.rda")
cat("Done!\n")
