########################################
# Analysis of gene modules and ploting #
########################################


load("neftel_signatures.rda")

cluster_list_filtered = cluster_list[filtered_clusters]
signature_list_filtered = signature_list[filtered_clusters]


filter_signatures <- function(signature_list_filtered, cluster_list_filtered,
                              seurat_obj, min_exp_sig = 0.1) {
  # filter genes in signatures with expression less than min_exp_sig in 
  # the signature to avoid non-important genes
  
  # new filtered signature genes
  signature_genes_filtered = list()
  # find markers for all clusters
  cluster_names = names(signature_list_filtered)
  # iterate over each marker dataframe for each cluster 
  for(cluster_name in cluster_names) {
    genes_signature = signature_list_filtered[[cluster_name]]
    cells_signature = cluster_list_filtered[[cluster_name]]
    genes_data = as.matrix(seurat_obj@data[genes_signature, cells_signature])
    genes_data = gene_filter(as.matrix(genes_data), exp_level = min_exp_sig, n_samples = 10)
    
    #genes_mean = rowMeans(genes_data)
    # filter genes with low expression
    #genes = names(genes_mean[genes_mean > min_exp_sig])
    genes = rownames(genes_data)
    signature_genes_filtered[[cluster_name]] <- genes
  }
  
  return(signature_genes_filtered)
}

signature_list_filtered = filter_signatures(signature_list_filtered, cluster_list_filtered,
                                            kura_all, min_exp_sig = 0.1)
# keep signatures with more than 50 genes
signature_list_filtered = signature_list_filtered[lapply(signature_list_filtered, length) > 50]
# remove big clusters (k > 5)
remove_clusters = c("c1", "c2", "c3", "c4", "c5", "c6")
all_clusters = names(signature_list_filtered)
valid_clusters = subset(all_clusters, !(all_clusters %in% remove_clusters))

cluster_list_filtered = cluster_list_filtered[valid_clusters]
signature_list_filtered = signature_list_filtered[valid_clusters]


signature_overlap = sapply(signature_list_filtered, 
                           function(x) sapply(signature_list_filtered,
                                              function(y) jaccardSets(x,y)))

#signature_dist = as.dist(signature_overlap)
signature_dist = as.dist(1 - cor(signature_overlap))

#signature_dist[is.na(signature_dist)] <- 0
#signature_dist[is.nan(signature_dist)] <- 0
# hierarchical clustering
cluster_signatures = hclust(signature_dist, method="average")

plot(cluster_signatures)
cutree(cluster_signatures, k = 4)
dendextend::cutree(cluster_signatures2, k = 4)
#correlation 
cormM = cor(signature_overlap, method="pearson")


library(dendextend)
dend = as.dendrogram(hclust(signature_dist, method = "average")) # can be used as cluster_rows and columns arg instead of T
dend <- click_rotate(dend, continue = TRUE)

# keep clusters with > 3 signatures
remove_clusters = c("c202", "c21", "c203", "c383", "c8", "c18", "c17")
all_clusters = names(signature_list_filtered)
valid_clusters = subset(all_clusters, !(all_clusters %in% remove_clusters))

cluster_list_filtered = cluster_list_filtered[valid_clusters]
signature_list_filtered = signature_list_filtered[valid_clusters]

signature_overlap = sapply(signature_list_filtered, 
                           function(x) sapply(signature_list_filtered,
                                              function(y) jaccardSets(x,y)))

signature_dist = as.dist(1 - cor(signature_overlap))
cormM = cor(signature_overlap, method="pearson")
dend = as.dendrogram(hclust(signature_dist, method = "average"))
dend <- dendextend::click_rotate(dend, continue = TRUE)


pdf("heatmap_signatures_kura_all.pdf")
heatmap_signatures = Heatmap(cormM, show_column_names = F, show_row_dend = T, 
                             show_column_dend = F, show_row_names = F, 
                             name = "Corr", row_names_gp = gpar(fontsize = 8),
                             column_names_gp = gpar(fontsize = 8),
                             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(50),
                             cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
                             column_dend_reorder = F,
                             row_split = 6, column_split = 6, 
                             row_title = NULL, column_title = c("A", "B", "C", "D", "E", "F"),
                             column_title_gp = gpar(fontsize = 8),
                             row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
                             #heatmap_width = unit(10, "cm"), heatmap_height = unit(12, "cm"),
                             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                                         title = "Correlation",
                                                         labels_gp = gpar(fontsize = 8),
                                                         legend_height = unit(2, "cm"),
                                                         border = "black"))

draw(heatmap_signatures)

plot(dend)


################################################################################
# Creating modules
################################################################################


create_meta_modules <- function(seurat_obj, clusters_dend, cluster_list_filtered,
                                 signature_list_filtered, height, n_sigs = 3) {
  # create modules by determining which genes belong to which modules
  # seurat_obj: seurat object with data containing expression matrix
  # clusters: dendrogram object created for clustering
  # cluster_list_filtered = cluster / cells list
  # signature_list_filtered = cluster / signatures
  # height: Where to cut cluster to get clusters
  
  # get clusters with signature names
  clusters = dendextend::cutree(clusters_dend, h = height)
  # numbers of modules present
  modules = unique(clusters)
  avg_modules = list()
  exp_data = as.data.frame(as.matrix(seurat_obj@data))
  
  for(module in modules) {
    # get names of the signatures for the module
    signatures = names(clusters[clusters == module])
    module_len = length(signatures) # number of members in the module
    # get all genes in the module
    #genes_module = unique(Reduce(c, signature_list_filtered[signatures]))
    genes_module = table(Reduce(c, signature_list_filtered[signatures]))
    #genes_module = names(genes_module[genes_module > (0.3*module_len)])
    genes_module = names(genes_module[genes_module >= n_sigs])
    cells_module = unique(Reduce(c, cluster_list_filtered[signatures]))
    
    #genes_module = genes_module[genes_module > 0.1]
    # keep log_ratios for all signatures of a module
    
    all_cells = colnames(seurat_obj@data)
    rest_cells = subset(all_cells, !(all_cells %in% cells_module))
    # expression in module
    avg_exp_module = rowMeans(exp_data[genes_module, cells_module])
    # expression in rest of cells
    avg_exp_rest = rowMeans(exp_data[genes_module, rest_cells])
    # ratio expression
    #log2ratio = log2(avg_exp_module / avg_exp_rest)
    log2ratio = avg_exp_module - avg_exp_rest
    
    # append logs
    avg_modules[[module]] <- log2ratio
  }
  return(avg_modules)
}


meta_modules_h_0.95 = create_meta_modules(kura_all, dend, cluster_list_filtered,
                                          signature_list_filtered, h = 0.95)

modules = lapply(meta_modules_h_0.95, sort, decreasing = T)
modules = lapply(modules, names)
names(modules) <- c("B", "C", "A", "D", "E", "F")


save(cluster_list_filtered, signature_list_filtered, dend, signature_overlap,
     cormM, meta_modules, modules, heatmap_signatures,
     meta_modules_create3_h_0.95, file = "meta_modules_neftel.rda")

load("meta_modules_neftel.rda")


################################################################################
# Signature gene expression                                                    #
################################################################################

create_module_data <- function(seurat_obj, modules, signatures_order, 
                               cluster_list_filtered, top = 50) {
  # create dataframes for gene avg expression of modules across signatures
  data_all = list()
  i = 1
  for(module in modules) {
    genes_module = head(module, top)
    data_module = c()
    for(signature in signatures_order) {
      cells_signature = cluster_list_filtered[[signature]]
      gene_means = Matrix::rowMeans(seurat_obj@data[genes_module, cells_signature])
      data_module = cbind(data_module, gene_means)
    }
    colnames(data_module) <- signatures_order
    data_all[[i]] <- data_module
    i = i + 1
  }
  return(data_all)
}

x = lapply(meta_modules_create3_h_0.95, sort, decreasing = T)
y = lapply(x, names)
signatures_order = labels(dend) # gets the leaf order of the dendrogram 

modules_genes_exp = create_module_data(kura_all, y, signatures_order,
                                       cluster_list_filtered, top = 100)


################################################################################
# Plotting heatmap
################################################################################

# Calculate the frequency of cells from each signature for bottom annotation
# calculating frequencies for samples and states
metadata = kura_all@meta.data
# fix state names
#metadata[metadata$res.0.3 == "0",]
#head(metadata)

cell_freqs <- function(x, metadata, column = 3) {
  # Calculate frequencies of cells in each condition (sample = column 3)
  # Returns vector of frequencies of cells belonging to each condition
  
  freq = table(metadata[x, column])
  freq = (freq / length(x))*100
  return(freq)
}

# reordering to calculate in the order of the clusters for heatmap
cluster_cells_order = cluster_list_filtered[signatures_order]

cluster_freqs = t(sapply(cluster_cells_order, cell_freqs, metadata, 3))
reorder_cols = c("ctrlI", "1uM", "2.5uM", "5uM", "10uM", "20uM", "40uM2nd", 
                 "80uM", "160uM", "320uM")
cluster_freqs = cluster_freqs[, reorder_cols]
colnames(cluster_freqs) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
cluster_freqs

new_colors = c("#8DA3A6","#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
               "#D0B541", "#E67F33", "#CE2220", "#80271E")

column_ha = HeatmapAnnotation(foo = anno_barplot(cluster_freqs, border = F,
                                                 gp = gpar(fill = new_colors,
                                                           col = new_colors,
                                                           lwd = 0.1)),
                              show_annotation_name = c(foo = FALSE))


#draw(ha)
states_col = metadata$res.0.3
states_col[states_col == 0] <- "C"
states_col[states_col == 1] <- "S1"
states_col[states_col == 2] <- "S2"
states_col[states_col == 3] <- "S3"
states_col[states_col == 4] <- "S4"
states_col[states_col == 5] <- "S5"
states_col[states_col == 6] <- "S5"

# giving new names to samples for plotting
metadata$res.0.3 <- as.factor(states_col)

state_freqs = t(sapply(cluster_cells_order, cell_freqs, metadata, column = 6))
state_freqs

state_colors = c("#8DA3A6", "#a9d466", "#4cc687", "#00b8ab", "#247a98", "#3f3f69", "#3f3f69")

column_ha = HeatmapAnnotation(sample = anno_barplot(cluster_freqs, border = F,
                                                    gp = gpar(fill = new_colors,
                                                              col = new_colors,
                                                              lwd = 0.1,
                                                              fontsize = 7),
                                                    height = unit(0.7, "cm")),
                              state = anno_barplot(state_freqs, border = F,
                                                   gp = gpar(fill = state_colors,
                                                             col = state_colors,
                                                             lwd = 0.1,
                                                             fontsize = 7),
                                                   height = unit(0.7, "cm")),
                              show_annotation_name = c(sample = F, state = F))

# function to scale each module separately
scale_data <- function(x) {
  data_scale = t(scale(t(x), center = T, scale = T))
  data_scale[data_scale <= -3] = -3
  data_scale[data_scale >= 3] = 3
  return(data_scale)
}


modules_scaled = lapply(modules_genes_exp, scale_data)

## Function to show genes on heatmap
show_genes <- function(genes_to_show, data_scaled) {
  # return annotation for plotting heatmap
  gene_indices = match(genes_to_show, rownames(data_scaled))
  gene_annotation = rowAnnotation(genes = anno_mark(at = gene_indices, 
                                                    labels = genes_to_show,
                                                    labels_gp = gpar(fontsize = 6),
                                                    link_width = unit(3, "mm")))
  return(gene_annotation)  
}

# Labeling individual genes
genes_to_show = c("IFI6", "IFI27", "STAT1", "KRT8", "BCAM", "DSP", "CDKN2A", 
                  "SOX17", "PARP1", "CDH6", "PAX8")
gene_annotation = show_genes(genes_to_show, modules_scaled[[3]])
ht1 = Heatmap(modules_scaled[[3]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), show_column_names = F, border = T,
              show_heatmap_legend = F, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation)
draw(ht1)

genes_to_show = c("FN1", "KRT8", "CLDN4", "PRKCB", "KRT18", "IGF2", "KLK6", "CRABP2", "RBP1", "FYN")
gene_annotation = show_genes(genes_to_show, modules_scaled[[1]])
ht2 = Heatmap(modules_scaled[[1]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), show_column_names = F, border = T,
              show_heatmap_legend = F, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation)

genes_to_show = c("COL3A1", "MEST", "SERPINE2", "CYP1B1", "PKDCC", "FN1", "HES1", "SPARC",
                  "BCAT1", "CCND1")
gene_annotation = show_genes(genes_to_show, modules_scaled[[2]])
ht3 = Heatmap(modules_scaled[[2]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), show_column_names = F, border = T,
              show_heatmap_legend = F, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation)

genes_to_show = c("FN1", "PHLDA1", "CD44", "EGR1", "FOS", "VIM",
                  "DUSP1", "SERPINE1", "DAB2", "ZFP36", "PLK2", "PPP1R15A", "PPP1R3C",
                  "NT5E", "EGFR", "RPL13A", "RPL26", "RPL34", "RPL7", "RPS10", "RPS25",
                  "RPS9")
genes_to_show = c("FN1", "TGFBI", "PHLDA1", "CD44", "EGR1", "FOS", "VIM",
                  "DUSP1", "SERPINE1", "ZFP36", "PLK2", "DEPTOR",
                  "NT5E", "EGFR", "RPL26", "RPL34")
gene_annotation = show_genes(genes_to_show, modules_scaled[[4]])
ht4 = Heatmap(modules_scaled[[4]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), show_column_names = F, border = T,
              show_heatmap_legend = F, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation)
#draw(ht4)

genes_to_show = c("RPS21", "RPL28", "RPS5", "CD44", "DDIT4", "RPS26", "GAPDH",
                  "HMOX1", "RPL38", "VIM", "FTH1", "FTL", "NDRG1", "RACK1",
                  "GPX4", "PHLDA1", "NT5E", "PFKP", "SCD", "HIF1A", "ASNS",
                  "BNIP3", "TKT", "HIGD2A", "ALDOC", "SMAD3", "CDKN1A", "VEGFA")

genes_to_show = c("RPS21", "RPS5", "CD44", "DDIT4", "GAPDH",
                  "HMOX1", "VIM", "NDRG1", "RACK1",
                  "GPX4", "PHLDA1", "PFKP", "HIF1A",
                  "TKT", "VEGFA")
gene_annotation = show_genes(genes_to_show, modules_scaled[[5]])
ht5 = Heatmap(modules_scaled[[5]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), show_column_names = F, border = T,
              show_heatmap_legend = F, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation)
#draw(ht5)

genes_to_show = c("CYP1B1", "EPRS", "NQO1", "ASNS", "KDM3A",
                  "ATF6", "FBXO2", "CPSF1", "MXI1", "IMPDH2", "CDC25B",
                  "MKNK2", "PAK4", "ETHE1")
gene_annotation = show_genes(genes_to_show, modules_scaled[[6]])
ht6 = Heatmap(modules_scaled[[6]], cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 7), show_column_names = T,
              show_heatmap_legend = T, col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50), border = T,
              right_annotation = gene_annotation,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                          title = "Expression",
                                          labels_gp = gpar(fontsize = 7),
                                          legend_height = unit(2, "cm"),
                                          border = "black"))

ht_list = ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% column_ha
draw(ht_list, ht_gap = unit(0.06, "cm"))


################################################################################
# Scoring cells with gene modules                                              #
################################################################################

mtrnr_genes <- grep(pattern = "^MTRNR", x = rownames(kura_all@data), value = TRUE)

# get top 100 genes and filter out for ribo, mito and humanin genes
genes_score = lapply(modules, function(x) subset(x, !(x %in% c(ribo_genes, mt_genes, mtrnr_genes)))[1:100])
genes_score = lapply(modules, function(x) subset(x, !(x %in% c(ribo_genes, mt_genes, mtrnr_genes))))


kura_all = AddModuleScore(kura_all, genes.list = genes_score)
module_scores = kura_all@meta.data[7:12]
colnames(module_scores) <- c("B", "C", "A", "D", "E", "F")
module_scores_scaled = as.data.frame(scale(module_scores, scale = T))

# find the highes cluster score for all cells
high_scores = apply(module_scores_scaled, 1, which.max)
#high_scores = apply(kura_scores, 1, which.max)

reorder_cells <- function(module_cells, module_scores_scaled, column_index) {
  # reorder cells by highest module score
  cells = module_scores_scaled[module_cells, column_index]
  names(cells) <- module_cells
  cells_ord = names(sort(cells, decreasing = T))
  return(cells_ord)
}

# get cells with Module 3 (S1 genes) higher score
moduleA_cells = names(high_scores[high_scores == 3])
moduleA_cells = reorder_cells(moduleA_cells, module_scores_scaled, 3)
#module3_cells = reorder_cells(module3_cells, kura_scores, 3)

moduleB_cells = names(high_scores[high_scores == 1])
moduleB_cells = reorder_cells(moduleB_cells, module_scores_scaled, 1)
#module1_cells = reorder_cells(module1_cells, kura_scores, 1)

moduleC_cells = names(high_scores[high_scores == 2])
moduleC_cells = reorder_cells(moduleC_cells, module_scores_scaled, 2)
#module2_cells = reorder_cells(module2_cells, kura_scores, 2)

moduleD_cells = names(high_scores[high_scores == 4])
moduleD_cells = reorder_cells(moduleD_cells, module_scores_scaled, 4)
#module4_cells = reorder_cells(module2_cells, kura_scores, 4)

moduleE_cells = names(high_scores[high_scores == 5])
moduleE_cells = reorder_cells(moduleE_cells, module_scores_scaled, 5)
#module5_cells = reorder_cells(module5_cells, kura_scores, 5)

moduleF_cells = names(high_scores[high_scores == 6])
moduleF_cells = reorder_cells(moduleF_cells, module_scores_scaled, 6)
#module6_cells = reorder_cells(module6_cells, kura_scores, 6)

cells_reordered = c(moduleA_cells, moduleB_cells, moduleC_cells,
                    moduleD_cells, moduleE_cells, moduleF_cells)

test = as.matrix(t(module_scores_scaled[cells_reordered, c(3,1,2,4,5,6)]))
test[test <= -3] = -3
test[test >= 3] = 3


sample_annotation = as.character(kura_all@meta.data[cells_reordered, ]$orig.ident)
sample_annotation = stringr::str_replace_all(sample_annotation, 
                                             c("ctrlI" = "C", "1uM" = "T1",
                                               "2.5uM" = "T2.5", "5uM" = "T5",
                                               "10uM" = "T10", "20uM" = "T20",
                                               "40uM2nd" = "T40", "80uM" = "T80",
                                               "160uM" = "T160", "320uM" = "T320"))

sample_annotation[sample_annotation == "3T20"] <- "T320"
table(sample_annotation)

state_annotation = as.character(kura_all@meta.data[cells_reordered, ]$res.0.3)
state_annotation = stringr::str_replace_all(state_annotation, 
                                            c("0" = "C", "1" = "S1",
                                              "2" = "S2", "3" = "S3",
                                              "4" = "S4", "5" = "S5",
                                              "6" = "S5"))
table(state_annotation)

new_colors = c("#8DA3A6","#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
               "#D0B541", "#E67F33", "#CE2220", "#80271E")

state_colors = c("#8DA3A6", "#a9d466", "#4cc687", "#00b8ab", "#247a98", "#3f3f69", "#3f3f69")


sample_state_ann = HeatmapAnnotation(treatment = sample_annotation,
                                     state = state_annotation,
                                     #simple_anno_size = unit(2, "mm"),
                                     col = list(treatment = c("C" = "#8DA3A6",
                                                              "T1" = "#B997C7", "T2.5" = "#824D99",
                                                              "T5" = "#4E78C4", "T10" = "#57A2AC",
                                                              "T20" = "#7EB875", "T40" = "#D0B541", #"T40-W" = "#dbcc92",
                                                              "T80" = "#E67F33", "T160" = "#CE2220",
                                                              "T320" = "#80271E"),
                                                state = c("C" = "#8DA3A6", "S1" = "#a9d466",
                                                          "S2" = "#4cc687", "S3" = "#00b8ab",
                                                          "S4" = "#247a98", "S5" = "#3f3f69")),
                                     annotation_legend_param = list(title_gp = gpar(fontsize = 7),
                                                                    #title = "Treatment",
                                                                    labels_gp = gpar(fontsize = 7)),
                                     show_annotation_name = c(treatment = FALSE, state = FALSE))

pdf("heatmap_signature_scores.pdf")
cells_scores = Heatmap(test, cluster_rows = F, cluster_columns = F,
                       show_row_names = T, show_column_dend = F, 
                       show_column_names = F, show_heatmap_legend = T, 
                       col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
                       bottom_annotation = sample_state_ann,
                       row_names_gp = gpar(fontsize = 7),
                       heatmap_width = unit(12, "cm"), heatmap_height = unit(6, "cm"),
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                                   title = "z-score",
                                                   labels_gp = gpar(fontsize = 7),
                                                   legend_height = unit(2, "cm"),
                                                   border = "black"))

draw(cells_scores)
dev.off()


