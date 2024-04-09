library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("/path/pdx_talazo_analysis")

load("data_all_g1_raw_r2.rda")
load("data_all_g1_raw_r1.rda")

load("pdx_all_combined.rda") # created later
#################################
# Seurat object for Replicate 1 #
#################################

pdx_all_r1 = CreateSeuratObject(data_all_g1_raw_r1, names.delim = "_", names.field = 2)
pdx_all_r1[["percent.mt"]] <- PercentageFeatureSet(pdx_all_r1, pattern = "^MT-")
pdx_all_r1 = NormalizeData(pdx_all_r1, normalization.method = "LogNormalize", scale.factor = 10000)
pdx_all_r1 = FindVariableFeatures(pdx_all_r1, selection.method = "vst", nfeatures = 1000)
pdx_all_r1 = ScaleData(pdx_all_r1, vars.to.regress = c("nCount_RNA", "percent.mt"))

pdx_all_r1 <- RunPCA(pdx_all_r1, features = VariableFeatures(object = pdx_all_r1))
DimPlot(pdx_all_r1, reduction = "pca")
ElbowPlot(pdx_all_r1)
pdx_all_r1 <- FindNeighbors(pdx_all_r1, dims = 1:10)
pdx_all_r1 <- FindClusters(pdx_all_r1, resolution = 0.3)
pdx_all_r1 <- RunUMAP(pdx_all_r1, dims = 1:10)
DimPlot(pdx_all_r1, reduction = "umap")
DimPlot(pdx_all_r1, reduction = "umap", group.by = "orig.ident")
pdx_all_r1_markers = FindAllMarkers(pdx_all_r1, test.use = "MAST",
                                    only.pos = T, 
                                    min.pct = 0.2, min.cells.feature = 10)
pdx_all_r1_markers = pdx_all_r1_markers[pdx_all_r1_markers$p_val_adj < 0.01, ]
FeaturePlot(pdx_all_r1, features = c("NOV", "IFI27", "SOX17", "NDRG1", "KRT8"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

#################################
# Seurat object for Replicate 2 #
#################################

pdx_all_r2 = CreateSeuratObject(data_all_g1_raw_r2, names.delim = "_", names.field = 2)
pdx_all_r2[["percent.mt"]] <- PercentageFeatureSet(pdx_all_r2, pattern = "^MT-")
pdx_all_r2 = NormalizeData(pdx_all_r2, normalization.method = "LogNormalize", scale.factor = 10000)
pdx_all_r2 = FindVariableFeatures(pdx_all_r2, selection.method = "vst", nfeatures = 1000)
pdx_all_r2 = ScaleData(pdx_all_r2, vars.to.regress = c("nCount_RNA", "percent.mt"))

pdx_all_r2 <- RunPCA(pdx_all_r2, features = VariableFeatures(object = pdx_all_r2))
DimPlot(pdx_all_r2, reduction = "pca")
ElbowPlot(pdx_all_r2)
pdx_all_r2 <- FindNeighbors(pdx_all_r2, dims = 1:10)
pdx_all_r2 <- FindClusters(pdx_all_r2, resolution = 0.3)
pdx_all_r2 <- RunUMAP(pdx_all_r2, dims = 1:10)
DimPlot(pdx_all_r2, reduction = "umap")
DimPlot(pdx_all_r2, reduction = "umap", group.by = "orig.ident")
pdx_all_r2_markers = FindAllMarkers(pdx_all_r2, test.use = "MAST",
                                    only.pos = T, 
                                    min.pct = 0.2, min.cells.feature = 10)
pdx_all_r2_markers = pdx_all_r2_markers[pdx_all_r2_markers$p_val_adj < 0.01, ]
FeaturePlot(pdx_all_r2, features = c("RPL37", "VIM", "SOX2", "UXT", "KRT8", "LDHB",
                               "IFI6", "B2M", "STAT1", "SQSTM1", "QSOX1", "CD44", "CRABP1"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


#######################
# FastMNN with Seurat #
#######################

# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)


common_var_genes = unique(intersect(pdx_all_r1@assays$RNA@var.features, 
                                    pdx_all_r2@assays$RNA@var.features))
pdx_all_r1[["batch"]] <- "R1"
pdx_all_r2[["batch"]] <- "R2"

# merging both data
pdx_all_combined <- merge(pdx_all_r1, y = pdx_all_r2)

pdx_all_combined = NormalizeData(pdx_all_combined, normalization.method = "LogNormalize", scale.factor = 10000)
pdx_all_combined = FindVariableFeatures(pdx_all_combined, selection.method = "vst", nfeatures = 1000)
pdx_all_combined = ScaleData(pdx_all_combined, vars.to.regress = c("nCount_RNA", "percent.mt"))

pdx_all_combined <- RunPCA(pdx_all_combined, features = common_var_genes)

pdx_all_combined <- RunFastMNN(object.list = SplitObject(pdx_all_combined, split.by = "batch"),
                               features = common_var_genes)

pdx_all_combined <- RunUMAP(pdx_all_combined, reduction = "mnn", dims = 1:6)
pdx_all_combined <- FindNeighbors(pdx_all_combined, reduction = "mnn", dims = 1:6)
pdx_all_combined <- FindClusters(pdx_all_combined, resolution = 0.3)
DimPlot(pdx_all_combined, group.by = c("batch", "orig.ident", "RNA_snn_res.0.3"), ncol = 3)

pdx_all_markers = FindAllMarkers(pdx_all_combined, test.use = "MAST",
                                    only.pos = T, group.by = "RNA_snn_res.0.3",
                                    min.pct = 0.2, min.cells.feature = 10)

pdx_all_markers = pdx_all_markers[pdx_all_markers$p_val_adj < 0.01, ]

FeaturePlot(pdx_all_combined, features = c("VIM", "QSOX1", "CD44", "B2M", "LGALS3", "FTL",
                                     "NDRG1", "LAMTOR5", "NDUFA4", "PSMB4"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

save(pdx_all_combined, pdx_all_markers, file = "pdx_all_combined.rda")

# save seurat object and markers
new_meta_pdx_dose = new_meta
save(pdx_all_combined, pdx_all_markers, new_meta_pdx_dose, table_pdx_persister_markers,
     file = "pdx_all_combined_seurat.rda")
load("pdx_all_combined_seurat.rda")

DimPlot(pdx_all_combined, reduction = "umap", group.by = "orig.ident")


########################################################################################
# Analysis to show Dose dependent genes using pairwise diff. expressed genes and mfuzz #
########################################################################################

# subset rep1 and rep2 separately
pdx_all_combined_r1 = subset(x = pdx_all_combined, subset = batch == "R1")
pdx_all_combined_r2 = subset(x = pdx_all_combined, subset = batch == "R2")

labels_r1 = as.character(unique(pdx_all_combined_r1@meta.data$orig.ident))
labels_r2 = as.character(unique(pdx_all_combined_r2@meta.data$orig.ident))


pairwise_diff_expression <-function(seurat_obj, ctrl_id = "V4") {
  # Calculates the pairwise differential expression across cells of all ids
  # To compare "conserved" and "specific" genes differentially expressed among
  # the conditions
  # input: seurat object
  
  # get labels
  labels = as.character(unique(seurat_obj@meta.data$orig.ident))
  diff_expression = list()
  for(query_id in labels) {
    for(subject_id in labels) {
      # check if they are identical, then skip
      if(query_id == subject_id) {
        next
      } else {
        query_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident == query_id, ])
        subject_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident == subject_id, ])
        data_pair = FindMarkers(seurat_obj, ident.1 = query_cells, ident.2 = subject_cells,
                                min.cells.feature = 10, min.pct = 0.2,
                                test.use = "MAST", only.pos = T)
        # filter the significant ones
        data_pair = data_pair[data_pair$p_val_adj < 0.01, ]
        name = paste(query_id, subject_id, sep = "_") 
        diff_expression[[name]] <- data_pair
      }
    }
  }
  # Adding all conditions vs V and V vs all conditions
  v_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident == ctrl_id, ])
  t_cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data$orig.ident != ctrl_id, ])
  # test all vs V
  data_pair = FindMarkers(seurat_obj, ident.1 = t_cells, ident.2 = v_cells,
                          min.cells.feature = 10, min.pct = 0.2,
                          test.use = "MAST", only.pos = T)
  # filter the significant ones
  data_pair = data_pair[data_pair$p_val_adj < 0.01, ]
  name = paste("All", ctrl_id, sep = "_") 
  diff_expression[[name]] <- data_pair
  
  # test V vs all
  data_pair = FindMarkers(seurat_obj, ident.1 = v_cells, ident.2 = t_cells,
                          min.cells.feature = 10, min.pct = 0.2,
                          test.use = "MAST", only.pos = T)
  # filter the significant ones
  data_pair = data_pair[data_pair$p_val_adj < 0.01, ]
  name = paste(ctrl_id, "All", sep = "_") 
  diff_expression[[name]] <- data_pair
  
  return(diff_expression)
}

table_diff_expression_r1 = pairwise_diff_expression(pdx_all_combined_r1)
table_diff_expression_r2 = pairwise_diff_expression(pdx_all_combined_r2, ctrl_id = "V3")

# find genes shared across all treated conditions vs vehicle
shared_v_genes_r1 = Reduce(intersect, list(rownames(table_diff_expression_r1$V4_VLD4c),
                                          rownames(table_diff_expression_r1$V4_LD3b),
                                          rownames(table_diff_expression_r1$V4_ND4a),
                                          rownames(table_diff_expression_r1$V4_HD2d)))

shared_v_genes_r2 = Reduce(intersect, list(rownames(table_diff_expression_r2$V3_LD4b),
                                           rownames(table_diff_expression_r2$V3_LD4b),
                                           rownames(table_diff_expression_r2$V3_ND5a),
                                           rownames(table_diff_expression_r2$V3_HD4d)))

shared_all_genes_r2 = Reduce(intersect, list(rownames(table_diff_expression_r2$VLD3c_V3),
                                             rownames(table_diff_expression_r2$LD4b_V3),
                                             rownames(table_diff_expression_r2$ND5a_V3),
                                             rownames(table_diff_expression_r2$HD4d_V3)))



# these genes can be used for monocle and mfuzz for example
genes_diff_expressed_r1 = unique(Reduce(union, list(rownames(table_diff_expression_r1$VLD4c_V4),
                                                 rownames(table_diff_expression_r1$LD3b_V4),
                                                 rownames(table_diff_expression_r1$ND4a_V4),
                                                 rownames(table_diff_expression_r1$HD2d_V4),
                                                 rownames(table_diff_expression_r1$All_V4),
                                                 rownames(table_diff_expression_r1$V4_All))))

genes_diff_expressed_r2 = unique(Reduce(union, list(rownames(table_diff_expression_r2$VLD3c_V3),
                                                    rownames(table_diff_expression_r2$LD4b_V3),
                                                    rownames(table_diff_expression_r2$ND5a_V3),
                                                    rownames(table_diff_expression_r2$HD4d_V3),
                                                    rownames(table_diff_expression_r2$All_V3),
                                                    rownames(table_diff_expression_r2$V3_All))))

save(table_diff_expression_r1, table_diff_expression_r2,
     genes_diff_expressed_r1, genes_diff_expressed_r2, shared_all_genes_r1,
     shared_all_genes_r2, file = "diff_expressed_genes_dose.rda")



genes_diff_r1_r2 = intersect(genes_diff_expressed_r1, genes_diff_expressed_r2)



# these genes can be used for zavit or mfuzz: Union of diff. expressed genes
all_genes_diff_expressed_r1 = unique(Reduce(union, list(rownames(table_diff_expression_r1$VLD4c_V4),
                                                 rownames(table_diff_expression_r1$LD3b_V4),
                                                 rownames(table_diff_expression_r1$ND4a_V4),
                                                 rownames(table_diff_expression_r1$HD2d_V4),
                                                 #rownames(table_diff_expression_r1$All_V4),
                                                 rownames(table_diff_expression_r1$HD2d_VLD4c),
                                                 rownames(table_diff_expression_r1$HD2d_LD3b),
                                                 rownames(table_diff_expression_r1$HD2d_ND4a),
                                                 rownames(table_diff_expression_r1$ND4a_LD3b),
                                                 rownames(table_diff_expression_r1$ND4a_VLD4c),
                                                 rownames(table_diff_expression_r1$LD3b_VLD4c),
                                                 #rownames(table_diff_expression_r1$V4_All),
                                                 rownames(table_diff_expression_r1$V4_))))

# these genes can be used for zavit or mfuzz: Union of diff. expressed genes
all_genes_diff_expressed_r2 = unique(Reduce(union, list(rownames(table_diff_expression_r2$VLD3c_V3),
                                                        rownames(table_diff_expression_r2$LD4b_V3),
                                                        rownames(table_diff_expression_r2$ND5a_V3),
                                                        rownames(table_diff_expression_r2$HD4d_V3),
                                                        rownames(table_diff_expression_r2$All_V3),
                                                        rownames(table_diff_expression_r2$HD4d_VLD3c),
                                                        rownames(table_diff_expression_r2$HD4d_LD4b),
                                                        rownames(table_diff_expression_r2$HD4d_ND5a),
                                                        rownames(table_diff_expression_r2$ND5a_VLD3c),
                                                        rownames(table_diff_expression_r2$ND5a_LD4b),
                                                        rownames(table_diff_expression_r2$LD4b_VLD3c),
                                                        rownames(table_diff_expression_r2$V3_All))))

all_genes_diff_expressed = unique(intersect(all_genes_diff_expressed_r1, all_genes_diff_expressed_r2))

# average expression for each time point
labels = as.character(unique(pdx_all_combined_r1@meta.data$orig.ident))
list_means = list()
for(label in labels) {
  cells = rownames(pdx_all_combined_r1@meta.data[pdx_all_combined_r1@meta.data$orig.ident == label, ])
  avg_cells = rowMeans(as.matrix(pdx_all_combined_r1@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}

avg_data_r1 = do.call(cbind, list_means)

# Compute avg data for R2
labels = as.character(unique(pdx_all_combined_r2@meta.data$orig.ident))
list_means = list()
for(label in labels) {
  cells = rownames(pdx_all_combined_r2@meta.data[pdx_all_combined_r2@meta.data$orig.ident == label, ])
  avg_cells = rowMeans(as.matrix(pdx_all_combined_r2@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}

avg_data_r2 = do.call(cbind, list_means)


#### MFUZZ GENES CLUSTERING ####
library(Mfuzz)

avg_data_genes_r1 = as.matrix(avg_data_r1[all_genes_diff_expressed, ])
eset_r1 = new('ExpressionSet', exprs = avg_data_genes_r1)

avg_data_genes_r2 = as.matrix(avg_data_r2[all_genes_diff_expressed, ])
eset_r2 = new('ExpressionSet', exprs = avg_data_genes_r2)

eset.s_r1 = standardise(eset_r1)
mestimate(eset.s_r1)

eset.s_r2 = standardise(eset_r2)
mestimate(eset.s_r2)


eset_clust_r1 <- mfuzz(eset.s_r1, c = 7, m = 1.8)
eset_clust_r2 <- mfuzz(eset.s_r2, c = 7, m = 1.8)
#pdf(file = "mfuzz_filtered2.pdf", width = 13, height = 10)
mfuzz.plot2(eset.s_r2, cl = eset_clust_r1, mfrow=c(4,5), 
            time.labels = colnames(avg_data_genes_r1), centre = T)

mfuzz.plot2(eset.s_r2, cl = eset_clust_r2, mfrow=c(4,5), 
            time.labels = colnames(avg_data_genes_r2), centre = T)

dev.off()

# clusters of interest
#clusters = c(1, 2, 5, 7, 12, 14, 15)
clusters = c(1, 2, 3, 4, 5, 6, 7)

genes_clusters_r1 = names(eset_clust_r1$cluster[eset_clust_r1$cluster %in% clusters])
genes_clusters_r2 = names(eset_clust_r2$cluster[eset_clust_r2$cluster %in% clusters])

# get only genes in membership data that correspond to clusters I am interested
cluster_members_r1 = eset_clust_r1$membership[genes_clusters_r1, ]
cluster_members_r2 = eset_clust_r2$membership[genes_clusters_r2, ]

# membership threshold
m_th = 0.3
# filter genes in the selected clusters with high membership values
cluster_members_r1 = cluster_members_r1[rowSums(cluster_members_r1 > m_th) >= 1, ]
cluster_genes_data_r1 = avg_data_genes_r1[rownames(cluster_members_r1), ]

cluster_members_r2 = cluster_members_r2[rowSums(cluster_members_r2 > m_th) >= 1, ]
cluster_genes_data_r2 = avg_data_genes_r2[rownames(cluster_members_r2), ]

# geting the original expression as z-score for selected genes
cluster_genes_data_scaled_r1 = as.data.frame(t(scale(t(cluster_genes_data_r1))))
cluster_genes_data_scaled_r2 = as.data.frame(t(scale(t(cluster_genes_data_r2))))

# change colnames to be numbers
colnames(cluster_genes_data_scaled_r1) <- c("V_R1", "D1_R1", "D2_R1", "D3_R1", "D4_R1")
colnames(cluster_genes_data_scaled_r2) <- c("V_R2", "D1_R2", "D2_R2", "D3_R2", "D4_R2")

# add gene names
cluster_genes_data_scaled_r1$gene <- rownames(cluster_genes_data_scaled_r1)
cluster_genes_data_scaled_r2$gene <- rownames(cluster_genes_data_scaled_r2)

# adding the cluster information on genes
selected_clusters_r1 = as.character(eset_clust_r1$cluster[cluster_genes_data_scaled_r1$gene])
selected_clusters_r2 = as.character(eset_clust_r2$cluster[cluster_genes_data_scaled_r2$gene])

# add the clusters to data
cluster_genes_data_scaled_r1$cluster <- selected_clusters_r1
cluster_genes_data_scaled_r2$cluster <- selected_clusters_r2

# make it in long format to plot
cluster_genes_data_scaled_melted_r1 = reshape2::melt(cluster_genes_data_scaled_r1)
cluster_genes_data_scaled_melted_r2 = reshape2::melt(cluster_genes_data_scaled_r2)


ggplot(cluster_genes_data_scaled_melted_r2, aes(x = variable, y = value)) +
  geom_smooth(aes(group = gene), se=F, colour="gray", size=0.6, alpha = 0.4) +
  geom_smooth(aes(group = cluster), se=F, colour="black", size=1, method = "loess") +
  facet_wrap(~ cluster, ncol = 3) +
  xlab("") + ylab("Normalized expression") +
  theme_minimal() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=10, colour="black"),
        strip.background=element_rect(colour="white", fill="white"))

# getting genes that show consistent patterns with clusters in both replicates, 
# i.e up (dose) and down clusters (v), so show the intersection of them

# down genes are in clusters 1,6 for R1 and 1,3 for R2: Cluster numbers change every time!!
down_r1 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("1", "6"), ]$gene)
down_r2 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("1", "3"), ]$gene)
down_all = unique(intersect(down_r1, down_r2))

# down genes are in clusters 1,6 for R1 and 1,3 for R2: Cluster numbers change every time!!
up_r1 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("2", "3", "5"), ]$gene)
up_r2 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("2", "6", "7"), ]$gene)
up_all = unique(intersect(up_r1, up_r2))


save(cluster_genes_data_scaled_melted_r1, cluster_genes_data_scaled_melted_r2,
     eset_clust_r1, eset_clust_r2, file = "mfuzz_analysis_dose.rda")

# merge data frames
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

avg_data_genes_r1_scaled = t(scale(t(avg_data_genes_r1), center = T))
avg_data_genes_r2_scaled = t(scale(t(avg_data_genes_r2), center = T))

avg_data_genes_r1_r2_scaled = merge_exp(list(avg_data_genes_r1_scaled[c(down_all, up_all), ],
                                             avg_data_genes_r2_scaled[c(down_all, up_all), ]))

avg_data_genes_r1_r2 = merge_exp(list(avg_data_genes_r1[c(down_all, up_all), ],
                                      avg_data_genes_r2[c(down_all, up_all), ]))

colnames(avg_data_genes_r1_r2_scaled) <- c("V_1", "D1_1", "D2_1", "D3_1", "D4_1",
                                           "V_2", "D1_2", "D2_2", "D3_2", "D4_2")
colnames(avg_data_genes_r1_r2) <- c("V_1", "D1_1", "D2_1", "D3_1", "D4_1",
                                    "V_2", "D1_2", "D2_2", "D3_2", "D4_2")


## Gene order
down_r1_c1 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("1"), ]$gene)
down_r1_c6 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("6"), ]$gene)

up_r1_c3 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("3"), ]$gene)
up_r1_c5 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("5"), ]$gene)
up_r1_c2 = unique(cluster_genes_data_scaled_melted_r1[cluster_genes_data_scaled_melted_r1$cluster %in% c("2"), ]$gene)

down_r2_c3 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("3"), ]$gene)
down_r2_c1 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("1"), ]$gene)

up_r2_c6 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("6"), ]$gene)
up_r2_c7 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("7"), ]$gene)
up_r2_c2 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("2"), ]$gene)

up_r1_all = unique(c(up_r1_c3, up_r1_c5, up_r1_c2))
up_r2_all = unique(c(up_r2_c6, up_r2_c7, up_r2_c2))

up_all_test = intersect(up_r1_all, up_r2_all)

down_r1_all = unique(c(down_r1_c1, down_r1_c6))
down_r2_all = unique(c(down_r2_c3, down_r2_c1))

down_all_test = intersect(down_r1_all, down_r2_all)

avg_data_genes_r1_r2_scaled = merge_exp(list(avg_data_genes_r1_scaled[c(down_all_test, up_all_test), ],
                                             avg_data_genes_r2_scaled[c(down_all_test, up_all_test), ]))
colnames(avg_data_genes_r1_r2_scaled) <- c("V_1", "D1_1", "D2_1", "D3_1", "D4_1",
                                           "V_2", "D1_2", "D2_2", "D3_2", "D4_2")



up_r2 = unique(cluster_genes_data_scaled_melted_r2[cluster_genes_data_scaled_melted_r2$cluster %in% c("2", "6", "7"), ]$gene)

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

genes_to_show = c("LDHB", "ENO1", "KRT8", "GAPDH", "PARP1", "PEG10", "PCNA", "CRABP2",
                  "LDHA", "CYC1", "HMGCR", "CDK4", "XPO1", "XRCC6",
                  "STAT2", "IFI27", "MX1", "JAK1", "PARP14", "NFKBIA",
                  "SOX2", "NPC2", "SQSTM1", "QSOX1", "KLF4", "NFE2L2",
                  "PHLDA2", "PLK2", "FTL", "FTH1", "CDKN2A", "DCXR")


genes_order = c(unique(down_all_test), unique(up_all_test))
avg_data_genes_r1_r2_scaled_order = avg_data_genes_r1_r2_scaled[genes_order, c(1,6,7,2,3,8,4,9,5,10)]

gene_annotation = show_genes(genes_to_show, avg_data_genes_r1_r2_scaled_order)


Heatmap(avg_data_genes_r1_r2_scaled_order, cluster_rows = F, cluster_columns = F, 
              show_row_names = F, show_column_dend = F,
              row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
              show_column_names = T, border = T,
              show_row_dend = F, show_heatmap_legend = T, 
              col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
              right_annotation = gene_annotation,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                          title = "Expression",
                                          labels_gp = gpar(fontsize = 7),
                                          legend_height = unit(2, "cm"),
                                          border = "black"))
draw(ht1)

save(eset, eset.s, cluster_genes_data_scaled_melted, cluster_genes_data_scaled,
     avg_data_genes, avg_data_genes_scaled, file = "mfuzz_clusters.rda")


# rename clusters 
new_meta = pdx_all_combined_metadata
new_meta$new_id <- stringr::str_replace_all(new_meta$orig.ident, 
                                            c("V4" = "V", "V3" = "V",
                                              "VLD3c" = "D1", "VLD4c" = "D1",
                                              "LD3b" = "D2", "LD4b" = "D2",
                                              "ND5a" = "D3", "ND4a" = "D3",
                                              "HD2d" = "D4", "HD4d" = "D4"))

new_meta$new_id <- factor(new_meta$new_id,
                          levels = c("V", "D1", "D2", "D3", "D4"))

new_meta$all <- c("All")

### rename clusters to make them ordered by vehicle frequency
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$seurat_clusters, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$new_cluster_id, 
                                                    c("A" = "0", "B" = "4", "C" = "1",
                                                      "D" = "2", "E" = "3"))

table_pdx_persister_markers = pdx_all_markers

### rename clusters to make them ordered by vehicle frequency
table_pdx_persister_markers$cluster <- stringr::str_replace_all(table_pdx_persister_markers$cluster, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
table_pdx_persister_markers$cluster <- stringr::str_replace_all(table_pdx_persister_markers$cluster, 
                                                    c("A" = "0", "B" = "4", "C" = "1",
                                                      "D" = "2", "E" = "3"))

write.table(table_pdx_persister_markers, file = "pdx_persistence_markers.tsv", row.names = F,
            col.names = T, quote = F, sep = "\t")

#######################
# Bubble plot markers #
#######################


genes_to_show = c("SOX17", "PAX8", "WT1", "KRT8", "CDKN2A", "IFI6", "IFI27", "IFITM3",
                  "VIM", "CD44", "SMAD3", "DDIT4", "CITED2", "SGK1", "PGK1", "PFKP", "ACLY",
                  "SCD", "MBOAT7", "SREBF1", "CTPS1", "TYMS", "GPX4", "STEAP3", "FTL",
                  "SOD2", "GGH", "GSTM3",
                  "PRDX1", "CYP1B1", "TXNRD1", "NQO1", "SQSTM1","NFE2L2")

## dose
avg_expression_dose = DotPlot(pdx_all_combined, features = genes_to_show, scale = F)
avg_expression_dose = avg_expression_dose$data
avg_expression_dose$experiment = "Persistence"

# rename clusters to my new names
### rename clusters to make them ordered by vehicle frequency
avg_expression_dose$new_cluster_id <- stringr::str_replace_all(avg_expression_dose$id, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
avg_expression_dose$new_cluster_id <- stringr::str_replace_all(avg_expression_dose$new_cluster_id, 
                                                    c("A" = "0", "B" = "4", "C" = "1",
                                                      "D" = "2", "E" = "3"))

## Adapted
avg_expression_adapt = DotPlot(pdx_all_adapted_combined_no_batch, features = genes_to_show, scale = F)
avg_expression_adapt = avg_expression_adapt$data
avg_expression_adapt$experiment = "Resistance"

### rename clusters to make them ordered by vehicle frequency
avg_expression_adapt$new_cluster_id <- stringr::str_replace_all(avg_expression_adapt$id, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
avg_expression_adapt$new_cluster_id <- stringr::str_replace_all(avg_expression_adapt$new_cluster_id, 
                                                    c("A" = "6", "B" = "8", "C" = "7",
                                                      "D" = "9", "E" = "5"))

avg_expression_dose_adapt = rbind(avg_expression_dose, avg_expression_adapt)

# now have to select clusters and gene and scale by gene together
scale_genes_individual <- function(data_frame_long) {
  genes = unique(data_frame_long$features.plot)
  new_data_frame = data.frame()
  for(gene in genes) {
    d = data_frame_long[data_frame_long$features.plot == gene, ]
    scaled_exp_together = scale(d$avg.exp.scaled, center = T)
    d$scaled_exp_together = scaled_exp_together
    new_data_frame = rbind(new_data_frame, d)
  }
  return(new_data_frame)
}

avg_expression_dose_adapt_scaled = scale_genes_individual(avg_expression_dose_adapt)
avg_expression_dose_adapt_scaled$new_cluster_id <- factor(avg_expression_dose_adapt_scaled$new_cluster_id,
                                                          levels = c("9","8","7","6","5",
                                                                     "4","3","2","1","0"))

ggplot(avg_expression_dose_adapt_scaled,
       aes(x = features.plot, y = new_cluster_id, size = pct.exp, fill = scaled_exp_together)) + 
  #geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 21) +
  #scale_size_continuous(range = c(1, 4)) +
  #scale_fill_gradientn(colours = pals::coolwarm(100)) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F") +
  facet_grid(experiment ~ ., scales = "free_y") +
  ylab("") +
  xlab("") +
  theme_bw() +
  #coord_fixed(ratio=3) +
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")

ggplot(avg_expression_dose_adapt_scaled,
       aes(x = new_cluster_id, y = features.plot, size = pct.exp, fill = scaled_exp_together)) + 
  #geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 21) +
  #scale_size_continuous(range = c(1, 4)) +
  #scale_fill_gradientn(colours = pals::coolwarm(100)) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F") +
  facet_grid(. ~ experiment, scales = "free_x") +
  ylab("") +
  xlab("") +
  theme_bw() +
  #coord_fixed(ratio=3) +
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")


save_plot("pdx_markers_bubble.pdf", pdx_markers_bubble, base_height = 4, base_width = 10)


##########################
# Plot Cluster frequency #
##########################

new_meta = pdx_all_combined@meta.data
new_meta$new_id <- stringr::str_replace_all(new_meta$orig.ident, 
                                            c("V4" = "V", "V3" = "V",
                                              "VLD3c" = "D1", "VLD4c" = "D1",
                                              "LD3b" = "D2", "LD4b" = "D2",
                                              "ND5a" = "D3", "ND4a" = "D3",
                                              "HD2d" = "D4", "HD4d" = "D4"))
new_meta$new_id <- factor(new_meta$new_id,
                          levels = c("V", "D1", "D2", "D3", "D4"))

# to save for GEO with replicate ids
new_meta$new_id <- stringr::str_replace_all(new_meta$orig.ident, 
                                            c("V4" = "V1", "V3" = "V2",
                                              "VLD4c" = "D1.1", "VLD3c" = "D1.2",
                                              "LD3b" = "D2.1", "LD4b" = "D2.2",
                                              "ND4a" = "D3.1", "ND5a" = "D3.2",
                                              "HD2d" = "D4.1", "HD4d" = "D4.2"))



new_meta$all <- c("All")

### rename clusters to make them ordered by vehicle frequency
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$seurat_clusters, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$new_cluster_id, 
                                                    c("A" = "0", "B" = "4", "C" = "1",
                                                      "D" = "2", "E" = "3"))

# write metadata for GEO and replace the orig ident and cell suffix by defined D1, D2, etc.
write.table(new_meta, file = "pdx_persisters_metadata_g1.tsv", row.names = T, col.names = T,
            quote = F, sep = "\t")


####################################################
# Plot barplot of frequencies of cells in clusters #
####################################################

# The order really matters! first clusters then condition
freq_clusters = new_meta %>% group_by(new_cluster_id, new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_clusters) <- c("cluster", "condition", "n", "freq")

freq_all = new_meta %>% group_by(new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
freq_all$cluster <- "All"
freq_all = freq_all[, c(4,1,2,3)]
colnames(freq_all) <- c("cluster", "condition", "n", "freq")

freq_cluster_all = rbind(as.data.frame(freq_all), as.data.frame(freq_clusters))

#############
# New plot 
####

freq_cluster = new_meta %>% group_by(new_id, new_cluster_id) %>%
  dplyr::summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_cluster) <- c("condition", "cluster", "n", "freq")

freq_cluster$cluster <- factor(freq_cluster$cluster,
                               levels = c("0", "1", "2","3","4"))

freq_cluster$condition <- factor(freq_cluster$condition,
                                 levels = c("V","D1","D2","D3","D4"))

#freq_cluster$condition <- factor(freq_cluster_all$condition,
#                                  levels = rev(c("C", "T5", "T10", "T20", "T40")))

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

#freq_cluster_all$cluster <- factor(freq_cluster_all$cluster,
#                                   levels = c("4", "3", "2", "1", "0", "All"))
pdx_dose_freq_cluster_plot = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  scale_fill_manual(values = cluster_colors,
                    limits = c("0","1", "2", "3", "4")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white"))

############### freq plot ends here ############

######################
# Plot Umap sample   #
######################

pdx_all_combined_umap_sample = DimPlot(pdx_all_combined, reduction = "umap", group.by = "orig.ident")
# fix cluster ids as frequency plot
pdx_all_combined_umap_sample$data$new_sample_id <- pdx_all_combined_umap_sample$data$orig.ident
pdx_all_combined_umap_sample$data$new_sample_id <- stringr::str_replace_all(pdx_all_combined_umap_sample$data$new_sample_id,
                                                                            c("V4" = "V1", "V3" = "V2",
                                                                              "VLD4c" = "D1.1", "VLD3c" = "D1.2",
                                                                              "LD3b" = "D2.1", "LD4b" = "D2.2",
                                                                              "ND4a" = "D3.1", "ND5a" = "D3.2",
                                                                              "HD2d" = "D4.1", "HD4d" = "D4.2"))

pdx_all_combined_umap_sample$data$new_sample_id <- factor(pdx_all_combined_umap_sample$data$new_sample_id,
                                                       levels = c("V1", "V2", "D1.1", "D1.2",
                                                                  "D2.1", "D2.2", "D3.1", "D3.2",
                                                                  "D4.1", "D4.2"))

sample_colors = c("#8da3a6", "#6b7d7f", "#cae8c1", "#bbd3b4", "#81c07a", 
                  "#709e69", "#30893b", "#296d2f", "#0f4615", "#0a2d0d")

umap_plot_dose_sample = ggplot(pdx_all_combined_umap_sample$data, 
                               aes(y = UMAP_2, x = UMAP_1, color = new_sample_id)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = sample_colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Sample", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))


ggsave("umap_plot_dose_sample.pdf", plot = umap_plot_dose_sample, 
       width = 3.5, height = 3, units = c("cm"),
       dpi = 300, limitsize = T, scale = 3)

######################
# Plot Umap clusters #
######################

pdx_all_combined_umap = DimPlot(pdx_all_combined, reduction = "umap")
pdx_all_combined_umap$data$new_cluster_id <- new_meta$new_cluster_id

colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

umap_plot_dose = ggplot(pdx_all_combined_umap$data, aes(y = UMAP_2, x = UMAP_1, color = new_cluster_id)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Cluster", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))

ggsave("umap_plot_dose.pdf", plot = umap_plot_dose, 
       width = 3.5, height = 3, units = c("cm"),
       dpi = 300, limitsize = T, scale = 3)