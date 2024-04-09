

library(Seurat)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

setwd("/path/pc9_resistance/")

# get cells and CMOs table
cells_cmo = read.csv("pc9_R1/pc9_R1_08/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv")
# filter only those with one CMO
cells_cmo = cells_cmo[cells_cmo$num_features == 1, ]
barcodes = cells_cmo$cell_barcode
cmos = cells_cmo$feature_call

# CMO correspondence
cmo_doses = stringr::str_replace_all(cmos, c("CMO301" = "C", "CMO302" = "T0008", "CMO303" = "T0016",
                                             "CMO304" = "T0032", "CMO305" = "T0075", "CMO306" = "T0150",
                                             "CMO307" = "T0300", "CMO308" = "T0600", "CMO309" = "T1.2",
                                             "CMO310" = "P1.2"))

barcodes_cmo = paste(cells_cmo$cell_barcode, cmo_doses, sep = "_")

# Read 10X matrix
pc9_R1 = Read10X("pc9_R1/pc9_R1_08/outs/per_sample_outs/PC9_R1/count/sample_filtered_feature_bc_matrix")
# reorder the matrix based on the barcodes order
pc9_R1$`Gene Expression` = pc9_R1$`Gene Expression`[, barcodes]
# Rename cells using barcodes_cmo correspondence
colnames(pc9_R1$`Gene Expression`) <- barcodes_cmo

pc9_R1 = CreateSeuratObject(counts = pc9_R1$`Gene Expression`, project = "pc9_R1",
                             min.cells = 12, min.features = 500, names.delim = "_",
                             names.field = 2)
# percent mitochondrial genes
pc9_R1[["percent.mt"]] <- PercentageFeatureSet(pc9_R1, pattern = "^MT-")
VlnPlot(pc9_R1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# filtering data, by transcripts and mitocondrial
# threshold is mean + 2 standard deviations 
up_count_threshold = mean(pc9_R1@meta.data$nCount_RNA) + (2*(sd(pc9_R1@meta.data$nCount_RNA)))
pc9_R1 = subset(pc9_R1, subset = nCount_RNA > 700 & nCount_RNA < up_count_threshold & percent.mt < 20)

# Normalize data
pc9_R1 = NormalizeData(pc9_R1, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes
pc9_R1 = FindVariableFeatures(pc9_R1, selection.method = "vst", nfeatures = 2000)
# scale data 
pc9_R1 = ScaleData(pc9_R1, vars.to.regress = c("nCount_RNA", "percent.mt"))
# PCA
pc9_R1 = RunPCA(pc9_R1, features = VariableFeatures(object = pc9_R1))
DimPlot(pc9_R1, reduction = "pca")
ElbowPlot(pc9_R1)
DimPlot(pc9_R1, reduction = "pca", group.by = "orig.ident")


# Clustering
pc9_R1 = FindNeighbors(pc9_R1, dims = 1:10)
pc9_R1 = FindClusters(pc9_R1, resolution = 0.5)
pc9_R1 = RunUMAP(pc9_R1, dims = 1:10, reduction = "pca")
DimPlot(pc9_R1, reduction = "umap", group.by = "orig.ident")

#######################
# Subsetting G1 cells #
#######################

s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

pc9_R1 = CellCycleScoring(pc9_R1, s.features = s.genes, 
                            g2m.features = g2m.genes, set.ident = TRUE)
# check umap, pca on cell cyle
DimPlot(pc9_R1, reduction = "pca", group.by = "Phase")
DimPlot(pc9_R1, reduction = "umap", group.by = "Phase")

# plot frequency of G1, G2M, S cells by sample

# control cells were mostly cycling rapidly, while treated ones were mostly G1!
g1_cells = rownames(pc9_R1@meta.data[pc9_R1@meta.data$Phase == "G1", ])
# subset G1 cells only
pc9_R1_g1 = subset(pc9_R1, cells = g1_cells)

# Normalize data
pc9_R1_g1 = NormalizeData(pc9_R1_g1, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes
pc9_R1_g1 = FindVariableFeatures(pc9_R1_g1, selection.method = "vst", nfeatures = 2000)
# scale data 
pc9_R1_g1 = ScaleData(pc9_R1_g1, vars.to.regress = c("nCount_RNA", "percent.mt"))
# PCA
pc9_R1_g1 = RunPCA(pc9_R1_g1, features = VariableFeatures(object = pc9_R1_g1))
DimPlot(pc9_R1_g1, reduction = "pca")
DimPlot(pc9_R1_g1, reduction = "pca", group.by = "orig.ident")

ElbowPlot(pc9_R1_g1)

# Clustering
pc9_R1_g1 = FindNeighbors(pc9_R1_g1, dims = 1:10)
pc9_R1_g1 = FindClusters(pc9_R1_g1, resolution = 0.2)
pc9_R1_g1 = RunUMAP(pc9_R1_g1, dims = 1:10, reduction = "pca")
DimPlot(pc9_R1_g1, reduction = "umap")
DimPlot(pc9_R1_g1, reduction = "umap", group.by = "orig.ident")

# remove the annoying small cluster not meaningful
pc9_R1_g1 = subset(pc9_R1_g1, idents = c(0,1,2,3,4))
DimPlot(pc9_R1_g1, reduction = "umap")


######################
# Plot Umap Sample   #
######################

pc9_g1_umap_sample = DimPlot(pc9_R1_g1, reduction = "umap", group.by = "orig.ident")

# including persisters
colors = c("#8DA3A6", "#b74980", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
           "#D0B541", "#E67F33", "#CE2220")

pc9_umap_plot_sample = ggplot(pc9_g1_umap_sample$data, aes(y = UMAP_2, x = UMAP_1, color = orig.ident)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Treatment", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10))

######################
# Plot umap clusters #
######################

# rename clusters to put in order of more resistant
pc9_R1_g1 <- RenameIdents(object = pc9_R1_g1,
                              "0" = "2",
                              "1" = "3",
                              "2" = "4",
                              "3" = "1",
                              "4" = "0")
# assign new order
pc9_R1_g1@active.ident <- factor(x = pc9_R1_g1@active.ident, levels = c("0", "1", "2", "3", "4"))

pc9_R1_umap_cluster = DimPlot(pc9_R1_g1, reduction = "umap") 

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

pc9_umap_plot_cluster = ggplot(pc9_R1_umap_cluster$data, aes(y = UMAP_2, x = UMAP_1, color = ident)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = cluster_colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Cluster", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10))


##########################
# Plot Cluster frequency #
##########################

new_meta = pc9_R1_g1@meta.data
new_meta$new_clusters <- pc9_R1_umap_cluster$data$ident

########################################################
# Replot barplot with frequencies and renamed clusters #
########################################################

library(dplyr)
library(ggplot2)

freq_cluster = new_meta %>% group_by(orig.ident, new_clusters) %>%
  dplyr::summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_cluster) <- c("condition", "cluster", "n", "freq")

freq_cluster$cluster <- factor(freq_cluster$cluster,
                               levels = c("0", "1", "2", "3", "4"))

freq_cluster$condition <- factor(freq_cluster$condition,
                                 levels = c("C","P1.2","T0008","T0016","T0032",
                                            "T0075","T0150","T0300","T0600", "T1.2"))

colors = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
           "#D0B541", "#E67F33", "#CE2220")

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

pc9_freq_cluster_plot = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  scale_fill_manual(values = cluster_colors,
                    limits = c("0","1", "2", "3", "4")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pc9_umap_plot_sample + pc9_umap_plot_cluster + pc9_freq_cluster_plot


################
# Find markers #
################
pc9_R1_g1_markers = FindAllMarkers(pc9_R1_g1, test.use = "MAST",
                                   only.pos = T, 
                                   min.pct = 0.2, min.cells.feature = 10)
# filter pvalue
pc9_R1_g1_markers = pc9_R1_g1_markers[pc9_R1_g1_markers$p_val_adj < 0.01, ]


pc9_R1_g1_markers_Per = FindMarkers(pc9_R1_g1, test.use = "MAST", ident.1 = "2",
                                ident.2 = "3",
                                   only.pos = T, 
                                   min.pct = 0.2, min.cells.feature = 10)
pc9_R1_g1_markers_Per = pc9_R1_g1_markers_Per[pc9_R1_g1_markers_Per$p_val_adj < 0.01, ]

FeaturePlot(pc9_R1_g1, features = c("CD44", "ITGA6", "LAMC2"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

VlnPlot(pc9_R1_g1, features = c("CD44", "CD24", "GLS"), group.by = "seurat_clusters")
cluster_markers = list(c0_markers, c1_markers, c2_markers,
                       c3_markers, c4_markers)

##################################################################
# Assigning the new identity to meta data to do diff. expression #
##################################################################

pc9_R1_g1@meta.data$new_id <- pc9_R1_g1@active.ident
Idents(pc9_R1_g1) <- pc9_R1_g1@meta.data$'new_id'

# Find markers
pc9_R1_g1_markers_new = FindAllMarkers(pc9_R1_g1, test.use = "MAST",
                                       only.pos = T, group.by = "new_id",
                                       min.pct = 0.2, min.cells.feature = 10)
# filter pvalue
pc9_R1_g1_markers_new = pc9_R1_g1_markers_new[pc9_R1_g1_markers_new$p_val_adj < 0.01, ]

FeaturePlot(pc9_R1_g1, 
            features = c("CD24", "CD44", "KRT6A", "EGFR", "AXL", "CDKN3", "MKI67", "PTTG1", "TUBA1B"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

save(pc9_R1, pc9_R1_g1, pc9_R1_g1_markers, pc9_R1_g1_markers_new, file = "pc9_R1_seurat.rda")
write.table(pc9_R1_g1_markers_new, file = "pc9_markers.txt", quote = F, col.names = T,
            row.names = T, sep = "\t")
load("pc9_R1_seurat.rda")

# save table for GEO
geo_pc9 = as.matrix(pc9_R1@assays$RNA@counts)
write.table(geo_pc9, file = "pc9_all_C_to_T1.2.csv", sep = ",", col.names = T, row.names = T,
            quote = F)

# save metadata for GEO
write.table(pc9_R1_g1@meta.data, file = "pc9_metadata_g1.tsv", col.names = T, row.names = T,
            sep = "\t", quote = F)

####################################################
# Plot barplot of frequencies of cells in clusters #
####################################################

library(dplyr)
library(ggplot2)

# The order really matters, first clusters then condition
meta_data_pc9 = pc9_R1_g1@meta.data
freq_clusters = meta_data_pc9 %>% group_by(seurat_clusters, orig.ident) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_clusters) <- c("cluster", "condition", "n", "freq")

freq_clusters = new_meta %>% group_by(seurat_clusters, new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_clusters) <- c("cluster", "condition", "n", "freq")

freq_all = new_meta %>% group_by(new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
freq_all$cluster <- "All"
freq_all = freq_all[, c(4,1,2,3)]
colnames(freq_all) <- c("cluster", "condition", "n", "freq")

freq_cluster_all = rbind(as.data.frame(freq_all), as.data.frame(freq_clusters))

freq_cluster_all$condition <- factor(freq_cluster_all$condition)
freq_cluster_all$cluster <- factor(freq_cluster_all$cluster)
freq_cluster_all$cluster <- factor(freq_cluster_all$cluster,
                                   levels = c("All", "0", "1", "2", "3", "4"))

freq_cluster_all$condition <- factor(freq_cluster_all$condition,
                                     levels = rev(c("C", "T5", "T10", "T20", "T40")))
freq_clusters$condition <- factor(freq_clusters$condition,
                                 levels = rev(c("C", "T003", "T006", "T0125", "T0250",
                                                "T05", "T1", "T2", "T4")))

colors = c("#8DA3A6","#4E78C4","#57A2AC","#7EB875","#D0B541")

frequency_plot = ggplot(freq_clusters, aes(y = freq, x = cluster, fill = condition)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8)) +
  theme(strip.background =element_rect(fill="white"))
dev.off()


##########################
# Dotplot of the markers #
##########################

genes_to_show = c("KRT6A", "CD44", "CD24", "TMSB4X", "LCN2", "NFKBIA", "HLA-A", 
                  "FTH1", "AKR1B1", "AKR1C1", "SQSTM1", 
                  "TXNRD1", "TKT", "GCLM",
                  "CDKN3", "MKI67", "PTTG1", "TUBA1B",
                  "SCD", "G6PD", "SLC7A11", "GLS", "GSTK1", "MGST3", "TSTD1")


avg_expression = DotPlot(pc9_R1_g1, features = genes_to_show)
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))

avg_expression$id <- factor(avg_expression$id, levels = c(4,3,2,1,0))

ggplot(avg_expression,
       aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  geom_point(shape = 21) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F") +
  ylab("") +
  xlab("") +
  theme_minimal() +
  #coord_fixed(ratio=3) +
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")



#######################################
# Trying to analyze them individually #
#######################################

setwd("/path/pc9_resistance/")

load("pc9_R1_seurat.rda")
load("pc9_R1_g1_individual.rda")
source("plotting.R")
source("seurat_wrapper.R")

# split into original samples

pc9_R1_g1_split = SplitObject(pc9_R1_g1, split.by = "orig.ident")

# now cluster them individually and plot!

seurat_individual_analysis <- function(seurat_obj_ind, prefix = "pc9_R1_g1") {
  # For each seurat object in a list, do the pipeline individually to find 
  # individual clusters. Returns a list with processed object and markers
  
  # Find variable genes
  sample_id = unique(seurat_obj_ind@meta.data$orig.ident)
  sample_id = paste(prefix, sample_id, sep = "_")
  seurat_obj_ind = FindVariableFeatures(seurat_obj_ind, selection.method = "vst", nfeatures = 1500)
  # scale data 
  seurat_obj_ind = ScaleData(seurat_obj_ind, vars.to.regress = c("nCount_RNA", "percent.mt"))
  # PCA
  seurat_obj_ind = RunPCA(seurat_obj_ind, features = VariableFeatures(object = seurat_obj_ind))
  
  # Clustering
  seurat_obj_ind = FindNeighbors(seurat_obj_ind, dims = 1:6)
  seurat_obj_ind = FindClusters(seurat_obj_ind, resolution = 0.3)
  seurat_obj_ind = RunUMAP(seurat_obj_ind, dims = 1:6, reduction = "pca")
  
  # Find markers
  seurat_obj_ind_markers = FindAllMarkers(seurat_obj_ind, test.use = "MAST",
                                      only.pos = T, 
                                      min.pct = 0.2, min.cells.feature = 10)
  # filter pvalue
  seurat_obj_ind_markers = seurat_obj_ind_markers[seurat_obj_ind_markers$p_val_adj < 0.01, ]
  
  # add IDs to the output
  sample_id = unique(seurat_obj_ind@meta.data$orig.ident)
  sample_id = paste(prefix, sample_id, sep = "_")
  markers_id = paste(sample_id, "markers", sep = "_")
  seurat_out = list(seurat_obj_ind, seurat_obj_ind_markers)
  names(seurat_out) <- c(sample_id, markers_id)
  return(seurat_out)
}

# Run Seurat pipeline for all samples
pc9_R1_g1_individual = lapply(pc9_R1_g1_split, seurat_individual_analysis)
save(pc9_R1_g1_individual, pc9_R1_g1_split, file = "pc9_R1_g1_individual.rda")

# plot and collect individual umap plots: Colors list in colors.R already has the 
# cluster colors in order

plots_umap_list <- list()
plot_titles = names(pc9_R1_g1_individual)
load("plotting.R")
# excluding persisters here
for(i in 1:(length(pc9_R1_g1_individual) - 1)) {
#  print(cluster_colors_list_pc9_r1[[i]])
  umap_plot = plot_umap(pc9_R1_g1_individual[[i]][[1]],
                        cluster_colors_list_pc9_r1[[i]],
                        title = plot_titles[i])
  plots_umap_list[[i]] <- umap_plot
}

DimPlot(pc9_R1_g1_individual$T1.2$pc9_R1_g1_T1.2, reduction = "umap")
# saved as 6 x 8 in pdf in my local computer, file: pc9_R1_g1_umaps_cluster.pdf
# location: ../projects/new_lines_evolution/pc9_dabrafenib/pc9_R1
do.call("grid.arrange", c(plots_umap_list, ncol = 3))

#plot_umap(pc9_R1_g1_individual$T4$pc9_R1_g1_T4,
#          colors = cluster_colors_new[30:33])

FeaturePlot(pc9_R1_g1_individual$C$pc9_R1_g1_C, 
            features = c("CD44", "CD24", "AXL", "MKI67", "CDKN3"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

##########################################
# Bulk analysis                          #
##########################################

library(ComplexHeatmap)
library(RColorBrewer)

# get top pc1 and pc2 genes for each individual seurat object (sample)
top_pc1 = sapply(sapply(pc9_R1_g1_individual,"[[",1), get_top_pc_genes, pc = 1, n = 40)
top_pc2 = sapply(sapply(pc9_R1_g1_individual,"[[",1), get_top_pc_genes, pc = 2, n = 40)

top_pc_genes = unique(c(as.vector(top_pc1), as.vector(top_pc2)))
var_genes = pc9_R1_g1@assays$RNA@var.features

# do on all of them together
# Doing this helps with the fucking error:
Csparse_validate = "CsparseMatrix_validate"
pc9_g1_avg = AverageExpression(pc9_R1_g1, group.by = "orig.ident")
pc9_g1_avg = pc9_g1_avg$RNA

# pairwise expression correlation on top PC genes
cormM = cor(pc9_g1_avg[var_genes, ], method="spearman")
cormM = cor(pc9_g1_avg, method="spearman")
cormM = cor(pc9_g1_avg[, -c(2)], method="spearman")
cormM = cor(pc9_g1_avg[top_pc_genes, -c(2)], method="spearman")

library(dendextend)
dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T
dend <- click_rotate(dend, continue = TRUE)

pdf("heatmap_pc9_bulk.pdf")
h2 = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
             show_column_dend = F, show_row_names = T, 
             name = "Corr", row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
             cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
             column_dend_reorder = F,
             heatmap_width = unit(9, "cm"), heatmap_height = unit(8, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "Spearman",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))

draw(h2)
dev.off()

#############################################
# Finding States across individual clusters #
#############################################


# Creating an average matrix of all samples per cluster
avg_data_list <- list()
sample_id = names(pc9_R1_g1_individual)
# exclude persisters
for(i in 1:(length(pc9_R1_g1_individual) - 1)) {
  pc9_avg_data_ind = AverageExpression(pc9_R1_g1_individual[[i]][[1]],
                                        group.by = "seurat_clusters")
  pc9_avg_data_ind = pc9_avg_data_ind$RNA
  # add sample names to the cluster columns
  colnames(pc9_avg_data_ind) <- paste(sample_id[i], colnames(pc9_avg_data_ind), sep = "_")
  # add in a list
  avg_data_list[[i]] <- pc9_avg_data_ind
  # merge the dataframes 
  pc9_g1_avg_cluster = merge_exp(avg_data_list)
}

# try my own method for averages of clusters to see if there is a difference with Seurat's ?

# Get top diff. expressed genes for each condition/cluster
# get top pc1 and pc2 genes for each individual seurat object (sample)
top_diff = sapply(sapply(pc9_R1_g1_individual, "[", 2), get_top_diff_genes, n = 100)
top_diff = unique(unlist(top_diff))

# select and correlate clusters based on diff. expressed genes
# pairwise expression correlation on top PC genes
cormM = cor(pc9_g1_avg_cluster[top_diff, -c(20,21,22,23)], method="spearman")
cormM = cor(pc9_g1_avg_cluster[top_pc_genes, ], method="spearman")

cormM = cor(t(scale(t(pc9_g1_avg_cluster[var_genes, ]))), method="spearman")
cormM = cor(pc9_g1_avg_cluster[var_genes, ], method="spearman")


library(dendextend)
dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T
dend <- click_rotate(dend, continue = TRUE)

h2 = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
             show_column_dend = F, show_row_names = T, 
             name = "Corr", row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
             cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
             column_dend_reorder = F,
             heatmap_width = unit(10, "cm"), heatmap_height = unit(9, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "Corr",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))

draw(h2)
dev.off()


plots_markers_list <- list()
plot_titles = names(pc9_R1_g1_individual)
for(i in 1:length(pc9_R1_g1_individual)) {
  umap_plot_markers = FeaturePlot(pc9_R1_g1_individual[[i]][[1]],
                          features = c("HLA-C"), ncol = 3)
  plots_markers_list[[i]] <- umap_plot_markers
}

VlnPlot(pc9_R1_g1_individual$T4$pc9_R1_g1_T4, features = c("TENM2"), ncol = 1)

do.call("grid.arrange", c(plots_markers_list, ncol = 3))


