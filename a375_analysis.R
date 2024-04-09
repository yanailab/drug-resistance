

library(Seurat)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("/path/a375_resistance/")

# get cells and CMOs table
cells_cmo = read.csv("a375_R1_08/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv")
# filter only those with one CMO
cells_cmo = cells_cmo[cells_cmo$num_features == 1, ]
barcodes = cells_cmo$cell_barcode
cmos = cells_cmo$feature_call

# CMO correspondence
cmo_doses = stringr::str_replace_all(cmos, c("CMO301" = "C", "CMO302" = "T003", "CMO303" = "T006",
                                             "CMO304" = "T0125", "CMO305" = "T0250", "CMO306" = "T05",
                                             "CMO307" = "T1", "CMO308" = "T2", "CMO309" = "T4"))

barcodes_cmo = paste(cells_cmo$cell_barcode, cmo_doses, sep = "_")

# Read 10X matrix
a375_R1 = Read10X("a375_R1_08/outs/per_sample_outs/A375_R1/count/sample_filtered_feature_bc_matrix")
# reorder the matrix based on the barcodes order
a375_R1$`Gene Expression` = a375_R1$`Gene Expression`[, barcodes]
# Rename cells using barcodes_cmo correspondence
colnames(a375_R1$`Gene Expression`) <- barcodes_cmo

a375_R1 = CreateSeuratObject(counts = a375_R1$`Gene Expression`, project = "a375_R1",
                             min.cells = 9, min.features = 500, names.delim = "_",
                             names.field = 2)
# percent mitochondrial genes
a375_R1[["percent.mt"]] <- PercentageFeatureSet(a375_R1, pattern = "^MT-")
VlnPlot(a375_R1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# filtering data, by transcripts and mitocondrial
# threshold is mean + 2 standard deviations 
up_count_threshold = mean(a375_R1@meta.data$nCount_RNA) + (2*(sd(a375_R1@meta.data$nCount_RNA)))
a375_R1 = subset(a375_R1, subset = nCount_RNA > 700 & nCount_RNA < up_count_threshold & percent.mt < 20)

# Normalize data
a375_R1 = NormalizeData(a375_R1, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes
a375_R1 = FindVariableFeatures(a375_R1, selection.method = "vst", nfeatures = 2000)
# scale data 
a375_R1 = ScaleData(a375_R1, vars.to.regress = c("nCount_RNA", "percent.mt"))
# PCA
a375_R1 = RunPCA(a375_R1, features = VariableFeatures(object = a375_R1))
DimPlot(a375_R1, reduction = "pca")
ElbowPlot(a375_R1)

# Clustering
a375_R1 = FindNeighbors(a375_R1, dims = 1:10)
a375_R1 = FindClusters(a375_R1, resolution = 0.5)
a375_R1 = RunUMAP(a375_R1, dims = 1:10, reduction = "pca")
DimPlot(a375_R1, reduction = "umap", group.by = "orig.ident")

#######################
# Subsetting G1 cells #
#######################

s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

a375_R1 = CellCycleScoring(a375_R1, s.features = s.genes, 
                            g2m.features = g2m.genes, set.ident = TRUE)
# check umap, pca on cell cyle
DimPlot(a375_R1, reduction = "umap", group.by = "Phase")
# control cells were mostly cycling rapidly, while treated ones were mostly G1!
g1_cells = rownames(a375_R1@meta.data[a375_R1@meta.data$Phase == "G1", ])
# subset G1 cells only
a375_R1_g1 = subset(a375_R1, cells = g1_cells)

# Normalize data
a375_R1_g1 = NormalizeData(a375_R1_g1, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes
a375_R1_g1 = FindVariableFeatures(a375_R1_g1, selection.method = "vst", nfeatures = 2000)
# scale data 
a375_R1_g1 = ScaleData(a375_R1_g1, vars.to.regress = c("nCount_RNA", "percent.mt"))
# PCA
a375_R1_g1 = RunPCA(a375_R1_g1, features = VariableFeatures(object = a375_R1_g1))
DimPlot(a375_R1_g1, reduction = "pca")
ElbowPlot(a375_R1_g1)

# Clustering
a375_R1_g1 = FindNeighbors(a375_R1_g1, dims = 1:10)
a375_R1_g1 = FindClusters(a375_R1_g1, resolution = 0.4)
a375_R1_g1 = RunUMAP(a375_R1_g1, dims = 1:10, reduction = "pca")
DimPlot(a375_R1_g1, reduction = "umap")
DimPlot(a375_R1_g1, reduction = "umap", group.by = "orig.ident")

# Find markers
a375_R1_g1_markers_new = FindAllMarkers(a375_R1_g1, test.use = "MAST",
                                   only.pos = T, 
                                   min.pct = 0.2, min.cells.feature = 10)
# filter pvalue
a375_R1_g1_markers_new = a375_R1_g1_markers[a375_R1_g1_markers$p_val_adj < 0.01, ]
FeaturePlot(a375_R1_g1, features = c("ATF4"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

cluster_markers = list(c0_markers, c1_markers, c2_markers,
                       c3_markers, c4_markers, c5_markers,
                       c6_markers, c7_markers)
names(cluster_markers) <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7")
save(cluster_markers, a375_R1_g1_markers, file = "cluster_markers.rda", version = 2)

save(a375_R1, a375_R1_g1, a375_R1_g1_markers, a375_R1_g1_markers_new, file = "a375_R1_seurat.rda")
write.table(a375_R1_g1_markers_new, file = "a375_markers.txt", col.names = T, row.names = T,
            sep = "\t", quote = F)

load("a375_R1_seurat.rda")

# save table for GEO
geo_a375 = as.matrix(a375_R1@assays$RNA@counts)
write.table(geo_a375, file = "a375_all_C_to_T4.csv", sep = ",", col.names = T, row.names = T,
            quote = F)

# Save metadata for GEO
write.table(a375_R1_g1@meta.data, file = "a375_metadata_g1.tsv", sep = "\t", col.names = T,
            row.names = T, quote = F)


######################
# Plot Umap Sample   #
######################

a375_g1_umap_sample = DimPlot(a375_R1_g1, reduction = "umap", group.by = "orig.ident")

colors = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
           "#D0B541", "#E67F33", "#CE2220")

a375_umap_sample_plot = ggplot(a375_g1_umap_sample$data, aes(y = UMAP_2, x = UMAP_1, color = orig.ident)) +
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
a375_R1_g1 <- RenameIdents(object = a375_R1_g1,
                          "0" = "6",
                          "1" = "2",
                          "2" = "3",
                          "3" = "7",
                          "4" = "1",
                          "5" = "0",
                          "6" = "4",
                          "7" = "5")
# assign new order
a375_R1_g1@active.ident <- factor(x = a375_R1_g1@active.ident,
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7"))

a375_R1_umap_cluster = DimPlot(a375_R1_g1, reduction = "umap")

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd", "#df9c08",
                   "#946316", "#dd0177")


a375_umap_cluster_plot = ggplot(a375_R1_umap_cluster$data, aes(y = UMAP_2, x = UMAP_1, color = ident)) +
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

new_meta = a375_R1_g1@meta.data
new_meta$new_clusters <- a375_R1_umap_cluster$data$ident

########################################################
# Replot barplot with frequencies and renamed clusters #
########################################################

library(dplyr)
library(ggplot2)

# The order really matters! first clusters then condition
freq_cluster = new_meta %>% group_by(orig.ident, new_clusters) %>%
  dplyr::summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_cluster) <- c("condition", "cluster", "n", "freq")

freq_cluster$cluster <- factor(freq_cluster$cluster,
                               levels = c("0", "1", "2", "3", "4", "5", "6","7"))

freq_cluster$condition <- factor(freq_cluster$condition,
                               levels = c("C","T003","T006","T0125","T0250",
                                          "T05","T1","T2","T4"))

colors = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
           "#D0B541", "#E67F33", "#CE2220")

a375_freq_cluster_plot = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  scale_fill_manual(values = cluster_colors,
                    limits = c("0","1", "2", "3", "4", "5", "6", "7")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a375_umap_sample_plot + a375_umap_cluster_plot + a375_freq_cluster_plot

dev.off()


####################################################
# Plot barplot of frequencies of cells in clusters #
####################################################

library(dplyr)
library(ggplot2)

# The order really matters! first clusters then condition
meta_data_a375 = a375_R1_g1@meta.data
freq_clusters = meta_data_a375 %>% group_by(orig.ident, seurat_clusters) %>%
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


#pdf("state_frequency_plot_kura.pdf")
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
  #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  #scale_fill_manual(values = state_colors,
  #                  labels = c("C", "S1", "S2", "S3", "S4", "S5"),
  #                  name = "State") +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8)) +
  theme(strip.background =element_rect(fill="white"))
dev.off()


##################################################################
# Assigning the new identity to meta data to do diff. expression #
##################################################################

a375_R1_g1@meta.data$new_id <- a375_R1_g1@active.ident
Idents(a375_R1_g1) <- a375_R1_g1@meta.data$'new_id'

# Find markers
a375_R1_g1_markers_new = FindAllMarkers(a375_R1_g1, test.use = "MAST",
                                        only.pos = T, group.by = "new_id",
                                        min.pct = 0.2, min.cells.feature = 10)
# filter pvalue
a375_R1_g1_markers_new = a375_R1_g1_markers_new[a375_R1_g1_markers_new$p_val_adj < 0.01, ]

FeaturePlot(a375_R1_g1, 
            features = c("ROBO1", "BMP5", "SYT4", "NEFL"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


##########################
# Dotplot of the markers #
##########################

genes_to_show = c("SOX10", "CD74", "HLA-DRA", "HLA-DRB1", "FOS", "JUN",
                  "GAS7", "DCT","NGFR", "MITF", "RAP1B", "TMSB4X",
                  "ROBO1", "BMP5", "SYT4",
                  "TRPC4", "FOXP1", "EXT1", "SOX9", "FN1", "HIF1A", "PFKP",
                  "GPX4", "FTH1", "FTL", "IMPDH2",
                  "ATF4", "PHGDH", "SHMT2", "PSAT1", "GSTP1", 
                  "GSTK1", "MGST1", "MGST3", "SLC3A2",
                  "TKT", "TALDO1", "AKR1B1", "AKR1C3", "TXNRD1", "CYP1B1")


avg_expression = DotPlot(a375_R1_g1, features = genes_to_show)
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))
# rename clusters, already done with the active ident.
#avg_expression$id <- stringr::str_replace_all(avg_expression$id, 
#                                                c("1" = "A", "3" = "B", "0" = "C",
#                                                  "4" = "D", "2" = "E"))
#avg_expression$id <- stringr::str_replace_all(avg_expression$id, 
#                                                c("A" = "0", "B" = "1", "C" = "2",
#                                                  "D" = "3", "E" = "4"))

avg_expression$id <- factor(avg_expression$id, levels = c(7,6,5,4,3,2,1,0))

ggplot(avg_expression,
       aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  #geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_point(shape = 21) +
  #scale_size_continuous(range = c(1, 4)) +
  #scale_fill_gradientn(colours = pals::coolwarm(100)) +
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

setwd("/path/a375_resistance/")

load("a375_R1_seurat.rda")
load("a375_R1_g1_individual.rda")
source("plotting.R")
source("seurat_wrapper.R")

# split into original samples

a375_R1_g1_split = SplitObject(a375_R1_g1, split.by = "orig.ident")

# now cluster them individually and plot!

seurat_individual_analysis <- function(seurat_obj_ind, prefix = "a375_R1_g1") {
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
a375_R1_g1_individual = lapply(a375_R1_g1_split, seurat_individual_analysis)
save(a375_R1_g1_individual, a375_R1_g1_split, file = "a375_R1_g1_individual.rda")

# plot and collect individual umap plots: Colors list in colors.R already has the 
# cluster colors in order

plots_umap_list <- list()
plot_titles = names(a375_R1_g1_individual)
for(i in 1:length(a375_R1_g1_individual)) {
  umap_plot = plot_umap(a375_R1_g1_individual[[i]][[1]],
                        cluster_colors_list[[i]],
                        title = plot_titles[i])
  plots_umap_list[[i]] <- umap_plot
}

# saved as 6 x 8 in pdf in my local computer, file: a375_R1_g1_umaps_cluster.pdf
# location: ../projects/new_lines_evolution/A375_dabrafenib/A375_R1
do.call("grid.arrange", c(plots_umap_list, ncol = 3))

#plot_umap(a375_R1_g1_individual$T4$a375_R1_g1_T4,
#          colors = cluster_colors_new[30:33])

FeaturePlot(a375_R1_g1_individual$T4$a375_R1_g1_T4, 
            features = c("SOX10", "SOX9", "MITF", "AXL", "SOX5", "DDIT4", "DDIT3"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)

##########################################
# Bulk analysis                          #
##########################################

library(ComplexHeatmap)
library(RColorBrewer)

# get top pc1 and pc2 genes for each individual seurat object (sample)
top_pc1 = sapply(sapply(a375_R1_g1_individual,"[[",1), get_top_pc_genes, pc = 1, n = 40)
top_pc2 = sapply(sapply(a375_R1_g1_individual,"[[",1), get_top_pc_genes, pc = 2, n = 40)

top_pc_genes = unique(c(as.vector(top_pc1), as.vector(top_pc2)))
var_genes = a375_R1_g1@assays$RNA@var.features

# do on all of them together
a375_g1_avg = AverageExpression(a375_R1_g1, group.by = "orig.ident")
a375_g1_avg = a375_g1_avg$RNA

# pairwise expression correlation on top PC genes
cormM = cor(a375_g1_avg[var_genes, ], method="spearman")
cormM = cor(a375_g1_avg, method="spearman")
cormM = cor(a375_g1_avg[, -c(2)], method="spearman")
cormM = cor(a375_g1_avg[top_pc_genes, -c(2)], method="spearman")

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
             heatmap_width = unit(7, "cm"), heatmap_height = unit(6, "cm"),
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
sample_id = names(a375_R1_g1_individual)
for(i in 1:length(a375_R1_g1_individual)) {
  a375_avg_data_ind = AverageExpression(a375_R1_g1_individual[[i]][[1]],
                                        group.by = "seurat_clusters")
  a375_avg_data_ind = a375_avg_data_ind$RNA
  # add sample names to the cluster columns
  colnames(a375_avg_data_ind) <- paste(sample_id[i], colnames(a375_avg_data_ind), sep = "_")
  # add in a list
  avg_data_list[[i]] <- a375_avg_data_ind
  # merge the dataframes 
  a375_g1_avg_cluster = merge_exp(avg_data_list)
}

# try my own method for averages of clusters to see if there is a difference with Seurat's ?

# Get top diff. expressed genes for each condition/cluster
# get top pc1 and pc2 genes for each individual seurat object (sample)
top_diff = sapply(sapply(a375_R1_g1_individual, "[", 2), get_top_diff_genes, n = 100)
top_diff = unique(unlist(top_diff))

# select and correlate clusters based on diff. expressed genes
# pairwise expression correlation on top PC genes
cormM = cor(a375_g1_avg_cluster[top_diff, ], method="spearman")
cormM = cor(a375_g1_avg_cluster[top_pc_genes, ], method="spearman")

cormM = cor(t(scale(t(a375_g1_avg_cluster[var_genes, ]))), method="spearman")
cormM = cor(a375_g1_avg_cluster[var_genes, ], method="spearman")


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
             heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "Corr",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))

draw(h2)
dev.off()

FeaturePlot(a375_R1_g1_individual$T0250$a375_R1_g1_T006, 
            features = c("ROBO1", "CD74", "SOX10"), ncol = 3)


plots_markers_list <- list()
plot_titles = names(a375_R1_g1_individual)
for(i in 1:length(a375_R1_g1_individual)) {
  umap_plot_markers = FeaturePlot(a375_R1_g1_individual[[i]][[1]],
                          features = c("NRG1"), ncol = 3)
  plots_markers_list[[i]] <- umap_plot_markers
}

VlnPlot(a375_R1_g1_individual$T4$a375_R1_g1_T4, features = c("TENM2"), ncol = 1)

do.call("grid.arrange", c(plots_markers_list, ncol = 3))


