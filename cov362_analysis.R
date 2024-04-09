library(Seurat)

setwd("/path/cov362_10X")
load("cov362_g1_raw.rda")


#################################
# Seurat pipeline               #
#################################

ribo_genes = c("RPL3","RPL4","RPL5","RPL6","RPL7","RPL7A","RPL8","RPL9","RPL10",
               "RPL10A","RPL11","RPL12","RPL13","RPL13A","RPL14","RPL15","RPL17",
               "RPL18","RPL18A","RPL19","RPL21","RPL22","RPL23","RPL23A","RPL24",
               "RPL26","RPL27","RPL27A","RPL28","RPL29","RPL30","RPL31","RPL32",
               "RPL34","RPL35","RPL35A","RPL36","RPL36A","RPL37","RPL37A","RPL38",
               "RPL39","RPL40","RPL41","RPLP0","RPLP1","RPLP2","RPSA","RPS2","RPS3",
               "RPS3A","RPS4X","RPS4Y","RPS5","RPS6","RPS7","RPS8","RPS9","RPS10","RPS11",
               "RPS12","RPS13","RPS14","RPS15","RPS15A","RPS16","RPS17","RPS18","RPS19",
               "RPS20","RPS21","RPS23","RPS24","RPS25","RPS26","RPS27","RPS27A","RPS28",
               "RPS29","RPS30")

seurat_pipeline <- function(x, dimensions) {
  # select ribo genes that are found in the matrix
  ribosomal_genes = rownames(x[rownames(x) %in% ribo_genes, ])
  x = CreateSeuratObject(x, names.delim = "_", names.field = 1)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x[["percent.ribo"]] <- PercentageFeatureSet(x, features = ribosomal_genes)
  x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  x = ScaleData(x, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))
  x = RunPCA(x, features = VariableFeatures(object = x))
  x = FindNeighbors(x, dims = dimensions)
  x = FindClusters(x, resolution = 0.3)
  x = RunUMAP(x, dims = dimensions)
  return(x)
}

seurat_diff_expression <- function(x) {
  x = FindAllMarkers(x, test.use = "MAST", only.pos = T, 
                     min.pct = 0.2, min.cells.feature = 10)
  x = x[x$p_val_adj < 0.01, ]
  return(x)
}


cov362_g1 = seurat_pipeline(cov362_g1_raw_data, dimensions = 6)

ribosomal_genes = rownames(cov362_g1_raw_data[rownames(cov362_g1_raw_data) %in% ribo_genes, ])

cov362_g1 = CreateSeuratObject(cov362_g1_raw_data, names.delim = "_", names.field = 2)
cov362_g1[["percent.mt"]] <- PercentageFeatureSet(cov362_g1, pattern = "^MT-")
cov362_g1[["percent.ribo"]] <- PercentageFeatureSet(cov362_g1, features = ribosomal_genes)
cov362_g1 = NormalizeData(cov362_g1, normalization.method = "LogNormalize", scale.factor = 10000)
cov362_g1 = FindVariableFeatures(cov362_g1, selection.method = "vst", nfeatures = 1000)
cov362_g1 = ScaleData(cov362_g1, vars.to.regress = c("nCount_RNA", "percent.mt"))

cov362_g1 <- RunPCA(cov362_g1, features = VariableFeatures(object = cov362_g1))
DimPlot(cov362_g1, reduction = "pca")
ElbowPlot(cov362_g1)
cov362_g1 <- FindNeighbors(cov362_g1, dims = 1:6)
cov362_g1 <- FindClusters(cov362_g1, resolution = 0.3)
cov362_g1 <- RunUMAP(cov362_g1, dims = 1:6)
DimPlot(cov362_g1, reduction = "umap")
DimPlot(cov362_g1, reduction = "umap", group.by = "orig.ident")
cov362_g1_markers = FindAllMarkers(cov362_g1, test.use = "MAST",
                                    only.pos = T, 
                                    min.pct = 0.2, min.cells.feature = 10)
cov362_g1_markers = cov362_g1_markers[cov362_g1_markers$p_val_adj < 0.01, ]
FeaturePlot(cov362_g1, features = c("SOX17", "GLS"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


######################
# Plot umap clusters #
######################

# rename clusters to put in order of more resistant
cov362_g1 <- RenameIdents(object = cov362_g1,
                          "0" = "2",
                          "1" = "1",
                          "2" = "0")
# assign new order
cov362_g1@active.ident <- factor(x = cov362_g1@active.ident,
                                 levels = c("0", "1", "2"))

cov362_g1@meta.data$new_id <- cov362_g1@active.ident
Idents(cov362_g1) <- cov362_g1@meta.data$'new_id'

cov362_umap_cluster = DimPlot(cov362_g1, reduction = "umap")

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d")


cov362_g1_umap_cluster_plot = ggplot(cov362_umap_cluster$data, aes(y = UMAP_2, x = UMAP_1, color = ident)) +
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

new_meta = cov362_g1@meta.data
new_meta$orig.ident <- stringr::str_replace_all(new_meta$orig.ident, c("T0" = "C"))

# adding the fixed cluster names from previous plot
new_meta$new_clusters <- cov362_g1_umap_cluster_plot$data$ident

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
                               levels = c("0", "1", "2"))

freq_cluster$condition <- factor(freq_cluster$condition,
                                 levels = c("C","T5","T10","T20","T40"))

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d")

cov362_freq_cluster_plot = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  scale_fill_manual(values = cluster_colors,
                    limits = c("0","1", "2")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


######################
# Plot Umap Sample   #
######################

cov362_g1_umap_sample = DimPlot(cov362_g1, reduction = "umap", group.by = "orig.ident")

# rename T0 to C
cov362_g1_umap_sample$data$orig.ident <- stringr::str_replace_all(cov362_g1_umap_sample$data$orig.ident, c("T0" = "C"))
# reorder factors
cov362_g1_umap_sample$data$orig.ident <- factor(cov362_g1_umap_sample$data$orig.ident, 
                                                levels = c("C", "T5", "T10", "T20", "T40"))

colors = c("#8DA3A6","#4E78C4","#57A2AC","#7EB875","#D0B541")

cov362_g1_umap_sample_plot = ggplot(cov362_g1_umap_sample$data, aes(y = UMAP_2, x = UMAP_1, color = orig.ident)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Treatment", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))

cov362_g1_umap_sample_plot + cov362_g1_umap_cluster_plot + cov362_freq_cluster_plot

cov362_g1_new_metadata = new_meta
save(cov362_g1_markers, cov362_g1, cov362_frequency_plot, cov362_g1_umap_cluster_plot,
     cov362_g1_umap_sample_plot, cov362_g1_umap, cov362_g1_new_metadata, cov362_g1_markers_new_id,
     file = "cov362_g1_analysis.rda")

load("cov362_g1_analysis.rda")

####
# Plotting gene markers
###

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

# calculate average expression by State
new_metadata = cov362_g1@meta.data

new_metadata$new_id <- stringr::str_replace_all(new_metadata$seurat_clusters, 
                                                c("0" = "C", "1" = "B", "2" = "A"))
new_metadata$new_id <- stringr::str_replace_all(new_metadata$new_id, 
                                                c("A" = "0", "B" = "1", "C" = "2"))
## Saving diff. expressed genes in clusters
cov362_g1_markers_new_id = cov362_g1_markers
cov362_g1_markers_new_id$cluster = stringr::str_replace_all(cov362_g1_markers_new_id$cluster, 
                                                    c("0" = "C", "1" = "B", "2" = "A"))
cov362_g1_markers_new_id$cluster <- stringr::str_replace_all(cov362_g1_markers_new_id$cluster, 
                                                c("A" = "0", "B" = "1", "C" = "2"))

# Write table for paper suppl
write.table(cov362_g1_markers_new_id, file = "cov362_markers.tsv", row.names = F,
            col.names = T, quote = F, sep = "\t")

# Write metadata for GEO 
write.table(cov362_g1@meta.data, file = "cov362_metadata_g1.tsv", row.names = T,
            col.names = T, quote = F, sep = "\t")



labels = as.character(unique(new_metadata$new_id))
list_means = list()
for(label in labels) {
  cells = rownames(new_metadata[new_metadata$new_id == label, ])
  avg_cells = rowMeans(as.matrix(cov362_g1@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}
cov362_g1_avg_data = do.call(cbind, list_means)
cov362_g1_avg_data = cov362_g1_avg_data[, c(2,3,1)]

genes_to_show = c("SOX17", "PAX8", "KRT8", "IFI27", "IFI35", "FN1", "VIM", "CD44",
                  "SMAD3", "HIF1A", "CITED2", "SGK1", "PGK1", "PFKP", "ACLY", "MBOAT7", "CTPS1", "TYMS",
                  "GPX4", "FTH1", "FTL",
                  "PHGDH", "SHMT2", "ASNS", "GSTP1", "GDF15",  
                  "TKT", "SOD2", "AKR1B1", "GGH", "GSTM3", "PRDX1", "CYP1B1", "TXNRD1", "SQSTM1", "NFE2L2")


avg_expression = DotPlot(cov362_g1, features = genes_to_show)
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))

avg_expression$id <- factor(avg_expression$id, levels = c(2,1,0))

ggplot(avg_expression,
       aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  geom_point(shape = 21) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F",
                       breaks = c(-1,0,1)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")

save_plot("cov362_markers_bubble.pdf", cov362_markers_bubble, base_height = 3, base_width = 9.5)



####################
# Bulk correlation #
####################

new_metadata = cov362_g1@meta.data
labels = as.character(unique(new_metadata$orig.ident))
list_means = list()
for(label in labels) {
  cells = rownames(new_metadata[new_metadata$orig.ident == label, ])
  avg_cells = rowMeans(as.matrix(cov362_g1@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}
cov362_g1_avg_bulk = do.call(cbind, list_means)
cov362_g1_avg_bulk = cov362_g1_avg_bulk[, c(2,3,4,1,5)]

# try with variable genes
var_genes <- VariableFeatures(cov362_g1)
seurat_df <- GetAssayData(cluster3.seurat.obj)[var_genes,]

cormM = cor((cov362_g1_avg_bulk[var_genes, ]), method="spearman")


library(dendextend)
dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T
dend <- click_rotate(dend, continue = TRUE)

h2 = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
             show_column_dend = F, show_row_names = T, 
             name = "Corr", row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 7),
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
             cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
             column_dend_reorder = F,
             heatmap_width = unit(7, "cm"), heatmap_height = unit(6, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                         title = "Correlation",
                                         labels_gp = gpar(fontsize = 7),
                                         legend_height = unit(2, "cm")))
draw(h2)


