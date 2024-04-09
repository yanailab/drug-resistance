library(Seurat)
library(ggplot2)

setwd("/path/ovsaho_10X")
load("ovsaho_g1_raw.rda")


ovsaho_g1 = CreateSeuratObject(ovsaho_g1_raw_data, names.delim = "_", names.field = 2)
ovsaho_g1[["percent.mt"]] <- PercentageFeatureSet(ovsaho_g1, pattern = "^MT-")
ovsaho_g1 = NormalizeData(ovsaho_g1, normalization.method = "LogNormalize", scale.factor = 10000)
ovsaho_g1 = FindVariableFeatures(ovsaho_g1, selection.method = "vst", nfeatures = 1000)
ovsaho_g1 = ScaleData(ovsaho_g1, vars.to.regress = c("nCount_RNA", "percent.mt"))

ovsaho_g1 <- RunPCA(ovsaho_g1, features = VariableFeatures(object = ovsaho_g1))
DimPlot(ovsaho_g1, reduction = "pca")
ElbowPlot(ovsaho_g1)
ovsaho_g1 <- FindNeighbors(ovsaho_g1, dims = 1:10)
ovsaho_g1 <- FindClusters(ovsaho_g1, resolution = 0.3)
ovsaho_g1 <- RunUMAP(ovsaho_g1, dims = 1:10)
DimPlot(ovsaho_g1, reduction = "umap")
DimPlot(ovsaho_g1, reduction = "umap", group.by = "orig.ident")
ovsaho_g1_markers = FindAllMarkers(ovsaho_g1, test.use = "MAST",
                                    only.pos = T, 
                                    min.pct = 0.2, min.cells.feature = 10)
ovsaho_g1_markers = ovsaho_g1_markers[ovsaho_g1_markers$p_val_adj < 0.01, ]
FeaturePlot(ovsaho_g1, features = c("PGK1", "LDHA", "PHLDA1", "RPL30", "RPS25", "RACK1"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


######################
# Plot umap clusters #
######################

# rename clusters to put in order of more resistant
ovsaho_g1 <- RenameIdents(object = ovsaho_g1,
                          "0" = "2",
                          "1" = "0",
                          "2" = "4",
                          "3" = "1",
                          "4" = "3")
# assign new order
ovsaho_g1@active.ident <- factor(x = ovsaho_g1@active.ident,
                                  levels = c("0", "1", "2", "3", "4"))

ovsaho_g1@meta.data$new_id <- ovsaho_g1@active.ident
Idents(ovsaho_g1) <- ovsaho_g1@meta.data$'new_id'

ovsaho_umap_cluster = DimPlot(ovsaho_g1, reduction = "umap")

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")


ovsaho_g1_umap_cluster_plot = ggplot(ovsaho_umap_cluster$data, aes(y = UMAP_2, x = UMAP_1, color = ident)) +
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

new_meta = ovsaho_g1@meta.data
new_meta$orig.ident <- stringr::str_replace_all(new_meta$orig.ident, c("T0" = "C"))
# adding the fixed cluster names from previous plot
new_meta$new_clusters <- ovsaho_g1_umap_cluster_plot$data$ident

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
                               levels = c("0", "1", "2", "3", "4"))

freq_cluster$condition <- factor(freq_cluster$condition,
                                 levels = c("C","T5","T10","T20","T40"))

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

ovsaho_freq_cluster_plot = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
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


######################
# Plot Umap Sample   #
######################

ovsaho_g1_umap_sample = DimPlot(ovsaho_g1, reduction = "umap", group.by = "orig.ident")

# rename T0 to C
ovsaho_g1_umap_sample$data$orig.ident <- stringr::str_replace_all(ovsaho_g1_umap_sample$data$orig.ident, c("T0" = "C"))
# reorder factors
ovsaho_g1_umap_sample$data$orig.ident <- factor(ovsaho_g1_umap_sample$data$orig.ident, 
                                                levels = c("C", "T5", "T10", "T20", "T40"))

colors = c("#8DA3A6","#4E78C4","#57A2AC","#7EB875","#D0B541")

ovsaho_g1_umap_sample_plot = ggplot(ovsaho_g1_umap_sample$data, aes(y = UMAP_2, x = UMAP_1, color = orig.ident)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Treatment", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))


ovsaho_g1_umap_sample_plot + ovsaho_g1_umap_cluster_plot + ovsaho_freq_cluster_plot

ovsaho_g1_new_metadata = new_meta

save(ovsaho_g1_markers, ovsaho_g1, ovsaho_freq_cluster_plot, ovsaho_g1_umap_cluster_plot,
     ovsaho_g1_umap_sample_plot, ovsaho_g1_umap, ovsaho_g1_new_metadata, ovsaho_g1_markers_new_id,
     file = "ovsaho_g1_analysis.rda")

load("ovsaho_g1_analysis.rda")

library(gridExtra)

ovsaho_g1_umap_all = grid.arrange(ovsaho_g1_umap_sample_plot, ovsaho_g1_umap_cluster_plot, ovsaho_frequency_plot,
             nrow = 1)

save_plot("ovsaho_g1_umap_all.pdf", ovsaho_g1_umap_all, base_height = 3, base_width = 10)



###########################
# Bubble plot for markers #
###########################

genes_to_show = c("SOX17", "PAX8", "WT1", "KRT8", "IFI27", "ISG15", "VIM", "CD24", "CD44",
                  "HIF1A", "DDIT4", "DUSP1", "PGK1", "LDHA", "IMPDH2",
                  "GPX4", "FTH1", "FTL",
                  "ATF4", "PHGDH", "ASNS", "GDF15", "TKT","SOD2", "AKR1B1", "NQO1", "SQSTM1", "NFE2L2")


avg_expression = DotPlot(ovsaho_g1, features = genes_to_show)
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
  labs(size="Percent", fill="Avg. Exp") +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.3, 'cm'), legend.position="bottom")

save_plot("ovsaho_markers_bubble.pdf", ovsaho_markers_bubble, base_height = 3, base_width = 8.5)

