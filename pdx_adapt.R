library(Seurat)
library(ggplot2)
library(dplyr)

setwd("/path/pdx_talazo_analysis")

load("pdx_adapted_raw_g1.rda")

#################################
# Seurat object for Replicate 1 #
#################################

# REGRESS RIBO GENES TOO

raw_data_list = list(white1_g1_raw, brown_g1_raw, res1_g1_raw, res2_g1_raw,
                     o174v_g1_raw, o174_t2_g1_raw, o174_t1_g1_raw)

names(raw_data_list) <- c("v1", "v2", "res1_1", "res1_2", "res3_v",
                          "res3_1", "res3_2")

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

### Fix Names field to 1, and try other parameters for integration such as dimensions, var genes, etc.
seurat_pipeline <- function(x) {
    # select ribo genes that are found in the matrix
    ribosomal_genes = rownames(x[rownames(x) %in% ribo_genes, ])
    x = CreateSeuratObject(x, names.delim = "_", names.field = 1)
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]] <- PercentageFeatureSet(x, features = ribosomal_genes)
    x = NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
    x = ScaleData(x, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))
    x = RunPCA(x, features = VariableFeatures(object = x))
    x = FindNeighbors(x, dims = 1:10)
    x = FindClusters(x, resolution = 0.3)
    x = RunUMAP(x, dims = 1:10)
    return(x)
}

seurat_diff_expression <- function(x) {
  x = FindAllMarkers(x, test.use = "MAST", only.pos = T, 
                 min.pct = 0.2, min.cells.feature = 10)
  x = x[x$p_val_adj < 0.01, ]
  return(x)
}

seurat_objects = lapply(raw_data_list, seurat_pipeline)
# set new ids
seurat_objects$v1$new_id <- "V1"
seurat_objects$v2$new_id <- "V2"
seurat_objects$res1_1$new_id <- "R1_1"
seurat_objects$res1_2$new_id <- "R1_2"
seurat_objects$res3_v$new_id <- "R3_V"
seurat_objects$res3_1$new_id <- "R3_1"
seurat_objects$res3_2$new_id <- "R3_2"


seurat_markers = lapply(seurat_objects, seurat_diff_expression)

save(seurat_objects, file = "pdx_adapt_seurat_list_g1.rda")


FeaturePlot(seurat_objects$res3_1, features = c("SOX17", "CD44", "VIM", "GAPDH", "MYC",
                                                "DDIT4", "PHLDA1", "KRT8", "KRT18", "WT1", "PAX8"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


################
# All together #
################

# merging multiple data
pdx_all_adapted_combined_no_batch <- merge(x = seurat_objects$v1, 
                                           y = c(seurat_objects$v2, seurat_objects$res1_1,
                                                 seurat_objects$res1_2, seurat_objects$res3_v,
                                                 seurat_objects$res3_1, seurat_objects$res3_2))

pdx_all_adapted_combined_no_batch = NormalizeData(pdx_all_adapted_combined_no_batch, 
                                                  normalization.method = "LogNormalize", scale.factor = 10000)
pdx_all_adapted_combined_no_batch = FindVariableFeatures(pdx_all_adapted_combined_no_batch, 
                                                         selection.method = "vst", nfeatures = 1000)
pdx_all_adapted_combined_no_batch = ScaleData(pdx_all_adapted_combined_no_batch, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))

# remove mito and ribo genes from var genes
mito_genes =  grep(pattern = "^MT-", x = rownames(pdx_all_adapted_combined_no_batch@assays$RNA@data), value = TRUE)
genes_to_remove = c(ribo_genes, mito_genes)

var_genes = VariableFeatures(pdx_all_adapted_combined_no_batch)
# get top 1500 var genes excluding ribo and mt
var_genes = var_genes[!(var_genes %in% genes_to_remove)][1:1500]


pdx_all_adapted_combined_no_batch = RunPCA(pdx_all_adapted_combined_no_batch, features = var_genes)
ElbowPlot(pdx_all_adapted_combined_no_batch)
pdx_all_adapted_combined_no_batch = FindNeighbors(pdx_all_adapted_combined_no_batch, dims = 1:10)
pdx_all_adapted_combined_no_batch = FindClusters(pdx_all_adapted_combined_no_batch, resolution = 0.2)
pdx_all_adapted_combined_no_batch = RunUMAP(pdx_all_adapted_combined_no_batch, dims = 1:10)

DimPlot(pdx_all_adapted_combined_no_batch, group.by = c("new_id", "RNA_snn_res.0.2"), ncol = 2)

pdx_all_markers_no_batch = FindAllMarkers(pdx_all_adapted_combined_no_batch, test.use = "MAST",
                                 only.pos = T, #group.by = "RNA_snn_res.0.2",
                                 min.pct = 0.2, min.cells.feature = 10)

pdx_all_markers_no_batch = pdx_all_markers_no_batch[pdx_all_markers_no_batch$p_val_adj < 0.01, ]


save(pdx_all_adapted_combined_no_batch, pdx_all_markers_no_batch, new_meta_pdx_resistant,
     table_pdx_adapted_markers,
     file = "pdx_all_adapted_seurat_g1.rda")


load("pdx_all_adapted_seurat_g1.rda")

####################################################
# Plot barplot of frequencies of cells in clusters #
####################################################

pdx_all_adapted_combined_meta = pdx_all_adapted_combined@meta.data
pdx_all_adapted_combined_meta$new_id <- stringr::str_replace_all(pdx_all_adapted_combined_meta$orig.ident, 
                                            c("white1" = "V1", "brown" = "V2",
                                              "res1" = "R1_1", "res2" = "R1_2",
                                              "o174v" = "R3_V", "o174t1" = "R3_2",
                                              "o174t2" = "R3_1"))

pdx_all_adapted_combined_meta$new_id <- factor(pdx_all_adapted_combined_meta$new_id,
                          levels = c("V1", "V2", "R1_1", "R1_2", "R3_V", "R3_1", "R3_2"))
pdx_all_adapted_combined_meta$all <- c("All")

pdx_all_adapted_combined_meta_no_batch = pdx_all_adapted_combined_no_batch@meta.data
pdx_all_adapted_combined_meta_no_batch$new_id <- stringr::str_replace_all(pdx_all_adapted_combined_meta_no_batch$orig.ident, 
                                                                 c("white1" = "V1", "brown" = "V2",
                                                                   "res1" = "R1_1", "res2" = "R1_2",
                                                                   "o174v" = "R3_V", "o174t1" = "R3_2",
                                                                   "o174t2" = "R3_1"))

pdx_all_adapted_combined_meta_no_batch$new_id <- factor(pdx_all_adapted_combined_meta_no_batch$new_id,
                                               levels = c("V1", "V2", "R1_1", "R1_2", "R3_V", "R3_1", "R3_2"))

# The order really matters! first clusters then condition
freq_clusters = pdx_all_adapted_combined_meta %>% group_by(seurat_clusters, new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_clusters) <- c("cluster", "condition", "n", "freq")

freq_clusters = pdx_all_adapted_combined_meta_no_batch %>% group_by(seurat_clusters, new_id) %>%
  summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_clusters) <- c("cluster", "condition", "n", "freq")

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


freq_cluster_all = rbind(as.data.frame(freq_all), as.data.frame(freq_clusters))

freq_cluster_all$condition <- factor(freq_cluster_all$condition)
freq_cluster_all$cluster <- factor(freq_cluster_all$cluster)
freq_cluster_all$cluster <- factor(freq_cluster_all$cluster,
                                   levels = c("All", "0", "1", "2", "3", "4","5", "6", "7", "8"))
freq_cluster_all$condition <- factor(freq_cluster_all$condition,
                                     levels = c("R3_2", "R3_1", "R3_V", "R1_2", "R1_1", "V1", "V2"))


colors = c("#8DA3A6","#6b7d7f", "#21773c", "#1b4414", 
           "#36a6bb", "#2674b0", "#2f327d")
freq_cluster_all$cluster <- factor(freq_cluster_all$cluster,
                                   levels = c("8", "7", "6", "5", "4", "3", "2", "1", "0", "All"))
ggplot(freq_cluster_all, aes(y = freq, x = cluster, fill = condition)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = colors, limits = c("V1", "V2", "R1_1", "R1_2", 
                               "R3_V", "R3_1", "R3_2")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=8)) +
  theme(strip.background =element_rect(fill="white"))
dev.off()


DimPlot(pdx_all_adapted_combined_no_batch, reduction = "umap", group.by = "new_id")

##########
# New plot frequency
#########

# renaming clusters 
new_meta = pdx_all_adapted_combined_meta_no_batch

new_meta$new_id <- factor(new_meta$new_id,
                          levels = c("V1", "V2", "R1_1", "R1_2", "R3_V", "R3_1", "R3_2"))
new_meta$all <- c("All")

# changing new ids to replicates
### rename clusters to make them ordered by vehicle frequency
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$seurat_clusters, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$new_cluster_id, 
                                                    c("A" = "6", "B" = "8", "C" = "7",
                                                      "D" = "9", "E" = "5"))

new_meta$new_cell <- paste(rownames(new_meta), new_meta$new_id, sep = "_")
new_meta_pdx_resistant = new_meta

# Save metadata for GEO submission
write.table(new_meta, file = "pdx_resistance_metadata_g1.tsv", row.names = T, col.names = T,
            quote = F, sep = "\t")

###
# Save markers table for submission
###

table_pdx_adapted_markers = pdx_all_markers_no_batch
table_pdx_adapted_markers$cluster <- stringr::str_replace_all(table_pdx_adapted_markers$cluster, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
table_pdx_adapted_markers$cluster <- stringr::str_replace_all(table_pdx_adapted_markers$cluster, 
                                                    c("A" = "6", "B" = "8", "C" = "7",
                                                      "D" = "9", "E" = "5"))

write.table(table_pdx_adapted_markers, file = "pdx_adapted_markers.tsv", row.names = F,
            col.names = T, quote = F, sep = "\t")

FeaturePlot(pdx_all_adapted_combined_no_batch, features = c("LCN2", "GSPT1", "ATF4", "SHMT2"),
            reduction = "umap", cols = c("lightgrey", "blue"), 
            pt.size = 0.75)


##########################
# Plot Cluster frequency #
##########################

new_meta = pdx_all_adapted_combined_meta_no_batch

new_meta$new_id <- factor(new_meta$new_id,
                          levels = c("V1", "V2", "R1_1", "R1_2", "R3_V", "R3_1", "R3_2"))

# changing new ids to replicates
### rename clusters to make them ordered by vehicle frequency
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$seurat_clusters, 
                                                    c("0" = "A", "1" = "B", "2" = "C",
                                                      "3" = "D", "4" = "E"))
new_meta$new_cluster_id <- stringr::str_replace_all(new_meta$new_cluster_id, 
                                                    c("A" = "6", "B" = "8", "C" = "7",
                                                      "D" = "9", "E" = "5"))
new_meta$new_id_norep <- stringr::str_replace_all(new_meta$new_id, 
                                                  c("V1" = "V", "V2" = "V", "R1_1" = "R1",
                                                    "R1_2" = "R1", "R3_V" = "R3", 
                                                    "R3_1" = "R3", "R3_2" = "R3"))

freq_cluster = new_meta %>% group_by(new_id, new_cluster_id) %>%
  dplyr::summarise(n = n()) %>% dplyr::mutate(freq = (n / sum(n)) * 100)
colnames(freq_cluster) <- c("condition", "cluster", "n", "freq")

freq_cluster$cluster <- factor(freq_cluster$cluster,
                               levels = c("5", "6", "7","8","9"))

freq_cluster$condition <- factor(freq_cluster$condition,
                                 levels = c("V1","V2","R1_1","R1_2","R3_V",
                                            "R3_1", "R3_2"))

cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")

pdx_adapt_frequency_cluster = ggplot(freq_cluster, aes(y = freq, x = condition, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("") +
  ylab("% of cells") +
  labs(fill = "") +
  theme_classic() +
  scale_fill_manual(values = cluster_colors,
                    limits = c("5","6", "7", "8", "9")) +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white"))

############### freq plot ends here ############

#################################
# Plot Umap Sample and clusters #
#################################

pdx_all_adapted_combined_umap_no_batch = DimPlot(pdx_all_adapted_combined_no_batch, reduction = "umap")
# fix cluster ids as frequency plot
pdx_all_adapted_combined_umap_no_batch$data$new_cluster_id <- new_meta$new_cluster_id
pdx_all_adapted_combined_umap_no_batch$data$new_id <- new_meta$new_id


sample_colors = c("#8DA3A6","#6b7d7f", "#c6497d", "#8f3152", 
                  "#36a6bb", "#2674b0", "#2f327d")

ggplot(pdx_all_adapted_combined_umap_no_batch$data, aes(y = UMAP_2, x = UMAP_1, 
                                                        color = new_id)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = sample_colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Sample", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))

# cluster umap
cluster_colors = c("#2777b4", "#f67f11", "#2ea02d", "#d62727", "#9467bd")
umap_cluster_plot = ggplot(pdx_all_adapted_combined_umap_no_batch$data, aes(y = UMAP_2, x = UMAP_1, 
                                                                            color = new_cluster_id)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = cluster_colors) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() +
  guides(color=guide_legend("Cluster", override.aes = list(size = 2))) +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))

ggsave("umap_plot_adapt.pdf", plot = umap_cluster_plot, 
       width = 3.5, height = 3, units = c("cm"),
       dpi = 300, limitsize = T, scale = 3)

# all together
(umap_sample_plot | umap_cluster_plot)

