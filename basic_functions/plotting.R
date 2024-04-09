
# To plot individual umaps for each sample
plot_umap <- function(seurat_obj, colors, title) {
  # Plot umap figure for individual samples
  # First use the DimPlot from seurat to get the basis plot
  # then modify it
  plot_umap = DimPlot(seurat_obj, reduction = "umap")
  plot_umap = plot_umap + scale_color_manual(values = colors, name = "") +
    ggtitle(title) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(text = element_text(size = 8),
          axis.line=element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
  return(plot_umap)  
} 
