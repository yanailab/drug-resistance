###############################
# Single cell quality control 
# Author: GSF 
# Date: 03-24-17
###############################

library(scater)
library(scran)
library(ggplot2)
library(magrittr)
library(RColorBrewer)
# my set of functions to filter gene expression
source("filter_gene_expression.R")

# things to improve! Condition variables using annotation file
# histogram with density or bars (?), how to visualize for multiple conditions?

# TODO! # FUNCTION FOR SUBSETING GENES BASED ON A LIST OR FILE
# FOR HEATMAP PLOTTING FOR EX

#########################################################################
# My Own functions for doing basically the same things that Scater does #
#########################################################################


mt_genes = c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1",
             "MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6","MT-RNR1",
             "MT-RNR2","MT-TA","MT-TC","MT-TD","MT-TE","MT-TF","MT-TG","MT-TH",
             "MT-TI","MT-TK","MT-TL1","MT-TL2","MT-TM","MT-TN","MT-TP","MT-TQ",
             "MT-TR","MT-TS1","MT-TS2","MT-TT","MT-TV","MT-TW","MT-TY")

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

## still work on this... 
plotQC_scater <- function(annotation_file, exp_matrix) {
  # comment
  # Still have to figure out this Annotated dataframe (deals with multiple
  # samples)
  
  pd <- new("AnnotatedDataFrame", data = annotation_file)
  rownames(pd) <- pd$cells
  
  # gene counts dataframe
  gene_df <- data.frame(Gene = rownames(exp_matrix))
  rownames(gene_df) <- gene_df$Gene
  fd = new("AnnotatedDataFrame", data = gene_df)
  example_sceset = newSCESet(countData = exp_matrix, phenoData = pd,
                             featureData = fd)
  example_sceset <- calculateQCMetrics(example_sceset)
}


###########################################
# 1) Functions for filtering genes (rows) #
###########################################

remove_mt_genes <- function(exp_matrix, mt_genes = mt_genes) {
  # Removes mitochondrial genes from the expression matrix
  
  exp_matrix = exp_matrix[!(rownames(exp_matrix) %in% mt_genes), ]
  return(exp_matrix)
}


remove_ribo_genes <- function(exp_matrix, ribo_genes = ribo_genes) {
  # Removes mitochondrial genes from the expression matrix
  
  exp_matrix = exp_matrix[!(rownames(exp_matrix) %in% ribo_genes), ]
  return(exp_matrix)
}

# Filtering Genes all at once
filter_genes <- function(exp_matrix,
                         min_exp = 1,
                         min_samples = 10,
                         mito_genes = TRUE) {
  # Using pipes from magrittr to filter all genes at once:
  # min_exp: Minimum expression level for genes
  # min_samples: Minimum number of samples (cells) for the gene to be expressed
  # mito_genes: Remove or not mitochondrial genes. Default = True
  exp_matrix = exp_matrix %>% gene_filter(exp_level = min_exp,
                                          n_samples = min_samples) %>%
    remove_ribo_genes(ribo_genes = ribo_genes)
  # If we want to remove mitochodrial genes (Default = True)
  if(mito_genes == TRUE) {
    exp_matrix = exp_matrix %>% remove_mt_genes(mt_genes = mt_genes)
  }
  return(exp_matrix)
}


########################################################
# 2) Functions for filtering by gene expression levels #
########################################################

# call for example tpm() and filtering functions: Ex:
# TPM
# exp_matrix_tpm = tpm(exp_matrix)
# get rid of very lowly expressed genes with default parameters (exp 1, s = 1)
# exp_matrix_tpm = log2(gene_filter(exp_matrix_tpm, exp_level = 1, n_samples = 1) + 1)

# 2) filtering by dynamic genes 
# call for example gene_order_by() and reorder the dataframe based on a statistic. Ex:
# gene_order_by(exp_matrix, FUN). FUN can be fano_factor, cv, std, var, mean, median

# 3) filtering genes by a function. Get those above a threshold based on mean, median, std
# call for example gene_filter_by(). Ex:
# gene_filter_by(exp_matrix, FUN, n=1). FUN usually a function like mean, median; 
#                                       n (multiplying factor, ex: 2 x std) 

# 4) We can also convert back TPM to Raw counts with tpm_to_raw()
# See more from literature

##############################################
# 3) Functions for filtering cells (columns) #
##############################################

# Filter cells with sum of expression < than X. Can be either Raw counts or TPM. Ex:
# sample_filter(exp_matrix, exp_level = n)

# filter by cells with > % of mitocondrial genes
# filter by > number of X transcripts (average from diagnostics ?)
# filter by > number of genes expressed (we don't want very low N of genes)

# See more from literature

mitochondrial_filter <- function(exp_matrix, proportion = 10) {
  # subsets gene expression matrix to keep samples with proportion of expressed
  # mitochondrial genes < proportion
  # Like removing the "bad columns" or bad cells
  # Args:
  #   exp_matrix: the original gene expression matrix read in
  #   proportion: threshold for the percentage of mitochondrial genes accounting
  #              for the total expression. (usually < 10 or 15%)
  
  exp_matrix = exp_matrix[ , (colSums(exp_matrix[rownames(exp_matrix) %in% mt_genes, ]) /
                              colSums(exp_matrix)*100) < proportion, drop=F]
  
  if (ncol(exp_matrix) == 0) {
    print("Empty Columns! No columns matched the criteria.")
  }
  return(exp_matrix)
}

ribosomal_filter <- function(exp_matrix, proportion = 30) {
  # subsets gene expression matrix to keep samples with proportion of expressed
  # ribosomal genes < proportion
  # Like removing the "bad columns" or bad cells
  # Args:
  #   exp_matrix: the original gene expression matrix read in
  #   proportion: threshold for the percentage of mitochondrial genes accounting
  #              for the total expression. (usually < 10 or 15%)
  
  exp_matrix = exp_matrix[ , (colSums(exp_matrix[rownames(exp_matrix) %in% ribo_genes, ]) /
                                colSums(exp_matrix)*100) < proportion, drop=F]
  
  if (ncol(exp_matrix) == 0) {
    print("Empty Columns! No columns matched the criteria.")
  }
  return(exp_matrix)
}


transcript_filter <- function(exp_matrix, min, max = 1000000) {
  # filter out columns (cells) with < transcripts (UMIs) than the threshold
  
  exp_matrix = exp_matrix[ ,colSums(exp_matrix) > min]
  
  if(max != 1000000) {
    exp_matrix = exp_matrix[ ,colSums(exp_matrix) < max]
  }
  return(exp_matrix)
}


feature_filter <- function(exp_matrix, min, max = 500000, exp_level = 0) {
  # filter out columns (cells) with < number of genes expressed at "exp_level"
  # min = minimum number of genes expressed
  # max = maxiumum number of genes expressed
  # exp_level threshold default > 0 count, which means at least 1 count
  
  exp_matrix = exp_matrix[ ,apply(exp_matrix, 2, 
                            function(x) length(x[x > exp_level])) > min]
  if(max != 50000) {
    exp_matrix = exp_matrix[ ,apply(exp_matrix, 2, 
                                  function(x) length(x[x > exp_level])) < max]
  }
  return(exp_matrix)  
} 


filter_cells_by_genes <- function(exp_matrix, gene_ids, th = 0) {
  # filter out cells that have sum of expression for the selected genes (gene_ids)
  # equal or higher than the threshold.
  # Returns the expression matrix removing the cells that have the sum of the 
  # expression of selected genes higher than the threshold
  
  exp_matrix = exp_matrix[gene_ids, ]
  sum_exp = colSums(exp_matrix)
  cell_ids = which(sum_exp >= th)
  return(exp_matrix[, -cell_ids])
}



##############################################
# 4) Filtering QC cells and rows with Scater #
##############################################


scater_QC_metrics <- function(sce) {
  # Calculate QC metrics for sce object
  # get mitochondrial genes
  
  is.mito = grepl("^MT-|^mt-", rownames(sce))
  is.ribo = grepl("^RPL\\d+$|^RPL\\d+[A-Z]$|^RPS\\d+$|^RPS\\d+[A-Z]$|^RPS[A-Z]$",
                  rownames(sce))
  sce = calculateQCMetrics(sce, feature_controls = list(Mt = is.mito, Rb = is.ribo))
  return(sce)
}


scater_hist_QC <- function(sce) {
  # Histograms of QC metrics
  data = data.frame(counts = sce@colData$total_counts, genes = sce@colData$total_features,
                    mito = sce@colData$pct_counts_Mt,
                    ribo = sce@colData$pct_counts_Rb)
  # UMIs (Library size)
  hist_umi = ggplot(data, aes(x = counts)) +
               geom_histogram(col="black", fill="white") +
               xlab("Number of UMIs") +
               ylab("Number of cells") +
               theme_bw() +
               theme(axis.text = element_text(color="black"))
  
  # Gene number (N expressed genes)
  hist_gene = ggplot(data, aes(x = genes)) +
    geom_histogram(col="black", fill="white") +
    xlab("Number of genes expressed") +
    ylab("Number of cells") +
    theme_bw() +
    theme(axis.text = element_text(color="black"))
  
  # Gene number (N expressed genes)
  hist_mito = ggplot(data, aes(x = mito)) +
                geom_histogram(col="black", fill="white") +
                xlab("Mitochondrial genes proportion (%)") +
                ylab("Number of cells") +
                theme_bw() +
                theme(axis.text = element_text(color="black"))

  # Gene number (N expressed genes)
  hist_ribo = ggplot(data, aes(x = ribo)) +
                geom_histogram(col="black", fill="white") +
                xlab("Ribosomal genes proportion (%)") +
                ylab("Number of cells") +
                theme_bw() +
                theme(axis.text = element_text(color="black"))
 return(list(hist_umi, hist_gene, hist_mito, hist_ribo))
}


scater_plot_QC_mt_rb_counts <- function(sce, mt_lines = c()) {
  # Scatter plot for transcript number vs mitochondrial proportion and ribosomal
  # to show the thresholds for selecting cells
  
  # density mito
  colors = colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))
  d_mito = densCols(sce@colData$total_counts, sce@colData$pct_counts_Mt,
                    colramp = colors)
  # density ribo
  colors = colorRampPalette(rev(brewer.pal(n = 10, name = "Spectral")))
  d_ribo = densCols(sce@colData$total_counts, sce@colData$pct_counts_Rb,
                    colramp = colors)
  
  data = data.frame(counts = sce@colData$total_counts,
                    mito = sce@colData$pct_counts_Mt,
                    ribo = sce@colData$pct_counts_Rb,
                    d_mito = d_mito,
                    d_ribo = d_ribo)
  # scatter_plot mito vs transcripts
  scatter_mito = ggplot(data) +
                   geom_point(aes(counts, mito, col = d_mito), size = 1) +
                   scale_color_identity() +
                   xlab("Number of UMIs") +
                   ylab("Mitochondrial proportion (%)") +
                   ylim(c(0,100)) +
                   theme_bw() +
                   theme(axis.text = element_text(color="black"))
  # adding dashed lines for the chosen filtering threshold 
  if(length(mt_lines) > 0) {
    scatter_mito = scatter_mito +
      geom_vline(xintercept=c(mt_lines[1], mt_lines[2]), linetype="dotted") +
      geom_hline(yintercept=c(mt_lines[3]), linetype="dotted")
  }
  
  # scatter_plot ribo vs transcripts
  scatter_ribo = ggplot(data) +
                   geom_point(aes(counts, ribo, col = d_ribo), size = 1) +
                   scale_color_identity() +
                   xlab("Number of UMIs") +
                   ylab("Ribosomal proportion (%)") +
                   ylim(c(0,100)) +
                   theme_bw() +
                   theme(axis.text = element_text(color="black"))  
  return(list(scatter_mito, scatter_ribo))
}


scater_plot_QC <- function(sce, N = 50) {
  # Top N expressed genes
  
  fontsize <- theme(axis.text=element_text(size=10), axis.title=element_text(size=12))
  plotQC(sce, type = "highest-expression", n=N) + fontsize
}


scater_pca_QC <- function(sce) {
  # plot pca based on QC metrics. Study better this!
  fontsize <- theme(axis.text=element_text(size=10), axis.title=element_text(size=12))
  plotPCA(sce, pca_data_input="pdata") + fontsize
}


scater_summary_QC <- function(sce) {
  # From a given sce object, just get the summary of QC metrics and arrange
  # in a table to look prettier.
  sum_counts = as.vector(summary(sce@colData$total_counts))
  sum_feat = as.vector(summary(sce@colData$total_features))
  sum_mt = as.vector(summary(sce@colData$pct_counts_Mt))
  sum_rb = as.vector(summary(sce@colData$pct_counts_Rb))
  sum_sd = c(sd(sce@colData$total_counts), sd(sce@colData$total_features),
             sd(sce@colData$pct_counts_Mt),
             sd(sce@colData$pct_counts_Rb))
  sum_dat = data.frame(stat = c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "max"),
                       counts = sum_counts, feature = sum_feat,
                       mito = sum_mt, ribo = sum_rb)
  # add standard deviation row
  add_sd = data.frame(stat = "sd", counts = sum_sd[1], feature = sum_sd[2],
                      mito = sum_sd[3], ribo = sum_sd[4])
  sum_dat = rbind(sum_dat, add_sd)
  return(sum_dat)
}

###################################
# Scater filtering cells: Columns #
###################################

scater_filter_cells <- function(sce, threshold_cf=3, threshold_mt=1.5) {
  # Plot library sizes (transcripts) and number of genes per cell
  
  libsize.drop <- isOutlier(sce$total_counts, nmads=threshold_cf, type="lower", log=TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads=threshold_cf, type="lower", log=TRUE)
  mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=threshold_mt, type="higher")
  # filtered dataset
  sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
  # summary of filtered
  print(data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
                           ByMito=sum(mito.drop), Remaining=ncol(sce)))
  return(sce)
}

# filter cells and genes that are too highly expressed?? representing most of the expression, like > 20%
scater_filter_genes <- function(sce, threshold = 1, N = 10) {
  # Filter out lowly expressed genes (threshold = 1 count(s) in at least N cells)
  numcells <- nexprs(sce, byrow=TRUE, threshold = threshold)
  keep <- numcells >= N
  sce = sce[keep, ]
  return(sce)
}

###########################################
# 5) Functions for merging dataframes     #
###########################################


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


#############################
# Impute data with DrImpute #
#############################

dr_impute <- function(exp_matrix, k1=5, k2=10) {
  # impute data with DrImpute.
  # exp_matrix: Filtered expression matrix, removing bad cells, and bad genes?
  # k1 and k2 are parameters for K-neighbors used by DrImpute.
  library(DrImpute)
  
  cell_ids = colnames(exp_matrix)
  exp_imp <- DrImpute(as.matrix(exp_filter), ks=k1:k2)
  exp_imp <- as.data.frame(exp_imp)
  colnames(exp_imp) <- cell_ids
  return(exp_imp)
}

#######################################################################
# Functions to work with average expression level of single cell data #
#######################################################################

summary_exp <- function(norm_exp_matrix, raw_exp_matrix, threshold = 0, norm = F,
                        suffix = "S1") {
  # for an expression matrix, return the average normalized expression 
  # across all cells, the sum of raw counts, the frequency of cells in which
  # the gene is expressed above threshold for raw or normalized expression.
  # also store the standard deviation.
  
  # get average expression of normalized data
  avg_norm = round(rowMeans(norm_exp_matrix), 4)
  # get standard deviation of avg_norm
  std_norm = round(apply(norm_exp_matrix, 1, sd), 4)
  # avg count
  avg_count = round(rowMeans(raw_exp_matrix), 4)
  # sum of raw counts
  sum_count = rowSums(raw_exp_matrix[, colnames(norm_exp_matrix)])
  # frequency of expression
  if(norm == T) {
    freq = apply(norm_exp_matrix, 1,
                 function(x) (sum(x > threshold) / ncol(norm_exp_matrix)))
  } else {
    freq = apply(raw_exp_matrix, 1,
                 function(x) (sum(x > threshold) / ncol(raw_exp_matrix)))
  }
  freq = round(freq, 3)
  
  data = data.frame(avg_norm = avg_norm, std_norm = std_norm, avg_count = avg_count,
                    sum_count = sum_count, freq = freq)
  rownames(data) <- rownames(norm_exp_matrix)
  
  if(norm == T) {
    freq_col = "freq_norm"
  } else {
    freq_col = "freq_count"
  }
  # fixing colnames
  colnames(data)[5] <- freq_col
  colnames(data) <- paste(colnames(data), suffix, sep = "_")
  
  return(data)
}


