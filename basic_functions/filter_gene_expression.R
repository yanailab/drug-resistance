##################################################
# Functions for filtering gene expression matrix #
##################################################
                                                
# filter_gene_expression.R                         
#                                 
# Author: Gustavo S. Franca                      

tpm <- function(x, fact = 10^6) {
  # normalize gene expression matrix by TPM: 
  # (gene exp / sum of expression of all genes in sample)):
  # x/colsum(x) * 10^6 (counts per million)
  # Args:
  #   x: gene expression matrix as dataframe
  norm_x = sweep(x, 2, colSums(x), '/') * (fact)
  return(norm_x)
}

tpm_median <- function(x) {
  # normalize gene expression matrix by TPM: 
  # (gene exp / sum of expression of all genes in sample)):
  # x/colsum(x) * 10^6 (counts per million)
  # Args:
  #   x: gene expression matrix as dataframe
  norm_x = sweep(x, 2, colSums(x), '/') * median(colSums(x))
  return(norm_x)
}

median_norm <- function(x) {
  # normalize an expression matrix by median normalization: same method
  # of Flo's smoothing. 
  # Args:
  #   x: expression matrix (raw counts)
  scale_factor = median(colSums(x)) / colSums(x)
  norm_x = sweep(x, 2, scale_factor, "*")
  return(norm_x)
}

anscombe_transform <- function(x) {
  # use anscombe transformation on expression data
  ansc_x = 2*(sqrt(x + 3/8.0))
  return(ansc_x)
}

freeman_transform <- function(x) {
  # use freeman-tukey tranformation on expression data
  freeman_x = sqrt(x) + sqrt(x + 1)
  return(freeman_x)
}

mean_scale <- function(x, row = T) {
  # from a dataframe, scale each value by taking value - mean(xi), where xi
  # is the row (dim = 1) or column (dim = 2). Purpose is to scale by mean.
  if(row == T) {
    norm_x = sweep(x, 1, rowMeans(x), "-")
  } else {
    norm_x = sweep(x, 2, colMeans(x), "-")
  }
  return(norm_x)
}

mean_rm_zero <- function(x) {
  # create the mean function giving NaN for expression = 0 in 2 element vector
  # Args:
  #   x: numeric vector
  if ((length(x) == 2) && (any(x == 0))) {
    return(0)
  }
  else {
    return(mean(x))
  }
}

std <- function(x) {
  # create the sd as function (sd is not a function in R?))
  # Args:
  #   x: numeric vector
  if ((length(x) == 2) && (any(x == 0))) {
    return(0)
  }
  else {
    return(sd(x))
  }
}

fano_factor <- function(x) {
  # calculate the fano factor (variance/mean)
  # Args:
  #   x: numeric vector
  f = var(x) / mean(x)
  return(f)
}

cv <- function(x) {
  # calculate the coefficient of variation (sd/mean)
  # Args:
  #   x: numeric vector
  if ((length(x) == 2) && (any(x == 0))) {
    return(NaN)
  }
  else {
    return((std(x) / mean(x))*100)
  }
}

gene_filter <- function(exp_matrix, exp_level=1, n_samples=1) {
  # subsets the raw data matrix to keep only expressed genes
  # Like removing the "bad rows"
  # Args:
  #   dataset: the original gene expression matrix read in
  #   exp_level: threshold for minimum expression level (default = 1)
  #   n_samples: threshold for minimum n of samples with exp_level (default = 1)
  
  exp_matrix = exp_matrix[rowSums(exp_matrix > exp_level) >= n_samples, ]
  # just an alert if no rows matches the criteria
  if (nrow(exp_matrix) == 0) {
    print("Empty Rows! No rows matched the criteria.")
  }
  return(exp_matrix)
}

gene_filter_by <- function(exp_matrix, FUN, n=1) {
  # subsets the raw data matrix to keep only genes expressed above a threshold
  # Like removing the "bad rows" based on a statistic (ex: mean, median, std)
  # For example, calculates the mean expression for all genes in the matrix and
  # just keep those that are above the mean.
  # Args:
  #   exp_matrix: the original gene expression matrix read in
  #   FUN: the statistic function to define the threshold for minimum expression
  #        level.
  #   n = Number of times * threshold (ex: for n * Standard Deviations)
  
  # gives the stat (ex: mean, median, sd) for all genes across samples
  threshold = FUN(as.matrix(exp_matrix))
  # calculate stat for each gene across all samples (rows = 1)
  genes = apply(exp_matrix, 1, FUN)
  # filter the vector of genes above the stat
  genes_filtered = genes[genes > threshold * n, drop=F]
  # subset the original data set for those genes above the stats by gene names
  exp_matrix = exp_matrix[names(genes_filtered), ]
  return(exp_matrix)
}

sample_filter <- function(exp_matrix, exp_level=1) {
  # subsets gene expression matrix to keep samples with reliable expression
  # levels (threshold will depend if it is in counts, tpm, rpkm, etc.)
  # Like removing the "bad columns"
  # Args:
  #   exp_matrix: the original gene expression matrix read in
  #   exp_level: threshold for the sum of expression level in the sample
  #              (default = 1). If the values are counts, the threshold
  #              will be usually a high value (> 100,000 or > 30,000)

  exp_matrix = as.data.frame(exp_matrix[ ,colSums(exp_matrix) > exp_level, drop=F])
  
  if (ncol(exp_matrix) == 0) {
    print("Empty Columns! No columns matched the criteria.")
  }
  return(exp_matrix)
}

gene_order_by <- function(exp_matrix, FUN) {
  # calculates the specific statistic (Fano Factor, Coefficient Variation)
  # and reorder the dataframe based on these values (highest to lowest)
  # Args:
  #   exp_matrix: the original gene expression matrix read in
  #   FUN: a function to be applied in genes (ex. fano factor, CV, mean, median, std)
  # calculate fano_factor and sort them
  
  # return a vector of genes sorted by the result of FUN 
  exp_sort = sort(apply(exp_matrix, 1, FUN), decreasing=T)
  # reorder the expression matrix based on the order of genes above
  exp_ordered = exp_matrix[names(exp_sort), ]
  return(exp_ordered)  
}

tpm_to_raw <- function(exp, rownames_exp_tpm) {
  # Function to get back the raw expression matrix from 
  # the expression matrix in tpm that was filtered
  # Args:
  #   exp: raw expression matrix with all entries
  #   rownames_exp_tpm: rownames (genes) from the tpm exp matrix
  #   returns expression matrix raw with the same genes in tpm.
  
  exp_raw = exp[row.names(exp) %in% rownames_exp_tpm, ]
  return(exp_raw)
}

summarize_replicates <- function(exp_matrix, FUN=mean, N=2) {
  # calculates a function (mean, media, sd, var, etc) for summarizing 
  # replicates to reduce the expression matrix columns to a single
  # observation for each group of replicates. Example: Given 3 reps for
  # a sample, stay with the mean between them.
  # Args:
  #   exp_matrix: the expression matrix for genes (rows) and samples(cols)
  #   FUN: a function to be applied in genes (ex. mean, median, sd, var)
  #        for the respective replicates for each sample.
  #     N: number of adjacent columns to get. Default 2 (will get 2 columns)  
  
  summ_exp_matrix = c()
  for (i in seq(1, ncol(exp_matrix), by=N)) {
    # build column indices
    cols = seq(i, (i + N) - 1, by=1)
    # apply a function for each row for specified columns
    res = apply(exp_matrix[, cols], 1, FUN)
    # bind the current column with previous ones
    summ_exp_matrix = cbind(summ_exp_matrix, res)
    # get column names to join as a single column name
    colname = paste(colnames(exp_matrix[, cols]), collapse="_", sep="_")
    # renaming last column name
    colnames(summ_exp_matrix)[ncol(summ_exp_matrix)] <- colname
  }
  return(summ_exp_matrix)
}

summarize_replicates_jump <- function(exp_matrix, FUN=mean,
                                           N_elements=3, by=1,
                                           N_cols=9) {
  # calculates a function (mean, media, sd, var, etc) for summarizing 
  # replicates to reduce the expression matrix columns to a single
  # observation for each group of replicates. Example: Given 3 reps for
  # a sample, stay with the mean between them.
  # Args:
  #   exp_matrix: the expression matrix for genes (rows) and samples(cols)
  #   FUN: a function to be applied in genes (ex. mean, median, sd, var)
  #        for the respective replicates for each sample.
  #   N_elements: number of columns to summarize (same as number of replicates)
  #   by: steps to skip, usually = 1
  #   N_cols: Number of columns to have in the end of the final matrix
  
  summ_exp_matrix = c()
  for (i in seq(1, ncol(exp_matrix) / N_elements, by=by)) {
    # build column indices
    cols = seq(i, ncol(exp_matrix), by=N_cols)
    # apply a function for each row for specified columns
    res = apply(exp_matrix[, cols], 1, FUN)
    # bind the current column with previous ones
    summ_exp_matrix = cbind(summ_exp_matrix, res)
    # get column names to join as a single column name
    colname = paste(colnames(exp_matrix[, cols]), collapse="_", sep="_")
    # renaming last column name
    colnames(summ_exp_matrix)[ncol(summ_exp_matrix)] <- colname
  }
  return(summ_exp_matrix)
}

zavit <- function(exp_matrix, angle = 360) {
  # Order genes according to expression profiles in samples
  # Itai's version of Zavit script in R. Usually useful for
  # ordering genes expressed according to each time point/sample
  # Have to study better the theory behind that (PCA), angles, etc.
  # Returns the ordered expression matrix to be ploted in a heatmap for ex.
  # Args:
  #   exp_matrix: the expression matrix for genes (rows) and samples(cols)
  #               normalized in tpm (log2 + 1)
  
  exp_tpm_scale = t(scale(t(exp_matrix))) # scale by row, thats why t(t())
  # principal component with scaled matrix
  exp_tpm_pca = princomp(exp_tpm_scale, cor=F)
  
  X = exp_tpm_pca$scores[, c("Comp.1")]
  Y = exp_tpm_pca$scores[, c("Comp.2")]
  
  library("iemisc")
  t = atan2d(X, Y) %% angle
  t_ordered = sort(t, decreasing = T)
  
  exp_tpm_ord = exp_matrix[names(t_ordered), ]
  return(exp_tpm_ord)
}


filter_gene_exp_time_course <- function(exp_matrix, exp_th = 0.1, fold_th = 2) {
  # filter out genes in an expression matrix that have avg expression lower than
  # exp_th and fold change expression in the highest expressed sample divided
  # by the lowest expressed sample < fold_th.
  
  # filter genes expressed < 0.1
  exp_matrix = gene_filter(exp_matrix = exp_matrix, exp_level = exp_th, n_samples = 1)
  
  # function to calculate fold
  fold_ch <- function(x) {
    x_sort = sort(x, decreasing = T)
    pass = (x_sort[1] / x_sort[length(x_sort)]) > fold_th
    return(pass)
  }
  # genes that pass filter are TRUE, use them to select from exp_matrix
  genes_select = apply(exp_matrix, 1, fold_ch)
  exp_matrix = exp_matrix[genes_select, ]
  return(exp_matrix)
}


filter_gene_var_time_course <- function(exp_matrix, var_th = 0) {
  # filter out genes with variance lower than var_th
  
  var_genes = apply(exp_matrix, 1, var)
  var_genes = var_genes[var_genes > var_th]
  exp_matrix = exp_matrix[names(var_genes), ]
  return(exp_matrix)
}

#################################
# Typical workflow              #
# Vignette for multiple plots   #
#################################

# Examples
# x =         
#            sample_0001 sample_0002 sample_0003
# 2L52.1            0           2           7
# 2RSSE.1           0           0           2
# 2RSSE.2           1           0           0
# 3R5.1             7           4           5

# -- Filtering data --
# by genes not expressed
# a = gene_filter(x, exp_level=X, n_samples=Y) # remove genes (rows) not expressed 
                          # at X level in at least Y samples (default = 1 for both)

# by bad samples
# a = sample_filter(x, exp_level=X) # remove columns whose sum is less than X

# by genes not expressed below a threshold from some statistic (mean, median, std)
# a = gene_filter_by(x, FUN (mean, median, std), n = 1) # get the mean, median or std
                                        # of all genes, and filter by this threshold. 
                                        # n is the number of times (ex: 2*std)

# -- Ordering data --
# by ordering the expression matrix based on Fano factor or CV, mean, median, std
# a = gene_order_by(x, fano_factor)

# -- Normalizing data --
# Normalize by TPM (not using gene length normalization - good for CEL-Seq protocol)
# a = tpm(x) # will return expression matrix in tpm

# -- Calculate mean expression across replicates --
# get the respective columns of replicates and returns the dataframe with the means
# default: Will get pairs of adjacent columns as replicates
# exp_tpm_mean = summarize_replicates(exp_tpm)

