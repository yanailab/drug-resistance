library("iemisc")

zavit <- function(exp_matrix) {
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
  
  t = atan2d(X, Y) %% 360
  t_ordered = sort(t, decreasing = F)
  
  exp_tpm_ord = exp_matrix[names(t_ordered), ]
  return(exp_tpm_ord)
}

