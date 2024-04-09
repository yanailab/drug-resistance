
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)

source("/home/gu/R/myfunctions/filter_gene_expression.R")
source("/home/gu/R/myfunctions/differential_expression.R")

setwd("/path/kuramochi_olaparib_evolution/atac-seq")

#load("atac_deseq2.rda")
load("atac_deseq2_downsample.rda")

# data for adapted, ctrl and persisters - downsampled to same read depth
# have to remove the starting "#" in header, otherwise colnames are messed up
count_table = read.table("consensus_peaks_no_summit_counts_down_all.tab", header = T, sep = "\t")
annotation_table = read.csv("consensus_peaks_homer_annotation_down_all.txt", header = T, sep = "\t")

colnames(annotation_table)[1] <- "PeakID"
# fix start because homer made the coords 1-based and we want 0-based
annotation_table$Start <- annotation_table$Start - 1

# fix chr name
count_table$chr <- paste("chr", count_table$chr, sep = "")
rownames(count_table) <- paste(count_table$chr, count_table$start,
                               count_table$end, sep = "_")

coords = paste(annotation_table$Chr, annotation_table$Start,
               annotation_table$End, sep = "_")
rownames(annotation_table) <- coords
peak_ids = annotation_table[, 1:4]
peak_ids = as.character(peak_ids[rownames(count_table), ]$PeakID)

annotation_table = annotation_table[rownames(count_table), ]

rownames(count_table) <- peak_ids

##########################
# Correlating replicates #
##########################

counts = count_table[, 4:ncol(count_table)]

# tpm 
counts_norm = log2(tpm(counts) + 1)

var_peaks = apply(counts_norm, 1, fano_factor)
var_peaks = sort(var_peaks, decreasing = T)[1:2000]

#correlation 
cormM = cor(counts_norm, method="spearman")

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
             heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                         title = "Spearman",
                                         labels_gp = gpar(fontsize = 6),
                                         legend_height = unit(2, "cm")))

draw(h2)

##############################
# DESEQ2 diff. accessibility #
##############################

library(DESeq2)

# counts table
counts = count_table[, 4:ncol(count_table)]
# metadata
columnData = data.frame(condition = c(rep(c("C", "T5", "T10", "T20", "T40", "T80", "T160", "T320"), 2),
                                      "P10", "P10", "P320", "P320"),
                        batch = c(rep("b1", 8), rep("b2", 8), rep("b3", 4)))

rownames(columnData) <- colnames(counts)

# DESEq obj
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = columnData,
                              design = ~ condition)

dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# including persisters
dds$condition <- factor(dds$condition, levels = c("C", "T5", "T10", "T20", "T40",
                                                  "T80", "T160", "T320", "P10", "P320"))
dds <- DESeq(dds)
res <- results(dds)
res

res_sample <- results(dds, contrast = c("condition", "C", "P320"), alpha = 0.05)
summary(res_sample)

res_sample = res_sample[(res_sample$padj < 0.01) & (res_sample$log2FoldChange <= -log2(3)), ]
dim(res_sample)

vsd <- vst(dds)
pca_plot = plotPCA(vsd, "condition")

plotPCA(vsd, "condition", ntop = 2000)

new_colors = c("#8DA3A6","#4E78C4","#57A2AC","#7EB875",
               "#D0B541", "#E67F33", "#CE2220", "#80271E")

new_colors_per = c("#8DA3A6","#4E78C4","#57A2AC","#7EB875",
                   "#D0B541", "#E67F33", "#CE2220", "#80271E")
pca_plot + 
  #scale_color_manual(values=new_colors) +
  theme_minimal()

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

save(dds, count_table, annotation_table, res, file = "atac_deseq2_downsample.rda")


########################################
# Saving bed file for motif enrichment #
########################################


select_diff_peaks <- function(dds, samples, ref_sample_index = 1, up_treat = T,
                              padj = 0.001, fold = 3) {
  # selects the up regulated genes across all conditions and returns a list
  # of dataframes with diff. upregulated sets.
  # dds: DESeq2 object with results 
  # samples: vector with sample names. Ex: c("C", "T5", "T10", ...)
  # ref_sample_index: Is the "Ctrl" or the sample as reference for diff. expression
  # up_treat: If the upregulation is in the treated or control: by default 
  # Deseq apparently reports the Negative fold changes in the treated and positive
  # in the controls, or query and subject
  # padj: adjusted p-val threshold
  # fold: minimum fold change 
  
  selected_peaks = list()
  for(i in samples) {
    # skipping the comparison with itself
    if(i == samples[ref_sample_index]) {
      next
    } else {
      condition = c("condition", samples[ref_sample_index], i)
      res_sample <- DESeq2::results(dds, contrast = condition, alpha = 0.05) # get only signficant
      # selection criteria
      if(up_treat == T) {
        # take the negative fold changes
        res_sample = res_sample[(res_sample$padj < padj) & (res_sample$log2FoldChange <= -log2(fold)), ]
      } else {
        # take the positive fold changes
        res_sample = res_sample[(res_sample$padj < padj) & (res_sample$log2FoldChange >= log2(fold)), ]        
      }
      # order by p-val
      res_sample = res_sample[order(res_sample$padj), ]
      res_sample = as.data.frame(res_sample)
      comparison = paste(samples[ref_sample_index], i, sep = "_")
      #print(head(res_sample))
      selected_peaks[[comparison]] <- res_sample
    }
  }
  return(selected_peaks)
}


create_bed_file_peaks <- function(selected_peaks, bed_peaks, up_treat = T) {
  # create bed file tables for downstream purposes (ex: HOMER motif enrichment)
  # selected_peaks: differentially regulated peaks list
  # bed_peaks: all peaks in bed file format for homer
  
  comparisons = names(selected_peaks)
  for(c in comparisons) {
    diff_peaks = rownames(selected_peaks[[c]])
    bed_out = bed_peaks[bed_peaks$V4 %in% diff_peaks, ]
    
    # invert the name for first T5 then C
    # fix the sample order just to be more meaningful 
    comp = stringr::str_split(c, "_", simplify = T)
    comp = paste(comp[2], comp[1], sep = "_")
    if(up_treat == T) {
      file_name = paste("up_peaks", comp, "all.bed", sep = "_")
    } else {
      file_name = paste("down_peaks", comp, "all.bed", sep = "_")      
    }
    write.table(bed_out, file = file_name, col.names = F, row.names = F,
                sep = "\t", quote = F)
  }
}

# selecting comparisons all vs ctrl
peaks_selected = select_diff_peaks(dds, c("C", "T5", "T10", "T20", "T40",
                                          "T80", "T160", "T320", "P10", "P320"),
                         up_treat = T, padj = 0.001)

# writing files for homer motif enrichment
bed_peaks = read.table("path/consensus_peaks_no_summit_center_homer_down.bed", head = F)
create_bed_file_peaks(peaks_selected, bed_peaks, up_treat = F)

create_bed_file_peaks(peaks_selected, bed_peaks, up_treat = T)


################################################################################
# Heatmap for differentially accessible regions                                #
################################################################################


# selecting comparisons all vs ctrl: Upregulated
peaks_selected_up = select_diff_peaks(dds, c("C", "T5", "T10", "T20", "T40",
                                          "T80", "T160", "T320", "P10", "P320"),
                                   up_treat = T, padj = 0.001)
# select upregulated peaks
all_up_peaks = lapply(peaks_selected_up, rownames)
# subset to exclude the persisters
#all_up_peaks = all_up_peaks[!(names(all_up_peaks) %in% c("C_P10", "C_P320"))]
# merge peak names
all_up_peaks = unique(unlist(all_up_peaks))
# counts table from DESEQ, get the normalized table
counts_norm = counts(dds, normalized=T)
# select only the up peaks and from adapted lines
counts_norm_up = counts_norm[all_up_peaks, 1:(ncol(counts_norm))]
# including adapted + persisters
counts_norm_up = counts_norm[all_up_peaks, ]
# reorder the columns 
# reorder the columns including persisters
counts_norm_up = counts_norm_up[, c(1, 9, 17, 18, 19, 20, 2, 10, 3, 11, 
                                    4, 12, 5, 13, 6, 14, 7, 15, 8, 16)]

counts_norm_up_avg = summarize_replicates(counts_norm_up, mean, 2)
colnames(counts_norm_up_avg) <- c("C", "P10", "P320", "T5", "T10", "T20",
                                  "T40", "T80", "T160", "T320")

# scale z-score
set.seed(1234)
counts_norm_up_avg_zavit = t(scale(t(zavit(counts_norm_up_avg))))
counts_norm_up_avg_zavit = t(scale(t(zavit(counts_norm_up_avg, angle = 270))))

h2 = Heatmap(counts_norm_up_avg_zavit, show_column_names = T, show_row_dend = F, 
             show_column_dend = F, show_row_names = F, 
             name = "",
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(30),
             cluster_rows = F, cluster_columns = F,
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))
draw(h2)

# bubble plot for N of up peaks 
n_peaks_up = sapply(all_up_peaks, length)
n_peaks_up = data.frame(n_peaks = n_peaks_up, Sample = names(n_peaks_up))
n_peaks_up$Sample <- c("T5", "T10", "T20", "T40", "T80", "T160", "T320", "P10", "P320")
n_peaks_up$Sample <- factor(n_peaks_up$Sample, levels = c("P10", "P320", "T5", 
                                                          "T10", "T20", "T40", "T80", "T160", "T320"))
n_peaks_up$Category <- "up"

ggplot(n_peaks_up, aes(x = Sample, y = Category, size = n_peaks_up)) +
  geom_point(aes(size = n_peaks), shape = 21) +
  #scale_size_continuous(range = c(1, 5)) +
  scale_size(range = c(2,12), breaks=c(1000,3000,5000,7000),
             #labels=c("0","60","120","180","240",">=300"),
             guide="legend") +
  ylab("") +
  xlab("") +
  
  theme_minimal() 

###########################
# For Downregulated peaks #
###########################

# selecting comparisons all vs ctrl: Dowregulated
peaks_selected_down = select_diff_peaks(dds, c("C", "T5", "T10", "T20", "T40",
                                          "T80", "T160", "T320", "P10", "P320"),
                                   up_treat = F, padj = 0.001)
# select upregulated peaks
all_down_peaks = lapply(peaks_selected_down, rownames)
# subset to exclude the persisters
#all_down_peaks = all_down_peaks[!(names(all_down_peaks) %in% c("C_P10", "C_P320"))]
# merge peak names
all_down_peaks = unique(unlist(all_down_peaks))
# counts table from DESEQ, get the normalized table
counts_norm = counts(dds, normalized=T)
# select only the down peaks and from adapted lines
counts_norm_down = counts_norm[all_down_peaks, 1:(ncol(counts_norm) - 4)]
# including persisters
counts_norm_down = counts_norm[all_down_peaks, ]
# reorder the columns 
counts_norm_down = counts_norm_down[, c(1, 9, 17, 18, 19, 20, 2, 10, 3, 11, 
                                        4, 12, 5, 13, 6, 14, 7, 15, 8, 16)]


counts_norm_down_avg = summarize_replicates(counts_norm_down, mean, 2)
colnames(counts_norm_down_avg) <- c("C", "P10", "P320", "T5", "T10", "T20",
                                    "T40", "T80", "T160", "T320")

# scale z-score
# counts_norm_down_avg = t(scale(t(counts_norm_down_avg))) for clustering
counts_norm_down_avg_zavit = t(scale(t(zavit(counts_norm_down_avg))))


h2 = Heatmap(counts_norm_down_avg_zavit, show_column_names = T, show_row_dend = F, 
             show_column_dend = F, show_row_names = F, 
             name = "",
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(30),
             cluster_rows = F, cluster_columns = F,
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))
draw(h2)

################################################################################
# Merge Zavit data from Down and upregulated peaks to plot in a single heatmap #
################################################################################

counts_norm_down_up_avg_zavit = rbind(counts_norm_down_avg_zavit,
                                      counts_norm_up_avg_zavit)

direction_ann = c(rep("closed", nrow(counts_norm_down_avg_zavit)), 
                  rep("open", nrow(counts_norm_up_avg_zavit)))

# top annotation for states
dir_ann = rowAnnotation(df = data.frame(direction = direction_ann), 
                        simple_anno_size = unit(1.5, "mm"),
                        col = list(direction = c("closed" = "#984fa3",
                                                 "open" = "#4daf4a")),
                        show_annotation_name = F,
                        gap = unit(0, "mm"),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 8),
                                                       title = "",
                                                       labels_gp = gpar(fontsize = 8)))

h2 = Heatmap(counts_norm_down_up_avg_zavit, show_column_names = T, show_row_dend = F, 
             show_column_dend = F, show_row_names = F, 
             name = "", left_annotation = dir_ann,
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(30),
             cluster_rows = F, cluster_columns = F, border = "black",
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))
draw(h2)

################################################
# bubble plot for total N of up and down peaks #
################################################

n_peaks_down = sapply(all_down_peaks, length)
n_peaks_down = data.frame(n_peaks = n_peaks_down, Sample = names(n_peaks_down))
n_peaks_down$Sample <- c("T5", "T10", "T20", "T40", "T80", "T160", "T320", "P10", "P320")
n_peaks_down$Sample <- factor(n_peaks_down$Sample, levels = c("P10", "P320", "T5", "T10", 
                                                              "T20", "T40", "T80", "T160", "T320"))
n_peaks_down$Category <- "down"

# merge up and down peaks
n_peaks_total = rbind(n_peaks_up, n_peaks_down)

ggplot(n_peaks_total, aes(x = Sample, y = Category, size = n_peaks_up)) +
  geom_point(aes(size = n_peaks), shape = 21) +
  #scale_size_continuous(range = c(1, 5)) +
  scale_size(range = c(2,12), breaks=c(1000,3000,5000,7000),
             #labels=c("0","60","120","180","240",">=300"),
             guide="legend") +
  ylab("") +
  xlab("") +
  
  theme_minimal() 


###########################################
# Upset plot for shared up and down peaks #
###########################################

library(UpSetR)

upset(fromList(all_up_peaks), #sets = c("T320_C", "T10_C", "P320_C", "P10_C"),
      keep.order = T, order.by = "freq")


################################################################################
# Processing HOMER output for TF enrichment                                    #
################################################################################

setwd("/path/kuramochi_olaparib_evolution/atac-seq/deseq2_peaks_all_downsampled/homer_enrichment_down/")
setwd("/path/kuramochi_olaparib_evolution/atac-seq/deseq2_peaks_all_downsampled/homer_enrichment_up/")

filenames <- list.files(pattern='knownResults.txt', recursive=TRUE, full.names = TRUE)

list_data<- lapply(filenames, read.csv, head = T, sep = "\t")

filter_tfs <- function(data, pval_th = 1e-5) {
  # filter TFs with lower pval
  return(data[data$P.value <= pval_th, ])
}

set_rownames <- function(x) {
  rownames(x) <- make.names(x[,1], unique = T)
  return(x)
}

list_data_rownames = lapply(list_data, set_rownames)
all_tfs = unique(unlist(lapply(list_data_rownames, rownames)))

get_values <- function(x, all_tfs, column) {
  x = x[all_tfs, column]
  return(x)
}

#########################
# get log P-val columns #
#########################

tf_data = lapply(list_data_rownames, get_values, all_tfs, 4)

tf_data = as.data.frame(do.call(cbind, tf_data))
col_names = stringr::str_split(filenames, pattern = "/", simplify = T)[, 2]
colnames(tf_data) <- col_names
rownames(tf_data) <- all_tfs
tf_data_pval = tf_data
tf_data_pval = abs(tf_data_pval[, c(1,2,8,3,5,7,9,4,6)])

# adding the CNV overlap/no-overlap ones
# reorder columns and transform to abs values
tf_data_pval = abs(tf_data_pval[, c(8,1,3,7,9,2,6,4,5)])
head(tf_data_pval)

# filter the low high pvals and make the log pval positive: log(pval = 1e-5 = 11.5)
tf_data_pval = tf_data_pval[rowSums(tf_data_pval > 23.02) >= 1,] # 1e-10 in at least one sample
# get the significant ones for percentage subsetting
significant_tfs = rownames(tf_data_pval)

############################
# Get the original p-value #
############################
tf_data = lapply(list_data_rownames, get_values, all_tfs, 3)
# get rid of the fucking %
tf_data = lapply(tf_data, function(x) as.numeric(sub("%", "", x)))

tf_data = as.data.frame(do.call(cbind, tf_data))
col_names = stringr::str_split(filenames, pattern = "/", simplify = T)[, 2]
colnames(tf_data) <- col_names
rownames(tf_data) <- all_tfs
tf_data_orig_pval = tf_data
#tf_data_orig_pval = abs(tf_data_orig_pval[, c(6,1,3,5,7,2,4)])
# adding persisters
tf_data_orig_pval = abs(tf_data_orig_pval[, c(1,2,8,3,5,7,9,4,6)])

# adding the CNV overlap/no-overlap ones
# reorder columns and transform to abs values
tf_data_orig_pval = abs(tf_data_orig_pval[, c(8,1,3,7,9,2,6,4,5)])

tf_data_orig_pval = tf_data_orig_pval[significant_tfs, ]

###############################
# Get the N of TFs with motif #
###############################
tf_data = lapply(list_data_rownames, get_values, all_tfs, 6)
# get rid of the fucking %
tf_data = lapply(tf_data, function(x) as.numeric(sub("%", "", x)))

tf_data = as.data.frame(do.call(cbind, tf_data))
col_names = stringr::str_split(filenames, pattern = "/", simplify = T)[, 2]
colnames(tf_data) <- col_names
rownames(tf_data) <- all_tfs
tf_data_number = tf_data
#tf_data_number = abs(tf_data_number[, c(6,1,3,5,7,2,4)])
# add persisters
tf_data_number = abs(tf_data_number[, c(1,2,8,3,5,7,9,4,6)])

# adding the CNV overlap/no-overlap ones
# reorder columns and transform to abs values
tf_data_number = abs(tf_data_number[, c(8,1,3,7,9,2,6,4,5)])

tf_data_number = tf_data_number[significant_tfs, ]

##########################################
# Get the TFs now with Percentage values #
##########################################
tf_data = lapply(list_data_rownames, get_values, all_tfs, 7)
# get rid of the fucking %
tf_data = lapply(tf_data, function(x) as.numeric(sub("%", "", x)))

tf_data = as.data.frame(do.call(cbind, tf_data))
col_names = stringr::str_split(filenames, pattern = "/", simplify = T)[, 2]
colnames(tf_data) <- col_names
rownames(tf_data) <- all_tfs

# reorder columns and transform to abs values
#tf_data_percent = abs(tf_data[, c(6,1,3,5,7,2,4)])
# add persisters
tf_data_percent = abs(tf_data[, c(1,2,8,3,5,7,9,4,6)])


# adding the CNV overlap/no-overlap ones
# reorder columns and transform to abs values
tf_data_percent = abs(tf_data[, c(8,1,3,7,9,2,6,4,5)])
# filter to keep just the significant ones
tf_data_percent = tf_data_percent[significant_tfs, ]

# saving objects for later to save
tf_data_percent_up = tf_data_percent
tf_data_pval_up = tf_data_pval
tf_data_number_up = tf_data_number
tf_data_orig_pval_up = tf_data_orig_pval

tf_data_percent_down = tf_data_percent
tf_data_pval_down = tf_data_pval
tf_data_number_down = tf_data_number
tf_data_orig_pval_down = tf_data_orig_pval


################################################################################
# Clustering TFs based on patterns of gain / pvalue?
################################################################################

fix_tf_names <- function(tf_data) {
  tf_names = strsplit(rownames(tf_data), split = "\\.\\.")
  tf_names = sapply(tf_names, function(x) x[[1]])
  # deal with duplicated names if any
  tf_names = make.unique(tf_names, sep = "_")
  rownames(tf_data) <- tf_names
  return(tf_data)
}

tf_data_percent_up_fix = fix_tf_names(tf_data_percent_up)
tf_data_percent_up_fix = t(scale(t(tf_data_percent_up_fix), center = T, scale = T))
colnames(tf_data_percent_up_fix) <- c("T5", "T10", "T20", "T40", "T80", "T160", "T320", "T320_no_CNV", "T320_CNV")

tf_data_percent_down_fix = fix_tf_names(tf_data_percent_down)
tf_data_percent_down_fix = t(scale(t(tf_data_percent_down_fix), center = T, scale = T))
colnames(tf_data_percent_down_fix) <- c("T5", "T10", "T20", "T40", "T80", "T160", "T320")

# reorder cluster columns
test = dist(t(tf_data_percent_up_fix))
dend = as.dendrogram(hclust(test, method = "average"))
dend <- dendextend::click_rotate(dend, continue = TRUE)

pdf("")
h2 = Heatmap(tf_data_percent_up_fix, show_column_names = T, show_row_dend = T, 
             show_column_dend = T, show_row_names = T, 
             name = "", row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(30),
             cluster_rows = T, cluster_columns = dend,
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                         title = "",
                                         labels_gp = gpar(fontsize = 6),
                                         legend_height = unit(2, "cm")))

draw(h2)


tf_data_center = t(scale(t(tf_data), center = T, scale = F))
tf_data_center


#########################################
# Transcription factor enrichment plots #
#########################################

create_tf_data <- function(tf_data) {
  tf_data$tf <- rownames(tf_data)
  tf_names = strsplit(rownames(tf_data), split = "\\.\\.")
  tf_names = sapply(tf_names, function(x) x[[1]])
  # deal with duplicated names if any
  tf_names = make.unique(tf_names, sep = "_")
  tf_data$tf <- tf_names
  return(reshape2::melt(tf_data))
}

# creating dataframes for each variable: percentage of tf, pval and Number
tf_data1 = create_tf_data(tf_data_percent_up)
tf_data2 = create_tf_data(tf_data_pval_up)
tf_data3 = create_tf_data(tf_data_number_up)
tf_data4 = create_tf_data(tf_data_orig_pval_up)

tf_data5 = create_tf_data(tf_data_percent_down)
tf_data6 = create_tf_data(tf_data_pval_down)
tf_data7 = create_tf_data(tf_data_number_down)
tf_data8 = create_tf_data(tf_data_orig_pval_down)

# merge percent and pval
tf_data_merged_up = merge(merge(tf_data1, tf_data2, by = c("tf", "variable")), 
                          tf_data3, by = c("tf", "variable"))
tf_data_merged_up = merge(tf_data_merged_up, tf_data4, by = c("tf", "variable"))

colnames(tf_data_merged_up) <- c("tf", "sample", "percent", "pval", "n_tf", "orig_pval")
tf_data_merged_up$direction <- "open"

## downregulated
tf_data_merged_down = merge(merge(tf_data5, tf_data6, by = c("tf", "variable")), 
                                  tf_data7, by = c("tf", "variable"))
tf_data_merged_down = merge(tf_data_merged_down, tf_data8, by = c("tf", "variable"))

colnames(tf_data_merged_down) <- c("tf", "sample", "percent", "pval", "n_tf", "orig_pval")
tf_data_merged_down$direction <- "closed"

# merging both up and down peaks
tf_data_merged_all = merge(tf_data_merged_up, tf_data_merged_down, all = T)
# when not merging up and up
tf_data_merged_all = tf_data_merged_up
# down
tf_data_merged_all = tf_data_merged_down

# fix the original p-values for values that are so low that took 0. Add the minimum
# as possible.
x = sort(tf_data_merged_all$orig_pval)
x = x[x != 0]
min_val = x[1]

tf_data_merged_all$orig_pval[tf_data_merged_all$orig_pval == 0] <- min_val

# picking representative TFs for groups - paper can use TEAD, RUNX, ETS, generic
subset_tfs_direction <- function(data_all_tf, tfs_up, tfs_down) {
  # subsets open and closed TFs for plotting
  data_sub = data_all_tf[((data_all_tf$tf %in% tfs_up) & (data_all_tf$direction == "open")) |
                         ((data_all_tf$tf %in% tfs_down) & (data_all_tf$direction == "closed")), ]
  # reordering the same way
  data_sub$tf <- factor(data_sub$tf, levels = c(tfs_up, tfs_down))
  data_sub$direction <- factor(data_sub$direction, levels = c("open", "closed"))
  return(data_sub)
}

test = subset_tfs_direction(tf_data_merged_all, c("Atf4.bZIP", "Nrf2.bZIP", "AP.1.bZIP",
                                                  "ETS1.ETS", "RUNX.Runt"),
                                                c("Pax8.Paired.Homeobox",
                                                  "WT1.Zf", "Sox17.HMG",
                                                  "TEAD.TEA", "NF1.CTF"))
# fixing annoying TF names
x = as.character(test$tf)
x = stringr::str_replace_all(x, c("Atf4.bZIP" = "ATF4", "Nrf2.bZIP" = "NRF2",
                                  "AP.1.bZIP" = "AP1", "Pax8.Paired.Homeobox" = "PAX8",
                                  "WT1.Zf" = "WT1", "Sox17.HMG" = "SOX17",
                                  "ETS1.ETS" = "ETS1", "RUNX.Runt" = "RUNX",
                                  "TEAD.TEA" = "TEAD", "NF1.CTF" = "NF1"))
test$tf <- x
test$tf <- factor(test$tf, levels = c("ATF4", "NRF2", "AP1", "PAX8", "WT1", "SOX17"))
test$tf <- factor(test$tf, levels = c("AP1", "ETS1", "RUNX", "ATF4", "NRF2",
                                      "SOX17", "TEAD", "WT1", "NF1", "PAX8"))
# fixing sample names
y = as.character(test$sample)
y = stringr::str_replace_all(y, c("up_T5_C" = "T5", "up_T10_C" = "T10",
                                  "up_T20_C" = "T20", "up_T40_C" = "T40",
                                  "up_T80_C" = "T80", "up_T160_C" = "T160",
                                  "up_T320_C" = "T320", "up_P10_C" = "P10",
                                  "up_P320_C" = "P320",
                                  "down_T5_C" = "T5", "down_T10_C" = "T10",
                                  "down_T20_C" = "T20", "down_T40_C" = "T40",
                                  "down_T80_C" = "T80", "down_T160_C" = "T160",
                                  "down_T320_C" = "T320", "down_P10_C" = "P10",
                                  "down_P320_C" = "P320"))

test$sample <- y
test$sample <- factor(test$sample, levels = c("P10", "P320", "T5", "T10", "T20", "T40", "T80",
                                              "T160", "T320"))

# subsetting for quick plot
tf_data_merged_up_select = tf_data_merged_up[tf_data_merged_up$tf %in% c("AP.1.bZIP",
                                                                         "Nrf2.bZIP",
                                                                         "Atf4.bZIP"), ]
# subsetting for quick plot
tf_data_merged_down_select = tf_data_merged_down[tf_data_merged_down$tf %in% c("Sox17.HMG",
                                                                               "Pax8.Paired.Homeobox",
                                                                               "WT1.Zf"), ]

test_t320 = test[test$sample == "T320", ]
ggplot(test_t320, aes(y = percent, x = reorder(tf, percent))) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge", colour="black",
           aes(fill = log10(pval))) +
  geom_text(aes(label = n_tf), nudge_y = 5) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  #facet_wrap(~tf, scales = "free_y", nrow = 3) +
  xlab("") +
  ylab("Peaks with motif (%)") +
  labs(fill = "-log(log(p-value))") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_colourbar(frame.colour = "black", ticks.colour = "black",
                                ticks.linewidth = 1, frame.linewidth = 1))

# This will be to show the proportion of TFs for all peaks once I run again
# should that include persisters? I don't think so

test_t320 = test[test$sample == "T320" & test$direction == "open", ]
test_t320 = test[test$sample == "T320" & test$direction == "closed", ]
test_t320 = test[test$sample == "T320", ]
test_t320$tf <- factor(test_t320$tf, levels = c("SOX17", "TEAD", "WT1", "NF1", "PAX8",
                                      "AP1", "ETS1", "RUNX", "ATF4", "NRF2"))

direction_colors = c("#4daf4a","#984fa3")
ggplot(test_t320, aes(y = percent, x = tf, fill = direction)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge", colour="black") +
  scale_fill_manual(values = direction_colors) +
  #geom_text(aes(label = n_tf), nudge_y = 5) +
  #scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  #coord_flip() +
  #facet_wrap(~tf, scales = "free_y", nrow = 3) +
  xlab("") +
  ylab("Peaks with motif (%)") +
  labs(fill = "direction") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 

# as barplot
#tf_data_merged_up_select
# exclude persisters
test_adapt = test[(!test$sample %in% c("P10", "P320")) & test$tf %in% c("AP1", "ATF4", "NRF2"), ]
test_adapt = test[(!test$sample %in% c("P10", "P320")) & test$tf %in% c("SOX17", "WT1", "PAX8"), ]
ggplot(test_adapt, aes(y = -percent, x = sample)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge", colour="black",
           aes(fill = log10(pval))) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  facet_wrap(~tf, scales = "free_y", nrow = 3) +
  
  xlab("") +
  ylab("Peaks with motif (%)") +
  labs(fill = "-log(log(p-value))") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.3, 'cm')) +
  guides(fill = guide_colourbar(frame.colour = "black", ticks.colour = "black",
                                ticks.linewidth = 1, frame.linewidth = 1))

############
# Barplots #
############

save(tf_data_merged_all, tf_data_merged_down, tf_data_merged_up, file = "tf_enrichment_downsample.rda")
load("tf_enrichment_downsample.rda")

tf_plots_up = list()
tfs = c("AP.1.bZIP", "Nrf2.bZIP", "Atf4.bZIP")
for(trans_fact in tfs) {
  temp_data = tf_data_merged_up[tf_data_merged_up$tf == trans_fact, ]
  tf_plot = ggplot(temp_data, aes(y = percent, x = sample)) +
    geom_bar(stat = "identity", width = 0.8, position = "dodge", colour="black",
             aes(fill = (pval))) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    xlab("") +
    ylab("Peaks with motif (%)") +
    labs(fill = "-log(p-value)") +
    ggtitle(trans_fact) +
    theme_classic() +
    theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.key.size = unit(0.3, 'cm')) +
    guides(fill = guide_colourbar(frame.colour = "black", ticks.colour = "black",
                                  ticks.linewidth = 1, frame.linewidth = 1))
  tf_plots_up[[trans_fact]] = tf_plot
}

library(gridExtra)
do.call("grid.arrange", c(tf_plots_up, ncol = 1, nrow = 3))


tf_plots_down = list()
tfs = c("Sox17.HMG","Pax8.Paired.Homeobox","WT1.Zf")
for(trans_fact in tfs) {
  temp_data = tf_data_merged_down[tf_data_merged_down$tf == trans_fact, ]
  tf_plot = ggplot(temp_data, aes(y = percent, x = sample)) +
    geom_bar(stat = "identity", width = 0.8, position = "dodge", colour="black",
             aes(fill = (pval))) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    xlab("") +
    ylab("Peaks with motif (%)") +
    labs(fill = "-log(p-value)") +
    ggtitle(trans_fact) +
    theme_classic() +
    theme(axis.text = element_text(color="black"), text = element_text(size=10)) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.key.size = unit(0.3, 'cm')) +
    guides(fill = guide_colourbar(frame.colour = "black", ticks.colour = "black",
                                  ticks.linewidth = 1, frame.linewidth = 1))
  tf_plots_down[[trans_fact]] = tf_plot
}

library(gridExtra)
do.call("grid.arrange", c(tf_plots_down, ncol = 1, nrow = 3))


#########################################
# bubble plot for N of peaks in each TF #
#########################################

tf_plots_up_n = list()
tfs = c("AP.1.bZIP", "Nrf2.bZIP", "Atf4.bZIP")
for(trans_fact in tfs) {
  temp_data = tf_data_merged_up[tf_data_merged_up$tf == trans_fact, ]
  tf_plot = ggplot(temp_data, aes(x = sample, y = direction, size = n_tf)) +
    geom_point(aes(size = n_tf), shape = 21) +
    #scale_size_continuous(range = c(1, 5)) +
    scale_size(range = c(2,12),
               #labels=c("0","60","120","180","240",">=300"),
               guide="legend") +
    ylab("") +
    xlab("") +
    ggtitle(trans_fact) +
    theme_minimal() 
  tf_plots_up_n[[trans_fact]] = tf_plot
}

library(gridExtra)
do.call("grid.arrange", c(tf_plots_up_n, ncol = 1, nrow = 3))


tf_plots_down_n = list()
tfs = c("Sox17.HMG","Pax8.Paired.Homeobox","WT1.Zf")
for(trans_fact in tfs) {
  temp_data = tf_data_merged_down[tf_data_merged_down$tf == trans_fact, ]
  tf_plot = ggplot(temp_data, aes(x = sample, y = direction, size = n_tf)) +
    geom_point(aes(size = n_tf), shape = 21) +
    #scale_size_continuous(range = c(1, 5)) +
    scale_size(range = c(2,12),
               #labels=c("0","60","120","180","240",">=300"),
               guide="legend") +
    ylab("") +
    xlab("") +
    ggtitle(trans_fact) +
    theme_minimal() 
  tf_plots_down_n[[trans_fact]] = tf_plot
}

do.call("grid.arrange", c(tf_plots_down_n, ncol = 1, nrow = 3))


###############
# Bubble plot #
###############

# as bubble plot
ggplot(test[test$direction == "open", ], aes(x = sample, y = tf, size = -log10(orig_pval), fill = percent)) + 
  geom_point(shape = 21) +
  #scale_colour_gradient(low = "gray", high = "blue", name = "-log10(q-val)") + 
  #scale_color_viridis_c(name = "% of peaks") + 
  scale_fill_distiller(palette = "RdBu") +
  #scale_size(range = c(2,12), breaks = seq(0, max(test_plot$percent), by = 10)) +
  scale_size(range = c(1,12), breaks=c(0,60,120,180,240,300),
             labels=c("0","60","120","180","240",">=300"),
             guide="legend") +
  facet_wrap(~direction, scales = "free_y", ncol = 1, nrow = 2) +
#  coord_fixed(ratio=1, expand = T) +
  ylab("") +
  xlab("") +
  labs(fill = "peaks with motif (%)", size = "-log10(p-value)") +
  theme_minimal() +
  theme(axis.text = element_text(color="black"), text = element_text(size=10),
        #axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.4,
                                frame.colour = "black", ticks.colour = "black",
                                ticks.linewidth = 1, frame.linewidth = 1,
                                barheight = 0.8),
         size = guide_legend(title.position="top", title.hjust = 0.5,
                             label.position = "bottom"))


############################
# Peak numbers and overlap #
############################

setwd("/home/gu/Dropbox/projects/kuramochi_olaparib_evolution/atac-seq/deseq2_peaks")


filenames <- list.files(pattern='.bed', recursive=TRUE, full.names = TRUE)

bed_data<- lapply(filenames, read.csv, head = F, sep = "\t")

# getting the names of the samples
col_names = stringr::str_split(filenames, pattern = "/", simplify = T)[, 2]
col_names = stringr::str_split(col_names, pattern = "_", simplify = T)[, c(3,1)]

# adding names to the list
names_bed_data = paste(col_names[,1], col_names[, 2], sep = "_")
names(bed_data) <- names_bed_data

# create a table with N of peaks in each file
n_peaks = (sapply(bed_data, nrow))
n_peaks = data.frame(sample = col_names[,1], group = col_names[,2], n_peaks = n_peaks)

n_peaks$sample <- factor(n_peaks$sample, levels = c("T5", "T10", "T20", "T40",
                                                    "T80", "T160", "T320"))

ggplot(n_peaks, aes(y = n_peaks, x = sample, fill = group)) +
  geom_bar(stat = "identity", width = 0.95, position = "dodge") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#6baed6", "#a54d69"), name = "") +
  scale_y_continuous(limits = c(0, max(n_peaks$n_peaks) + 100), 
                     breaks = seq(0, max(n_peaks$n_peaks) + 100, by = 1000)) +
  xlab("") +
  ylab("Number of peaks") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8)) +
  theme(strip.background =element_rect(fill="white"))

####################################
# Testing for overlap across peaks #
####################################

# get peak names only
peaks_data = lapply(bed_data, function(x) as.character(x[, 4]))

peaks_up = peaks_data[c("T5_up", "T10_up", "T20_up", "T40_up", "T80_up", "T160_up", "T320_up")]
peaks_down = peaks_data[c("T5_down", "T10_down", "T20_down", "T40_down", "T80_down", "T160_down", "T320_down")]

# for pairwise jaccard indices
peaks_up_overlap = sapply(peaks_up, 
                         function(x) sapply(peaks_up,
                                            function(y)bayesbio::jaccardSets(x,y)))
# for pairwise jaccard indices
peaks_down_overlap = sapply(peaks_down, 
                          function(x) sapply(peaks_down,
                                             function(y)bayesbio::jaccardSets(x,y)))


peaks_dist = as.dist(1 - cor(peaks_down_overlap))
cormM = cor(peaks_down_overlap, method="pearson")
dend = as.dendrogram(hclust(peaks_dist, method = "average"))
dend <- dendextend::click_rotate(dend, continue = TRUE)


pdf("")
h2 = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
             show_column_dend = T, show_row_names = T, 
             name = "", row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
             cluster_rows = dend, cluster_columns = dend,
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "Correlation",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))

draw(h2)


#############################
# Plot PCA for or all peaks #
#############################

vsd <- vst(dds)
all_up_peaks = sort(unique(c(unlist(peaks_up))))
all_down_peaks = sort(unique(c(unlist(peaks_down))))
all_diff_peaks = sort(unique(c(unlist(peaks_down), unlist(peaks_up))))

pca_plot = plotPCA(vsd[all_up_peaks, ], "condition")
pca_plot = plotPCA(vsd[all_down_peaks, ], "condition")
pca_plot = plotPCA(vsd[all_diff_peaks, ], "condition", ntop = 3000)
pca_plot = plotPCA(vsd, "condition", ntop = 2000)

# reorder groups
pca_plot$data$group <- factor(pca_plot$data$group, levels = c("C", "P10", "P320",
                                                              "T5", "T10", "T20",
                                                              "T40", "T80", "T160", "T320"))
new_colors = c("#8DA3A6", "#b74980", "#cc8801", "#4E78C4","#57A2AC","#7EB875",
               "#D0B541", "#E67F33", "#CE2220", "#80271E")

ggplot(data = pca_plot$data, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = group), size = 2) +
    scale_colour_manual(values = new_colors) +
    xlab(pca_plot$labels$x) + ylab(pca_plot$labels$y) +
    theme_bw() +
    theme(axis.text = element_text(color="black"), text = element_text(size=10),
          legend.position="bottom", legend.title= element_blank())


##########################################################################
# Correlation between N of diff expressed peaks and total peaks detected #
# to rule out the claim that more peaks are diff. exp because more peaks 
# detected
##########################################################################

n_peaks_up = sapply(all_up_peaks, length)
n_peaks_down = sapply(all_down_peaks, length)

# information from server, files of consensus peaks for both reps
total_peaks_detected = c(#27024, control
                         45280, 41319, 28624, 20744, 31826, 38373, 36247,
                         79198, 81795)
total_peaks_diff = n_peaks_up + n_peaks_down
cor.test(total_peaks_detected, total_peaks_diff)

data_n_peaks_corr = data.frame(sample = c("T5", "T10", "T20", "T40", "T80",
                                          "T160", "T320", "P10", "P320"),
                               total_peaks = total_peaks_detected,
                               total_diff_peaks = total_peaks_diff)

data_n_peaks_corr$sample <- factor(data_n_peaks_corr$sample,
                                   levels = c("P10", "P320", 
                                              "T5", "T10", "T20", "T40", "T80",
                                              "T160", "T320"))
data_n_peaks_corr$total_peaks <- as.numeric(data_n_peaks_corr$total_peaks)
data_n_peaks_corr$total_diff_peaks <- as.numeric(data_n_peaks_corr$total_diff_peaks)


library(ggplot2)

colors = c("#b74980", "#cc8801", "#4E78C4","#57A2AC","#7EB875",
           "#D0B541", "#E67F33", "#CE2220", "#80271E")

ggplot(data_n_peaks_corr, aes(x = total_diff_peaks, y = total_peaks, color = sample)) +
  geom_point() +
  
  geom_smooth(method = "lm", se = F) +
  ylab("Detected peaks") +
  xlab("Differentially accessible peaks") +
  #stat_smooth(method = "lm",
  #            formula = y ~ x,
  #            geom = "smooth") +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks=seq(0, 15000, 2000)) +
  scale_y_continuous(breaks=seq(15000, 85000, 10000)) +
  theme_bw() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))
  #scale_color_manual(values = colors)


ggplot(data_n_peaks_corr, aes(x = total_diff_peaks, y = total_peaks) ) +
  geom_point(aes(color = sample)) +
  geom_smooth(method = "lm", se = FALSE, color = "#bbbbbb", size = 0.5) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks=seq(0, 15000, 2000)) +
  scale_y_continuous(breaks=seq(15000, 85000, 10000)) +
  ylab("Detected peaks") +
  xlab("Differentially accessible peaks") +
  theme_bw() +
  theme(axis.text = element_text(color="black"), text = element_text(size=8))

