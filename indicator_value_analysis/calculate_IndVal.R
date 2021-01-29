#!/usr/bin/Rscript
library(labdsv)
library(tidyverse)
library(vegan)

calculate_indival <- function (feature_count, metadata, category, numitr = 10, p_adjust_method = "fdr", random_seed = 19941119) {
  # feature_count: a matrix or data.frame with samples as columns and features as rows
  # metadata: a data.frame with categories as columns and samples as rows
  # category: categories that are included
  # numitr: the number of randomizations to iterate to calculate probabilities
  # p_adjust_method: method for adjustin p-values in multiple comparisons in the p.adjust function within the stats package
  # random_seed
  transposed_feature_count <- t(feature_count)
  group <- as.character(metadata[match(row.names(transposed_feature_count), row.names(metadata)), category])
  IndVal <- indval(transposed_feature_count, group, numitr)
  
  indval_stat <- data.frame(features = row.names(feature_count))
  indval_stat$indval <- IndVal$indcls[match(names(IndVal$indcls), indval_stat$features)]
  indval_stat$indicated_subtype <- apply(IndVal$indval, 1, function(x) colnames(IndVal$indval)[which(x == max(x))])
  indval_stat$p_value <- IndVal$pval[match(names(IndVal$pval), indval_stat$features)]
  adjusted_p_value <- p.adjust(IndVal$pval, method = p_adjust_method)
  indval_stat$adjusted_p_value <- adjusted_p_value[match(names(adjusted_p_value), indval_stat$features)]
  
  return(indval_stat)
}

