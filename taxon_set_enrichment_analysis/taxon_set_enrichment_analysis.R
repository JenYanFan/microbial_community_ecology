#!/usr/bin/Rscript
source("./indicator_value_analysis/calculate_IndVal.R")
library(clusterProfiler)

conduct_TSEA <- function(feature_count, metadata, category, subtype, taxonomy, rank, alpha = 0.05, random_seed = 19941119, p_adjust_method = "fdr") {
  IndVal <- calculate_indival(feature_count, metadata, "country")
  IndVal$local_statistics <- ifelse(IndVal$indicated_subtype == subtype, IndVal$indval, IndVal$indval * -1)
  
  taxon_list <- IndVal$local_statistics
  names(taxon_list) <- IndVal$features
  taxon_list <- sort(taxon_list, decreasing = TRUE)
  
  TERM2GENE <- data.frame(term = taxonomy[, rank][match(names(taxon_list), row.names(taxonomy))],
                          taxon = names(taxon_list))
  TSEA_stat <- GSEA(taxon_list, TERM2GENE = TERM2GENE, pvalueCutoff =  alpha, seed = random_seed, pAdjustMethod = p_adjust_method )
  return(TSEA_stat)
}

