#!/usr/bin/Rscript
source("./indicator_value_analysis/calculate_IndVal.R")
library(clusterProfiler)
library(enrichplot)
# read files
otu_count <- read.table("./data/2020-12-31_P1-E1_diabimmune_otu-count.txt")
metadata <- read.table("./data/2020-12-31_P1-E1_diabimmune_metadata.txt")
taxonomy <- read.table("./data/2020-12-31_P1-E1_diabimmune_taxonomy.txt")
rarified_otu_count <- t(rrarefy(t(otu_count), sample = min(colSums(otu_count))))
rarified_otu_count <- rarified_otu_count[rowSums(rarified_otu_count) > 0, ]

IndVal <- calculate_indival(rarified_otu_count, metadata, "country")
IndVal$local_statistics <- ifelse(IndVal$indicated_subtype == "FIN", IndVal$indval, IndVal$indval * -1)

taxon_list <- IndVal$local_statistics
names(taxon_list) <- IndVal$features
taxon_list <- sort(taxon_list, decreasing = TRUE)

TERM2GENE <- data.frame(term = taxonomy$Class[match(names(taxon_list), row.names(taxonomy))],
                        taxon = names(taxon_list))
a <- GSEA(taxon_list, TERM2GENE = TERM2GENE, pvalueCutoff =  0.05, seed = 19941119, pAdjustMethod = "fdr" )
gseaplot2(a, geneSetID = 1:2)

dotplot(a)
barplot(a)
ridgeplot(a)
heatplot(a)
