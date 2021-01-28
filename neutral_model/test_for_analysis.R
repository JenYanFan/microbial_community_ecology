#!/usr/bin/Rscript
source("./neutral_model/create_null_model.R")

# read files
otu_count <- read.table("./data/2020-12-31_P1-E1_diabimmune_otu-count.txt")
metadata <- read.table("./data/2020-12-31_P1-E1_diabimmune_metadata.txt")
rarified_otu_count <- t(rrarefy(t(otu_count), sample = min(colSums(otu_count))))
rarified_otu_count <- rarified_otu_count[rowSums(rarified_otu_count) > 0, ]

# neutral distance ----------
set.seed(19941119)
neutral_distance <- calculate_neutral_distance(rarified_otu_count, replication = 20)
neutral_distance <- lapply(neutral_distance, as.matrix)
observed_distance <- as.matrix(vegdist(t(rarified_otu_count), method = "bray"))
raup_crick <- calculate_raup_crick(observed_distance, neutral_distance)
deterministic_ratio <- calculate_deterministic_ratio(observed_distance, neutral_distance)

raup_crick_dist <- distance_to_dataframe(raup_crick)
DR_dist <- distance_to_dataframe(deterministic_ratio)

jotaro <- merge(raup_crick_dist , metadata, by.x = c("sample_1"), by.y = c("sampleid"))

ggplot(data = jotaro, aes(x = gender, y = distance, group = gender)) +
  geom_jitter()


