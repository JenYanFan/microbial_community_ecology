#!/usr/bin/Rscript
library(vegan)
library(tidyverse)
library(parallel)

calculate_betaNTI <- function(feature_count, phylogenetic_tree, replication = 10, cores = 1L) {
  observed_betaMNTD <- as.matrix(comdistnt(t(feature_count), cophenetic(phylogenetic_tree), abundance.weighted = TRUE))
  neutral_betaMNTD <- mclapply(1:replication,
                               mc.cores = getOption("mc.cores", cores),
                               function(i) {
                                 neutral_betaMNTD <-  as.matrix(comdistnt(t(feature_count),taxaShuffle(cophenetic(phylogenetic_tree)), abundance.weighted = TRUE, exclude.conspecifics = FALSE))
                                 return(neutral_betaMNTD)
                               }
  )
  betaNTI <- (observed_betaMNTD - Reduce(neutral_betaMNTD, f = "+")/length(neutral_betaMNTD))/apply(simplify2array(neutral_betaMNTD), 1:2, sd)
  return(betaNTI)
}
calculate_deterministic_ratio <- function(observed_distance, neutral_distance) {
  neutral_distance <- as.matrix(neutral_distance)
  mean_neutral_distance <- Reduce(neutral_distance, f = "+")/length(neutral_distance)
  deterministic_ratio <- (observed_distance - mean_neutral_distance)/observed_distance
  return(deterministic_ratio)
}
calculate_extreme_similarity <- function(observed_similarity, mean_neutral_similarity) {ifelse(observed_similarity >= mean_neutral_similarity, 1, 0)}
calculate_neutral_distance <- function (feature_count, replication = 10, cores = 1L, method = "bray") {
  neutral_distance <- 
    mclapply(1:replication,
             mc.cores = getOption("mc.cores", cores),
             function(i) {
               neutral_count <- create_neutral_count(feature_count)
               neutral_distance <- vegdist(t(neutral_count), method = method)
               return(neutral_distance)
             }
    )
  return(neutral_distance)
}
calculate_normalized_selection_strength <- function(observed_distance, neutral_distance) {
  observed_similarity <- 1 - observed_distance
  neutral_similarity <- lapply(neutral_distance, function(x) {1 - x})
  mean_neutral_similarity <- Reduce(neutral_similarity, f = "+")/length(neutral_similarity)
  extreme_similarity  <- calculate_extreme_similarity(observed_similarity, mean_neutral_similarity)
  observed_selection_strength <- sum(calculate_selection_strength(observed_similarity, mean_neutral_similarity))
  minimal_selection_strength <- min(sapply(neutral_similarity, function(x) {sum(calculate_selection_strength(x, mean_neutral_similarity))}))
  maximal_selection_strength <- sum(calculate_selection_strength(extreme_similarity, mean_neutral_similarity))
  normalized_selection_strength <- (observed_selection_strength - minimal_selection_strength)/(maximal_selection_strength - minimal_selection_strength)
  return(normalized_selection_strength)
}
calculate_raup_crick <- function(observed_distance, neutral_distance, split_ties = TRUE, classic_matric = FALSE) {
  neutral_distance <- as.matrix(neutral_distance)
  rc <- Reduce(lapply(neutral_distance, 
                      function(x) {
                        ifelse(split_ties == TRUE, 
                               return((observed_distance > x) + (observed_distance == x)/2),
                               return(observed_distance == x))}), f = "+")/length(neutral_distance)
  ifelse(classic_matric == TRUE, return(rc), return((rc - 0.5)*2))
}
calculate_selection_strength <- function(observed_similarity, mean_neutral_similarity) {
  x <- observed_similarity
  y <- mean_neutral_similarity
  return(ifelse(x >= y, (x - y)/x, (x - y)/(x - 1)))
}
create_neutral_community <- function (count_distribution, incidence_distribution, richness, community_size) {
  features <- names(count_distribution)
  selected_features <- sample(features, size = richness, replace = FALSE, prob = incidence_distribution)
  selected_count_distribution <- count_distribution[names(count_distribution) %in% selected_features]
  selected_count_distribution <- selected_count_distribution[match(selected_features, names(selected_count_distribution))]
  selected_features <- c(sample(selected_features, size = community_size - length(selected_features), prob = selected_count_distribution, replace = TRUE), selected_features)
  neutral_community <- table(selected_features)
  return(neutral_community)
}
create_neutral_count <- function (feature_count) {
  feature_distribution <- get_feature_distribution(feature_count)
  neutral_count <- 
    apply(feature_count,
          MARGIN = 2, 
          function(x) {
            neutral_community <- create_neutral_community(feature_distribution[["count"]], feature_distribution[["incidence"]], sum(x > 0), sum(x))
            x <- x*0
            x[match(names(neutral_community), names(x))] <- neutral_community
            return(x)
            }
          )
  return(neutral_count)
}
distance_to_dataframe <- function(distance) {
  distance <- as.matrix(distance, labels = TRUE)
  upper_tri <- upper.tri(distance)
  distance_dataframe <- data.frame(
    sample_1 = rownames(distance)[col(distance)[upper_tri]],
    sample_2 = rownames(distance)[row(distance)[upper_tri]],
    distance = distance[upper_tri]
  )
  return(distance_dataframe)
}
get_feature_distribution <- function(feature_count) {
  feature_incidence <- ifelse(feature_count > 0, 1, 0)
  feature_distribution <- list()
  feature_distribution[["count"]] <- rowSums(feature_count)
  feature_distribution[["incidence"]] <- rowSums(feature_incidence)
  return(feature_distribution)
}



