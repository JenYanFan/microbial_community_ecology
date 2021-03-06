#!/usr/bin/Rscript
source("./optimize_filtering_parameters/optimize_filtering_parameters.R")
# Parameters
required_len <- 253 + 20
read_len <- 175
min_len <- required_len - read_len
inspected_maxEE <- seq(1, 10, 2)
# Read files
path="./data"
fqF <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fqR <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Get post-filtering stats
filter_stat <- get_post_filtering_stat_for_multiple_files(fqF,fqR,inspected_maxEE, required_len, read_len)
# Extract optimal length cutoff and parameter candidate
candidate  <- prioritize_parameters(filter_stat)
# Visualize results
plot_data <- pivot_longer(candidate, 
                          cols = c("fo_diff", "pass_ratio"),
                          names_to = "names",
                          values_to = "values"
)
plot_data$names <- 
  factor(plot_data$names, levels = c("pass_ratio", "fo_diff"),
         labels = c("Remained read ratio", "First-order difference"))
p <- ggplot(data = plot_data,
            aes(x = EE_index, y = values)) +
  geom_vline(aes(xintercept = candidate[which.max(candidate$fo_diff), "EE_index"]),
             linetype = "dashed", color = "red", size = 1) +
  labs(x = "Index of expected error threshold for forward and reverse reads") + 
  geom_line(size = 1) +
  facet_wrap(~names, nrow = 2, scales = "free_y", strip.position = "left", shrink = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", color = "black"),
        strip.placement = "outside")

ggsave(p, filename = "./j.jpg",
       width = 8, height = 8, dpi = 300)

# --------------

filter_stat <- expand.grid(maxEEf = 1:10, 
                           lenF = 98:175, 
                           maxEEr = 1:10, 
                           lenR = 98:175) 
mtx_filt_stat <- as.matrix(filter_stat)
