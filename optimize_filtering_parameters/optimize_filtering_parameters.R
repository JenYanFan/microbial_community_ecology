#!/usr/bin/Rscript
packages <- c("ShortRead", "tidyverse")
install.packages(setdiff(packages, rownames(installed.packages())))
invisible(sapply(packages, function(x) require(x, character.only = TRUE, quietly = TRUE)))
# Get post-filtering stats
get_post_filtering_stat <- function (errProbF, errProbR, inspected_maxEE, required_len, read_len) {
  min_len <- required_len - read_len
  # Inspected parameter set
  filt_stat <- expand.grid(maxEEf = inspected_maxEE,
                             lenF = min_len:read_len,
                             maxEEr = inspected_maxEE,
                             lenR = min_len:read_len)
  filt_stat <- filt_stat[(filt_stat$lenF + filt_stat$lenR) == required_len, ]
  # Expected errors for reads with various lengths
  eeF <- lapply(lapply(min_len:ncol(errProbF), sequence),
                       function (x) { return(rowSums(errProbF[, x])) })
  eeR <- lapply(lapply(min_len:ncol(errProbR), sequence),
                       function (x) { return(rowSums(errProbR[, x])) })
  # Clear unuse variables
  errProbF <- NULL
  errProbR <- NULL
  # Count reads with expected errors <= MaxEE
  filt_stat <- 
    Reduce(lapply(split(filt_stat, row.names(filt_stat)),
           function (x) {
             x$pass_reads <- sum((eeF[[x$lenF + 1 - min_len]] <= x$maxEEf) * (eeR[[x$lenR + 1 - min_len]] <= x$maxEEr), na.rm = TRUE)
             x$total_reads <- length(eeF[[1]])
             return(x)
             }
           ), f = rbind)
  return(filt_stat)
}

# Integrate post-filtering stats from multiple samples
get_post_filtering_stat_for_multiple_files <- function(fileF, fileR, inspected_maxEE, required_len, read_len) {
  errProbF <- 10^-(as(quality(readFastq(fileF[1])), "matrix")/10)
  errProbR <- 10^-(as(quality(readFastq(fileR[1])), "matrix")/10)
  filt_stat <- get_post_filtering_stat(errProbF, errProbR, inspected_maxEE, required_len, read_len)
  for(i in 2:length(fileF)) {
    errProbF <- 10^-(as(quality(readFastq(fileF[i])), "matrix")/10)
    errProbR <- 10^-(as(quality(readFastq(fileR[i])), "matrix")/10)
    tmp <- get_post_filtering_stat(errProbF, errProbR, inspected_maxEE, required_len, read_len)
    filt_stat[, c("pass_reads", "total_reads")] <- filt_stat[, c("pass_reads", "total_reads")] + tmp[, c("pass_reads", "total_reads")]
  }
  filt_stat$pass_ratio <- filt_stat$pass_reads/filt_stat$total_reads
  filt_stat$EE_index <- (filt_stat$maxEEf^2 + filt_stat$maxEEr^2)^0.5
  return(filt_stat)
}

# Extract optimal length cutoff and parameter candidate
prioritize_parameters <- function(filt_stat) {
  optimal_len <- by(filt_stat, filt_stat[, c("maxEEr", "maxEEf")], function(x) rownames(x)[which.max(x$pass_reads)])
  filt_stat <- filt_stat[unlist(optimal_len), ]
  candidate <- filt_stat[order(filt_stat$EE_index, -filt_stat$pass_ratio), ]
  for (i in 2:nrow(candidate)) {
    if (candidate[i, "pass_ratio"] <= max(candidate[1:i-1, "pass_ratio"])) {
      candidate[i, "EE_index"] <- NA
    }
  }
  candidate <- candidate[!is.na(candidate$EE_index), ]
  fo_diff <- sapply(1:nrow(candidate),
                    function(pos){
                      left <- (candidate[pos, "pass_ratio"] - candidate[1, "pass_ratio"]) / (candidate[pos, "EE_index"])
                      right <- (candidate[nrow(candidate), "pass_ratio"] - candidate[pos, "pass_ratio"]) / (candidate[nrow(candidate) , "EE_index"]- candidate[pos, "EE_index"])
                      return(left - right)
                    }
  )
  candidate$fo_diff <- fo_diff
  candidate <- candidate[order(-candidate$fo_diff), ]
  return(candidate)
}







