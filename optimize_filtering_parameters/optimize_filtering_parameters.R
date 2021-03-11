#!/usr/bin/Rscript
packages <- c("ShortRead", "tidyverse")
install.packages(setdiff(packages, rownames(installed.packages())))
invisible(sapply(packages, function(x) require(x, character.only = TRUE, quietly = TRUE)))
# Get post-filtering stats
get_post_filtering_stat <- function (errProbF, errProbR, inspected_maxEE, required_len, read_len) {
  min_len <- required_len - read_len
  filter_stat <- expand.grid(maxEEf = inspected_maxEE, 
                             lenF = min_len:read_len, 
                             maxEEr = inspected_maxEE, 
                             lenR = min_len:read_len)
  filter_stat <- filter_stat[(filter_stat$lenF + filter_stat$lenR) == required_len, ]
  e <- new.env()
  e$prm <- list()
  filter_stat <- lapply(split(filter_stat, row.names(filter_stat)),
                        function (x) {
                          maxEEf <- as.character(x$maxEEf)
                          maxEEr <- as.character(x$maxEEr)
                          lenF <- as.character(x$lenF)
                          lenR <- as.character(x$lenR)
                          # Inspect hash table for calculated values
                          if(is.null(e$prm[["F"]][[maxEEf]][[lenF]])){
                            e$prm[["F"]][[maxEEf]][[lenF]] <- rowSums(errProbF[, 1:x$lenF]) <= x$maxEEf}
                          if(is.null(e$prm[["R"]][[maxEEr]][[lenR]])){
                            e$prm[["R"]][[maxEEr]][[lenR]] <- rowSums(errProbR[, 1:x$lenR]) <= x$maxEEr}
                          # Calculate reads with EE below maxEEf and below maxEEr
                          x$pass_reads <- sum(e$prm[["R"]][[maxEEr]][[lenR]] * e$prm[["F"]][[maxEEf]][[lenF]], na.rm = TRUE)
                          x$total_reads <- nrow(errProbF)
                          return(x)
                        }
  )
  filter_stat <- Reduce(filter_stat, f = rbind.data.frame)
  return(filter_stat)
}
# Integrate post-filtering stats from multiple samples
get_post_filtering_stat_for_multiple_files <- function(fileF, fileR, inspected_maxEE, required_len, read_len) {
  errProbF <- 10^-(as(quality(readFastq(fileF[1])), "matrix")/10)
  errProbR <- 10^-(as(quality(readFastq(fileR[1])), "matrix")/10)
  filter_stat <- get_post_filtering_stat(errProbF, errProbR, inspected_maxEE, required_len, read_len)
  for(i in 2:length(fileF)) {
    errProbF <- 10^-(as(quality(readFastq(fileF[i])), "matrix")/10)
    errProbR <- 10^-(as(quality(readFastq(fileR[i])), "matrix")/10)
    tmp <- get_post_filtering_stat(errProbF, errProbR, inspected_maxEE, required_len, read_len)
    filter_stat[, c("pass_reads", "total_reads")] <- filter_stat[, c("pass_reads", "total_reads")] + tmp[, c("pass_reads", "total_reads")]
  }
  filter_stat$pass_ratio <- filter_stat$pass_reads/filter_stat$total_reads
  filter_stat$EE_index <- (filter_stat$maxEEf^2 + filter_stat$maxEEr^2)^0.5
  return(filter_stat)
}
# Extract optimal length cutoff and parameter candidate
prioritize_parameters <- function(filter_stat) {
  optimal_len <- by(filter_stat, filter_stat[, c("maxEEr", "maxEEf")], function(x) rownames(x)[which.max(x$pass_reads)])
  filter_stat <- filter_stat[unlist(optimal_len), ]
  candidate <- filter_stat[order(filter_stat$EE_index, -filter_stat$pass_ratio), ]
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







