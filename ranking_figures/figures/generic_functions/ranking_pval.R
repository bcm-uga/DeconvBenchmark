library(foreach)

# Permutation test
agg_diff <- function(candidate_i, candidate_j, scores_Z, agg_process,
                     scores_to_keep = scores_to_keep) {
  source("../generic_functions/ranking_process.R")

  final_score <- agg_process(scores = scores_Z, scores_to_keep = scores_to_keep)

  final_score_i <- final_score$overall[final_score$candidate == candidate_i]
  final_score_j <- final_score$overall[final_score$candidate == candidate_j]
  return(as.numeric(final_score_i - final_score_j))
}

paired_permutation_test <- function(scores, agg_process,
                                    #top_n=5,
                                    n_permutations = 1000, scores_to_keep = NULL,
                                    candidates = NULL) {
  if (is.null(candidates)) {
    candidates <- agg_process(scores = scores, scores_to_keep = scores_to_keep) %>%
      #top_n(top_n, overall) %>%
      arrange(desc(overall)) %>%
      pull(candidate)
  }

  scores <- scores %>%
        arrange(candidate, dataset, name_score)

  n_pairs <- sum(1:(length(candidates) - 1))
  done <- 1
  pairs <- rep("", n_pairs)
  pvalues <- rep(NA, n_pairs)
  distributions <- list()
  for (i in 1:(length(candidates) - 1)) {
    print(paste(done, n_pairs, sep = " / "))
    candidate_i <- candidates[i]

    for (j in (i + 1):length(candidates)) {
      candidate_j <- candidates[j]
      pairs[done] <- paste(candidate_i, candidate_j, sep = ">")

      scores_Z <- scores %>%
        filter(candidate %in% c(candidate_i, candidate_j))
      n <- nrow(scores_Z) / 2

      s0 <- agg_diff(candidate_i, candidate_j, scores_Z, agg_process, scores_to_keep = scores_to_keep)
      s_k = pbapply::pbsapply(seq(n_permutations), function(p) {
        exchange <- sample(c(T,F), n, replace = TRUE)
        idx_X <- sapply(seq_along(exchange), function(idx)
          ifelse(exchange[idx],idx+n,idx))
        idx_Y <- sapply(idx_X, function(x) ifelse(x>n,x-n,x+n))
        scores_XY <- scores_Z
        scores_XY$trendval <- scores_Z$trendval[c(idx_X, idx_Y)]
        agg_diff(candidate_i, candidate_j, scores_XY, agg_process, scores_to_keep = scores_to_keep)
      })
      s_k <- c(s0, s_k)
      pvalues[done] <- mean(s_k >= s0)
      distributions[[pairs[done]]] <- s_k
      done <- done + 1
    }
  }
  names(pvalues) <- pairs
  return(list(pvalues = pvalues, distributions = distributions))
}