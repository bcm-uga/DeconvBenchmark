library(tidyr)

# reorder score matrix to get a candidate x judge matrix instead of a data.frame
make_pavao_matrix <- function(df) {
  df_wide <- df %>%
    pivot_wider(id_cols = candidate, names_from = judge, values_from = val)
  mat <- as.matrix(df_wide %>% select(-candidate))
  rownames(mat) <- df_wide %>% pull(candidate)
  return(mat)
}

# theoretical criteria
maj_criterion <- function(matrix,winner) {
  get_best_candidate_per_judge <- apply(matrix, 2, function(x) names(which.max(x)))
  most_voted <- sort(table(get_best_candidate_per_judge), decreasing = TRUE)[1]
  get_best_candidate <- ifelse(most_voted>ncol(matrix)/2,names(most_voted),NA)
  if (!is.na(get_best_candidate)) {
    return(get_best_candidate==winner)
  }
  else {
    return(NA)
  }
}

condorcet_criterion <- function(matrix,winner) {
  winner_idx <- which(rownames(matrix)==winner)
  is_cond_winner <- rep(TRUE, nrow(matrix))
  for (i in seq(nrow(matrix) - 1)) {
    for (j in (i+1):(nrow(matrix))) {
      sub_matrix <- matrix[c(i,j),]
      get_best_candidate_per_judge <- unlist(apply(sub_matrix, 2, function(z) names(which.max(z))))
      n_voted <- table(get_best_candidate_per_judge)
      most_voted <- max(n_voted)
      if (most_voted > ncol(sub_matrix) / 2) {
        pair_loser_idx <- ifelse(i == which.max(n_voted), j, i)
        is_cond_winner[pair_loser_idx] <- FALSE
      } else {
        is_cond_winner[i] <- FALSE
        is_cond_winner[j] <- FALSE
      }
    }
  }
  cond_winner_idx <- is_cond_winner[is_cond_winner == TRUE]
  return(ifelse(!is.null(cond_winner_idx), cond_winner_idx == winner_idx, NA))
}

consistency_criterion <- function(score_df, winner, f_M, n_parts=3, n_bootstrap=100) {
  judge_df <- score_df %>%
    select(c(dataset, name_score)) %>%
    unique()
  n_judge <- nrow(judge_df)
  winner_bootstrap <- rep(NA, n_bootstrap)
  for (k in seq(n_bootstrap)) {
    parts <- split(sample(seq(n_judge)), sort(seq(n_judge)%%n_parts))
    winner_part <- rep(NA, n_parts)
    # determine parts
    for (i in seq(n_parts)) {
      part <- parts[[i]]
      sub_score_df <- score_df %>%
        right_join(judge_df[part,])
      final_score <- f_M(scores = sub_score_df)
      winner_part[i] <- final_score %>%
        filter(!(grepl("CIBERSORT4", candidate) | grepl("svr4", candidate))) %>%
        arrange(desc(values)) %>%
        top_n(n=1, wt=values) %>%
        pull(candidate)
    }
    winner_part <- NA
    if (length(unique(winner_part)) == 1) {
      winner_part <- unique(winner_part)
    }
    if (!is.na(winner_part) & winner_part != winner) return(FALSE)
  }
  winner_bootstrap <- winner_bootstrap[!is.na(winner_bootstrap)]
  if (length(unique(winner_bootstrap)) == 1) {
    score <- unique(winner_bootstrap) == winner
  } else {
    score <- NA
  }
  return(score)
}

participation_criterion <- function(score_df, final_score, f_M, n_bootstrap=100) {
  judge_df <- score_df %>%
    select(c(dataset, name_score)) %>%
    unique()
  n_judge <- nrow(judge_df)
  candidates <- final_score %>%
    pull(candidate)
  for (k in seq(n_bootstrap)) {
    candidates_boot <- sample(candidates, 2)
    sub_final_score <- final_score %>%
      right_join(data.frame(candidate = candidates_boot))
    candidates_boot <- sub_final_score %>%
      arrange(desc(values)) %>%
      pull(candidate)
    diff_candidates <- sub_final_score[sub_final_score$candidate == candidates_boot[1], "values"] -
      sub_final_score[sub_final_score$candidate == candidates_boot[2], "values"]
    sub_score_df <- score_df %>%
        right_join(judge_df[-sample(n_judge, 1), ])
    final_score_boot <- f_M(scores = sub_score_df)
    diff_candidates_boot <- final_score_boot[final_score_boot$candidate == candidates_boot[1], "values"] -
      final_score_boot[final_score_boot$candidate == candidates_boot[2], "values"]
    cond <- diff_candidates <= diff_candidates_boot
    if (!is.na(cond) & cond) return(FALSE)
  }
  return(TRUE)
}

iia_criterion <- function(score_df, final_score, f_M, n_bootstrap=100) {
  candidates <- final_score %>%
    pull(candidate)
  for (k in seq(n_bootstrap)) {
    candidates_boot <- sample(candidates, 2)
    best_candidate <- final_score %>%
      right_join(data.frame(candidate = candidates_boot)) %>%
      arrange(desc(values)) %>%
      top_n(1) %>%
      pull(candidate)
    sub_score_df <- score_df %>%
        right_join(data.frame(candidate = candidates_boot))
    final_score_boot <- f_M(scores = sub_score_df)
    best_candidate_boot <- final_score_boot %>%
      arrange(desc(values)) %>%
      top_n(1) %>%
      pull(candidate)
    if (best_candidate != best_candidate_boot) return(FALSE)
  }
  return(TRUE)
}

# empirical criteria
avg_rank <- function(matrix,winner) {
  winner_idx <- which(rownames(matrix)==winner)
  rank_matrix <- apply(matrix,2,function(x) rank(1-x))
  avgrank <- mean(rank_matrix[winner_idx,])
  return(1-(avgrank-1)/(ncol(matrix)-1))
}

condorcet_rate <- function(matrix, winner) {
  winner_idx <- which(rownames(matrix)==winner)
  condorcet_win <- sapply(seq(nrow(matrix)), function(x) {
      maj_winner <- NA
      if (winner_idx != x) {
        sub_matrix <- matrix[c(winner_idx,x),]
        get_best_candidate_per_judge <- unlist(apply(sub_matrix, 2, function(z) names(which.max(z))))
        nb_votes_winner <- sort(table(get_best_candidate_per_judge), decreasing = TRUE)[winner]
        maj_winner <- ifelse(nb_votes_winner>ncol(sub_matrix)/2, winner, 'other')
      }
      return(maj_winner)
  })[-winner_idx]
  return(mean(condorcet_win==winner))
}

generalization_criterion <- function(score_df, f_M, scores_to_keep=NULL, n_bootstrap=100, p_excluded_judges=.1) {
  judge_df <- score_df %>% ungroup() %>%
    filter(name_score %in% scores_to_keep) %>%
    select(name_score) %>%
    unique()
  n_judge <- nrow(judge_df)
  n_excluded_judges <- max(round(n_judge * p_excluded_judges),1)
  if (n_excluded_judges==1) {n_bootstrap = min(n_judge,n_bootstrap)}
  generalization_res <- rep(NA, n_bootstrap)
  for (k in seq(n_bootstrap)) {
    excluded_judges_idx <- sample(n_judge, n_excluded_judges)
    if (n_excluded_judges==1) {excluded_judges_idx=k}
    excluded_judges = judge_df$name_score[excluded_judges_idx]
    final_score_train <- f_M(scores = score_df, scores_to_keep = scores_to_keep[!(scores_to_keep %in% excluded_judges)]) %>%
      mutate(rank1 = rank(1 - overall)) %>%
      select(candidate,rank1)
    final_score_test = NULL
    for (i in seq(n_excluded_judges)) {
      rank_1judge = f_M(scores = score_df, scores_to_keep = excluded_judges[i]) %>%
        ungroup() %>%
        mutate(rank2 = rank(1 - overall)) %>%
        select(candidate,rank2)
      if (is.null(final_score_test)) {
        final_score_test=rank_1judge
      } else {final_score_test=inner_join(final_score_test,rank_1judge,by="candidate")}
    }
    all_ranks = inner_join(final_score_test,final_score_train,by="candidate")
    cols_test = grep("rank2",colnames(all_ranks),value=T)
    generalization_res[k] <- sum(sapply(cols_test, function(x)
      cor(all_ranks$rank1,all_ranks[,x], method="spearman")/n_excluded_judges))
  }
  return(mean(generalization_res, na.rm = T))
}
