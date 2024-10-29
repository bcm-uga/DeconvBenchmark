weights_dic_values=list("consensus"=.5, "raw_perf"=1, "stab_perf"=.5, "time"=.5)

#######
####### BASIC FUNCTIONS ####### 
#######
require(dplyr, quietly = T)
library(ggplot2)
library(see)
options(dplyr.summarise.inform = FALSE)

geomMean <- function(x) {
  exp(mean(log(x), na.rm=T))
}
weighMean <- function(x,weights) {
  weighted.mean(x, weights)
}
weighgeomMean <- function(x,w) {
  prod(mapply(function(a,b) a^b, a=x, b=w), na.rm=T)^(1/sum(w))
}

CenterScaleNorm <-function(x) {
  #center 
  tr1 = x - mean(x, na.rm=T)
  #scale
  tr2 = tr1/sd(x, na.rm=T)
  return(pnorm(tr2))
}

convert_to_cat <- function(dict, vect) {
  tmp = rep(names(dict), lengths(dict))
  cat <- sapply(vect, function(x) {
    ifelse(x %in% unlist(dict), tmp[unlist(dict)==x], NA)
  })
  return(unname(unlist(cat)))
}

convert_to_weights <- function(dict, vect) {
  cat <- sapply(vect, function(x) {
    ifelse(x %in% names(dict), dict[x], NA)
  })
  return(unname(unlist(cat)))
}

euclidean_vector_dist <- function(vec1,vec2) {
  vec1 <- vec1[(!is.na(vec1)) & (!is.na(vec2))]
  vec2 <- vec2[(!is.na(vec1)) & (!is.na(vec2))]
  sqrt(sum((vec1 - vec2)^2))
}

tbl_groups_nb <- function(tbl) {
  g <- groups(tbl)
  g <- unlist(lapply(g, as.character))
  nrow(unique(tbl[, g, drop = FALSE]))
}

#######
####### RANK FUNCTION ####### 
#######
# input scores1 is the df with all scores except time
# input scores2 is the df with time scores
coerce_pearson <- function(df_median_centerscale_transfo) {
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           val = NA,
           cat_score = "stab_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  return(df_median_centerscale_transfo)
}

# normalization of primary metrics
ranking_norm <- function(scores1, scores2, cat_score=NULL) {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("time time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s","mae consensus")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time time","time sd"),
                      "consensus"=c("mae consensus","consensus sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # log10 transfo of time for odg
  score_df$score[score_df$name_score=="time time"] = log10(1+score_df$score[score_df$name_score=="time time"])
  
  # normalize for each dataset and score (step 1)
  df_centerscale <- score_df %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(score))
  rm(score_df)
  
  # transform scores such that 1 is the best score (step 2)
  df_centerscale_transfo <- df_centerscale
  df_centerscale_transfo$trendval <- 1 - df_centerscale_transfo$normval
  df_centerscale_transfo$trendval[grep("pearson med",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson med",df_centerscale_transfo$name_score)]
  df_centerscale_transfo$trendval[grep("pearson perf",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson perf",df_centerscale_transfo$name_score)]
  rm(df_centerscale)
  
  return(df_centerscale_transfo)
}

# normalization of primary metrics + compute secondary metrics
ranking_step1 <- function(scores1, scores2, cat_score=NULL) {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("time time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s","mae consensus")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time time","time sd"),
                      "consensus"=c("mae consensus","consensus sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # log10 transfo of time for odg
  score_df$score[score_df$name_score=="time time"] = log10(1+score_df$score[score_df$name_score=="time time"])
  
  # compute stability sd_g and time across sims (step 0)
  tmp1 <- score_df %>%
    filter(name_score=="mae perf_g") %>%
    group_by(dataset, candidate)
  tmp2 <- score_df %>%
    filter(name_score=="pearson perf_g") %>%
    group_by(dataset, candidate)
  tmp3 <- score_df %>%
    filter(name_score=="rmse perf_g") %>%
    group_by(dataset, candidate)
  tmp4 <- score_df %>%
    filter(name_score=="time time") %>%
    group_by(dataset, candidate)
  tmp5 <- score_df %>%
    filter(name_score=="mae consensus") %>%
    group_by(dataset, candidate)
  std_df <- data.frame("val"= c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                                tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                                tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                                tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                                tmp5 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',tbl_groups_nb(tmp1)),
                                      rep('pearson sd_g',tbl_groups_nb(tmp2)),
                                      rep('rmse sd_g',tbl_groups_nb(tmp3)),
                                      rep('time sd',tbl_groups_nb(tmp4)),
                                      rep('consensus sd',tbl_groups_nb(tmp5))),
                       "dataset"=c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(dataset),
                                   tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(dataset),
                                   tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(dataset),
                                   tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(dataset),
                                   tmp5 %>% summarise(std=sd(score, na.rm=T)) %>% pull(dataset)),
                       "candidate"=c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                     tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                     tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                     tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                     tmp5 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate)))
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_df %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score, na.rm=T))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw perf and time and consensus categories (step 0)
  score_df_median <- score_df %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score, na.rm=T))
  time_df_median <- score_df %>%
    filter(cat_score=="time") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score, na.rm=T))
  consensus_df_median <- score_df %>%
    filter(cat_score=="consensus") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score, na.rm=T))
  df_median <- bind_rows(score_df_median, std_df, time_df_median, consensus_df_median)
  rm(std_df, std_df2, score_df, time_df_median, consensus_df_median)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  return(df_median_centerscale_transfo)
}

# global aggregation procedure
ranking_step3 <- function(scores, col_name, weights_dic_param1=weights_dic_values) {
  score_inter <- scores %>% filter(!is.na(get(col_name))) %>%
    group_by(candidate,dataset,cat_score) %>%
    summarise(score_g1=geomMean(get(col_name))) %>% #aggregate over scores
    group_by(candidate,dataset) %>%
    reframe(weights_cat=convert_to_weights(weights_dic_param1,cat_score),
            score_g2=weighgeomMean(score_g1, weights_cat)) %>% #aggregate over categories
    group_by(candidate) %>%
    summarise(overall=mean(score_g2, na.rm=F)) #aggregate over datasets
  score_final <- stack(score_inter[,2:ncol(score_inter)])
  score_final$candidate <- score_inter$candidate
  score_final = score_final %>% select(candidate, values) %>% rename("overall"=values)
  return(score_final)
}

# Sraw
ranking_raw_end <- function(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  if (!is.null(scores_to_keep)) {
    scores <- scores %>% filter(name_score %in% scores_to_keep)
  }
  return(ranking_step3(scores, 'trendval', weights_dic_param1=weights_dic_param))
}
ranking_raw <- function(scores1, scores2, cat_score=NULL, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  df_median_centerscale_transfo <- ranking_step1(scores1=scores1, scores2=scores2, cat_score=cat_score)
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  score_final <- ranking_raw_end(df_median_centerscale_transfo, scores_to_keep, weights_dic_param)
  return(score_final)
}

# Stopsis
ranking_topsis_end <- function(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  if (!is.null(scores_to_keep)) {
    scores <- scores %>% filter(name_score %in% scores_to_keep)
  }
  # Compute best and worst candidate + distances
  scores = scores %>% group_by(name_score,dataset) %>%
    mutate(archetype_best=max(trendval, na.rm=T),
           archetype_worst=min(trendval, na.rm=T)) %>% filter(!is.na(trendval)) %>%
    group_by(candidate,dataset,name_score,cat_score) %>%
    summarise(d_best=euclidean_vector_dist(trendval,archetype_best),
              d_worst=euclidean_vector_dist(trendval,archetype_worst),
              value=d_worst/(d_worst+d_best))
  # Compute topsis score
  return(ranking_step3(scores, 'value', weights_dic_param1=weights_dic_param))
}
ranking_topsis <- function(scores1, scores2, cat_score=NULL, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  df_median_centerscale_transfo <- ranking_step1(scores1=scores1, scores2=scores2, cat_score=cat_score)
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  score_final = ranking_topsis_end(df_median_centerscale_transfo, scores_to_keep, weights_dic_param)
  return(score_final)
}

# Srank
ranking_avgrank_end <- function(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  if (!is.null(scores_to_keep)) {
    scores <- scores %>% filter(name_score %in% scores_to_keep)
  }
  # Get the rank for each judge (dataset x name_score) per category and the avg rank
  ranks <- scores %>%
    group_by(name_score,dataset,cat_score) %>%
    mutate(rank=rank(trendval, ties.method = "average", na.last="keep"))
  ranks$rank = ranks$rank/max(ranks$rank, na.rm=T)
  return(ranking_step3(ranks,'rank',weights_dic_param1=weights_dic_param))
}
ranking_avgrank <- function(scores1, scores2, cat_score=NULL, scores_to_keep=NULL, weights_dic_param=weights_dic_values) {
  df_median_centerscale_transfo <- ranking_step1(scores1=scores1, scores2=scores2, cat_score=cat_score)
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  score_final = ranking_avgrank_end(df_median_centerscale_transfo, scores_to_keep, weights_dic_param)
  return(score_final)
}

# Sconsensus
ranking_consensus_end <- function(scores, scores_to_keep=NULL, weights_dic_param3=weights_dic_values) {
  rank1 <- ranking_raw_end(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_param3)
  rank2 <- ranking_topsis_end(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_param3)
  rank3 <- ranking_avgrank_end(scores, scores_to_keep=NULL, weights_dic_param=weights_dic_param3)
  rank <- bind_rows(rank1,rank2,rank3) %>% group_by(candidate) %>% summarise(overall=mean(overall, na.rm = T))
  score_final <- stack(rank[,2:ncol(rank)])
  score_final$candidate <- rank$candidate
  score_final = score_final %>% select(candidate, values) %>% rename("overall"=values)
  return(score_final)
}
ranking_consensus <- function(scores1, scores2, cat_score=NULL, scores_to_keep=NULL, weights_dic_param3=weights_dic_values) {
  rank1 <- ranking_raw(scores1, scores2, cat_score=cat_score, scores_to_keep=scores_to_keep, weights_dic_param=weights_dic_param3)
  rank2 <- ranking_topsis(scores1, scores2, cat_score=cat_score, scores_to_keep=scores_to_keep, weights_dic_param=weights_dic_param3)
  rank3 <- ranking_avgrank(scores1, scores2, cat_score=cat_score, scores_to_keep=scores_to_keep, weights_dic_param=weights_dic_param3)
  rank <- bind_rows(rank1,rank2,rank3) %>% group_by(candidate) %>% summarise(overall=mean(overall, na.rm = T))
  score_final <- stack(rank[,2:ncol(rank)])
  score_final$candidate <- rank$candidate
  score_final = score_final %>% select(candidate, values) %>% rename("overall"=values)
  return(score_final)
}