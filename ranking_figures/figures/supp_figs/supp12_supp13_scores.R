# Generate intermediate scores for supp figures 12 and 13
## ----
## Set parameters, put your own
## ----
date = "241025"
score_path = '../../compute_metrics/scores/'
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/load_scores_SB_invitro.R")
source("../generic_functions/load_scores_SB_invivo.R")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

## ----
## Load libraries
## ----
library(dplyr)

## ----
## Functions
## ----
convert_to_cat <- function(dict, vect) {
  tmp = rep(names(dict), lengths(dict))
  cat <- sapply(vect, function(x) {
    ifelse(x %in% unlist(dict), tmp[unlist(dict)==x], NA)
  })
  return(unname(unlist(cat)))
}
tbl_groups_nb <- function(tbl) {
  g <- groups(tbl)
  g <- unlist(lapply(g, as.character))
  nrow(unique(tbl[, g, drop = FALSE]))
}
CenterScaleNorm <-function(x) {
  #center 
  tr1 = x - mean(x, na.rm=T)
  #scale
  tr2 = tr1/sd(x, na.rm=T)
  return(pnorm(tr2))
}
coerce_pearson <- function(df_median_centerscale_transfo) {
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(Dataset,candidate,Source) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(Dataset,candidate,Source) %>%
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
euclidean_vector_dist <- function(vec1,vec2) {
  vec1 <- vec1[(!is.na(vec1)) & (!is.na(vec2))]
  vec2 <- vec2[(!is.na(vec1)) & (!is.na(vec2))]
  sqrt(sum((vec1 - vec2)^2))
}
geomMean <- function(x) {
  exp(mean(log(x), na.rm=T))
}
convert_to_weights <- function(dict, vect) {
  cat <- sapply(vect, function(x) {
    ifelse(x %in% names(dict), dict[x], NA)
  })
  return(unname(unlist(cat)))
}
weighgeomMean <- function(x,w) {
  prod(mapply(function(a,b) a^b, a=x, b=w), na.rm=T)^(1/sum(w))
}
weights_dic_values=list("consensus"=.5, "raw_perf"=1, "stab_perf"=.5, "time"=.5)

## ----
## Load scores
## ----
res = load_data(date,score_path)
scores_silico = res$scores
time_silico = res$time
rm(res)

res = load_data_vitro(score_path)
scores_vitro = res$scores
time_vitro = res$time
rm(res)

res = load_data_vivo(score_path)
scores_vivo = res$scores
time_vivo = res$time
rm(res)

scores = rbind(cbind(scores_silico,data.frame(Source='In silico')),
               cbind(scores_vitro,data.frame(Source='In vitro')),
               cbind(scores_vivo,data.frame(Source='In vivo')))
time = rbind(cbind(time_silico,data.frame(Source='In silico')),
             cbind(time_vitro,data.frame(Source='In vitro')),
             cbind(time_vivo,data.frame(Source='In vivo')))
rm(scores_silico,scores_vitro,scores_vivo,
   time_silico,time_vitro,time_vivo)

## ----
## Load bestFS and filter scores
## ----
bestFS = unlist(lapply(readRDS("../main_figs/figure2CD_figure3CD/df_res.rds"), function(x) unique(x$candidate)))

scores = scores %>% filter(candidate %in% bestFS)
time = time %>% filter(candidate %in% bestFS)

## ----
## Renormalization of the scores just using the max and min for all datasets/all methods
## ----
score_df <- bind_rows(scores,time)
colnames(score_df) <- c("score","name_score","simulation","Dataset","candidate","Source")

# implement the column cat_score
cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                  "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                  "time"=c("time time","time sd"))
score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)

# log10 transfo of time for odg
score_df$score[score_df$name_score=="time time"] = log10(1+score_df$score[score_df$name_score=="time time"])

# compute stability sd_g and time across sims (step 0)
tmp1 <- score_df %>%
  filter(name_score=="mae perf_g") %>%
  group_by(Dataset, candidate, Source)
tmp2 <- score_df %>%
  filter(name_score=="pearson perf_g") %>%
  group_by(Dataset, candidate, Source)
tmp3 <- score_df %>%
  filter(name_score=="rmse perf_g") %>%
  group_by(Dataset, candidate, Source)
tmp4 <- score_df %>%
  filter(name_score=="time time") %>%
  group_by(Dataset, candidate, Source)
std_df <- data.frame("val"= c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                              tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                              tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std),
                              tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(std)),
                     "name_score"=c(rep('mae sd_g',tbl_groups_nb(tmp1)),
                                    rep('pearson sd_g',tbl_groups_nb(tmp2)),
                                    rep('rmse sd_g',tbl_groups_nb(tmp3)),
                                    rep('time sd',tbl_groups_nb(tmp4))),
                     "Dataset"=c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Dataset),
                                 tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Dataset),
                                 tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Dataset),
                                 tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Dataset)),
                     "candidate"=c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                   tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                   tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate),
                                   tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(candidate)),
                     "Source"=c(tmp1 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Source),
                                tmp2 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Source),
                                tmp3 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Source),
                                tmp4 %>% summarise(std=sd(score, na.rm=T)) %>% pull(Source))) %>%
  filter(Source == 'In silico')

# compute stability sd_c and sd_s across sims (step 0)
std_df2 <- score_df %>%
  filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
  group_by(name_score, Dataset, candidate, Source) %>%
  summarise(val=median(score, na.rm=T)) %>%
  filter(Source == 'In silico')

# merge resulting sd dataframes (step 0)
std_df <- bind_rows(std_df, std_df2)

# compute all medians across sims for raw perf and time and consensus categories (step 0)
score_df_median <- score_df %>%
  filter(cat_score=="raw_perf") %>%
  group_by(name_score, Dataset, candidate, Source) %>%
  summarise(val=median(score, na.rm=T))
time_df_median <- score_df %>%
  filter(cat_score=="time") %>%
  group_by(name_score, Dataset, candidate, Source) %>%
  summarise(val=median(score, na.rm=T))
consensus_df_median <- score_df %>%
  filter(cat_score=="consensus") %>%
  group_by(name_score, Dataset, candidate, Source) %>%
  summarise(val=median(score, na.rm=T))
df_median <- bind_rows(score_df_median, std_df, time_df_median, consensus_df_median)
rm(std_df, std_df2, score_df, time_df_median, consensus_df_median)

# add the column cat_score
df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)

# normalize for each dataset and score (step 1)
df_median_centerscale <- df_median %>%
  group_by(name_score) %>%
  mutate(normval=CenterScaleNorm(val))
rm(df_median)

# transform scores such that 1 is the best score (step 2)
df_median_centerscale$trendval <- 1 - df_median_centerscale$normval
df_median_centerscale$trendval[grep("pearson med",df_median_centerscale$name_score)] <- 1 - df_median_centerscale$trendval[grep("pearson med",df_median_centerscale$name_score)]
df_median_centerscale$trendval[grep("pearson perf",df_median_centerscale$name_score)] <- 1 - df_median_centerscale$trendval[grep("pearson perf",df_median_centerscale$name_score)]

## ----
## Aggregate
## ----
# inter scoring
score_inter <- df_median_centerscale %>% 
  coerce_pearson() %>%
  group_by(name_score,Dataset) %>%
  mutate(archetype_best=max(trendval, na.rm=T),
         archetype_worst=min(trendval, na.rm=T)) %>%
  group_by(name_score,Dataset,cat_score) %>%
  mutate(rank=rank(trendval, ties.method = "average", na.last="keep")) %>%
  ungroup() %>%
  filter(!is.na(trendval)) %>%
  mutate(rank=rank/max(rank, na.rm=T)) %>%
  group_by(candidate,Dataset,name_score,cat_score) %>%
  mutate(d_best=euclidean_vector_dist(trendval,archetype_best),
         d_worst=euclidean_vector_dist(trendval,archetype_worst),
         topsis=d_worst/(d_worst+d_best)) %>%
  group_by(candidate,Dataset,cat_score,Source) %>%
  summarise(score_inter1_raw=geomMean(trendval),
            score_inter1_topsis=geomMean(topsis),
            score_inter1_rank=geomMean(rank),
            score_inter1=mean(score_inter1_raw,score_inter1_topsis,score_inter1_rank, trim=0)) %>%
  group_by(candidate,Dataset,Source) %>%
  summarise(weights_cat=convert_to_weights(weights_dic_values,cat_score),
         score_inter2_raw=weighgeomMean(score_inter1_raw, weights_cat),
         score_inter2_topsis=weighgeomMean(score_inter1_topsis, weights_cat),
         score_inter2_rank=weighgeomMean(score_inter1_rank, weights_cat),
         score_inter2=mean(score_inter2_raw,score_inter2_topsis,score_inter2_rank, trim=0)) %>%
  ungroup() %>%
  select(candidate,Dataset,score_inter2,Source) %>%
  filter(!duplicated(paste(Dataset,candidate))) %>%
  group_by(Dataset) %>%
  mutate(dataset_ease_decon = mean(score_inter2)) %>% ungroup() %>%
  mutate(Block=sapply(Dataset, function(x) strsplit(x,'-')[[1]][2]))
saveRDS(score_inter, file='supp12_supp13/scores_inter.rds')
