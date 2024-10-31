## ----
## Set parameters, put your own
## ----
date = "241025"
score_path = "../../compute_metrics/scores/"
source("../generic_functions/ranking_process.R")
source("../generic_functions/load_scores_SB_silico.R")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_dnam_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_dnam_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

## ----
## Load libraries
## ----
library(dplyr)
library(funkyheatmap)
library(tibble)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(see)

## ----
## Load scores
## ----
res = load_data(date, score_path)
scores = res$scores
time = res$time
rm(res)

## ----
## Compute intermediate and overall scores
## ----
ranks = ranking_consensus(scores1=scores, scores2=time)

# inter scoring
score_inter <- ranking_step1(scores1=scores, scores2=time) %>%
  coerce_pearson() %>%
  group_by(name_score,dataset) %>%
  mutate(archetype_best=max(trendval, na.rm=T),
         archetype_worst=min(trendval, na.rm=T)) %>%
  group_by(name_score,dataset,cat_score) %>%
  mutate(rank=rank(trendval, ties.method = "average", na.last="keep")) %>%
  ungroup() %>%
  filter(!is.na(trendval)) %>%
  mutate(rank=rank/max(rank, na.rm=T)) %>%
  group_by(candidate,dataset,name_score,cat_score) %>%
  mutate(d_best=euclidean_vector_dist(trendval,archetype_best),
         d_worst=euclidean_vector_dist(trendval,archetype_worst),
         topsis=d_worst/(d_worst+d_best)) %>%
  group_by(candidate,dataset,cat_score) %>%
  summarise(score_inter1_raw=geomMean(trendval),
            score_inter1_topsis=geomMean(topsis),
            score_inter1_rank=geomMean(rank),
            score_inter1=mean(score_inter1_raw,score_inter1_topsis,score_inter1_rank, trim=0)) %>%
  group_by(candidate,dataset) %>%
  mutate(weights_cat=convert_to_weights(weights_dic_values,cat_score),
         score_inter2_raw=weighgeomMean(score_inter1_raw, weights_cat),
         score_inter2_topsis=weighgeomMean(score_inter1_topsis, weights_cat),
         score_inter2_rank=weighgeomMean(score_inter1_rank, weights_cat),
         score_inter2=mean(score_inter2_raw,score_inter2_topsis,score_inter2_rank, trim=0)) %>%
  ungroup() %>%
  select(candidate,dataset,cat_score,score_inter1,score_inter2)

score_all = left_join(score_inter, ranks, by='candidate')

## ----
## Retrieve best FS per DeconvTool
## ----
best_FS_per_tool = ranks %>%
  mutate(setting=sapply(candidate, function(x)
    paste0(head(strsplit(x,"-")[[1]],2),collapse="-")),
         fs=sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
         block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1])) %>%
  group_by(setting) %>%
  filter(overall==max(overall)) %>% ungroup() %>%
  select(setting,fs,candidate,block)

## ----
## Plot figures 2 and 3 panels C,D
## ----
settings = list(meth_dnam_sup,meth_dnam_unsup,
                meth_rna_sup,meth_rna_unsup)
blocks = rep(c("dnam","rna"),each=2)
df_fig2 = lapply(seq_along(settings), function(x) {
  res = score_all %>%
    mutate(DeconvTool=sapply(candidate,function(y) strsplit(y,"-")[[1]][2])) %>%
    filter(candidate %in% best_FS_per_tool$candidate[best_FS_per_tool$block==blocks[x]],
           DeconvTool %in% settings[[x]])
  full_join(res %>% filter(!duplicated(paste(candidate,dataset))) %>%
              rename("aggregated":=score_inter2) %>%
              select(candidate,DeconvTool,dataset,overall,aggregated),
            plyr::join_all(lapply(sort(unique(res$cat_score)), function(y) {
              res %>% filter(cat_score==y) %>%
                rename(!!y:=score_inter1) %>%
                select(candidate,DeconvTool,dataset,overall,y)
            }), by=c('candidate','DeconvTool','dataset','overall'), type='full'),
            by=c('candidate','DeconvTool','dataset','overall'))})
names(df_fig2) = c("dnam-sup","dnam-unsup","rna-sup","rna-unsup")
saveRDS(df_fig2, paste0(folder,"/df_res.rds"))

for (x in seq_along(df_fig2)) {
  df_fig2[[x]] = df_fig2[[x]] %>%
    ungroup() %>%
    mutate(dataset=sapply(dataset,function(y) strsplit(y,"-")[[1]][1])) %>%
    group_by(candidate) %>%
    arrange(desc(overall)) %>%
    ungroup() %>%
    mutate(id = as.character(seq_len(n())),
           id0 = seq_len(n()),
           DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool))
  df_fig2[[x]]$rank = 0
  for (i in df_fig2[[x]]$candidate) {
    idx = which(df_fig2[[x]]$candidate==i)
    if (idx[1]==1) {df_fig2[[x]]$rank[idx] = 1}
    else {df_fig2[[x]]$rank[idx] = df_fig2[[x]]$rank[idx[1]-1]+1}
    df_fig2[[x]]$DeconvTool[idx[-2]] = ""
  }
}

column_info <- tribble(
  ~ id,             ~ name,                 ~ geom,  ~ group, ~ palette, ~ width, ~ legend, ~ hjust, ~ overlay,
  "DeconvTool",     "DeconvTool",           "text",  NA,      NA,        4,       FALSE,    0,       FALSE,
  "dataset",        "Datasets",             "text",  NA,      NA,        2,       FALSE,    NA,       FALSE,
  "overall",        "Overall score",        "bar",   NA,      "greys",   3,       FALSE,    NA,      FALSE,
  "aggregated",     "Aggregated scores",    "bar",   NA,      "greys",   3,       FALSE,    NA,      FALSE,
  "raw_perf",       "Raw performance",      "bar",   NA,      "blues",   3,       FALSE,    NA,      FALSE,
  "stab_perf",      "Stability performance","bar",   NA,      "reds",    3,       FALSE,    NA,      FALSE,
  "time",           "Time performance",     "bar",   NA,      "yellows", 3,       FALSE,    NA,      FALSE,
)

plot_list = lapply(df_fig2, function(x)
  funky_heatmap(
    x,
    column_info = column_info,
    row_info = data.frame(id = x$id, group = factor(x$rank)),
    scale_column = FALSE,
    palettes = list(
      blues = "#477CBD",
      reds = "#B54893",
      yellows = "#E8B658",
      greys = "#484748"),
    position_args = position_arguments(row_height = .3,
                                       row_space = 0.13,
                                       row_bigspace = 0.5)))
rects = lapply(plot_list, function(x) {
  idxs = which(!duplicated(x$layers[[3]]$data$xmax[x$layers[[3]]$data$label=="overall"]))
  data.frame(xmin = x$layers[[3]]$data$xmin[idxs],
             xmax = x$layers[[3]]$data$xmax[idxs],
             ymin = x$layers[[3]]$data$ymax[idxs],
             ymax = x$layers[[3]]$data$ymin[c(idxs-1,sum(x$layers[[3]]$data$label=="overall"))])
  })

heights = c(5.6,6.2,12.7,6.5)
for (i in seq_along(plot_list)) {
  print(plot_list[[i]]+geom_rect(data = rects[[i]], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                           fill = "#484748", linewidth = .3, color="deeppink4"))
  ggsave(paste0(folder,"/panel_",names(plot_list)[i],".pdf"), width = 10, height = heights[i])
}
