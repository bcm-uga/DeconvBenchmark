## ----
## Set parameters, put your own
## ----
date = "241025"
score_path = "../../compute_metrics/scores/"
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/ranking_process.R")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_dnam_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_dnam_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

## ----
## Load libraries
## ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(see)

## ----
## Load scores & keep only multi-omic datasets
## ----
res = load_data(date,score_path)
scores = res$scores
time = res$time
rm(res)
datasets = sapply(unique(scores$dataset), function(x) strsplit(x,"-")[[1]][1])
multi_omics = datasets[duplicated(datasets)]
scores = scores %>% filter(dataset %in% c(sapply(c("dnam","rna"), function(x) paste(multi_omics,x,sep="-"))))
time = time %>% filter(dataset %in% c(sapply(c("dnam","rna"), function(x) paste(multi_omics,x,sep="-"))))

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
best_fs = lapply(readRDS("figure2CD_figure3CD/df_res.rds"), function(x) unique(x$candidate))
df_fig = lapply(best_fs, function(x)
  score_all %>%
    filter(candidate %in% x))
ranks = ranks %>%
  filter(candidate %in% unlist(best_fs)) %>%
  mutate(DeconvTool=sapply(candidate,function(y) strsplit(y,"-")[[1]][2]),
         Block=sapply(candidate,function(y) strsplit(y,"-")[[1]][1]),
         Supervised=ifelse(DeconvTool%in%c(meth_dnam_sup,meth_rna_sup),'Supervised','Unsupervised')) %>%
  arrange(desc(overall))
ranks$DeconvTool[ranks$DeconvTool=="rlr"]='RLR'
ranks$DeconvTool[ranks$DeconvTool=="fardeep"]='FARDEEP'
ranks$DeconvTool[ranks$DeconvTool=="fardeep_sto"]='FARDEEP_sto'
ranks$DeconvTool[ranks$DeconvTool=="nnls"]='NNLS'
ranks$DeconvTool[ranks$DeconvTool=="ols"]='OLS'
ranks$DeconvTool[ranks$DeconvTool=="svr"]='SVR'
ranks$DeconvTool[ranks$DeconvTool=="elasticnet"]='Elastic net'
ranks$Block[ranks$Block=='dnam'] = 'DNAm'
ranks$Block[ranks$Block=='rna'] = 'RNA'
ranks$DeconvTool = factor(ranks$DeconvTool, levels=unique(ranks$DeconvTool))
  
## ----
## Plot figure 4A
## ----
ggplot(ranks, aes(y=overall,x=DeconvTool,color=Block)) +
  geom_point(shape=17,size=3) +
  facet_grid(~Supervised,scales='free_x') +
  scale_color_manual(values=c("gold2","royalblue3")) +
  ylab('Overall benchmark score') +
  xlab('') +
  theme_modern(axis.text.angle = 30, axis.text.size = 16) +
  geom_vline(xintercept = c(1:12), alpha=.1)
ggsave(paste0(folder,"/plot.pdf"),width=12,height=5)
