## ----
## Set parameters, put your own
## ----
set.seed(2024)
n_permutations = 1000
score_path = "../../compute_metrics/scores/"
date = "241025"
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/ranking_process.R")
source("../generic_functions/ranking_pval.R")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_dnam_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_dnam_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

## ----
## Load libraries
## ----
library(dplyr)
library(foreach)
library(funkyheatmap)
library(tibble)
library(ggplot2)
library(ggpubr)
library(see)

## ----
## Functions
## ----
scale_col_fun = function(x) {
  if (x%in%c(2,4)) {
    scale_color_metro(reverse=T) # colors for unsupervised panels
  }
  else if (x==1) {scale_color_manual(values=custom_palette[c(1,4,7,10,12)])} # colors for the supervised DNAm panel
  else if (x==3) {scale_color_manual(values=custom_palette)} # colors for the supervised RNA panel
}

## ----
## Compute p-values
## ----
res = load_data(date, score_path)
scores = res$scores
time = res$time
rm(res)
df_scores_inter = readRDS("figure2CD_figure3CD/df_res.rds")
df_scores_norm = lapply(df_scores_inter, function(x)
  ranking_step1(scores,time) %>%
  coerce_pearson() %>%
  filter(candidate %in% x$candidate))

# compute p-values
for (block in c("dnam","rna")) {
  print(paste("Running", block))
  for (class in c('sup','unsup')) {
    print(paste("Running", class))
    if (!file.exists(paste0(folder,"/pval_",block,"_",class,".rds"))) {
      pval_res <- paired_permutation_test(df_scores_norm[[paste(block,class,sep='-')]],
                                          ranking_consensus_end,
                                          n_permutations=n_permutations)
      saveRDS(pval_res,paste0(folder,"/pval_",block,"_",class,".rds"))
    }
  }
}

## ----
## Prepare df
## ----
df_pval <- lapply(df_scores_inter, function(y)
  y %>%
    mutate(Block = sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
           FS = sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
           DeconvTool = sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
           Supervised = ifelse(DeconvTool %in% c(meth_dnam_unsup,meth_rna_unsup),"Unsupervised","Supervised")) %>%
    filter(!duplicated(paste(DeconvTool,dataset))) %>%
    arrange(desc(overall)))
df_pval = lapply(df_pval, function(x) {
  df = data.frame(Dataset=c(rep("Overall",length(unique(x$candidate))),
                       x$dataset),
             Score=c(x$overall[!duplicated(x$candidate)],
                     x$aggregated),
             Block=unique(x$Block),
             DeconvTool=c(x$DeconvTool[!duplicated(x$candidate)],
                          x$DeconvTool),
             FS=c(x$FS[!duplicated(x$candidate)],
                  x$FS)) %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool))
  df %>%
    mutate(DeconvTool = factor(DeconvTool, levels=unique(df %>% filter(Dataset=="Overall") %>%
                                                           arrange(desc(Score)) %>%
                                                           pull(DeconvTool))),
           Dataset = factor(Dataset, levels=c(unique(Dataset[Dataset!='Overall']),'Overall')))
  })

signif_pvalues = lapply(c(t(outer(c("dnam","rna"), c("sup","unsup"), paste, sep="_"))), function(x) {
  pval_res = readRDS(paste0(folder,"/pval_",x,".rds"))$pvalues
  df = data.frame(pairs = names(pval_res), p = pval_res)
  rownames(df) = NULL
  df %>%
    mutate(group1 = sapply(pairs, function(x) strsplit(strsplit(x, ">")[[1]][1], "-")[[1]][2]),
           group2 = sapply(pairs, function(x) strsplit(strsplit(x, ">")[[1]][2], "-")[[1]][2]),
           p_lab = if_else(p < .001, "***",
                           if_else(p < .01, "**",
                                   if_else(p < .05, "*", ""))),
           group1 = case_when(group1=='nnls'~'NNLS',
                              group1=='rlr'~'RLR',
                              group1=='elasticnet'~'Elastic net',
                              group1=='fardeep'~'FARDEEP',
                              group1=='fardeepsto'~'FARDEEP_sto',
                              group1=='ols'~'OLS',
                              group1=='svr'~'SVR', .default = group1),
           group2 = case_when(group2=='nnls'~'NNLS',
                              group2=='rlr'~'RLR',
                              group2=='elasticnet'~'Elastic net',
                              group2=='fardeep'~'FARDEEP',
                              group2=='fardeepsto'~'FARDEEP_sto',
                              group2=='ols'~'OLS',
                              group2=='svr'~'SVR', .default = group2),
           group1 = factor(group1, levels = levels(df_pval[[sub('_','-',x)]]$DeconvTool)),
           group2 = factor(group2, levels = levels(df_pval[[sub('_','-',x)]]$DeconvTool)),
           Score=sapply(group1, function(y)
             df_pval[[sub('_','-',x)]] %>% filter(Dataset=='Overall',DeconvTool==y) %>% pull(Score))) %>%
    select(!pairs) %>%
    filter(p<=0.1)
  })
names(signif_pvalues) = c(t(outer(c("dnam","rna"), c("sup","unsup"), paste, sep="-")))

## ----
## Plot figure 2,3 panels A,B
## ----
for (i in seq_along(df_pval)) {
  ggplot(df_pval[[i]] %>%
           filter(Dataset=="Overall"), aes(x=DeconvTool, y=Score)) +
    geom_point(aes(color=DeconvTool), size=8, shape=17) +
    stat_pvalue_manual(data=signif_pvalues[[i]],
                       label = "p_lab", y = "Score", x = "group2", color = "group1", label.size = 9) +
    scale_col_fun(i) +
    ylim(c(0,1)) +
    ylab("Overall benchmark score") +
    theme_modern(axis.text.angle = 30, axis.text.size = 20) +
    guides(color="none", shape="none")
  ggsave(paste0(folder,"/scatterplot_",names(df_pval)[i],".pdf"), width = 8.5, height = 6.8)
}

## ----
## Plot supp figure 5
## ----
for (i in seq_along(df_pval)) {
  ggplot(df_pval[[i]] , aes(x=DeconvTool, y=Score)) +
    geom_point(aes(shape=Dataset, color=DeconvTool, size=Dataset)) +
    scale_col_fun(i) +
    scale_size_manual(values=c(rep(4.5,length(unique(df_pval[[i]]$Dataset))-1),6)) +
    scale_shape_manual(values = c(1:(length(unique(df_pval[[i]]$Dataset))-1),17)) +
    ylim(c(0,1)) +
    ylab("Benchmark score") +
    theme_modern(axis.text.angle = 30, axis.text.size = 22, axis.title.size = 24) +
    guides(color="none")
  ggsave(paste0("../supp_figs/supp5/",names(df_pval)[i],".pdf"), width = 10, height = 6.8)
}
