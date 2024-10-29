## ----
## Set parameters
## ----
library(dplyr)
library(funkyheatmap)
library(tibble)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(see)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
date = "231027"
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")
scale_fill_fun = function(x) {
  if (x%in%c(2,4)) {
    scale_fill_metro(reverse=T)
  }
  else if (x==3) {scale_fill_manual(values=custom_palette)}
  else if (x==1) {scale_fill_manual(values=custom_palette[c(1,4,7,10,12)])}
}

## ----
## Load df
## ----
df_fig2 = readRDS("../fig/fig2/df_res.rds")

## ----
## Plot supp = rank per dataset
## ----
df_fig2 = lapply(df_fig2, function(df) {
  x = df
  x$candidate = sapply(x$candidate, function(y)
    sub("-none"," - No feature selection",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("-hvf"," - Highly variable features",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("-toast"," - TOAST",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("met-","",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("rna-","",y))
  x$dataset = sapply(x$dataset, function(y)
    strsplit(y,'-')[[1]][1])
  x$candidate = factor(x$candidate, levels=x %>% arrange(desc(overall)) %>% pull(candidate) %>% unique())
  x = x %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool))
  x$DeconvTool = factor(x$DeconvTool, levels=x %>% arrange(desc(overall)) %>% pull(DeconvTool) %>% unique())
  x})

widths=c(7,7,8,7)
heights=c(4,4,4,4)
for (i in seq_along(df_fig2)) {
  print(ggplot(df_fig2[[i]], aes(x=dataset, y=aggregated, fill=DeconvTool)) +
          geom_bar(stat='identity', position='dodge') +
          scale_fill_fun(i) +
          ylab("Aggregated benchmark score") +
          xlab("") +
          theme_modern(axis.text.angle = 30, axis.title.size = 15))
  ggsave(paste0(folder,"/",names(df_fig2)[i],".pdf"), width = widths[i], height = heights[i])
}

## ----
## Plot supp 2 for all data types (vitro, vivo, etc)
## ----
scores_inter = readRDS("../fig/fig2_invitro_invivo/score_inter.rds")
scores_inter = lapply(scores_inter, function(x)
  x %>% arrange(Source,dataset) %>%
    mutate(dataset=factor(dataset, levels=unique(dataset))))
scores_inter = lapply(seq_along(scores_inter), function(x) {
  scores_inter[[x]]$candidate = sapply(scores_inter[[x]]$candidate, function(y)
    sub("-none"," - No feature selection",y))
  scores_inter[[x]]$candidate = sapply(scores_inter[[x]]$candidate, function(y)
    sub("-hvf"," - Highly variable features",y))
  scores_inter[[x]]$candidate = sapply(scores_inter[[x]]$candidate, function(y)
    sub("-toast"," - TOAST",y))
  scores_inter[[x]]$candidate = sapply(scores_inter[[x]]$candidate, function(y)
    sub("met-","",y))
  scores_inter[[x]]$candidate = sapply(scores_inter[[x]]$candidate, function(y)
    sub("rna-","",y))
  scores_inter[[x]]$candidate = factor(scores_inter[[x]]$candidate, levels=levels(df_fig2[[x]]$candidate))
  scores_inter[[x]]})
names(scores_inter) = names(df_fig2)
widths=c(10,10,15,10)
for (i in seq_along(scores_inter)) {
  print(ggplot(scores_inter[[i]], aes(x=dataset, y=aggregated, fill=candidate)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_fun(i) +
  ylab("Aggregated benchmark score") +
  xlab("") +
  theme_modern(axis.text.angle = 30))
  ggsave(paste0(folder,"/allsources_",names(scores_inter)[i],".pdf"), width = widths[i], height = heights[i])
}
