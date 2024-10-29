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
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

## ----
## Load scores & rank
## ----
source("../../src/2_ranking_ranking_procedure_functions.R")
df_scores = list()
df_time = list()
for (experiment in c('missinginT','addedinT')) {
  source("../../src/load_scores_SB_missing.R")
  res = load_data(experiment)
  scores = res$scores
  time = res$time
  rm(res)
  df_scores[[experiment]] = scores
  df_time[[experiment]] = time
}
source("../../src/load_scores_SB.R")
res = load_data("231027")
df_scores[["baseline"]] = res$scores %>% mutate(candidate=paste(candidate,"baseline",sep = "-"))
df_time[["baseline"]] = res$time %>% mutate(candidate=paste(candidate,"baseline",sep = "-"))
rm(res)
rm(scores,time)

# each setting is ranked independently
bestFS = unlist(lapply(readRDS("fig2/df_res.rds"), function(x) x %>% pull(candidate) %>% unique()))
winners = readRDS("fig2/df_res.rds")
winners = lapply(winners[grep('-sup',names(winners))], function(x) x %>% 
                   arrange(desc(overall)) %>% 
                   mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                                 DeconvTool=='rlr'~'RLR',
                                                 DeconvTool=='elasticnet'~'Elastic net',
                                                 DeconvTool=='fardeep'~'FARDEEP',
                                                 DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                                 DeconvTool=='ols'~'OLS',
                                                 DeconvTool=='svr'~'SVR', .default = DeconvTool)) %>%
                   pull(DeconvTool) %>% 
                   unique())

rank_consensus = lapply(seq_along(df_scores), function(x)
  ranking_consensus(df_scores[[x]],df_time[[x]]))
rank_consensus = do.call(rbind,rank_consensus)
rank_consensus = rank_consensus %>%
  mutate(Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
         FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
         Missing=sapply(candidate, function(x) strsplit(x,"-")[[1]][4]),
         Supervised=ifelse(DeconvTool %in% c(meth_met_sup,meth_rna_sup),"sup","unsup"),
         Setting=paste(Block,Supervised,sep="-"),
         candidate2=paste(Block,DeconvTool,FS,sep="-")) %>%
  filter(candidate2 %in% bestFS, Supervised=='sup')
rank_consensus$Missing[rank_consensus$Missing=="missinginT"] = "-1Type"
rank_consensus$Missing[rank_consensus$Missing=="baseline"] = "AllTypes"
rank_consensus$Missing[rank_consensus$Missing=="addedinT"] = "+1Type"
rank_consensus$Block[rank_consensus$Block=="met"] = "DNAm"
rank_consensus$Block[rank_consensus$Block=="rna"] = "RNA"

rank_consensus = rank_consensus %>%
  mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                DeconvTool=='rlr'~'RLR',
                                DeconvTool=='elasticnet'~'Elastic net',
                                DeconvTool=='fardeep'~'FARDEEP',
                                DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                DeconvTool=='ols'~'OLS',
                                DeconvTool=='svr'~'SVR', .default = DeconvTool))

## ----
## Plot scatter
## ----
rank_consensus$Missing = factor(rank_consensus$Missing, levels=c("-1Type","AllTypes","+1Type"))
rank_consensus$DeconvTool = factor(rank_consensus$DeconvTool, levels=unique(c(winners$`rna-sup`, winners$`met-sup`)))

ggplot(rank_consensus, aes(x=Missing, y=overall, color=DeconvTool)) +
  geom_point(shape=8) +
  scale_color_manual(values=custom_palette) +
  geom_line(aes(group=DeconvTool), linewidth=1.2) +
  xlab("") +
  facet_wrap(~Block) +
  theme_modern(legend.title.size = 16, legend.text.size = 14)
ggsave(paste0(folder,"/missing.pdf"), width=8,height=4)

rank_consensus %>%
  filter(Block=='RNA') %>%
  select(overall,DeconvTool,Missing) %>%
  group_by(DeconvTool) %>%
  mutate(alltype=max(overall),delta = alltype-overall) %>%
  filter(delta!=0,delta>0.25)
rank_consensus %>%
  filter(Block=='RNA') %>%
  select(overall,DeconvTool,Missing) %>%
  group_by(DeconvTool) %>%
  mutate(alltype=max(overall),delta = alltype-overall) %>%
  filter(delta!=0,delta<0.1)

