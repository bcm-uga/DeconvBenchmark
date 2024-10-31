## ----
## Set parameters, put your own
## ----
score_path = "../../compute_metrics/scores/"
source("../generic_functions/load_scores_SB_missing.R")
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/ranking_process.R")
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
library(ggplot2)
library(see)

## ----
## Load scores & rank
## ----
df_scores = list()
df_time = list()
for (experiment in c('241025_missing','241025_added')) {
  res = load_data_missing(experiment, score_path)
  scores = res$scores
  time = res$time
  rm(res)
  df_scores[[experiment]] = scores
  df_time[[experiment]] = time
}
res = load_data("241025", score_path)
df_scores[["baseline"]] = res$scores %>% mutate(candidate=paste(candidate,"baseline",sep = "-"))
df_time[["baseline"]] = res$time %>% mutate(candidate=paste(candidate,"baseline",sep = "-"))
rm(res,
   scores,time)

# each setting is ranked independently
bestFS = unlist(lapply(readRDS("figure2CD_figure3CD/df_res.rds"), function(x) x %>% pull(candidate) %>% unique()))
winners = readRDS("figure2CD_figure3CD/df_res.rds")
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
         NCellTypes=sapply(candidate, function(x) strsplit(x,"-")[[1]][4]),
         Supervised=ifelse(DeconvTool %in% c(meth_dnam_sup,meth_rna_sup),"sup","unsup"),
         Setting=paste(Block,Supervised,sep="-"),
         candidate2=paste(Block,DeconvTool,FS,sep="-")) %>%
  filter(candidate2 %in% bestFS, Supervised=='sup')
rank_consensus$NCellTypes[rank_consensus$NCellTypes=="241025_missing"] = "-1Type"
rank_consensus$NCellTypes[rank_consensus$NCellTypes=="baseline"] = "AllTypes"
rank_consensus$NCellTypes[rank_consensus$NCellTypes=="241025_added"] = "+1Type"
rank_consensus$Block[rank_consensus$Block=="dnam"] = "DNAm"
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
## Plot scatter figure 5 panel C
## ----
rank_consensus$NCellTypes = factor(rank_consensus$NCellTypes, levels=c("-1Type","AllTypes","+1Type"))
rank_consensus$DeconvTool = factor(rank_consensus$DeconvTool, levels=unique(c(winners$`rna-sup`, winners$`dnam-sup`)))

ggplot(rank_consensus, aes(x=NCellTypes, y=overall, color=DeconvTool)) +
  geom_point(shape=8) +
  scale_color_manual(values=custom_palette) +
  geom_line(aes(group=DeconvTool), linewidth=1.2) +
  xlab("") +
  facet_wrap(~Block) +
  theme_modern(legend.title.size = 16, legend.text.size = 14)
ggsave(paste0(folder,"/plot.pdf"), width=8,height=4)

# check the differences in overall score
rank_consensus %>%
  filter(Block=='RNA') %>%
  select(overall,DeconvTool,NCellTypes) %>%
  group_by(DeconvTool) %>%
  mutate(alltype=max(overall),delta = alltype-overall) %>%
  filter(delta!=0,delta>0.25)
rank_consensus %>%
  filter(Block=='RNA') %>%
  select(overall,DeconvTool,NCellTypes) %>%
  group_by(DeconvTool) %>%
  mutate(alltype=max(overall),delta = alltype-overall) %>%
  filter(delta!=0,delta<0.1)

