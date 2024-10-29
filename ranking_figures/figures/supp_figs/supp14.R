source("../../src/1_generate_score.R")
source("../../src/1_generate_time.R")
source("../../src/2_ranking_ranking_procedure_functions.R")

## ----
## Set parameters
## ----
library(ggplot2)
library(see)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
date = "231027"
n_candidates=2
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

## ----
## Prep for consensus deconv scores
## ----
source("../../src/load_scores_SB.R")
# res from fig2
res_fig2 = readRDS("../fig/fig2/df_res.rds")
deconv_methods = do.call(rbind,lapply(res_fig2, function(x)
  x %>% ungroup() %>% arrange(desc(overall)) %>%
    filter(!duplicated(candidate)) %>%
    top_n(n_candidates,overall) %>%
    mutate(FS=sapply(candidate, function(y) strsplit(y,'-')[[1]][3]),
           Block=sapply(candidate, function(y) strsplit(y,'-')[[1]][1]),
           Supervised=ifelse(DeconvTool%in%c(meth_met_sup,meth_rna_sup),'sup','unsup')) %>%
    select(DeconvTool,FS,Block,Supervised)))

## ----
## Consensus deconv scores and time
## ----
input_path <- "../0simu/simulations/"
deconv_path <- "../1SB/deconv/"
timing_path <- "../1SB/timing/"
date <- "231027"
score_methods <- c("rmse", "mae", "pearson")
if (!file.exists(paste0(folder,"/scores.rds"))) {
  scores = generate_score_SB_mean(input_path, deconv_path, date, score_methods, deconv_methods)
  saveRDS(scores, paste0(folder,"/scores.rds"))
} else {scores = readRDS(paste0(folder,"/scores.rds"))}
if (!file.exists(paste0(folder,"/time.rds"))) {
  time = generate_time_SB_mean(input_path, timing_path, date, deconv_methods)
  saveRDS(time, paste0(folder,"/time.rds"))
} else {time = readRDS(paste0(folder,"/time.rds"))}

## ----
## Rank
## ----
scores = scores %>%
  mutate(candidate=paste(block,candidate,class,sep="-"))
time = time %>%
  mutate(candidate=paste(block,candidate,class,sep="-"))
ranks = ranking_consensus(scores, time)
categories = c("raw_perf","stab_perf","time")
df_score_cat = lapply(categories, function(x) {
  df = ranking_step1(scores, time) %>%
    coerce_pearson() %>%
    filter(cat_score==x) %>%
    filter(!is.na(trendval)) %>%
    group_by(name_score,dataset) %>%
    mutate(archetype_best=max(trendval, na.rm=T),
           archetype_worst=min(trendval, na.rm=T),
           rank=rank(trendval, ties.method = "average", na.last="keep")) %>%
    ungroup() %>%
    mutate(rank=rank/max(rank, na.rm=T)) %>%
    group_by(candidate,dataset,name_score) %>%
    mutate(d_best=euclidean_vector_dist(trendval,archetype_best),
           d_worst=euclidean_vector_dist(trendval,archetype_worst),
           topsis=d_worst/(d_worst+d_best))
  if (x=='time') {df = df %>% filter(name_score=='time time')}
  df %>%
    group_by(candidate,dataset) %>%
    summarise(score_inter1_raw=geomMean(trendval),
           score_inter1_topsis=geomMean(topsis),
           score_inter1_rank=geomMean(rank),
           !!x:=mean(score_inter1_raw,score_inter1_topsis,score_inter1_rank, trim=0)) %>%
    mutate(DeconvTool = sapply(candidate, function(x) strsplit(x,'-')[[1]][2]),
           Block = sapply(candidate, function(x) strsplit(x,'-')[[1]][1]),
           sup = sapply(candidate, function(x) strsplit(x,'-')[[1]][3])) %>%
    select(dataset,candidate,x,DeconvTool,Block,sup)
})
names(df_score_cat) = categories

## ----
## Plot each setting with a delta barplot (overall score)
## ----
df_delta = ranks %>%
  mutate(Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         Block = ifelse(Block=="met","DNAm","RNA"),
         Class=sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
         Class = ifelse(Class=="sup","Supervised","Unsupervised"),
         Setting=paste(Block,Class,sep="\n"),
         DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2])) %>%
  group_by(Setting) %>%
  mutate(Consensus=overall[grep("Consensus",DeconvTool)],
         Delta=Consensus-overall) %>%
  filter(Delta!=0)
df_delta$DeconvTool = factor(df_delta$DeconvTool, levels=unique(unlist(lapply(res_fig2, function(x)
  x %>% filter(!duplicated(candidate)) %>% arrange(desc(overall)) %>%
    ungroup() %>% top_n(n_candidates,overall) %>% pull(DeconvTool)))))
df_delta$Setting = factor(df_delta$Setting, levels=c("DNAm\nSupervised","DNAm\nUnsupervised",
                                                     "RNA\nSupervised","RNA\nUnsupervised"))

ggplot(df_delta, aes(x=Setting, y=Delta, fill=DeconvTool)) +
  geom_col(position = 'dodge', color='black') +
  scale_fill_manual(values=c(custom_palette[c(1,7)],palette_metro()(5)[1:2],
                              custom_palette[12],palette_metro()(5)[3:4])) +
  ylab("\u0394 (Consensus-DeconvTool)") +
  geom_hline(yintercept = 0, col='red', linetype='dashed') +
  xlab("") +
  ylim(c(-.26,.26)) +
  theme_modern()
ggsave(paste0(folder,"/delta_barplot.pdf"), width=6.8, height=3.5, device = cairo_pdf)

## ----
## Plot each setting with a delta barplot (raw perf only)
## ----
df_delta_cat = lapply(names(df_score_cat), function(y) {
  df = df_score_cat[[y]] %>%
    mutate(Block = ifelse(Block=="met","DNAm","RNA"),
           Class = ifelse(sup=="sup","Supervised","Unsupervised"),
           Setting=paste(Block,Class,sep="\n")) %>%
    group_by(Setting,dataset) %>%
    mutate(n=n()) %>% ungroup() %>% filter(n==max(n)) %>% group_by(Setting,dataset) %>%
    mutate(Consensus=get(y)[grep("Consensus",DeconvTool)],
           Delta=Consensus-get(y)) %>%
    filter(Delta!=0)
  df$DeconvTool = factor(df$DeconvTool, levels=unique(unlist(lapply(res_fig2, function(x)
    x %>% filter(!duplicated(candidate)) %>% arrange(desc(overall)) %>%
      ungroup() %>% top_n(n_candidates,overall) %>% pull(DeconvTool)))))
  df$Setting = factor(df$Setting, levels=c("DNAm\nSupervised","DNAm\nUnsupervised",
                                           "RNA\nSupervised","RNA\nUnsupervised"))
  df})
names(df_delta_cat) = categories

for (i in categories) {
  ggplot(df_delta_cat[[i]], aes(x=Setting, y=Delta, fill=DeconvTool)) +
    geom_boxplot() +
    geom_point(color='grey20',position=position_dodge(width=0.75),aes(group=DeconvTool)) +
    geom_hline(yintercept = 0, col='red', linetype='dashed') +
    scale_fill_manual(values=c(custom_palette[c(1,7)],palette_metro()(5)[1:2],
                               custom_palette[12],palette_metro()(5)[3:4])) +
    ylab("\u0394 (Consensus-DeconvTool)") +
    xlab("") +
    ylim(c(-.75,.75)) +
    guides(color='none') +
    theme_modern()
  ggsave(paste0(folder,"/delta_barplot_",i,".pdf"), width=6.8, height=3.5, device = cairo_pdf)
}
