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

## ----
## Load scores
## ----
source("../../src/load_scores_SB.R")
res = load_data(date)
scores = res$scores
time = res$time
rm(res)

## ----
## Ranks and per dataset
## ----
source("../../src/2_ranking_ranking_procedure_functions.R")
ranks = ranking_consensus(scores1=scores, scores2=time)

# dataset scoring
score_datasets <- coerce_pearson(ranking_step1(scores1=scores, scores2=time)) %>%
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
         topsis=d_worst/(d_worst+d_best)) %>%
  group_by(candidate,dataset,cat_score) %>%
  summarise(score_inter1_raw=geomMean(trendval),
            score_inter1_topsis=geomMean(topsis),
            score_inter1_rank=geomMean(rank),
            score_inter1=mean(score_inter1_raw,score_inter1_topsis,score_inter1_rank, trim=0)) %>%
  group_by(candidate,dataset) %>%
  mutate(weights_cat=convert_to_weights(weights_dic_values,cat_score),
         aggregated_raw=weighgeomMean(score_inter1_raw, weights_cat),
         aggregated_topsis=weighgeomMean(score_inter1_topsis, weights_cat),
         aggregated_rank=weighgeomMean(score_inter1_rank, weights_cat),
         aggregated=mean(aggregated_raw,aggregated_topsis,aggregated_rank, trim=0)) %>%
  select(candidate,dataset,aggregated) %>%
  filter(!duplicated(candidate,dataset)) %>%
  mutate(Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         Setting=paste(Block,dataset))

## ----
## Retrieve best FS per DeconvTool
## ----
best_FS_per_tool = ranks %>%
  mutate(Setting=sapply(candidate, function(x)
    paste0(head(strsplit(x,"-")[[1]],2),collapse="-")),
    FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
    Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1])) %>%
  group_by(Setting) %>%
  filter(overall==max(overall)) %>% ungroup() %>%
  arrange(desc(Block),desc(overall)) %>%
  select(Setting,FS,candidate,Block)

best_FS_per_tool_datasets = score_datasets %>%
  mutate(Setting=sapply(candidate, function(x)
                    paste0(head(strsplit(x,"-")[[1]],2),collapse="-")),
         FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3])) %>%
  ungroup() %>% group_by(Setting,dataset) %>%
  filter(aggregated==max(aggregated)) %>% ungroup() %>%
  select(Setting,FS,dataset,candidate)

levelz = unique(sapply(best_FS_per_tool$Setting, function(x) strsplit(x,"-")[[1]][2]))
levelz = c(levelz[seq(which(levelz=='NMF')-1)],
           "MeDeCom",
           levelz[setdiff(seq(which(levelz=='NMF'),length(levelz)),which(levelz=="MeDeCom"))])

## ----
## Plot barplot
## ----
plotting_function = function(df, y, ylab) {
  ggplot(df, aes(x=x.label, y=get(y), fill=FS)) +
    facet_grid(~facet, labeller = labeller(facet=facet.labs), scales = "free_x", space = "free_x") +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_viridis_d() +
    ylab(ylab) +
    xlab("") +
    geom_text(aes(label=Best_FS, group=FS, y=get(y) + .005), col='red', position=position_dodge(.85)) +
    theme_classic(base_size = 11, base_family = "") +
    theme(plot.title = element_text(size = 15, face = "plain", margin = margin(0, 0, 20, 0)),
          plot.title.position = "plot", 
          legend.position = "right",
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 13),
          legend.key = element_blank(), 
          legend.spacing.x = unit(2, "pt"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.text.x = element_markdown(colour="red",
                                         angle = 45, 
                                         hjust = 1),
          axis.text = element_text(size = 12),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 13, face = "plain"),
          plot.tag = element_text(size = 15, face = "bold"), 
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

df <- ranks %>%
  mutate(Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
         FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3]),
         Supervised = ifelse(DeconvTool %in% c(meth_met_unsup,meth_rna_unsup),"No","Yes"),
         facet = paste(Block, Supervised)) %>%
  group_by(paste(Block,DeconvTool)) %>%
  arrange(desc(overall)) %>%
  ungroup()
df_datasets <- lapply(unique(score_datasets$Setting), function(dataXblock)
  score_datasets %>% filter(Setting==dataXblock) %>% group_by(candidate,Setting) %>%
    reframe(aggregated = unique(aggregated)) %>%
    mutate(Dataset = sapply(Setting,function(x) strsplit(x," ")[[1]][2]),
           Block = sapply(Setting,function(x) strsplit(x," ")[[1]][1]),
           FS = sapply(candidate,function(x) strsplit(x,"-")[[1]][3]),
           DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
           Supervised = ifelse(DeconvTool %in% c(meth_met_unsup,meth_rna_unsup),"No","Yes"),
           facet = paste(Block, Supervised),
           Setting = paste(Dataset,Block)) %>%
    select(candidate, aggregated, Dataset, Block, FS, DeconvTool, Supervised, Setting, facet) %>%
    group_by(DeconvTool) %>%
    arrange(desc(aggregated)) %>%
    ungroup())

df$FS[df$FS=="hvf"] = "Highly variable features"
df$FS[df$FS=="none"] = "No feature selection"
df$FS[df$FS=="toast"] = "TOAST"
df$FS = factor(df$FS, levels=c("No feature selection","Highly variable features","TOAST"))
df$DeconvTool = factor(df$DeconvTool, levels=levelz)
df_datasets = lapply(df_datasets, function(x) {
  x$FS[x$FS=="hvf"] = "Highly variable features"
  x$FS[x$FS=="none"] = "No feature selection"
  x$FS[x$FS=="toast"] = "TOAST"
  x$FS = factor(x$FS, levels=c("No feature selection","Highly variable features","TOAST"))
  x$DeconvTool = factor(x$DeconvTool, levels=levelz)
  x})

df = df %>% ungroup() %>%
  group_by(Block,DeconvTool) %>%
  mutate(Best_FS = sapply(candidate, function(x) ifelse(x %in% best_FS_per_tool$candidate,"*","")),
         x.label = paste("<span style = 'color: ",
                         ifelse(FS[Best_FS=="*"]=="No feature selection",
                                "#440154", ifelse(FS[Best_FS=="*"]=="Highly variable features",
                                                  "#21918c", "#d7c202")),
                         ";'>",
                         DeconvTool,
                         "</span>", sep = ""))
df_datasets = lapply(df_datasets, function(data)
  data %>% ungroup() %>%
    group_by(DeconvTool) %>%
    mutate(Best_FS = sapply(candidate, function(x) ifelse(x %in%
                                                            (best_FS_per_tool_datasets %>%
                                                               filter(dataset==unique(data$Dataset)) %>%
                                                               pull(candidate)),"*","")),
           x.label = paste("<span style = 'color: ",
                           ifelse(FS[Best_FS=="*"]=="No feature selection",
                                  "#440154", ifelse(FS[Best_FS=="*"]=="Highly variable features",
                                                    "#21918c", "#d7c202")),
                           ";'>",
                           DeconvTool,
                           "</span>", sep = "")))
facet.labs <- c("DNAm \n Unsupervised",
                "DNAm \n Supervised",
                "RNA \n Unsupervised",
                "RNA \n Supervised")
names(facet.labs) <- unique(sort(df$facet))

df$x.label = factor(df$x.label,
                    levels=c(sapply(levels(df$DeconvTool), function(x)
                      sapply(c("#440154","#21918c","#d7c202"), function(y)
                        paste0("<span style = 'color: ",y,";'>",x,"</span>")))))
df$facet = factor(df$facet, levels=c('met Yes','met No','rna Yes','rna No'))

plotting_function(df, 'overall', "Overall benchmark score")
ggsave(paste0(folder,"/all.pdf"), width = 12, height = 3.5)
plotting_function(df %>% filter(Block=='rna'), 'overall', "Overall benchmark score")
ggsave(paste0(folder,"/rna.pdf"), width = 10, height = 3.5)
plotting_function(df %>% filter(Block=='met'), 'overall', "Overall benchmark score")
ggsave(paste0(folder,"/met.pdf"), width = 8, height = 3.7)

df_datasets_rna = do.call(rbind,df_datasets) %>% filter(Block=='rna')
df_datasets_met = do.call(rbind,df_datasets) %>% filter(Block=='met')
df_datasets_rna$x.label = factor(df_datasets_rna$x.label,
                                 levels=c(sapply(levels(df_datasets_rna$DeconvTool), function(x)
                                   sapply(c("#440154","#21918c","#d7c202"), function(y)
                                     paste0("<span style = 'color: ",y,";'>",x,"</span>")))))
df_datasets_met$x.label = factor(df_datasets_met$x.label,
                                 levels=c(sapply(levels(df_datasets_met$DeconvTool), function(x)
                                   sapply(c("#440154","#21918c","#d7c202"), function(y)
                                     paste0("<span style = 'color: ",y,";'>",x,"</span>")))))
df_datasets_rna$facet = factor(df_datasets_rna$facet, levels=c('met Yes','met No','rna Yes','rna No'))
df_datasets_met$facet = factor(df_datasets_met$facet, levels=c('met Yes','met No','rna Yes','rna No'))

plot_datasets_rna = list()
for (x in seq_along(unique(df_datasets_rna$Dataset))) {
  if (x%in%c(1,4)) {
    plot_datasets_rna[[unique(df_datasets_rna$Dataset)[x]]] = plotting_function(df_datasets_rna %>% filter(Dataset==unique(df_datasets_rna$Dataset)[x]),
                                                                                  'aggregated',
                                                                                  "Aggregated benchmark score")
  } else {plot_datasets_rna[[unique(df_datasets_rna$Dataset)[x]]] = plotting_function(df_datasets_rna %>% filter(Dataset==unique(df_datasets_rna$Dataset)[x]),
                                                                                        'aggregated',
                                                                                        "")}
}
plot_datasets_rna[["Overall"]] = plotting_function(df %>% filter(Block=='rna'), 'overall', "Overall benchmark score")

plot_datasets_met = list()
for (x in seq_along(unique(df_datasets_met$Dataset))) {
  if (x%in%c(1,4)) {
    plot_datasets_met[[unique(df_datasets_met$Dataset)[x]]] = plotting_function(df_datasets_met %>% filter(Dataset==unique(df_datasets_met$Dataset)[x]),
                                                                                  'aggregated',
                                                                                  "Aggregated benchmark score")
  } else {plot_datasets_met[[unique(df_datasets_met$Dataset)[x]]] = plotting_function(df_datasets_met %>% filter(Dataset==unique(df_datasets_met$Dataset)[x]),
                                                                                        'aggregated',
                                                                                        "")}
}
plot_datasets_met[["Overall"]] = plotting_function(df %>% filter(Block=='met'), 'overall', "Overall benchmark score")

ggarrange(plotlist = plot_datasets_rna,
          common.legend = T,
          labels = c(sapply(unique(df_datasets_rna$Dataset), function(x) strsplit(x,"-")[[1]][1]),"Overall"))
ggsave(paste0(folder,"/datasets_rna.pdf"), width = 13, height = 8)
ggarrange(plotlist = plot_datasets_met,
          common.legend = T,
          labels = c(sapply(unique(df_datasets_met$Dataset), function(x) strsplit(x,"-")[[1]][1]),"Overall"))
ggsave(paste0(folder,"/datasets_met.pdf"), width = 13, height = 8)

## ----
## Plot categories
## ----
settings = list(meth_met_sup,meth_met_unsup,
                meth_rna_sup,meth_rna_unsup)
names(settings) = rep(c('met','rna'),each=2)
df_cat = lapply(seq_along(settings), function(x)
  ranking_step1(scores1=scores, scores2=time) %>% ungroup() %>%
    coerce_pearson() %>%
    filter(!is.na(trendval)) %>%
    arrange(candidate,dataset,name_score) %>%
    mutate(group_fs=sapply(candidate, function(y) paste(strsplit(y,'-')[[1]][-3],collapse='-'))) %>%
    group_by(group_fs,dataset,name_score) %>%
    mutate(trendval_none = trendval[grep("none",candidate)],
           trendval_diff = trendval-trendval_none,
           DeconvTool=sapply(candidate,function(y) strsplit(y,"-")[[1]][2]),
           FS=sapply(candidate,function(y) strsplit(y,"-")[[1]][3]),
           Block=sapply(candidate,function(y) strsplit(y,"-")[[1]][1])) %>% ungroup() %>%
    select(name_score,dataset,candidate,trendval_diff,DeconvTool,FS,Block) %>%
    filter(DeconvTool %in% settings[[x]],
           Block %in% names(settings)[x],
           FS!='none') %>%
    filter(candidate %in% best_FS_per_tool$candidate[-grep("none",best_FS_per_tool$candidate)])# %>%
    #group_by(candidate,name_score,DeconvTool,FS) %>%
    #summarise(trendval_diff_mean = mean(trendval_diff, na.rm = T))
  )
names(df_cat) = c("met-sup","met-unsup","rna-sup","rna-unsup")
df_cat = lapply(df_cat, function(x) {
  x$FS = factor(gsub("hvf","Highly variable features",
              gsub("toast","TOAST",x$FS)),levels=c("TOAST","Highly variable features"))
  x$name_score = factor(x$name_score, levels=c(c("mae perf_g", "rmse perf_g", "pearson perf_mean"),
                                               c("mae sd_g", "rmse sd_g", "pearson sd_mean"),
                                               c("time sd","time time")))
  x
})

heights = c(5,5,5,5)
for (i in seq_along(df_cat)) {
  ggplot(df_cat[[i]], aes(x=name_score, y=trendval_diff)) +
    geom_boxplot(position='dodge', width=.5, alpha=.7, aes(color=name_score, fill=FS)) +
    geom_line(aes(group=dataset), alpha=.1) +
    geom_hline(yintercept = 0, col='red', linetype='dotted', linewidth=.3) +
    scale_fill_manual(values = rev(viridis::viridis(3)[-1])) +
    scale_color_manual(values = c(rep("black",length(unique(df_cat[[i]]$name_score))-1),'orangered2')) +
    facet_wrap(~DeconvTool) +
    ylim(c(-1,1)) +
    ylab("\u0394 (FS-noFS)") +
    xlab("") +
    guides(color='none') +
    theme_modern(axis.text.angle = 80)
  ggsave(paste0(folder,"/mean_atomic_scores_",names(df_cat)[i],".pdf"), width=8, height=heights[i], device = cairo_pdf)
}

