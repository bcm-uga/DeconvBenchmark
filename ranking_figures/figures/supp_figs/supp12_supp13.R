## ----
## Set parameters, put your own
## Run script supp12_supp13_scores.R first
## ----
folder = "supp12_supp13"

## ----
## Load libraries
## ----
library(dplyr)
library(ggplot2)
library(see)
library(ggpubr)

## ----
## Functions
## ----
scale_col_source = function(source) {
  if (source=="In silico") {
    scale_color_manual(values=RColorBrewer::brewer.pal(name="Reds",n=7)[2:7])
  }
  else if (source=="In vitro") {
    scale_color_manual(values=RColorBrewer::brewer.pal(name="Blues",n=4)[2:4])
  }
  else if (source=="In vivo") {
    scale_color_manual(values=RColorBrewer::brewer.pal(name="Greens",n=5)[2:5])
  }
}

gplot_fun = function(data, fill_group) {
  ggplot(data, 
         aes(x=score_inter2, y=y.lab)) +
    geom_boxplot(aes(fill=.data[[fill_group]]), alpha=.6) +
    geom_point(position=position_dodge(width=0.75),size=2,aes(group=.data[[fill_group]], color=.data[[fill_group]])) +
    xlab('Renormalized aggregated benchmark score') +
    ylab("DeconvTool") +
    theme_classic(base_size = 11, base_family = "") +
    theme(plot.title = element_text(size = 15, face = "plain", margin = margin(0, 0, 20, 0)),
          plot.title.position = "plot", 
          legend.position = "right",
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 13),
          legend.key = element_blank(), 
          legend.spacing.x = unit(2, "pt"),
          axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.text.y = ggtext::element_markdown(colour="red",
                                                 hjust = 1),
          axis.text = element_text(size = 12),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 13, face = "plain"),
          plot.tag = element_text(size = 15, face = "bold"), 
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

## ----
## Load data
## ----
charac_rna = readRDS("supp11/charac_rna.rds")
charac_dnam = readRDS("supp11/charac_dnam.rds")
scores_inter = readRDS("supp12_supp13/scores_inter.rds")
scores_inter = scores_inter %>%
  mutate(Dataset = sapply(Dataset, function(x) strsplit(strsplit(x,"_")[[1]][1],"-")[[1]][1]),
         FS = sapply(candidate, function(x) strsplit(x,"-")[[1]][3]))
winners = readRDS("../main_figs/figure2CD_figure3CD/df_res.rds")

scores_inter$nFeatures = NA
scores_inter$nFeatures[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"nFeatures"])
scores_inter$nFeatures[scores_inter$Block=='dnam'] = 3e4
scores_inter$nFeatures[scores_inter$FS!='none'] = 1e3
scores_inter$nCellTypes = NA
scores_inter$nCellTypes[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"nCellTypes"])
scores_inter$nCellTypes[scores_inter$Block=='dnam'] = sapply(scores_inter$Dataset[scores_inter$Block=='dnam'], function(x) charac_dnam[x,"nCellTypes"])
scores_inter$Pearson = NA
scores_inter$Pearson[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Mean.R2"])
scores_inter$Pearson[scores_inter$Block=='dnam'] = sapply(scores_inter$Dataset[scores_inter$Block=='dnam'], function(x) charac_dnam[x,"Mean.R2"])
scores_inter$PhenoVol = NA
scores_inter$PhenoVol[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Phenotypic.volume"])
scores_inter$PhenoVol[scores_inter$Block=='dnam'] = sapply(scores_inter$Dataset[scores_inter$Block=='dnam'], function(x) charac_dnam[x,"Phenotypic.volume"])
scores_inter$Sparsity = NA
scores_inter$Sparsity[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Sparsity"])
scores_inter$Kurtosis = NA
scores_inter$Kurtosis[scores_inter$Block=='dnam'] = sapply(scores_inter$Dataset[scores_inter$Block=='dnam'], function(x) charac_dnam[x,"Kurtosis"])
scores_inter$Technology = NA
scores_inter$Technology[scores_inter$Block=='dnam'] = sapply(scores_inter$Dataset[scores_inter$Block=='dnam'], function(x) charac_dnam[x,"Technology"])

scores_inter = scores_inter %>% 
  mutate(HG=PhenoVol/Pearson) %>%
  ungroup()

scores_inter = scores_inter %>%
  mutate(Supervised = sapply(candidate, function(x) ifelse(x %in% c(winners$`dnam-sup`$candidate,
                                                                    winners$`rna-sup`$candidate),
                                                           "Supervised","Unsupervised")),
         DeconvTool = sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
         DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                DeconvTool=='rlr'~'RLR',
                                DeconvTool=='elasticnet'~'Elastic net',
                                DeconvTool=='fardeep'~'FARDEEP',
                                DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                DeconvTool=='ols'~'OLS',
                                DeconvTool=='svr'~'SVR', .default = DeconvTool),
         Block = case_when(Block=='dnam'~'DNAm',
                           Block=='rna'~'RNA'))
scores_inter$DeconvTool = factor(scores_inter$DeconvTool,
                                 levels=rev(do.call(rbind,winners) %>%
                                   filter(!duplicated(candidate)) %>%
                                   arrange(desc(overall)) %>%
                                   mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                                                 DeconvTool=='rlr'~'RLR',
                                                                 DeconvTool=='elasticnet'~'Elastic net',
                                                                 DeconvTool=='fardeep'~'FARDEEP',
                                                                 DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                                                 DeconvTool=='ols'~'OLS',
                                                                 DeconvTool=='svr'~'SVR', .default = DeconvTool)) %>%
                                   pull(DeconvTool) %>% unique()))

## ----
## Prep plot dfs
## ----
levels_deconv = levels(scores_inter$DeconvTool)
scores_inter = scores_inter %>% ungroup() %>%
  mutate(y.lab = sapply(seq_along(DeconvTool), function(idx)
    paste("<span style = 'color: ",
          ifelse(Supervised[idx]=="Supervised",
                 "grey40", "#E51400"),
          ";'>",
          DeconvTool[idx],
          "</span>", sep = "")),
    y.lab = factor(y.lab,
                   levels=c(sapply(levels_deconv, function(x)
                     sapply(c("grey40", "#E51400"), function(y)
                       paste0("<span style = 'color: ",y,";'>",x,"</span>")))))) %>%
  mutate(Composite_feature = HG/Sparsity)

scores_inter_silico = scores_inter %>%
  group_by(Block) %>%
  mutate(HG_groups=cut_interval(HG,n=2,labels=c("low","high"))) %>%
  filter(Source == 'In silico') %>%
  group_by(Block,DeconvTool,HG_groups) %>%
  mutate(median_perf_hg = median(score_inter2)) %>%
  group_by(Block,DeconvTool) %>%
  mutate(xmin=0.02,
         xmax=0.85,
         y_rect_pearson=ifelse(unique(median_perf_hg[HG_groups=="high"]) > unique(median_perf_hg[HG_groups=="low"]),
                               'good_trend','bad_trend')) %>%
  group_by(Block) %>%
  mutate(y_lab_num = as.numeric(factor(y.lab)))

scores_inter_silico_rna = scores_inter %>% filter(Block=='RNA') %>% ungroup() %>%
  mutate(Sparsity_groups=cut_interval(Sparsity,n=2,labels=c("low","high")),
         HG_groups=cut_interval(HG,n=2,labels=c("low","high")),
         Composite_groups=cut_interval(Composite_feature,n=2,labels=c("low","high"))) %>%
  filter(Source == 'In silico') %>%
  group_by(Block,DeconvTool,Sparsity_groups) %>%
  mutate(median_perf_sparsity = median(score_inter2)) %>%
  group_by(Block,DeconvTool,HG_groups) %>%
  mutate(median_perf_hg = median(score_inter2)) %>%
  group_by(Block,DeconvTool,Composite_groups) %>%
  mutate(median_perf_composite = median(score_inter2)) %>%
  group_by(Block,DeconvTool) %>%
  mutate(xmin=0.02,
         xmax=0.85,
         y_rect_sparsity=ifelse(unique(median_perf_sparsity[Sparsity_groups=="high"]) > unique(median_perf_sparsity[Sparsity_groups=="low"]),
                                'bad_trend','good_trend'),
         y_rect_hg=ifelse(unique(median_perf_hg[HG_groups=="high"]) > unique(median_perf_hg[HG_groups=="low"]),
                                'good_trend','bad_trend'),
         y_rect_composite=ifelse(unique(median_perf_composite[Composite_groups=="high"]) > unique(median_perf_composite[Composite_groups=="low"]),
                                 'good_trend','bad_trend')) %>%
  ungroup() %>%
  mutate(y_lab_num = as.numeric(factor(y.lab)))

## ----
## Plots supp13 panelA
## ----
ggarrange(plotlist=lapply(sort(unique(scores_inter$Source)), function(source)
  ggplot(scores_inter %>% filter(Source==source), aes(x=HG, y=dataset_ease_decon)) +
    geom_point(aes(color=Dataset, shape=Block), size=4) +
    scale_col_source(source) +
    geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
    xlab("Heterogeneity composite proxy") +
    ylab("Mean dataset performance") +
    theme_modern(base_size = 15)), nrow = 1)
ggsave(paste0(folder,'/supp13_panelA.pdf'), width=15, height=4)

## ----
## Plots silico only for HG + sensitivity (supp12)
## ----
ggplot(scores_inter, aes(x=score_inter2, y=y.lab)) +
  geom_boxplot(aes(fill=Supervised), outliers = F) +
  geom_point(aes(shape=Source), size=3) +
  facet_wrap(~Block) +
  scale_fill_manual(values=c("grey40", "#E51400")) +
  scale_shape_manual(values=c(16,2,3)) +
  theme_classic(base_size = 11, base_family = "") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(0, 0, 20, 0)),
        plot.title.position = "plot", 
        legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13),
        legend.key = element_blank(), 
        legend.spacing.x = unit(2, "pt"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 12),
        axis.text.y = ggtext::element_markdown(colour="red",
                                               hjust = 1, size = 15),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 13, face = "plain"),
        plot.tag = element_text(size = 15, face = "bold"), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  xlab('Renormalized aggregated benchmark score')
ggsave(paste0(folder,'/supp12.pdf'), width=10, height=6)

gplot_fun(scores_inter_silico, 'HG_groups') +
  facet_wrap(~Block, scales = 'free_y') +
  scale_fill_manual(values=c("firebrick1","seagreen1")) +
  scale_color_manual(values=c("firebrick3","seagreen3")) #+
  #geom_rect(data=scores_inter_silico %>% filter(Block=="RNA",y_rect_pearson=='good_trend'),
  #          aes(xmin=xmin, xmax=xmax,
  #              ymin=y_lab_num-.4, ymax=y_lab_num+.4), 
  #          col='grey50', fill=NA, linewidth=.2) +
  #geom_rect(data=scores_inter_silico %>% filter(Block=="DNAm",y_rect_pearson=='good_trend'),
  #          aes(xmin=xmin, xmax=xmax,
  #              ymin=y_lab_num-.4, ymax=y_lab_num+.4), 
  #          col='grey50', fill=NA, linewidth=.2)
ggsave(paste0(folder,'/supp13_panelB.pdf'), width=15, height=6)
