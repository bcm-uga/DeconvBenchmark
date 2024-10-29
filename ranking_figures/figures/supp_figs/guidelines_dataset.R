## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
library(ggpubr)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")
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

## ----
## Plot function
## ----
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
## Load data (methods perf and charac)
## ----
charac_rna = readRDS("characteristics_dataset/charac_rna.rds")
charac_met = readRDS("characteristics_dataset/charac_met.rds")
scores_inter = readRDS("guidelines_dataset/scores_inter.rds")
scores_inter = scores_inter %>%
  mutate(Dataset = sapply(Dataset, function(x) strsplit(strsplit(x,"_")[[1]][1],"-")[[1]][1]),
         FS = sapply(candidate, function(x) strsplit(x,"-")[[1]][3]))
winners = readRDS("../fig/fig2/df_res.rds")

scores_inter$nFeatures = NA
scores_inter$nFeatures[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"nFeatures"])
scores_inter$nFeatures[scores_inter$Block=='met'] = 3e4
scores_inter$nFeatures[scores_inter$FS!='none'] = 1e3
scores_inter$nCellTypes = NA
scores_inter$nCellTypes[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"nCellTypes"])
scores_inter$nCellTypes[scores_inter$Block=='met'] = sapply(scores_inter$Dataset[scores_inter$Block=='met'], function(x) charac_met[x,"nCellTypes"])
scores_inter$Pearson = NA
scores_inter$Pearson[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Mean.R2"])
scores_inter$Pearson[scores_inter$Block=='met'] = sapply(scores_inter$Dataset[scores_inter$Block=='met'], function(x) charac_met[x,"Mean.R2"])
scores_inter$PhenoVol = NA
scores_inter$PhenoVol[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Phenotypic.volume"])
scores_inter$PhenoVol[scores_inter$Block=='met'] = sapply(scores_inter$Dataset[scores_inter$Block=='met'], function(x) charac_met[x,"Phenotypic.volume"])
scores_inter$Sparsity = NA
scores_inter$Sparsity[scores_inter$Block=='rna'] = sapply(scores_inter$Dataset[scores_inter$Block=='rna'], function(x) charac_rna[x,"Sparsity"])
scores_inter$Kurtosis = NA
scores_inter$Kurtosis[scores_inter$Block=='met'] = sapply(scores_inter$Dataset[scores_inter$Block=='met'], function(x) charac_met[x,"Kurtosis"])
scores_inter$Technology = NA
scores_inter$Technology[scores_inter$Block=='met'] = sapply(scores_inter$Dataset[scores_inter$Block=='met'], function(x) charac_met[x,"Technology"])

scores_inter = scores_inter %>% 
  mutate(HG=PhenoVol/Pearson) %>%
  ungroup()

## ----
## Explore : Link each charac with global perf (all meth) -> focus on HG
## ----
colnames(scores_inter)
# Cannot look at global perf for nFeatures because it varies from methods with FS to methods w/o FS
ggplot(scores_inter, aes(x=nCellTypes, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  theme_modern() # nothing
ggplot(scores_inter, aes(x=1/Pearson, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  theme_modern() # interesting
ggplot(scores_inter, aes(x=PhenoVol, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  theme_modern() # interesting
ggplot(scores_inter, aes(x=HG, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  theme_modern() # interesting
ggplot(scores_inter %>% filter(Block=="rna"), aes(x=Sparsity, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  theme_modern() # nothing per source
ggplot(scores_inter %>% filter(Block=="met"), aes(x=Kurtosis, y=dataset_ease_decon)) +
  geom_point(aes(color=Dataset, shape=Source, size=Source)) +
  geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
  facet_wrap(~Block, scales = "free_x") +
  theme_modern() # nothing

## ----
## Explore : Link each charac with perf for each method
## ----
scores_inter = scores_inter %>%
  mutate(Supervised = sapply(candidate, function(x) ifelse(x %in% c(winners$`met-sup`$candidate,
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
         Block = case_when(Block=='met'~'DNAm',
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
# Check if methods are sensitive to dataset (only in silico)
ggplot(scores_inter %>% filter(Source=='In silico'), aes(x=score_inter2, y=DeconvTool)) +
  geom_boxplot(aes(fill=Supervised)) +
  geom_point() +
  scale_fill_manual(values=palette_metro()(2)) +
  facet_wrap(~Block) +
  guides(fill='none') +
  theme_modern()
# Check if it varies a lot across sources /!\ I don't have the same nb of datasets...
ggplot(scores_inter, aes(x=score_inter2, y=DeconvTool)) +
  geom_boxplot(aes(fill=Source), alpha=.6) +
  geom_point(position=position_dodge(width=0.75),aes(group=Source, color=Source)) +
  facet_wrap(~Block, scales="free_y") +
  theme_modern()
# Check if it varies a lot across Pearson cat
ggplot(scores_inter %>%
         group_by(Block) %>%
         mutate(Pearson_groups=cut_interval(1/Pearson,n=2)) %>%
         ungroup() %>% filter(Source=='In silico'), 
       aes(x=score_inter2, y=DeconvTool)) +
  geom_boxplot(aes(fill=Pearson_groups), alpha=.6) +
  geom_point(position=position_dodge(width=0.75),aes(group=Pearson_groups, color=Pearson_groups)) +
  facet_wrap(~Block, scales="free_y") +
  theme_modern()
# Check if it varies a lot across Sparsity cat
ggplot(scores_inter %>%
         filter(Block=='RNA') %>%
         mutate(Sparsity_groups=cut_interval(1/Sparsity,n=2)) %>% filter(Source=='In silico'), 
       aes(x=score_inter2, y=DeconvTool)) +
  geom_boxplot(aes(fill=Sparsity_groups), alpha=.6) +
  geom_point(position=position_dodge(width=0.75),aes(group=Sparsity_groups, color=Sparsity_groups)) +
  theme_modern()
# Check if it varies a lot across Techno
ggplot(scores_inter %>%
         filter(Block=='DNAm',Source=='In silico'), 
       aes(x=score_inter2, y=DeconvTool)) +
  geom_boxplot(aes(fill=Technology), alpha=.6) +
  geom_point(position=position_dodge(width=0.75),aes(group=Technology, color=Technology)) +
  theme_modern()

ggplot(scores_inter %>% filter(Source=='In silico',Block=="RNA"), aes(x=1/Pearson, y=score_inter2)) +
  geom_point(aes(color=DeconvTool)) + #, shape=Source)) +
  geom_line(aes(group=paste(DeconvTool,Source), color=DeconvTool)) +
  facet_wrap(~Supervised, ncol=1) +
  theme_modern()
ggplot(scores_inter %>% filter(Source=='In silico') %>% filter(Block=="DNAm"), aes(x=1/Pearson, y=score_inter2)) +
  geom_point(aes(color=DeconvTool)) + #, shape=Source)) +
  geom_line(aes(group=paste(DeconvTool,Source), color=DeconvTool)) +
  facet_wrap(~Supervised, ncol=1) +
  theme_modern()
ggplot(scores_inter %>% filter(Source=='In silico') %>% filter(Block=="RNA"), aes(x=1/Sparsity, y=score_inter2)) +
  geom_point(aes(color=DeconvTool)) + #, shape=Source)) +
  geom_line(aes(group=paste(DeconvTool,Source), color=DeconvTool)) +
  facet_wrap(~Supervised, ncol=1) +
  theme_modern()

## ----
## Prep plot dfs
## ----
#scores_inter$dataset = factor(scores_inter$dataset, levels=c(sort(grep("CL",unique(scores_inter$dataset),value=T)),
#                                                             sort(grep("MIX",unique(scores_inter$dataset),value=T)),
#                                                             sort(grep("REAL",unique(scores_inter$dataset),value=T))))
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
## Plots global
## ----
ggarrange(plotlist=lapply(sort(unique(scores_inter$Source)), function(source)
  ggplot(scores_inter %>% filter(Source==source), aes(x=HG, y=dataset_ease_decon)) +
    geom_point(aes(color=Dataset, shape=Block), size=4) +
    scale_col_source(source) +
    geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
    xlab("Heterogeneity composite proxy") +
    ylab("Mean dataset performance") +
    theme_modern(base_size = 15)), nrow = 1)
ggsave(paste0(folder,'/hg_global.pdf'), width=15, height=4)

ggarrange(plotlist=lapply(sort(unique(scores_inter$Source)), function(source)
  ggplot(scores_inter %>% filter(Source==source, Block=='RNA'), aes(x=Sparsity, y=dataset_ease_decon)) +
    geom_point(aes(color=Dataset, shape=Block), size=4) +
    scale_col_source(source) +
    geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
    xlab("Sparsity") +
    ylab("Mean dataset performance") +
    ylim(c(0.4,0.58)) +
    theme_modern()), nrow = 1)

ggarrange(plotlist=lapply(sort(unique(scores_inter$Source)), function(source)
  ggplot(scores_inter %>% filter(Source==source, Block=='RNA'), aes(x=Composite_feature, y=dataset_ease_decon)) +
    geom_point(aes(color=Dataset)) +
    geom_smooth(method='lm', se=T, color='black', linewidth=.4) +
    xlab("HG / Sparsity") +
    #ylim(c(0.4,0.65)) +
    #xlim(c(0,.2)) +
    ylab("Mean dataset performance") +
    theme_modern()))

## ----
## Plots silico only for HG + sensitivity
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
ggsave(paste0(folder,'/allsources_sensitivity.pdf'), width=10, height=6)

ggplot(scores_inter_silico, aes(x=score_inter2, y=y.lab)) +
  geom_boxplot(aes(fill=Supervised)) +
  geom_point() +
  facet_wrap(~Block) +
  scale_fill_manual(values=c("grey40", "#E51400")) +
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
        strip.text = element_text(face = "bold")) +
  xlab('Renormalized aggregated benchmark score')
ggsave(paste0(folder,'/silico_sensitivity.pdf'), width=10, height=6)

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
ggsave(paste0(folder,'/hg_silico.pdf'), width=15, height=6)

gplot_fun(scores_inter_silico_rna, 'Composite_groups') +
  geom_rect(data=scores_inter_silico_rna %>% filter(y_rect_composite=='good_trend'),
            aes(xmin=xmin, xmax=xmax,
                ymin=y_lab_num-.4, ymax=y_lab_num+.4), 
            col='grey50', fill=NA, linewidth=.2)
#ggsave(paste0(folder,'/composite_silico.pdf'), width=8, height=8)

scores_inter_silico$Technology = factor(scores_inter_silico$Technology, levels=c('800k','450k','EPIC'))
gplot_fun(scores_inter_silico %>% filter(Block=="DNAm"), 'Technology') +
  scale_color_bluebrown() +
  scale_fill_bluebrown()
#ggsave(paste0(folder,'/techno_silico.pdf'), width=12, height=8)
 