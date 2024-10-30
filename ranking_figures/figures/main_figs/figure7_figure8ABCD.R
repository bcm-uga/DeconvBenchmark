## ----
## Set parameters, put your own
## ----
score_path = "../../compute_metrics/scores/"
date = '241025'
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
source("../generic_functions/ranking_process.R")
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/load_scores_SB_invitro.R")
source("../generic_functions/load_scores_SB_invivo.R")

meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_dnam_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_dnam_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

scores_to_keep = c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean",
                   "time time")

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
library(fmsb)

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
## Load scores
## ----
res = load_data_vitro(score_path)
scores_vitro = res$scores
time_vitro = res$time
rm(res)
res = load_data_vivo(score_path)
scores_vivo = res$scores
time_vivo = res$time
rm(res)

## ----
## Scores inter and ranks
## ----
ranks_vitro = ranking_consensus(scores1=scores_vitro, scores2=time_vitro, scores_to_keep=scores_to_keep)
ranks_vivo = ranking_consensus(scores1=scores_vivo, scores2=time_vivo, scores_to_keep=scores_to_keep)
res_fig2 = readRDS("figure2_CD_figure3_CD/df_res.rds")
ranks = lapply(res_fig2, function(x)
  rbind(ranks_vitro %>% filter(candidate %in% unique(x$candidate)) %>%
    mutate(DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
           Source="In vitro"),
    ranks_vivo %>% filter(candidate %in% unique(x$candidate)) %>%
      mutate(DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
             Source="In vivo"),
    x %>% select(candidate,overall,DeconvTool) %>%
      filter(!duplicated(candidate)) %>%
      mutate(Source="In silico")))

# intermediate scoring
score_inter_vitro <- ranking_step1(scores1=scores_vitro, scores2=time_vitro) %>%
  filter(name_score %in% scores_to_keep) %>%
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
         aggregated=mean(score_inter2_raw,score_inter2_topsis,score_inter2_rank, trim=0)) %>%
  ungroup() %>%
  filter(!duplicated(paste(candidate,dataset))) %>%
  select(candidate,dataset,aggregated)
score_inter_vitro = left_join(score_inter_vitro, ranks_vitro, by='candidate')

score_inter_vivo <- ranking_step1(scores1=scores_vivo, scores2=time_vivo) %>%
  filter(name_score %in% scores_to_keep) %>%
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
         aggregated=mean(score_inter2_raw,score_inter2_topsis,score_inter2_rank, trim=0)) %>%
  ungroup() %>%
  filter(!duplicated(paste(candidate,dataset))) %>%
  select(candidate,dataset,aggregated)
score_inter_vivo = left_join(score_inter_vivo, ranks_vivo, by='candidate')

score_inter = lapply(res_fig2, function(x)
  rbind(score_inter_vitro %>% filter(candidate %in% unique(x$candidate)) %>%
          mutate(DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
                 Source="In vitro",
                 dataset=sapply(dataset, function(x) strsplit(x,"_")[[1]][1])) %>%
          select(candidate,DeconvTool,dataset,overall,aggregated,Source),
        score_inter_vivo %>% filter(candidate %in% unique(x$candidate)) %>%
          mutate(DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
                 Source="In vivo",
                 dataset=sapply(dataset, function(x) strsplit(x,"_")[[1]][1])) %>%
          select(candidate,DeconvTool,dataset,overall,aggregated,Source),
        x %>% select(candidate,DeconvTool,dataset,overall,aggregated) %>%
          mutate(Source="In silico",
                 dataset=sapply(dataset, function(x) strsplit(x,"-")[[1]][1]))))
saveRDS(score_inter, paste0(folder,'/score_inter.rds'))


ranks = mapply(function(x,y) {
  x %>% mutate(Setting=paste0(sub("dnam","DNAm",
                                  sub("rna","RNA",
                                      sub("unSupervised","Unsupervised",
                                          sub("sup","Supervised",strsplit(y,'-')[[1]])))),collapse='\n'))
}, x=ranks,y=names(ranks),SIMPLIFY = F)
score_inter = mapply(function(x,y) {
  x %>% mutate(Setting=paste0(sub("dnam","DNAm",
                                  sub("rna","RNA",
                                      sub("unSupervised","Unsupervised",
                                          sub("sup","Supervised",strsplit(y,'-')[[1]])))),collapse='\n'))
}, x=score_inter,y=names(score_inter),SIMPLIFY = F)

## ----
## Keep only methods that ran on ALL datasets for BOTH in vitro/in vivo sources
## ----
ran_all = lapply(score_inter, function(x)
  x %>% arrange(DeconvTool) %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool)) %>% 
    filter(Source!='In silico') %>%
    group_by(DeconvTool,Source) %>%
    summarise(n_dataset=n()) %>% arrange(Source) %>%
    group_by(Source) %>%
    mutate(max = max(n_dataset)) %>%
    filter(n_dataset == max) %>%
    group_by(DeconvTool) %>%
    filter(duplicated(DeconvTool)) %>% pull(DeconvTool))
ranks_filter = lapply(seq_along(ranks), function(x) {
  df = ranks[[x]] %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool)) %>%
    filter(DeconvTool %in% ran_all[[x]]) %>%
    group_by(Source) %>%
    mutate(scaled_overall=(overall-min(overall))/(max(overall)-min(overall)))
  left_join(ranks[[x]] %>%
              mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                            DeconvTool=='rlr'~'RLR',
                                            DeconvTool=='elasticnet'~'Elastic net',
                                            DeconvTool=='fardeep'~'FARDEEP',
                                            DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                            DeconvTool=='ols'~'OLS',
                                            DeconvTool=='svr'~'SVR', .default = DeconvTool)),df)})

mean_overall = lapply(ranks_filter, function(x) {
  df1 = data.frame(x %>%
                    group_by(Source) %>%
                    mutate(Source_rank=rank(1-overall)) %>%
                    group_by(DeconvTool) %>%
                    mutate(mean_overall=mean(scaled_overall)) %>%
                    ungroup() %>%
                    mutate(Global_rank=rank(1-mean_overall, na.last='keep'),
                           Global_rank=ceiling(Global_rank/3)))
  df = inner_join(left_join(left_join(df1 %>% filter(Source=='In silico') %>% select(DeconvTool,Source_rank) %>% rename("In silico"=Source_rank),
                                        df1 %>% filter(Source=='In vitro') %>% select(DeconvTool,Source_rank) %>% rename("In vitro"=Source_rank)),
                             df1 %>% filter(Source=='In vivo') %>% select(DeconvTool,Source_rank) %>% rename("In vivo"=Source_rank)),
                  df1 %>% filter(Source=='In silico') %>% select(DeconvTool,Global_rank) %>% rename("Global"=Global_rank))
  rownames(df) = df$DeconvTool
  df  %>% select(Global,`In silico`,`In vitro`,`In vivo`) %>% arrange(Global)
})

## ----
## Plot rank tables for methods that ran on ALL datasets for BOTH in vitro/in vivo sources (figure 7)
## ----
max_show = sapply(mean_overall, function(x) max(x$`In silico`, na.rm=T))
plot_table = lapply(mean_overall, function(data)
  ggtexttable(data[,c(2:4,1)],
              theme = ttheme(base_size = 8, base_style = "default",
                             colnames.style=colnames_style(color = "grey15", fill = "white", size=10),
                             rownames.style=rownames_style(color = "grey15", fill = "white", size=10, face="bold"))))
for (k in seq_along(plot_table)) {
  for (i in seq(2,tab_nrow(plot_table[[k]]))) {
    for (j in seq(2,tab_ncol(plot_table[[k]]))) {
      plot_table[[k]] <- plot_table[[k]] %>%
        table_cell_bg(row = i, column = j,
                      fill = viridis::rocket(max_show[k], direction=-1)[mean_overall[[k]][,c(2:4,1)][i-1,j-1]])
      if (!is.na(mean_overall[[k]][,c(2:4,1)][i-1,j-1])) { if (mean_overall[[k]][,c(2:4,1)][i-1,j-1] >= round(3/4*nrow(mean_overall[[k]]))) {
        plot_table[[k]] <- plot_table[[k]] %>%
          table_cell_font(row = i, column = j, size=8,
                          color = 'grey90')
      }}
    }
  }
}
names(plot_table) = names(ranks)
heights = c(2,2,4,2)
for (p in seq_along(plot_table)) {
  print(plot_table[[p]])
  ggsave(paste0(folder,"/fig7_table_global_ranks_",names(plot_table)[p],".pdf"), height=heights[p], width=5)
}

## ----
## Plot scatterplot for methods that ran on ALL datasets for BOTH in vitro/in vivo sources (figure 7)
## Plot scatterplot for all methods (supp figure 9)
## ----
ranks_filter = lapply(seq_along(ranks_filter), function(x) {
  ranks_filter[[x]]$DeconvTool = factor(ranks_filter[[x]]$DeconvTool, levels=rownames(mean_overall[[x]]))
  ranks_filter[[x]]})
score_inter = lapply(seq_along(score_inter), function(x) {
  score_inter[[x]] = score_inter[[x]] %>% mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                DeconvTool=='rlr'~'RLR',
                                DeconvTool=='elasticnet'~'Elastic net',
                                DeconvTool=='fardeep'~'FARDEEP',
                                DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                DeconvTool=='ols'~'OLS',
                                DeconvTool=='svr'~'SVR', .default = DeconvTool))
  score_inter[[x]]$DeconvTool = factor(score_inter[[x]]$DeconvTool,
                                       levels=unique(c(rownames(mean_overall[[x]]),
                                                     score_inter[[x]]$DeconvTool)))
  score_inter[[x]]$Source = factor(score_inter[[x]]$Source, levels=unique(sort(score_inter[[x]]$Source)))
  score_inter[[x]]})

ggarrange(plotlist = mapply(function(x,xbis) {
  p = ggplot(x %>% filter(DeconvTool %in% ran_all[[xbis]]), aes(x=Source,y=overall)) +
    geom_point(aes(color=DeconvTool), size=4) +
    geom_line(aes(group=DeconvTool, color=DeconvTool), linetype='dotted') +
    scale_col_fun(xbis) +
    xlab("") +
    ylab(ifelse(xbis%in%c(2,4),"","Overall\n benchmark score")) +
    theme_modern(axis.title.size = 18,
                 legend.title.size = 18,
                 legend.text.size = 18)
  if (xbis==4) {p = p + scale_color_manual(values=palette_metro(reverse = T)(5)[1:2])}
  p},
  x=ranks_filter,xbis=seq_along(ranks_filter), SIMPLIFY = F))
ggsave(paste0(folder,"/fig7_scatterplot.pdf"), width=12,height=8)

jitterVal <- lapply(score_inter, function(x)
  seq(from=-.4,to=.4,length.out=length(levels(x$DeconvTool)))[as.integer(x$DeconvTool)])
ggarrange(plotlist = lapply(seq_along(score_inter), function(x)
  ggplot(score_inter[[x]], aes(x=as.numeric(Source) + jitterVal[[x]])) +
    geom_point(aes(y=aggregated, color=DeconvTool),
               shape=4, size=2) +
    geom_line(aes(y=aggregated,
                  group=dataset), size=.1) +
    geom_point(aes(y=overall,color=DeconvTool), shape=1, size=3) +
    scale_col_fun(x) +
    xlab("") +
    scale_x_continuous(breaks=unique(sort(as.numeric(score_inter[[x]]$Source))),
                     labels=levels(score_inter[[x]]$Source)) +
    ylab("Overall benchmark score") +
    theme_modern(axis.title.size = 18,
                 legend.title.size = 18,
                 legend.text.size = 18)))
ggsave("../supp_figs/supp9/scatterplot.pdf", width=12,height=8)

## ----
## Prep spiderplots for figure 8
## ----
scores_norm_silico = ranking_norm(load_data(date, score_path)$scores,load_data(date, score_path)$time)
scores_norm_vitro = ranking_norm(scores_vitro,time_vitro)
scores_norm_vivo = ranking_norm(scores_vivo,time_vivo)
scores_norm = rbind(scores_norm_silico,scores_norm_vitro,scores_norm_vivo)
rm(scores_norm_silico,scores_norm_vitro,scores_norm_vivo)
scores_norm = lapply(res_fig2, function(x)
  scores_norm %>% filter(candidate %in% x$candidate))

spiderlist1 = lapply(scores_norm, function(setting) {
  res = lapply(unique(setting$candidate), function(tool) {
    tmp = setting %>% filter(candidate==tool)
    df = plyr::join_all(lapply(unique(tmp$name_score), function(y)
      data.frame(tmp %>% filter(name_score==y) %>%
                   group_by(dataset) %>%
                   summarise(trendval=median(trendval, na.rm=T)) %>%
                   rename(!!y:=trendval) %>% 
                   ungroup() %>% arrange(dataset))),
      by=c('dataset'), type='full')
    rownames(df) = df$dataset
    rbind(rep(1,ncol(df)-1),
          rep(0,ncol(df)-1),
          df[c(grep("perf_g",colnames(df), value=T),
               grep("med_",colnames(df), value=T),
               grep("sd_",colnames(df), value=T),
               grep("time",colnames(df), value=T))])
  })
  names(res)=unique(setting$candidate)
  res})

dataset_list = list('dnam'=rownames(spiderlist1$`dnam-sup`[[1]]),
                    'rna'=rownames(spiderlist1$`rna-sup`[[1]]))
dataset_silico = lapply(dataset_list, function(x) grep("CL",x)-2)
dataset_vitro = lapply(dataset_list, function(x) grep("MIX",x)-2)
dataset_vivo = lapply(dataset_list, function(x) grep("REAL",x)-2)
colors_in = lapply(dataset_list, function(x) rep(NA,length(x)-2))
for (i in seq_along(colors_in)) {
  colors_in[[i]][dataset_silico[[i]]] = RColorBrewer::brewer.pal(name="Reds",n=length(dataset_silico[[i]])+1)[seq(from=2,length.out=length(dataset_silico[[i]]))]
  colors_in[[i]][dataset_vitro[[i]]] = RColorBrewer::brewer.pal(name="Blues",n=length(dataset_vitro[[i]]))[if (length(dataset_vitro[[i]])<3) {seq(from=2,length.out=length(dataset_vitro[[i]]))}]
  colors_in[[i]][dataset_vivo[[i]]] = RColorBrewer::brewer.pal(name="Greens",n=length(dataset_vivo[[i]]))[if (length(dataset_vivo[[i]])<3) {seq(from=2,length.out=length(dataset_vivo[[i]]))}]
}


## ----
## Plot spiderplots (figure 8 panels A,B,C,D)
## ----
dev.off()
lapply(seq_along(spiderlist1), function(setting)
  lapply(seq_along(spiderlist1[[setting]]), function(method) {
    block = strsplit(names(spiderlist1)[setting],'-')[[1]][1]
    data = spiderlist1[[setting]][[method]][,c(1,rev(seq(2,ncol(spiderlist1[[setting]][[method]]))))]
    colnames(data) = c('RMSE','Time','Pearson_s (sd)','Pearson_c (sd)', 'Pearson_s', 'Pearson_c','Pearson_m','MAE')
    legends = sapply(rownames(spiderlist1[[setting]][[method]][-c(1,2),]), function(x)
      strsplit(strsplit(x,"_")[[1]][1],"-")[[1]][1])
    pdf(file = paste0(folder,"/fig8/",names(spiderlist1)[setting],"/",names(spiderlist1[[setting]])[method],".pdf"),
        width = 8)
    par(mar = c(1, 1, 1, 1))
    radarchart(data,
               axistype=0, seg=5, pcol=colors_in[[block]] , plwd=8 , plty=1,
               cglcol="grey", cglty=1, axislabcol="grey",
               caxislabels=seq(0,20,5), cglwd=0.8,
               vlcex=0.8 )
    legend(x=1.1, y=1.1,
           legend = legends[c(grep("CL",legends),grep("MIX",legends),grep("REAL",legends))],
           bty = "n", pch=20 , col=colors_in[[block]][c(grep("CL",legends),grep("MIX",legends),grep("REAL",legends))], text.col = "black",
           cex=1, pt.cex=1, x.intersp=.8, y.intersp=.8)
    dev.off()
  }))
