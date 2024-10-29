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
## Load df and scores
## ----
source("../../src/load_scores_SB.R")
res = load_data(date)
scores = res$scores
time = res$time
rm(res)
source("../../src/2_ranking_ranking_procedure_functions.R")
df_fig2 = readRDS("../fig/fig2/df_res.rds")
scores_norm = ranking_norm(scores,time)
scores_step1 = ranking_step1(scores,time)
scores_norm = lapply(df_fig2, function(data)
  scores_norm %>% filter(candidate %in% unique(data$candidate)) %>%
    select(score,name_score,dataset,candidate,trendval))
scores_step1 = lapply(df_fig2, function(data)
  scores_step1 %>% filter(candidate %in% unique(data$candidate)) %>%
    select(name_score,dataset,candidate,trendval))

## ----
## Plot supp = individual metrics
## ----
plots = lapply(seq_along(scores_norm), function(idx)
  lapply(unique(scores_norm[[idx]]$name_score), function(score_type) {
    data = scores_norm[[idx]] %>% ungroup() %>% filter(name_score==score_type) %>%
      mutate(DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
             DeconvTool=factor(DeconvTool, levels=df_fig2[[idx]] %>%
                                 arrange(desc(overall)) %>%
                                 filter(!duplicated(DeconvTool)) %>% pull(DeconvTool)))
    p = ggplot(data, aes(x=DeconvTool, y=score, color=dataset)) +
      geom_boxplot() +
      scale_color_metro_d() +
      ylab(score_type) +
      xlab("") +
      theme(axis.text.x = element_blank()) +
      theme_modern(axis.text.angle = 30)
    if (score_type=='time time') {p = p + scale_y_log10()}
    p
  }))
for (i in seq_along(plots)) {
  ggarrange(plotlist = plots[[i]], common.legend = T)
  ggsave(paste0(folder,"/boxplot_",names(scores_norm)[i],".pdf"), width=16, height=15)
}


## ----
## Plot supp = spider graph
## ----
library(fmsb)
spiderlist1 = lapply(scores_step1, function(setting) {
  res = lapply(unique(setting$candidate), function(tool) {
    tmp = setting %>% filter(candidate==tool)
    df = plyr::join_all(lapply(unique(tmp$name_score), function(y) {
      tmp %>% filter(name_score==y) %>%
        rename(!!y:=trendval) %>% ungroup() %>%
        select(dataset,y)
    }), by=c('dataset'), type='full')
    rownames(df) = df$dataset
    rbind(rep(1,length(unique(tmp$name_score))),
          rep(0,length(unique(tmp$name_score))),
          df[c(grep("perf_g",colnames(df), value=T),
               grep("med_",colnames(df), value=T),
               grep("sd_",colnames(df), value=T),
               grep("time",colnames(df), value=T))])
  })
  names(res)=unique(setting$candidate)
  res})
spiderlist2 = lapply(seq_along(scores_step1), function(setting) {
  res = lapply(unique(scores_step1[[setting]]$name_score), function(metric) {
    tmp = scores_step1[[setting]] %>% filter(name_score==metric)
    df = plyr::join_all(lapply(unique(tmp$candidate), function(y) {
      tmp %>% filter(candidate==y) %>%
        rename(!!y:=trendval) %>% ungroup() %>%
        select(dataset,y)
    }), by=c('dataset'), type='full')
    rownames(df) = df$dataset
    rbind(rep(1,length(unique(tmp$candidate))),
          rep(0,length(unique(tmp$candidate))),
          df[,-1][,df_fig2[[setting]] %>% filter(!duplicated(DeconvTool)) %>% arrange(desc(overall)) %>% pull(candidate)])
  })
  names(res)=unique(scores_step1[[setting]]$name_score)
  res})

colors_in=palette_metro()(5)
dev.off()
lapply(seq_along(spiderlist1), function(setting)
  lapply(seq_along(spiderlist1[[setting]]), function(method) {
    png(filename = paste0(folder,"/spider_method/",names(spiderlist1)[setting],"/",names(spiderlist1[[setting]])[method],".png"),
        width = 580)
    par(mar = c(2, 3, 3, 2))
    radarchart(spiderlist1[[setting]][[method]][,c(1,rev(seq(2,ncol(spiderlist1[[setting]][[method]]))))],
               axistype=0, seg=5, pcol=colors_in , plwd=4 , plty=1,
               cglcol="grey", cglty=1, axislabcol="grey",
               caxislabels=seq(0,20,5), cglwd=0.8,
               vlcex=0.8 )
    legend(x=1.1, y=1.1,
           legend = rownames(spiderlist1[[setting]][[method]][-c(1,2),]),
           bty = "n", pch=20 , col=colors_in , text.col = "black",
           cex=1, pt.cex=1, x.intersp=.8, y.intersp=.8)
    dev.off()
    }))
lapply(seq_along(spiderlist2), function(setting)
  lapply(seq_along(spiderlist2[[setting]]), function(score) {
    png(filename = paste0(folder,"/spider_metric/",names(df_fig2)[setting],"/",names(spiderlist2[[setting]])[score],".png"),
        width = 580)
    par(mar = c(2, 3, 3, 2))
    radarchart(spiderlist2[[setting]][[score]],
               axistype=0, seg=5, pcol=colors_in , plwd=2 , plty=1,
               cglcol="grey", cglty=1, axislabcol="grey",
               caxislabels=seq(0,20,5), cglwd=0.8,
               vlcex=0.8 )
    legend(x=1.1, y=1.1,
           legend = rownames(spiderlist2[[setting]][[score]][-c(1,2),]),
           bty = "n", pch=20 , col=colors_in , text.col = "black",
           cex=1, pt.cex=1, x.intersp=.8, y.intersp=.8)
    dev.off()
  }))
