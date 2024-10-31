## ----
## Set parameters, put your own
## ----
source("../generic_functions/ranking_process.R")
source("../generic_functions/ranking_criteria.R")
rank_process = c("Sconsensus","Sraw","Stopsis","Srank")
ranking_process=c(ranking_consensus,ranking_raw,ranking_topsis,ranking_avgrank)
ranking_process_end=c(ranking_consensus_end,ranking_raw_end,ranking_topsis_end,ranking_avgrank_end)
names(ranking_process_end) = rank_process
ranking_score = lapply(c("a","b"), function(x) {
  tmp=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean",
        "mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean",
        "time time","time sd")
  if (length(grep("b",x))==1) {tmp=tmp[!(tmp %in% c("mae sd_g","pearson sd_g","rmse sd_g","time sd","consensus sd"))]}
  tmp
})
names(ranking_score) = c("a","b")

date = "241025"
score_path = "../../compute_metrics/scores/"
source("../generic_functions/load_scores_SB_silico.R")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_dnam_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_dnam_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

## ----
## Load libraries
## ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(see)

## ----
## Load scores
## ----
res = load_data(date, score_path)
scores = res$scores
time = res$time
rm(res)
settings = list(meth_dnam_sup,meth_dnam_unsup,
              meth_rna_sup,meth_rna_unsup)
names(settings) = c('dnam-sup','dnam-unsup','rna-sup','rna-unsup')

scores_merged_norm = mapply(function(setting,block) {
  ranking_step1(scores, time) %>% coerce_pearson() %>%
    mutate(DeconvTool = sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
           Block = sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
           Setting=paste(Block,DeconvTool,sep="-")) %>%
    filter(Setting %in% paste(block,setting,sep='-'))},
  settings, rep(c("dnam","rna"),each=2), SIMPLIFY=F)

bestFS = unlist(lapply(readRDS("../main_figs/figure2CD_figure3CD/df_res.rds"), function(df) df$candidate))
scores_merged_norm = lapply(scores_merged_norm, function(x)
  x %>% filter(candidate %in% bestFS))

## ----
## Rank
## ----
ranks = lapply(seq_along(settings), function(blockXclass)
  lapply(ranking_process, function(process)
    lapply(ranking_score, function(scoreskeep)
      process(scores, time, scores_to_keep=scoreskeep) %>%
        mutate(DeconvTool = sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
               Block = sapply(candidate, function(x) strsplit(x,"-")[[1]][1])) %>%
        filter(Block==strsplit(names(settings)[blockXclass],"-")[[1]][1],
               DeconvTool%in%settings[[blockXclass]]))))

ranks = do.call(rbind,lapply(ranks, function(rk_setting)
  do.call(rbind,lapply(seq_along(rk_setting), function(rk_process)
    do.call(rbind,lapply(seq_along(rk_setting[[rk_process]]), function(rk)
      rk_setting[[rk_process]][[rk]] %>%
        mutate(Scores_selection = names(ranking_score)[rk],
               Block = sapply(candidate, function(y) strsplit(y, "-")[[1]][1]),
               FS = sapply(candidate, function(y) strsplit(y, "-")[[1]][3]),
               DeconvTool = sapply(candidate, function(y) strsplit(y, "-")[[1]][2]),
               Candidate = paste(DeconvTool,FS,sep='-'),
               Process = rank_process[rk_process],
               Supervised = ifelse(DeconvTool %in% c(meth_dnam_sup,meth_rna_sup),'Supervised','Unsupervised'),
               sup = ifelse(Supervised=='Supervised','sup','unsup'))))
    ))
  )) %>%
  filter(candidate %in% bestFS)

rank_order = lapply(names(settings), function(x) {
  tmp = ranks %>%
    filter(Block==strsplit(x,'-')[[1]][1],
           sup==strsplit(x,'-')[[1]][2]) %>%
    group_by(Process,Scores_selection) %>%
    mutate(Ranking=rank(1-overall, na.last="keep")) %>%
    select(Candidate,Ranking) %>% ungroup()
  df = plyr::join_all(lapply(rank_process, function(z)
    plyr::join_all(lapply(names(ranking_score), function(y)
      tmp %>% filter(Scores_selection==y,
                   Process==z) %>%
        rename(!!paste0(z,'_',toupper(y)):=Ranking) %>%
        select(Candidate,paste0(z,'_',toupper(y)))),
      by=c('Candidate'), type='full')),
    by=c('Candidate'), type='full')
  data = df[order(df$Sconsensus_A),]
  rownames(data) = data$Candidate
  rownames(data) = gsub("-none"," - No feature selection",
                        gsub("-toast"," - TOAST",
                             gsub("-hvf"," - Highly variable features",rownames(data))))
  data = data[,grep("_",colnames(data))]
  apply(data,2,floor)})
names(rank_order) = names(settings)

rank_order = lapply(rank_order, function(x) {
  rownames(x) = sapply(rownames(x), function(y) strsplit(y,' -')[[1]][1])
  rownames(x) = sapply(rownames(x), function(y) case_when(y=='nnls'~'NNLS',
                                                          y=='rlr'~'RLR',
                                                          y=='elasticnet'~'Elastic net',
                                                          y=='fardeep'~'FARDEEP',
                                                          y=='fardeepsto'~'FARDEEP_sto',
                                                          y=='ols'~'OLS',
                                                          y=='svr'~'SVR', .default = y))
  x
})

## ----
## Plot table to compare rankings (supp figure 2 panels A,B,C,D)
## ----
plot_table = lapply(rank_order, function(data)
  ggtexttable(data,
              theme = ttheme(base_size = 12,
                             colnames.style=colnames_style(color = "grey15", fill = "white", size=14),
                             rownames.style=rownames_style(color = "grey15", fill = "white", face="bold", size=14))))
max_show = sapply(rank_order, function(x) max(x, na.rm=T))
for (k in seq_along(plot_table)) {
  for (i in seq(2,tab_nrow(plot_table[[k]]))) {
    for (j in seq(2,tab_ncol(plot_table[[k]]))) {
      plot_table[[k]] <- plot_table[[k]] %>%
        table_cell_bg(row = i, column = j,
                      fill = viridis::rocket(max_show[k], direction=-1)[rank_order[[k]][i-1,j-1]])
      if (!is.na(rank_order[[k]][i-1,j-1])) { if (rank_order[[k]][i-1,j-1] >= round(3/4*nrow(rank_order[[k]]))) {
        plot_table[[k]] <- plot_table[[k]] %>%
          table_cell_font(row = i, column = j, size=12,
                          color = 'grey90')
      }}
    }
  }
}

heights = c(3,3,5,3)
for (p in seq_along(plot_table)) {
  print(plot_table[[p]])
  ggsave(paste0(folder,"/table_compa_ranks_",names(plot_table)[p],".pdf"), height=heights[p], width=11)
}
rm(i,j,k,p,bestFS,heights,max_show,
   meth_dnam_sup,meth_dnam_unsup,meth_rna_sup,meth_rna_unsup)

## ----
## Prepare Pavao
## ----
mat_trendval <- lapply(scores_merged_norm, function(df)
  df %>%
    ungroup() %>%
    mutate(judge=paste(name_score,dataset,sep="/")) %>%
    select(judge,candidate,trendval) %>%
    rename(val=trendval) %>%
    make_pavao_matrix())

df_pavao <- list()
for (setting in names(scores_merged_norm)) {
  print(paste0("Running setting ",setting))
  df_pavao[[setting]] = NULL
  for (process in seq_along(ranking_process)) {
    for (idx in seq_along(ranking_score)) {
      scores_to_keep = ranking_score[[idx]]
      print(paste0("Running process ", rank_process[process],"_",toupper(names(ranking_score)[idx])))
      winner <- grep(sapply(names(sort(rank_order[[setting]][,2*(process-1)+idx])[1]), function(y)
        case_when(y=='NNLS'~'nnls',
                  y=='RLR'~'rlr',
                  y=='Elastic net'~'elasticnet',
                  y=='FARDEEP'~'fardeep',
                  y=='FARDEEP_sto'~'fardeepsto',
                  y=='OLS'~'ols',
                  y=='SVR'~'svr', .default = y)),
                     unique(scores_merged_norm[[setting]]$candidate), value=T)
      avg_rank_crit <- avg_rank(mat_trendval[[setting]], winner)
      cond_rate <- condorcet_rate(mat_trendval[[setting]], winner)
      generalization_crit <- generalization_criterion(scores_merged_norm[[setting]],
                                                      ranking_process_end[[process]],
                                                      scores_to_keep=scores_to_keep)
      df <- data.frame("Value"=c(avg_rank_crit, cond_rate, generalization_crit),
                       "Criterion"=c("Average rank", "Condorcet rate", "Generalization criterion"),
                       "Process"=rank_process[process],
                       "Scores_selec"=toupper(names(ranking_score)[idx]))
      if (is.null(df_pavao[[setting]])) {
        df_pavao[[setting]] <- df
      } else {
        df_pavao[[setting]] <- rbind(df_pavao[[setting]], df)
      }
    }
  }
}

## ----
## Plot Pavao criteria (supp figure 2 panels E,F)
## ----
df_pavao = do.call(rbind,lapply(names(df_pavao), function(x)
  df_pavao[[x]] %>%
    mutate("Ranking process"=paste(Process,Scores_selec,sep="_"),
           Setting=x)))
df_pavao$Setting <- gsub("dnam-sup","DNAm - Supervised",
                         gsub("dnam-unsup","DNAm - Unsupervised",
                              gsub("rna-sup","RNA - Supervised",
                                   gsub("rna-unsup","RNA - Unsupervised",
                                        df_pavao$Setting))))
df_pavao$`Ranking process`=factor(df_pavao$`Ranking process`,
                                  levels=c(sapply(rank_process,function(x) paste0(x,"_",c("A","B")))))

ggplot(df_pavao, aes(x=`Ranking process`, y=Value)) +
  geom_point(aes(color=Setting)) +
  geom_boxplot(alpha=0) +
  facet_wrap(~Criterion) +
  scale_color_manual(values=c("#F39C12","#E74C3C","#2980B9","#B86F7F")) +
  ylab("Criterion value") +
  theme_modern(axis.text.angle = 30)
ggsave(paste0(folder,"/panelE.pdf"), width=12, height=4)

dot = lapply(c("A","B"), function(x) {
  tmp=c("RMSE", "MAE", "Pearson",
        "RMSE standard deviation", "MAE standard deviation", "Pearson standard deviation",
        "Time","Time standard deviation")
  if (x=="B") {tmp=tmp[-grep("standard deviation",tmp)]}
  tmp
})
names(dot)=c("A","B")
df_dot = data.frame(X=c(rep("A",length(dot$A)),
                        rep("B",length(dot$B))),
                    Y=unlist(dot),
                    Category=c(rep("Performance",3),
                               rep("Stability",3),
                               rep("Time",2),
                               rep("Performance",3),
                               rep("Time",1)))
df_dot$Y <- factor(df_dot$Y, levels=rev(c("RMSE", "MAE", "Pearson",
                                          "RMSE standard deviation", "MAE standard deviation", "Pearson standard deviation",
                                          "Time","Time standard deviation")))
df_dot2 = data.frame(X=rep("B",length(dot$A)-length(dot$B)),
                     Y=setdiff(dot$A,dot$B),
                     Category=c(rep("Stability",3),
                                rep("Time",1)))
ggplot(df_dot, aes(x=X,y=Y)) +
  geom_point(size=3, aes(color=Category)) +
  geom_point(data=df_dot2, size=3, alpha=.15, aes(color=Category)) +
  geom_line(aes(group=X), alpha=.7) +
  ylab("") +
  xlab("") +
  theme_modern() +
  scale_color_manual(values=c("#477CBD","#B54893","#E8B658"))
ggsave(paste0(folder,"/panelF.pdf"), width=5, height=3)
  