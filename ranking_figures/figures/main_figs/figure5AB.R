## ----
## Set parameters, put your own
## ----
dates = c("241025LessDisp","241025","241025MoreDisp") # simulations done with different levels of dispersion
name_file = 'noise'
#dates = c("241025LessSamples","241025") # simulations done with different number of samples
#name_file = 'sample'

score_path = "../../compute_metrics/scores/"
n_permutations <- 1000
source("../generic_functions/load_scores_SB_silico.R")
source("../generic_functions/ranking_process.R")
source("../generic_functions/ranking_pval.R")
folder = "figure5AB"

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
library(funkyheatmap)
library(tibble)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(see)

## ----
## Load scores & rank
## ----
df_scores = list()
df_time = list()
for (idx in seq_along(dates)) {
  date = dates[idx]
  res = load_data(date, score_path)
  scores = res$scores[grep("none",res$scores$candidate),]
  time = res$time[grep("none",res$time$candidate),]
  rm(res)
  scores$candidate = paste(scores$candidate,date,sep="-")
  time$candidate = paste(time$candidate,date,sep="-")
  df_scores[[idx]] = scores
  df_time[[idx]] = time
}
scores = do.call(rbind,df_scores)
time = do.call(rbind,df_time)
rank_consensus = ranking_consensus(scores,time)
rank_consensus = rank_consensus %>%
  mutate(Block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         DeconvTool=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]),
         Date=sapply(candidate, function(x) strsplit(x,"-")[[1]][4]),
         Supervised=ifelse(DeconvTool %in% c(meth_dnam_sup,meth_rna_sup),"sup","unsup"),
         Setting=paste(Block,Supervised,sep="-"))
rank_consensus$Date[rank_consensus$Date=="241025"] = "Baseline"
rank_consensus$Date[rank_consensus$Date=="241025LessSamples"] = "LessSamples"
rank_consensus$Date[rank_consensus$Date=="241025LessDisp"] = "LessDisp"
rank_consensus$Date[rank_consensus$Date=="241025MoreDisp"] = "MoreDisp"

## ----
## Delta dataframe
## ----
rank_consensus$Date = factor(rank_consensus$Date, levels=c("LessSamples","LessDisp","Baseline","MoreDisp"))
rank_consensus$Date = droplevels(rank_consensus$Date)
rank_consensus$DeconvTool2 = paste0(rank_consensus$DeconvTool," (",rank_consensus$Block,")")

delta_df = rank_consensus %>%
  group_by(DeconvTool,Block) %>%
  mutate(Delta=overall-overall[Date=="Baseline"]) %>%
  filter(Date!="Baseline") %>%
  mutate(Block=sub('dnam','DNAm',
                   sub('rna','RNA',Block)))

## ----
## Pval for each method between simu params
## ----
df_scores_norm = ranking_step1(scores,time) %>%
  coerce_pearson() %>%
  mutate(block=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
         class=ifelse(sapply(candidate, function(x) strsplit(x,"-")[[1]][2]) %in% c(meth_dnam_sup,meth_rna_sup),'sup','unsup'),
         simuparam=sapply(candidate, function(x) strsplit(x,"-")[[1]][4]),
         method=sapply(candidate, function(x) strsplit(x,"-")[[1]][2]))

# compute p-values
for (Block in unique(df_scores_norm$block)) {
  print(paste("Running", Block))
  for (Class in unique(df_scores_norm$class)) {
    print(paste("Running", Class))
    if (!file.exists(paste0(folder,"/pval_",name_file,"_",Block,"_",Class,".rds"))) {
      df_filter = df_scores_norm %>%
        filter(block==Block,
               class==Class)
      pval_list = list()
      for (Method in unique(df_filter$method)) {
        df_filter2 = df_filter %>% filter(method==Method)
        candidates = paste(paste(Block,Method,"none",sep='-'),rev(dates),sep='-')
        pval_res <- paired_permutation_test(df_filter2,
                                            ranking_consensus_end,
                                            n_permutations=n_permutations,
                                            candidates=candidates)
        pval_list[[Method]] = pval_res
      }
      saveRDS(pval_list,paste0(folder,"/pval_",name_file,"_",Block,"_",Class,".rds"))
    }
  }
}
pval_list = lapply(list.files(folder, pattern=paste0('pval_',name_file),full.names=T),readRDS)
names(pval_list) = sapply(sapply(list.files(folder, pattern=paste0('pval_',name_file)), function(x)
  paste0(strsplit(x,"_")[[1]][c(3,4)],collapse='-')), function(x)
    strsplit(x,".rds")[[1]][1])
pval_df = do.call(rbind,mapply(function(x,y) {
  do.call(rbind,lapply(x, function(z) {
    data.frame(comparison=names(z$pvalues),
               pval=unname(z$pvalues)) %>%
      mutate(DeconvTool=sapply(comparison, function(compa) strsplit(compa,"-")[[1]][2]),
             Block=sapply(y, function(setting) strsplit(setting,"-")[[1]][1]),
             Supervised=sapply(y, function(setting) strsplit(setting,"-")[[1]][2]),
             Date1=sapply(comparison, function(compa) strsplit(strsplit(compa,">")[[1]][1],"-")[[1]][4]),
             Date2=sapply(comparison, function(compa) strsplit(compa,"-")[[1]][7])) %>%
      filter(Date1=="241025" | Date2=="241025") %>%
      mutate(Date=ifelse(Date1=="241025",Date2,Date1)) %>%
      select(pval,DeconvTool,Block,Date,Supervised)}))
},x=pval_list,y=names(pval_list),SIMPLIFY=F))
pval_df$Date[pval_df$Date=="241025"] = "Baseline"
pval_df$Date[pval_df$Date=="241025LessSamples"] = "LessSamples"
pval_df$Date[pval_df$Date=="241025LessDisp"] = "LessDisp"
pval_df$Date[pval_df$Date=="241025MoreDisp"] = "MoreDisp"
pval_df$Block[pval_df$Block=="dnam"] = "DNAm"
pval_df$Block[pval_df$Block=="rna"] = "RNA"

pval_df = pval_df %>%
  group_by(DeconvTool,Block,Date) %>%
  mutate(p_lab = if_else(pval < .001, "***",
                         if_else(pval < .01, "**",
                                 if_else(pval < .05, "*", ""))),
         y1 = delta_df$Delta[delta_df$DeconvTool==DeconvTool &
                            delta_df$Block==Block &
                            delta_df$Date==Date],
         y = ifelse(y1<0, .02, y1))

## ----
## Plot figure 5 panels A,B
## ----
if (length(dates)>2) {my_comparisons <- list( c(levels(rank_consensus$Date)[1], levels(rank_consensus$Date)[2]),
                                              c(levels(rank_consensus$Date)[2], levels(rank_consensus$Date)[3]))
} else {my_comparisons <- list( c(levels(rank_consensus$Date)[1], levels(rank_consensus$Date)[2]))}


ggplot(delta_df %>% filter(Supervised=='unsup'), aes(x=Block, y=Delta, fill=DeconvTool)) +
  geom_col(width = 0.8, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, col='red', linetype='dashed') +
  ylab("\u0394 (Variation-Baseline)") +
  xlab("") +
  facet_wrap(~Date) +
  geom_text(data=pval_df %>% filter(Supervised=='unsup'),
            aes(x=Block, y=y, color=DeconvTool, label=p_lab),
            position = position_dodge(width = .9),
            size=6) +
  scale_color_metro(reverse=T) +
  scale_fill_metro(reverse=T) +
  guides(colour="none") +
  theme_modern(axis.title.size = 16, legend.title.size = 16, legend.text.size = 14)
ggsave(paste0(folder,"/delta_barplot_",name_file,".pdf"), width=10, height=4, device = cairo_pdf)

