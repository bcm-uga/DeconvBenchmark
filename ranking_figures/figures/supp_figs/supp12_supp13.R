## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

## ----
## Load data (methods perf and charac) - only IN SILICO
## ----
source('../../src/load_scores_SB.R')
charac_method = readRDS("characteristics_method/charac_meth.rds")
scores_inter = do.call(rbind,readRDS("../fig/fig2/df_res.rds")) %>%
  select(candidate,DeconvTool,overall) %>%
  filter(!duplicated(candidate)) %>%
  mutate(Block = sapply(candidate, function(x) strsplit(x,'-')[[1]][1]),
         Block = case_when(Block=='met'~'DNAm',
                           Block=='rna'~'RNA'),
         Supervised = ifelse(DeconvTool %in% c(meth_met_unsup,meth_rna_unsup),'Unsupervised','Supervised'),
         DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                DeconvTool=='rlr'~'RLR',
                                DeconvTool=='elasticnet'~'Elastic net',
                                DeconvTool=='fardeep'~'FARDEEP',
                                DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                DeconvTool=='ols'~'OLS',
                                DeconvTool=='svr'~'SVR', .default = DeconvTool))

scores_inter$Approach = sapply(scores_inter$DeconvTool, function(x) charac_method[x,'Approach'])
scores_inter$STO = sapply(scores_inter$DeconvTool, function(x) charac_method[x,'STO'])

scores_inter = scores_inter %>% ungroup() %>%
mutate(x.lab.approach = sapply(seq_along(DeconvTool), function(idx)
    paste("<span style = 'color: ",
          ifelse(Supervised[idx]=="Supervised",
                "grey40", "#E51400"),
          ";'>",
          Approach[idx],
          "</span>", sep = "")))

order_factor_approach = scores_inter %>%
  group_by(Approach) %>%
  summarise(perf_approach = mean(overall)) %>%
  arrange(desc(perf_approach)) %>% pull(Approach)
scores_inter$x.lab.approach = factor(scores_inter$x.lab.approach, levels=c(sapply(order_factor_approach, function(x)
  sapply(c("grey40", "#E51400"), function(y)
    paste0("<span style = 'color: ",y,";'>",x,"</span>")))))

## ----
## Plot : Link each charac with overall score
## ----
colnames(scores_inter)
scores_inter$DeconvTool = factor(scores_inter$DeconvTool,
                                 levels=c(sort(c(unique(c(meth_rna_sup,meth_met_sup)),"NNLS","RLR","Elastic net","FARDEEP","FARDEEP_sto","OLS","SVR")),
                                          sort(unique(c(meth_met_unsup,meth_rna_unsup)))))
scores_inter$DeconvTool = droplevels(scores_inter$DeconvTool)

ggplot(scores_inter, aes(x=x.lab.approach, y=overall)) +
  geom_point(aes(color=DeconvTool, shape=Block), size=4) +
  #geom_boxplot() +
  scale_color_manual(values=c(custom_palette,palette_metro()(length(unique(c(meth_met_unsup,meth_rna_unsup)))))) +
  ylab("Overall benchmark score") +
  theme_classic(base_size = 11, base_family = "") +
  theme(plot.title = element_text(size = 15, face = "plain", margin = margin(0, 0, 20, 0)),
        plot.title.position = "plot", 
        legend.position = "right",
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 16),
        legend.key = element_blank(), 
        legend.spacing.x = unit(2, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(colour="red",
                                               hjust = 1, angle=30),
        axis.text = element_text(size = 16),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16, face = "plain"),
        plot.tag = element_text(size = 15, face = "bold"), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(paste0(folder,'/approach_sensitivity.pdf'), width=12, height=7)

ggplot(scores_inter, aes(x=STO, y=overall)) +
  geom_point(aes(color=DeconvTool, shape=Block)) +
  geom_boxplot() +
  theme_modern()
