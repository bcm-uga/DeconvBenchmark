## ----
## Set parameters, put your own
## ----
folder = "supp6"

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
## Functions
## ----
scale_fill_fun = function(x) {
  if (x%in%c(2,4)) {
    scale_fill_metro(reverse=T)
  }
  else if (x==1) {scale_fill_manual(values=custom_palette[c(1,4,7,10,12)])}
  else if (x==3) {scale_fill_manual(values=custom_palette)}
}

## ----
## Load df
## ----
df_fig2 = readRDS("../main_figs/figure2CD_figure3CD/df_res.rds")

## ----
## Plot supp figure 6
## ----
df_fig2 = lapply(df_fig2, function(df) {
  x = df
  x$candidate = sapply(x$candidate, function(y)
    sub("-none"," - No feature selection",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("-hvf"," - Highly variable features",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("-toast"," - TOAST",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("dnam-","",y))
  x$candidate = sapply(x$candidate, function(y)
    sub("rna-","",y))
  x$dataset = sapply(x$dataset, function(y)
    strsplit(y,'-')[[1]][1])
  x$candidate = factor(x$candidate, levels=x %>% arrange(desc(overall)) %>% pull(candidate) %>% unique())
  x = x %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool))
  x$DeconvTool = factor(x$DeconvTool, levels=x %>% arrange(desc(overall)) %>% pull(DeconvTool) %>% unique())
  x})

widths=c(7,7,8,7)
heights=c(4,4,4,4)
for (i in seq_along(df_fig2)) {
  print(ggplot(df_fig2[[i]], aes(x=dataset, y=aggregated, fill=DeconvTool)) +
          geom_bar(stat='identity', position='dodge') +
          scale_fill_fun(i) +
          ylab("Aggregated benchmark score") +
          xlab("") +
          theme_modern(axis.text.angle = 30, axis.title.size = 15))
  ggsave(paste0(folder,"/",names(df_fig2)[i],".pdf"), width = widths[i], height = heights[i])
}
