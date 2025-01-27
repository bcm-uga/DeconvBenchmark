## ----
## Set parameters, put your own
## ----
score_path = "../../compute_metrics/scores/"
data_path = "../../../data/simulations/"
deconv_path = "../../../deconvolution/results/prediction/"
date = "241025"
n_sim = 10
folder = "figure6"
source("../../compute_metrics/generic_functions_metrics.R")

rare_cell_types <- list(
  BlCL = c("mDC", "Neut"),
  BrCL2 = c("Jurkat", "Thp1"),
  LuCL = c("CD19B", "CD4T", "CD8T", "Monocyte", "Neutrophil", "NKcell", "WBC"),
  PaCL2 = c("B cells", "CD4 T cells", "CD8 T cells", "Macrophages", "Neutrophils"))
rare_datasets = data.frame("Dataset"=c('PaCL2','PaCL2','LuCL','BrCL2','BlCL'),
                           "Omic"=c('rna','dnam','dnam','rna','rna'),
                           "height"=c(8,7,7,8,8))

custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

## ----
## Load libraries
## ----
library(dplyr)
library(ggplot2)
library(see)
library(ggtext)
library(grid)

## ----
## Functions
## ----
calculate_pearson_scores <- function(A_true, A_pred, rare_types) {
  if (all(is.na(A_pred))) {
    return(NA)
  } else {
    A_pred = A_pred[rownames(A_true),]
    pearson_scores <- sapply(rare_types, function(celltype) {
      if (celltype %in% rownames(A_true)) {
        return(metrics(A_true[celltype,], A_pred[celltype,], "pearson"))
      } else {warning("error in cell type")}
    })
    return(pearson_scores)
  }
}
generate_plot <- function(df_all, dataset, block) {
  df_filtered <- df_all %>%
    mutate(DeconvTool = case_when(DeconvTool=='nnls'~'NNLS',
                                  DeconvTool=='rlr'~'RLR',
                                  DeconvTool=='elasticnet'~'Elastic net',
                                  DeconvTool=='fardeep'~'FARDEEP',
                                  DeconvTool=='fardeepsto'~'FARDEEP_sto',
                                  DeconvTool=='ols'~'OLS',
                                  DeconvTool=='svr'~'SVR', .default = DeconvTool)) %>%
    filter(Dataset == dataset, Block == block) %>%
    mutate(y.label = paste("<span style = 'color: ",
                           ifelse(Supervised == "Supervised", "grey40", "#E51400"),
                           ";'>",
                           DeconvTool,
                           "</span>", sep = "")) %>%
    arrange(desc(Overall))
  df_filtered$y.label <- factor(df_filtered$y.label, levels = rev(unique(df_filtered$y.label)))
  
  custom_fill <- c("Common" = "black", "Rare" = "seagreen3")
  custom_color <- c("Common" = "black", "Rare" = "seagreen")
  
  print(ggplot(df_filtered, aes(x = Pearson, y = y.label)) +
          geom_point(aes(color = Type), position=position_dodge(.8)) +
          geom_boxplot(alpha=.5, aes(fill = Type)) +
          labs(title = paste(dataset, "-", block),
               x = "Pearson correlation",
               y = "Deconvolution Tool",
               color = "Cell Type Frequency") +
          scale_color_manual(values = custom_color) +
          scale_fill_manual(values = custom_fill) +
          guides(fill='none') +
          theme_classic(base_size = 14, base_family = "") +
          theme(plot.title = element_text(size = 15, face = "plain", margin = margin(0, 0, 20, 0)),
                plot.title.position = "plot",
                legend.position = "right",
                legend.text = element_blank(),
                legend.title = element_blank(),
                legend.key = element_blank(),
                legend.spacing.x = unit(2, "pt"),
                axis.line = element_blank(),
                axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                axis.title.y = element_blank(),
                axis.text.y = element_markdown(colour = "red",
                                               hjust = 1),
                axis.text = element_text(size = 20),
                axis.ticks.y = element_blank(),
                axis.title = element_text(size = 20, face = "plain"),
                plot.tag = element_text(size = 15, face = "bold"),
                panel.grid = element_blank()) +
          geom_segment(x=min(df_filtered$Pearson, na.rm=T)-.05, xend=min(df_filtered$Pearson, na.rm=T)-.05,
                       y=.5,yend=length(unique(df_filtered$DeconvTool))+.5) +
          geom_segment(x=max(df_filtered$Pearson, na.rm=T)+.05, xend=max(df_filtered$Pearson, na.rm=T)+.05,
                       y=.5,yend=length(unique(df_filtered$DeconvTool))+.5)) +
          geom_segment(x=min(df_filtered$Pearson, na.rm=T)-.05, xend=max(df_filtered$Pearson, na.rm=T)+.05,
                       y=.5,yend=.5) +
          geom_segment(y=seq(1,nrow(df_filtered))+.5,yend=seq(1,nrow(df_filtered))+.5,
                       x=min(df_filtered$Pearson, na.rm=T)-.05,xend=max(df_filtered$Pearson, na.rm=T)+.05)
}

## ----
## Load best FS
## ----
df_res <- readRDS("figure2CD_figure3CD/df_res.rds")
meth_unsup <- unique(unlist(lapply(df_res[grep("unsup", names(df_res))], function(x) unique(x$DeconvTool))))
df_res <- as.data.frame(bind_rows(df_res))

## ----
## Compute Pearson for all cell types + retrieve overall
## ----
df_all = do.call(rbind,lapply(names(rare_cell_types), function(dataset) {
  print(dataset)
  blocks = sapply(c("dnam","rna"), function(x)
    ifelse(length(list.files(paste0(data_path,x),pattern = dataset))!=0,x,""))
  blocks = blocks[blocks!=""]
  do.call(rbind,lapply(blocks, function(block) {
    candidates = df_res %>%
      mutate(Modality=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
             FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3])) %>%
      filter(Modality==block, !duplicated(candidate))
    Nsim = if (n_sim>9 & n_sim<20) {c(paste0("0",seq(9)),as.character(seq(10,n_sim)))} else if (n_sim<10) {paste0("0",seq(n_sim))} else {stop("Modify the code to define Nsim")}
    A_true <- lapply(Nsim, function(sim)
      readRDS(paste0(data_path,block,"/",date,"_",dataset,"_sim",sim,".rds"))$A_ref)
    do.call(rbind,lapply(seq(nrow(candidates)), function(candidate) {
      meth_class = ifelse(candidates$DeconvTool[candidate] %in% meth_unsup,'unsup','sup')
      A_pred <- lapply(Nsim, function(sim)
        readRDS(paste0(deconv_path,block,"/",meth_class,"/",date,"_",dataset,"_Apred_",candidates$FS[candidate],"_",candidates$DeconvTool[candidate],"_sim",sim,".rds"))) 
      all_types <- rownames(A_true[[1]])
      pearsons_sim <- sapply(seq_along(A_true), function(sim)
        calculate_pearson_scores(A_true[[sim]], A_pred[[sim]], all_types))
      if (all(is.na(pearsons_sim))) {pearson=NA} else {pearson=rowMeans(pearsons_sim, na.rm=T)}
      data.frame(Celltype=all_types,
                 Pearson=pearson,
                 DeconvTool=candidates$DeconvTool[candidate],
                 Dataset=dataset,
                 Overall=candidates$overall[candidate],
                 Type=ifelse(all_types %in% rare_cell_types[[dataset]], "Rare", "Common"),
                 Block=block)
    }))
  }))
}))

df_all = df_all %>%
  mutate(Supervised=ifelse(DeconvTool %in% meth_unsup,'Unsupervised','Supervised'))

## ----
## Plot Pearson for all cell types / all methods / 1 dataset for figure 6
## ----
for (i in seq(nrow(rare_datasets))) {
  generate_plot(df_all, rare_datasets$Dataset[i], rare_datasets$Omic[i]) + guides(fill='none',color='none')
  ggsave(paste0(folder,"/pearson_",rare_datasets$Dataset[i],"_",rare_datasets$Omic[i],".pdf"),
         width=8, height=rare_datasets$height[i])
}

## ----
## Compute the mean diff between Pearson of rare/common types
## ----
df_all %>%
  group_by(Dataset,Block,DeconvTool,Type) %>%
  filter(Block=='rna') %>%
  summarize(Mean_Pearson = mean(Pearson)) %>%
  mutate(Mean_Pearson_common = Mean_Pearson[Type=='Common'],
         Mean_Pearson_rare = Mean_Pearson[Type=='Rare'],
         Delta = Mean_Pearson_common-Mean_Pearson_rare) %>%
  filter(!duplicated(Delta)) %>%
  filter(Delta<.1)
