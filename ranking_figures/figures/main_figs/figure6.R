## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
library(ggtext)
library(grid)
custom_palette <- c("#D11141", "#F47F43", "#FFA500", "#FFD700", 
                    "#C0C0C0", "#808080", "darkolivegreen3", "darkolivegreen", 
                    "#76958F", "lightblue3", "#4682B4", "#005F6B")

folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
rare_cell_types <- list(
  Hoek = c("mDC", "Neut"),
  Cobos = c("Jurkat", "Thp1"),
  He = c("CD19B", "CD4T", "CD8T", "Monocyte", "Neutrophil", "NKcell", "WBC"),
  lot1 = c("B cells", "CD4 T cells", "CD8 T cells", "Macrophages", "Neutrophils"))

## ----
## Functions
## ----
calculate_pearson_scores <- function(A_true, A_pred, rare_types) {
  if (all(is.na(A_pred))) {
    return(NA)
  } else {
    source("../../src/1_generate_score.R")
    A_pred = A_pred[rownames(A_true),]
    pearson_scores <- sapply(rare_types, function(celltype) {
      if (celltype %in% rownames(A_true)) {
        return(score_perf(A_true[celltype,], A_pred[celltype,], "pearson"))
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
violin_data_celltype <- function(dataset, cell_type) {
  blocks = sapply(c("met","rna"), function(x)
    ifelse(length(list.files(paste0("../0simu/simulations/",x),pattern = dataset))!=0,x,""))
  blocks = blocks[blocks!=""]
  df_violin <- do.call(rbind, lapply(blocks, function(block) {
    candidates <- df_res %>%
      mutate(Modality = sapply(candidate, function(x) strsplit(x, "-")[[1]][1]),
             FS = sapply(candidate, function(x) strsplit(x, "-")[[1]][3])) %>%
      filter(Modality == block, !duplicated(candidate))
    A_true <- do.call(cbind, lapply(c(paste0("0",seq(9)),"10"), function(sim)
      readRDS(paste0("../0simu/simulations/",block,"/231027_",dataset,"_sim",sim,".rds"))$A_ref))
    do.call(rbind, lapply(seq(nrow(candidates)), function(candidate) {
      meth_class <- ifelse(candidates$DeconvTool[candidate] %in% meth_unsup, 'unsup', 'sup')
      A_pred <- do.call(cbind, lapply(c(paste0("0", seq(9)), "10"), function(sim) {
        readRDS(paste0("../1SB/deconv/", block, "/", meth_class, "/231027/231027_", dataset, "_Apred", candidates$FS[candidate], "_", candidates$DeconvTool[candidate], "_sim", sim, ".rds"))
      }))
      if (!all(is.na(A_pred))) {
        if (cell_type %in% rownames(A_true)) {
          data.frame(Proportion = c(A_true[cell_type, ], A_pred[cell_type, ]),
                     Type = rep(c("Truth", "Prediction"), each = ncol(A_true)),
                     DeconvTool = candidates$DeconvTool[candidate],
                     Block = block,
                     DeconvTool2 = paste(candidates$DeconvTool[candidate],block,sep='-'))
        }
      }
    }))
  }))
  return(df_violin)
}
violin_data_method <- function(dataset, method) {
  blocks = sapply(c("met","rna"), function(x)
    ifelse(length(list.files(paste0("../0simu/simulations/",x),pattern = dataset))!=0,x,""))
  blocks = blocks[blocks!=""]
  df_violin <- do.call(rbind, lapply(blocks, function(block) {
    candidates <- df_res %>%
      mutate(Modality = sapply(candidate, function(x) strsplit(x, "-")[[1]][1]),
             FS = sapply(candidate, function(x) strsplit(x, "-")[[1]][3])) %>%
      filter(Modality == block, !duplicated(candidate),
             DeconvTool == method)
    if (nrow(candidates)>0) {
      A_true <- do.call(cbind, lapply(c(paste0("0",seq(9)),"10"), function(sim)
        readRDS(paste0("../0simu/simulations/",block,"/231027_",dataset,"_sim",sim,".rds"))$A_ref))
      meth_class <- ifelse(candidates$DeconvTool %in% meth_unsup, 'unsup', 'sup')
      A_pred <- do.call(cbind, lapply(c(paste0("0", seq(9)), "10"), function(sim) {
          readRDS(paste0("../1SB/deconv/", block, "/", meth_class, "/231027/231027_", dataset, "_Apred", candidates$FS, "_", candidates$DeconvTool, "_sim", sim, ".rds"))
        }))
      data.frame(Proportion = c(A_true[rare_cell_types[[dataset]],], A_pred[rare_cell_types[[dataset]],]),
               Cell = c(rep(rare_cell_types[[dataset]],ncol(A_true)),
                        rep(rare_cell_types[[dataset]],ncol(A_pred))),
               Type = rep(c("Truth", "Prediction"), each = ncol(A_true)*length(rare_cell_types[[dataset]])),
               DeconvTool = candidates$DeconvTool,
               Block = block,
               DeconvTool2 = paste(candidates$DeconvTool,block,sep='-'))
    }
  }))
  return(df_violin)
}

## ----
## Load best FS
## ----
df_res <- readRDS("fig2/df_res.rds")
meth_unsup <- unique(unlist(lapply(df_res[grep("unsup", names(df_res))], function(x) unique(x$DeconvTool))))
df_res <- as.data.frame(bind_rows(df_res))

## ----
## Compute Pearson for all cell types + retrieve overall
## ----
df_all = do.call(rbind,lapply(names(rare_cell_types), function(dataset) {
  print(dataset)
  blocks = sapply(c("met","rna"), function(x)
    ifelse(length(list.files(paste0("../0simu/simulations/",x),pattern = dataset))!=0,x,""))
  blocks = blocks[blocks!=""]
  do.call(rbind,lapply(blocks, function(block) {
    candidates = df_res %>%
      mutate(Modality=sapply(candidate, function(x) strsplit(x,"-")[[1]][1]),
             FS=sapply(candidate, function(x) strsplit(x,"-")[[1]][3])) %>%
      filter(Modality==block, !duplicated(candidate))
    A_true <- lapply(c(paste0("0",seq(9)),"10"), function(sim)
      readRDS(paste0("../0simu/simulations/",block,"/231027_",dataset,"_sim",sim,".rds"))$A_ref)
    do.call(rbind,lapply(seq(nrow(candidates)), function(candidate) {
      meth_class = ifelse(candidates$DeconvTool[candidate] %in% meth_unsup,'unsup','sup')
      A_pred <- lapply(c(paste0("0",seq(9)),"10"), function(sim)
        readRDS(paste0("../1SB/deconv/",block,"/",meth_class,"/231027/231027_",dataset,"_Apred",candidates$FS[candidate],"_",candidates$DeconvTool[candidate],"_sim",sim,".rds"))) 
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
## Plot Pearson for all cell types / all methods / 1 dataset
## ----
rare_datasets = data.frame("Dataset"=c('lot1','lot1','He','Cobos','Hoek'),
                           "Omic"=c('rna','met','met','rna','rna'),
                           "height"=c(8,7,7,8,8))
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

## ----
## Violin Plots for Rare Cell Types
## ----

#Charge data for neutrophils
violin_neut1 <- violin_data_celltype('Hoek', 'Neut')
violin_neut2 <- violin_data_celltype('He', 'Neutrophil')
violin_neut3 <- violin_data_celltype('lot1', 'Neutrophils')

#Charge data for debCAM
violin_debcam1 <- violin_data_method('Cobos', 'debCAM')
violin_debcam2 <- violin_data_method('Hoek', 'debCAM')
violin_debcam3 <- violin_data_method('He', 'debCAM')
violin_debcam4 <- violin_data_method('lot1', 'debCAM')

df_atrue_neut <- subset(rbind(
  transform(violin_neut1, Dataset = 'Hoek'),
  transform(violin_neut2, Dataset = 'He'),
  transform(violin_neut3, Dataset = 'lot1')
), Type == "Truth")
df_atrue_debcam <- subset(rbind(
  transform(violin_debcam1, Dataset = 'Cobos'),
  transform(violin_debcam2, Dataset = 'Hoek'),
  transform(violin_debcam3, Dataset = 'He'),
  transform(violin_debcam4, Dataset = 'lot1')
), Type == "Truth")

df_apred_neut <- subset(rbind(
  transform(violin_neut1, Dataset = 'Hoek'),
  transform(violin_neut2, Dataset = 'He'),
  transform(violin_neut3, Dataset = 'lot1')
), Type == "Prediction")
df_apred_debcam <- subset(rbind(
  transform(violin_debcam1, Dataset = 'Cobos'),
  transform(violin_debcam2, Dataset = 'Hoek'),
  transform(violin_debcam3, Dataset = 'He'),
  transform(violin_debcam4, Dataset = 'lot1')
), Type == "Prediction")

df_atrue_neut$DeconvTool <- "Truth"
df_atrue_neut$DeconvTool2 <- "Truth"
df_atrue_debcam$DeconvTool <- "Truth"
df_atrue_debcam$DeconvTool2 <- "Truth"

df_combined_neut <- rbind(df_atrue_neut,df_apred_neut)
df_combined_debcam <- rbind(df_atrue_debcam,df_apred_debcam)

df_combined_neut$DeconvTool2 = factor(df_combined_neut$DeconvTool2, levels=c('Truth',sort(unique(df_combined_neut$DeconvTool2[df_combined_neut$DeconvTool2!='Truth']))))
df_combined_debcam$DeconvTool2 = factor(df_combined_debcam$DeconvTool2, levels=c('Truth',sort(unique(df_combined_debcam$DeconvTool2[df_combined_debcam$DeconvTool2!='Truth']))))

df_combined_neut = df_combined_neut %>% mutate(Supervised = ifelse(DeconvTool %in% meth_unsup, F, 
                                                         ifelse(DeconvTool2 == 'Truth', NA, T)))
df_combined_debcam = df_combined_debcam %>% mutate(Supervised = ifelse(DeconvTool %in% meth_unsup, F, 
                                                                   ifelse(DeconvTool2 == 'Truth', NA, T)))

df_combined_neut$Block[df_combined_neut$DeconvTool2=='Truth'] = NA
df_combined_debcam$Block[df_combined_debcam$DeconvTool2=='Truth'] = NA

ggplot(df_combined_neut %>% filter(!Supervised | is.na(Supervised)),
       aes(x = Dataset, y = Proportion, fill = DeconvTool2, color = Block)) +
  geom_violin(alpha = 0.9, position = position_dodge(width = 0.9), scale = 'width') +
  labs(
    title = "Neutrophils",
    x = "Dataset",
    y = "Proportion") +
  scale_color_flat() +
  scale_fill_metro() +
  theme_modern()
ggsave(paste0(folder,'/supp_violin_neutro_unsup.pdf'), width=10, height=5)

ggplot(df_combined_neut %>% filter(!Supervised | is.na(Supervised)),
       aes(x = DeconvTool2, y = Proportion, fill = Dataset, color = Block)) +
  geom_violin(alpha = 0.9, position = position_dodge(width = 0.9), scale = 'width') +
  labs(
    title = "Neutrophils",
    x = "Dataset",
    y = "Proportion") +
  scale_color_flat() +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,'Reds')) +
  theme_modern(axis.text.angle = 30)
ggsave(paste0(folder,'/supp_violin_neutro_unsup_v2.pdf'), width=10, height=5)

df_combined_debcam$Cell[df_combined_debcam$Cell=="CD19B"] = "B cells"
df_combined_debcam$Cell[df_combined_debcam$Cell=="CD4T"] = "CD4 T cells"
df_combined_debcam$Cell[df_combined_debcam$Cell=="CD8T"] = "CD8 T cells"
df_combined_debcam$Cell[df_combined_debcam$Cell=="Neut"] = "Neutrophils"
df_combined_debcam$Cell[df_combined_debcam$Cell=="Neutrophil"] = "Neutrophils"
ggplot(df_combined_debcam %>% filter(Block=='rna'),
       aes(x = Cell, y = Proportion, fill = Dataset)) +
  geom_violin(alpha = 0.9, position = position_dodge(width = 0.9), scale = 'width') +
  labs(
    title = "debCAM",
    x = "Cell Type",
    y = "Proportion") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,'Reds')) +
  theme_modern(axis.text.angle = 45)
ggsave(paste0(folder,'/supp_violin_debcam_rna.pdf'), width=10, height=5)
