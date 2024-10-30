## ----
## Set parameters, put your own
## ----
score_path = "../../compute_metrics/scores/"
data_path = "../../../data/simulations/"
deconv_path = "../../../deconvolution/results/prediction/"
time_path = "../../../deconvolution/results/timing/"
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

## ----
## Load libraries
## ----
library(dplyr)
library(ggplot2)
library(see)

## ----
## Load fig2/3 res
## ----
df_res <- readRDS("figure2_CD_figure3_CD/df_res.rds")

## ----
## Select 1 setting / 2 methods / 1 dataset with rare cell types: (you can put your own selection)
## Here we chose rna-sup / InstaPrism+CIBERSORT with no feature selection / 241025_PaCL2_sim01
## ----
res = df_res$`rna-sup` %>% filter(DeconvTool %in% c('InstaPrism','CIBERSORT'),
                                  dataset == 'PaCL2-rna')
proportions_insta = list('true'=readRDS(paste0(data_path,"rna/241025_PaCL2_sim01.rds"))$A_ref,
                         'pred'=readRDS(paste0(deconv_path,"/rna/sup/241025_PaCL2_Apred_none_InstaPrism_sim01.rds")))
proportions_ciber = list('true'=readRDS(paste0(data_path,"rna/241025_PaCL2_sim01.rds"))$A_ref,
                         'pred'=readRDS(paste0(deconv_path,"/rna/sup/241025_PaCL2_Apred_none_CIBERSORT_sim01.rds")))

time_insta = readRDS(paste0(time_path,"/rna/sup/241025_PaCL2_timing_none_InstaPrism_sim01.rds"))
time_ciber = readRDS(paste0(time_path,"/rna/sup/241025_PaCL2_timing_none_CIBERSORT_sim01.rds"))

df_insta = data.frame(Atrue = c(proportions_insta$true[sort(rownames(proportions_insta$true)),]),
                      Apred = c(proportions_insta$pred[sort(rownames(proportions_insta$pred)),]),
                      CellType = sort(rownames(proportions_insta$true)))
df_ciber = data.frame(Atrue = c(proportions_ciber$true[sort(rownames(proportions_ciber$true)),]),
                      Apred = c(proportions_ciber$pred[sort(rownames(proportions_ciber$pred)),]),
                      CellType = sort(rownames(proportions_ciber$true)))

# all immune cells are rare
rare_types = c("B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils")
df_insta$CellType = factor(df_insta$CellType, levels = c("Cancer basal","Cancer classical","Endothelial","Fibroblasts",
                                                         "B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils"))
df_ciber$CellType = factor(df_ciber$CellType, levels = c("Cancer basal","Cancer classical","Endothelial","Fibroblasts",
                                                         "B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils"))

## ----
## Plot (figure 8 panels E,F)
## ----
ggarrange(
  ggplot(df_insta, aes(x=Atrue, y=Apred)) +
    geom_point(aes(col=CellType)) +
    annotate(geom='text',
             x = 0.04,
             y = .72,
             hjust = 'inward',
             size = 4,
             label = paste0("RMSE = ",round(sqrt(mean((df_insta$Atrue-df_insta$Apred)^2)), 3))) +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(palette_social()(5)[-c(3)],
                                RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  ggplot(df_ciber, aes(x=Atrue, y=Apred)) +
    geom_point(aes(col=CellType)) +
    annotate(geom='text',
             x = 0.04,
             y = .76,
             hjust = 'inward',
             size = 4,
             label = paste0("RMSE = ",round(sqrt(mean((df_ciber$Atrue-df_ciber$Apred)^2)), 3))) +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(palette_social()(5)[-c(3)],
                                RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  common.legend = T
)
ggsave(paste0(folder,'/panelE.pdf'), width=10, height=5)

ggarrange(
  ggplot(df_insta %>% filter(CellType %in% rare_types), aes(x=Atrue, y=Apred, col=CellType)) +
    geom_point() +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.3)) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  ggplot(df_ciber %>% filter(CellType %in% rare_types), aes(x=Atrue, y=Apred, col=CellType)) +
    geom_point() +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.3)) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  common.legend = T
)
ggsave(paste0(folder,'/panelF.pdf'), width=10, height=5)

## ----
## Plot time (figure 8G)
## ----
time_insta = data.frame(cat = c("Time",'0'),
                        fraction = c(time_insta/200,(1-time_insta/200))) %>%
  mutate(ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1)))
time_ciber = data.frame(cat = c("Time",'0'),
                        fraction = c(time_ciber/200,(1-time_ciber/200))) %>%
  mutate(ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1)))

ggplot(time_insta, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=cat)) +
  geom_rect() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("white","orangered3")) +
  xlim(c(3.5, 4)) +
  guides(fill='none') +
  geom_segment(x=3.5,y=0,xend=4,yend=0,col='black') +
  geom_segment(x=3.5,y=time_insta$ymax[1],xend=4,yend=time_insta$ymax[1],col='black') +
  scale_y_continuous(breaks = seq(0,1,by=1/12)) +
  theme_bw(base_size = 0, base_line_size = 1) %+replace% 
  theme(panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "black"),
        legend.background = element_blank(), 
        legend.key = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), strip.background = element_blank(), 
        plot.background = element_blank())
ggsave(paste0(folder,'/panelG1.pdf'))
ggplot(time_ciber, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=cat)) +
  geom_rect() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("white","orangered3")) +
  xlim(c(3.5, 4)) +
  guides(fill='none') +
  geom_segment(x=3.5,y=0,xend=4,yend=0,col='black') +
  geom_segment(x=3.5,y=time_ciber$ymax[1],xend=4,yend=time_ciber$ymax[1],col='black') +
  scale_y_continuous(breaks = seq(0,1,by=1/12)) +
  theme_bw(base_size = 0, base_line_size = 1) %+replace% 
  theme(panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "black"),
        legend.background = element_blank(), 
        legend.key = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), strip.background = element_blank(), 
        plot.background = element_blank())
ggsave(paste0(folder,'/panelG2.pdf'))



