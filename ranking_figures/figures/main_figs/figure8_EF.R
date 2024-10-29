## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]
colors = c("#00916e", "#ffcf00", "#ee6123", "#fa003f")
colors = c("#230007", "#d7cf07", "#d98324", "#a40606")
colors = c("#ad343e", "#474747", "#f2af29", "#000000")
colors = c("#290F1F","#5f0f40","#9a031e","#cb793a")

## ----
## Load fig2 res
## ----
df_res <- readRDS("fig2/df_res.rds")

## ----
## Pick 2 methods : InstaPrism+CIBERSORT, 1 dataset : PaCL2 sim01, 1 setting : rna_sup
## ----
res = df_res$`rna-sup` %>% filter(DeconvTool %in% c('InstaPrism','CIBERSORT'),
                                  dataset == 'PaCL2-rna')
proportions_InstaPrism = list('true'=readRDS("../0simu/simulations/rna/231027_lot1_sim01.rds")$A_ref,
                       'pred'=readRDS("../1SB/deconv/rna/sup/231027/231027_lot1_Aprednone_InstaPrism_sim01.rds"))
proportions_ciber = list('true'=readRDS("../0simu/simulations/rna/231027_lot1_sim01.rds")$A_ref,
                       'pred'=readRDS("../1SB/deconv/rna/sup/231027/231027_lot1_Apredtoast_CIBERSORT_sim01.rds"))

df_InstaPrism = data.frame(Atrue = c(proportions_InstaPrism$true[sort(rownames(proportions_InstaPrism$true)),]),
                    Apred = c(proportions_InstaPrism$pred[sort(rownames(proportions_InstaPrism$pred)),]),
                    CellType = sort(rownames(proportions_InstaPrism$true)))
df_ciber = data.frame(Atrue = c(proportions_ciber$true[sort(rownames(proportions_ciber$true)),]),
                    Apred = c(proportions_ciber$pred[sort(rownames(proportions_ciber$pred)),]),
                    CellType = sort(rownames(proportions_ciber$true)))

# all immune cells are rare
immune = c("B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils")
df_InstaPrism$CellType = factor(df_InstaPrism$CellType, levels = c("Cancer basal","Cancer classical","Endothelial","Fibroblasts",
                                                                   "B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils"))
df_ciber$CellType = factor(df_ciber$CellType, levels = c("Cancer basal","Cancer classical","Endothelial","Fibroblasts",
                                                         "B cells","CD4 T cells","CD8 T cells","Macrophages","Neutrophils"))

## ----
## Plot
## ----
ggarrange(
  ggplot(df_InstaPrism, aes(x=Atrue, y=Apred)) +
    geom_point(aes(col=CellType)) +
    annotate(geom='text',
             x = 0.04,
             y = .72,
             hjust = 'inward',
             size = 4,
             label = paste0("RMSE = ",round(sqrt(mean((df_InstaPrism$Atrue-df_InstaPrism$Apred)^2)), 3))) +
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
ggsave(paste0(folder,'/scatter_global.pdf'), width=10, height=5)

ggarrange(
  ggplot(df_InstaPrism %>% filter(CellType %in% immune), aes(x=Atrue, y=Apred, col=CellType)) +
    geom_point() +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.3)) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  ggplot(df_ciber %>% filter(CellType %in% immune), aes(x=Atrue, y=Apred, col=CellType)) +
    geom_point() +
    ggpmisc::stat_poly_eq() +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(name="Greens",n=6)[2:6])) +
    theme_modern() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits=c(0,.3)) +
    geom_abline(intercept = 0, slope = 1, col='grey10', linewidth=.2),
  common.legend = T
)
ggsave(paste0(folder,'/scatter_rare.pdf'), width=10, height=5)

## ----
## Plot time
## ----
time_InstaPrism = readRDS("../1SB/timing/rna/sup/231027/231027_lot1_timingnone_InstaPrism_sim01.rds")
time_ciber = readRDS("../1SB/timing/rna/sup/231027/231027_lot1_timingtoast_CIBERSORT_sim01.rds")
time_InstaPrism = data.frame(cat = c("Time",'0'),
                             fraction = c(time_InstaPrism/200,(1-time_InstaPrism/200))) %>%
  mutate(ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1)))
time_ciber = data.frame(cat = c("Time",'0'),
                             fraction = c(time_ciber/200,(1-time_ciber/200))) %>%
  mutate(ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1)))

ggplot(time_InstaPrism, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3.5, fill=cat)) +
  geom_rect() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("white","orangered3")) +
  xlim(c(3.5, 4)) +
  guides(fill='none') +
  geom_segment(x=3.5,y=0,xend=4,yend=0,col='black') +
  geom_segment(x=3.5,y=time_InstaPrism$ymax[1],xend=4,yend=time_InstaPrism$ymax[1],col='black') +
  scale_y_continuous(breaks = seq(0,1,by=1/12)) +
  theme_bw(base_size = 0, base_line_size = 1) %+replace% 
  theme(panel.grid.major = element_line(colour = "black"),
        panel.grid.minor = element_line(colour = "black"),
        legend.background = element_blank(), 
        legend.key = element_blank(), panel.background = element_blank(), 
        panel.border = element_blank(), strip.background = element_blank(), 
        plot.background = element_blank())
ggsave(paste0(folder,'/time_InstaPrism.pdf'))
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
ggsave(paste0(folder,'/time_ciber.pdf'))



