## ----
## Add libs and params
## ----
library(ggplot2)
library(see)
library(dplyr)
library(ggpubr)
prop_ref0 = list(dBREAST=c(0.45, 0.1, 0.15, 0.3),
                 dPANCREAS=c(0.45, 0.1, 0.15, 0.3),
                 lot1=c(0.01, 0.04, 0.03, 0.15, 0.46, 0.01, 0.01, 0.29))
alpha0 = c("3","10","30")
alpha0_names = c("240130LessNoise","231027","240130MoreNoise")
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

## ----
## Load  sims
## ----
sims = lapply(c("dBREAST","dPANCREAS","lot1"), function(lot)
  lapply(alpha0_names, function(x)
  readRDS(paste0("../0simu/simulations/rna/",x,"_",lot,"_sim01.rds"))$A_ref))

## ----
## Merge tumor types for PDAC datasets -dPANCREAS&lot1
## ----
sims[[2]] = lapply(sims[[2]], function(x) {
  normal = x[-grep("TUM",rownames(x)),]
  TUM = colSums(x[grep("TUM",rownames(x)),])
  rbind(normal,TUM)
})
sims[[3]] = lapply(sims[[3]], function(x) {
  normal = x[-grep("Cancer",rownames(x)),]
  Cancer = colSums(x[grep("Cancer",rownames(x)),])
  rbind(normal,Cancer)
})

## ----
## Plot and conclude
## ----
df = do.call(rbind,mapply(function(x,y)
  data.frame("Celltype"=rep(rownames(x[[1]]),ncol(x[[1]])),
             "Proportions"=unlist(lapply(x,c)),
             "Alpha0"=c(sapply(seq_along(alpha0),function(z)
               rep(alpha0[z],length(x[[z]])))),
             "Dataset"=y),
  x=sims,y=names(prop_ref0), SIMPLIFY=F))
df$Alpha0 = factor(df$Alpha0, levels=rev(alpha0))
df$Celltype = factor(df$Celltype, levels=c(unique(df$Celltype[df$Dataset==names(prop_ref0)[1]])[order(prop_ref0[[1]])],
                                           unique(df$Celltype[df$Dataset==names(prop_ref0)[2]])[order(prop_ref0[[2]])],
                                           unique(df$Celltype[df$Dataset==names(prop_ref0)[3]])[order(prop_ref0[[3]])]))
prop_ref = lapply(prop_ref0,sort)

plot = list()
for (i in seq_along(unique(df$Dataset))) {
  plot[[i]] = ggplot(df %>% filter(Dataset==unique(df$Dataset)[i]), aes(x=Celltype, y=Proportions)) +
    geom_violin(aes(fill=Alpha0)) +
    scale_fill_viridis_d(option = 'plasma', name=expression(alpha[0])) +
    theme_modern(axis.text.angle = 20, axis.text.size = 15, legend.title.size = 20, legend.text.size = 18, axis.title.size = 18) +
    xlab("")
  for (j in seq_along(df %>% filter(Dataset==unique(df$Dataset)[i]) %>% pull(Celltype) %>% unique())) {
    plot[[i]] = plot[[i]] +
    geom_segment(x=j-.4,xend=j+.4,y=prop_ref[[i]][j],yend=prop_ref[[i]][j], color="grey30", linewidth=.3, linetype='dashed')
  }
}
ggarrange(plotlist=plot, ncol=1)
plot[[1]]
ggsave(paste0(folder,"/boxplot_dBREAST.pdf"), width=6, height=3.5)
plot[[2]]
ggsave(paste0(folder,"/boxplot_dPANCREAS.pdf"), width=6, height=3.5)
plot[[3]]
ggsave(paste0(folder,"/boxplot_lot1.pdf"), width=6, height=3.5)
