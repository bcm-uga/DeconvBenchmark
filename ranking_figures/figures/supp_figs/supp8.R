## ----
## Set parameters, put your own
## ----
# choose a dataset, here BrCL1_sim01
data_path = "../../../data/simulations/rna/"
dataset = "BrCL1_sim01"
folder = "supp8"
prop_ref0 = c(0.45, 0.1, 0.15, 0.3)
alpha0 = c("3","10","30")
alpha0_names = c("241025MoreDisp","241025","241025LessDisp")

## ----
## Load libraries
## ----
library(ggplot2)
library(see)
library(dplyr)

## ----
## Load  sims
## ----
sims = lapply(alpha0_names, function(x)
  readRDS(paste0(data_path,x,"_",dataset,".rds"))$A_ref)

## ----
## Plot and conclude
## ----
df = do.call(rbind,lapply(seq_along(sims), function(x)
  data.frame("Celltype"=rep(rownames(sims[[x]]),ncol(sims[[x]])),
             "Proportions"=c(sims[[x]]),
             "Alpha0"=alpha0[x])))
df$Alpha0 = factor(df$Alpha0, levels=rev(alpha0))
df$Celltype = factor(df$Celltype, levels=c(unique(df$Celltype)[order(prop_ref0)]))
prop_ref = sort(prop_ref0)

p = ggplot(df, aes(x=Celltype, y=Proportions)) +
    geom_violin(aes(fill=Alpha0)) +
    scale_fill_viridis_d(option = 'plasma', name=expression(alpha[0])) +
    theme_modern(axis.text.angle = 20, axis.text.size = 15, legend.title.size = 20, legend.text.size = 18, axis.title.size = 18) +
    xlab("")
for (j in seq_along(df %>% pull(Celltype) %>% unique())) {
  p = p +
    geom_segment(x=j-.4,xend=j+.4,y=prop_ref[j],yend=prop_ref[j], color="grey30", linewidth=.3, linetype='dashed')
}
p
ggsave(paste0(folder,"/boxplot_BrCL1.pdf"), width=6, height=3.5)
