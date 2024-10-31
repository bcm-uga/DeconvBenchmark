## ----
## Set parameters, put your own
## ----
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

## ----
## Load libraries
## ----
library(dplyr)
library(ComplexHeatmap)
library(circlize)

## ----
## Create methods charac, add your methods if needed
## ----
charac_meth = data.frame("Supervised"=c(rep("Supervised",4),
                                        rep("Unsupervised",3),
                                        rep("Supervised",1),
                                        rep("Unsupervised",3),
                                        rep("Supervised",7),
                                        rep("Unsupervised",2)),
                         "Omic"=c(rep("DNAm & RNA",7),
                                   rep("DNAm",4),
                                   rep("RNA",9)),
                         "Approach"=c("RPC-LS","SVR","Bayesian","CLS",
                                      "ICA","NMF","Convex",
                                      "CLS",
                                      "NMF","NMF","NMF",
                                      "SVR","LS",'LS','CLS','CLS','CLS','CLS',
                                      "NMF","LDA"),
                         "STO"=c(rep("Equality",7),"Inequality",
                                 rep("Equality",8),"Inequality",
                                 "Equality","Inequality","Equality"))
rownames(charac_meth) = c("RLR","CIBERSORT","InstaPrism","NNLS",
                          "ICA","NMF","debCAM",
                          "epidishCP",
                          "MeDeCom","RefFreeEWAS","EDec",
                          "SVR","Elastic net","OLS","WISP","DeconRNASeq","FARDEEP","FARDEEP_sto",
                          "PREDE","CDSeq")
charac_meth$Omic = factor(charac_meth$Omic, levels=unique(charac_meth$Omic))
saveRDS(charac_meth, file = paste0(folder,'/charac_methods.rds'))

## ----
## Plot heatmap supp figure 7
## ----
col_fun_sup = c(rep("white",length(charac_meth %>%
                                                 filter(Supervised=="Supervised") %>%
                                                 pull(Approach) %>% unique())),
                rep(c("olivedrab3","brown2"),2))
names(col_fun_sup) = c(charac_meth %>%
                         filter(Supervised=="Supervised") %>%
                         pull(Approach) %>% unique(),
                       rep(c("Equality","Inequality"),2))
col_fun_unsup = c(rep("white",length(charac_meth %>%
                                               filter(Supervised=="Unsupervised") %>%
                                               pull(Approach) %>% unique())),
                rep(c("olivedrab3","brown2"),2))
names(col_fun_unsup) = c(charac_meth %>%
                         filter(Supervised=="Unsupervised") %>%
                         pull(Approach) %>% unique(),
                       rep(c("Equality","Inequality"),2))

ht_sup = Heatmap(charac_meth %>%
                select(!Omic) %>%
                filter(Supervised=="Supervised") %>%
                select(!Supervised),
             col = col_fun_sup,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_order = seq(nrow(charac_meth %>% filter(Supervised=="Supervised"))),
             row_names_side = "left",
             row_split = charac_meth$Omic[charac_meth$Supervised=="Supervised"],
             row_title = unique(charac_meth$Omic[charac_meth$Supervised=="Supervised"]),
             row_title_gp = gpar(fontface='EUC', fontsize=10, fontfamily='HersheySerif',fill=c("palegreen1","lightgoldenrod1","lightblue2")),
             column_split = colnames(charac_meth %>% select(!Omic) %>% select(!Supervised)),
             column_title_gp = gpar(fontsize=0),
             show_heatmap_legend = F,
             border = T,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(charac_meth[charac_meth$Supervised=="Supervised",!colnames(charac_meth)%in%c("Omic","Supervised")][i,j], x, y)
             },)
ht_unsup = Heatmap(charac_meth %>%
                     select(!Omic) %>%
                     filter(Supervised=="Unsupervised") %>%
                     select(!Supervised),
                   col = col_fun_unsup,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_order = seq(nrow(charac_meth %>% filter(Supervised=="Unsupervised"))),
                   row_names_side = "left",
                   row_split = charac_meth$Omic[charac_meth$Supervised=="Unsupervised"],
                   row_title = unique(charac_meth$Omic[charac_meth$Supervised=="Unsupervised"]),
                   row_title_gp = gpar(fontface='EUC', fontsize=10, fontfamily='HersheySerif',fill=c("palegreen1","lightgoldenrod1","lightblue2")),
                   column_split = colnames(charac_meth %>% select(!Omic) %>% select(!Supervised)),
                   column_title_gp = gpar(fontsize=0),
                   show_heatmap_legend = F,
                   border = T,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(charac_meth[charac_meth$Supervised=="Unsupervised",!colnames(charac_meth)%in%c("Omic","Supervised")][i,j], x, y)
                   },)
pdf(paste0(folder,"/sup.pdf"))
draw(ht_sup,
     column_title = "Supervised methods",
     column_title_gp = gpar(fontsize = 16))
dev.off()
pdf(paste0(folder,"/unsup.pdf"))
draw(ht_unsup,
     column_title = "Unsupervised methods",
     column_title_gp = gpar(fontsize = 16))
dev.off()

