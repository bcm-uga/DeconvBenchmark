## ----
## Set parameters
## ----
library(dplyr)
library(ComplexHeatmap)
library(circlize)
folder = strsplit(basename(rstudioapi::getSourceEditorContext()$path),".R")[[1]]

## ----
## Load RNA reference profiles
## ----
# silico
lots_silico = list.files("../0simu/simulations/rna",
                         pattern=glob2rx(paste0("231027", "*ref.rds")))
names_silico = sapply(lots_silico, function(x)
  strsplit(x,"_")[[1]][2])
names_silico_CL = gsub("Cobos","BrCL2",
                       gsub("Hoek","BlCL",
                            gsub("dBREAST","BrCL1",
                                 gsub("lot1","PaCL2",
                                 gsub("dPANCREAS","PaCL1",names_silico)))))
data_silico = sapply(lots_silico, function(x) readRDS(paste0("../0simu/simulations/rna/",x)))

# vitro
lots_vitro = list.files("../1SB_invitro/data/", pattern='rna')
names_vitro = c("Cobos23","cometh_lot1/transcriptome")
data_vitro = sapply(names_vitro, function(x) readRDS(paste0("../../../datashare/",x,"/T_raw.rds")))

# vivo
lots_vivo = list.files("../1SB_invivo/data/", pattern='rna')
names_vivo = c("Deconer23/Linsley","Sturm19/Racle")
data_vivo = sapply(names_vivo, function(x) readRDS(paste0("../../../datashare/",x,"/T_rna_norm.rds")))
  
## ----
## Compute RNA charac
## ----
data_all_rna = c(data_silico,data_vitro,data_vivo)
charac_rna = data.frame("nFeatures" = sapply(data_all_rna, nrow),
                        "nCellTypes" = sapply(data_all_rna, ncol),
                        "Sparsity" = sapply(data_all_rna, function(x) mean(x==0)),
                        "Mean R2" = pbapply::pbsapply(data_all_rna, function(x)
                          mean(sapply(seq(ncol(x)-1), function(y)
                            mean(sapply(seq(y+1,ncol(x)), function(z) cor(x[,y],x[,z])))))),
                        "Phenotypic volume" = sapply(data_all_rna, function(x) {
                          eign = eigen(cov(t(x[TOAST::findRefinx(as.matrix(x), nmarker=ncol(x)-1),])), symmetric = T, only.values = T)$values
                          log(prod(eign[eign!=0])/(ncol(x)-1))
                          }),
                        "data_type" = c(rep("silico",length(data_silico)),
                                      rep("vitro",length(data_vitro)),
                                      rep("vivo",length(data_vivo))))
# Technology info is not critical as all in silico are Illumina

rownames(charac_rna) = c(names_silico_CL,
                         sapply(lots_vitro,function(x) strsplit(x,"_")[[1]][1]),
                         sapply(lots_vivo,function(x) strsplit(x,"_")[[1]][1]))
saveRDS(charac_rna,paste0(folder,"/charac_rna.rds"))

## ----
## Load DNAm reference profiles
## ----
# silico
lots_silico = list.files("../0simu/simulations/met",
                         pattern=glob2rx(paste0("231027", "*ref.rds")))
names_silico = sapply(lots_silico, function(x)
  strsplit(x,"_")[[1]][2])
names_silico_CL = gsub("He","LuCL",
                       gsub("dBREAST","BrCL1",
                            gsub("lot1","PaCL2",
                                 gsub("dPANCREAS","PaCL1",names_silico))))
data_silico = sapply(lots_silico, function(x) readRDS(paste0("../0simu/simulations/met/",x)))

# vitro
lots_vitro = list.files("../1SB_invitro/data/", pattern='met')
names_vitro = c("Koestler16","cometh_lot1/methylation")
data_vitro = sapply(names_vitro, function(x) readRDS(paste0("../../../datashare/",x,"/T_met_norm.rds")))

# vivo
lots_vivo = list.files("../1SB_invivo/data/", pattern='met')
names_vivo = c("Teschendorff17/Reinius","Teschendorff17/Liu")
data_vivo = sapply(names_vivo, function(x) readRDS(paste0("../../../datashare/",x,"/T_met_norm.rds")))

## ----
## Compute DNAm charac
## ----
data_all_met = c(data_silico,data_vitro,data_vivo)
charac_met = data.frame("nFeatures" = sapply(data_all_met, nrow),
                        "nCellTypes" = sapply(data_all_met, ncol),
                        "Mean R2" = pbapply::pbsapply(data_all_met, function(x) {
                          mean(sapply(seq(ncol(x)-1), function(y)
                            mean(sapply(seq(y+1,ncol(x)), function(z) cor(x[,y],x[,z])))))}),
                        "Phenotypic volume" = sapply(data_all_met, function(x) {
                          eign = eigen(cov(t(x[TOAST::findRefinx(as.matrix(x), nmarker=ncol(x)-1),])), symmetric = T, only.values = T)$values
                          log(prod(eign[eign!=0])/(ncol(x)-1))
                        }),
                        "Kurtosis" = sapply(data_all_met, function(x) mean(apply(x,2,moments::kurtosis))),
                        "data_type" = c(rep("silico",length(data_silico)),
                                        rep("vitro",length(data_vitro)),
                                        rep("vivo",length(data_vivo))),
                        "Technology"=c('EPIC','800k','450k','800k',
                                       '450k','800k',
                                       '450k','450k'))
rownames(charac_met) = c(names_silico_CL,
                         sapply(lots_vitro,function(x) strsplit(x,"_")[[1]][1]),
                         sapply(lots_vivo,function(x) strsplit(x,"_")[[1]][1]))
saveRDS(charac_met,paste0(folder,"/charac_met.rds"))

## ----
## Plot RNA heatmap
## ----
col_fun_rna = sapply(seq(ncol(charac_rna)-1),function(x)
  colorRamp2(seq(from=min(charac_rna[,x]),
                 to=max(charac_rna[,x]),
                 length.out=9), RColorBrewer::brewer.pal(name="Purples",n=9)))
names(col_fun_rna) = colnames(charac_rna)[seq(ncol(charac_rna)-1)]

names_ht_rna = c("nFeatures","nCellTypes","Sparsity","Mean R2","Phenotypic volume")
thd_rna = c(25000,8,.3,.8,100)
ht_rna = list(Heatmap(charac_rna %>% select(nFeatures),
              col = col_fun_rna$nFeatures,
              name = "nFeatures",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_side = "left",
              row_split = charac_rna$data_type,
              row_title = c("in silico"," in vitro","in vivo"),
              row_title_gp = gpar(fontface='EUC', fontsize=10, fontfamily='HersheySerif',fill=c("firebrick","navy","green4"),col=c("black","grey90","black")),
              row_order = c(order(rownames(charac_rna)[charac_rna$data_type=="silico"]),
                            order(rownames(charac_rna)[charac_rna$data_type=="vitro"])+sum(charac_rna$data_type=="silico"),
                            order(rownames(charac_rna)[charac_rna$data_type=="vivo"])+sum(charac_rna$data_type%in%c("silico","vitro"))),
              column_names_side = 'bottom',
              border = T,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(charac_rna[i,j],3), x, y,
                          gp=gpar(col=ifelse(round(charac_rna[i,j],3)>thd_rna[1],
                                             'white','black')))
              },))
ht_rna = c(ht_rna,lapply(seq(2,ncol(charac_rna)-1), function(column)
  Heatmap(charac_rna %>% select(colnames(charac_rna)[column]),
              col = col_fun_rna[[column]],
              name = names_ht_rna[column],
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_order = c(order(rownames(charac_rna)[charac_rna$data_type=="silico"]),
                            order(rownames(charac_rna)[charac_rna$data_type=="vitro"])+sum(charac_rna$data_type=="silico"),
                            order(rownames(charac_rna)[charac_rna$data_type=="vivo"])+sum(charac_rna$data_type%in%c("silico","vitro"))),
              row_names_gp = gpar(fontsize=0),
              column_names_side = 'bottom',
              row_split = charac_rna$data_type,
              border = T,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(charac_rna[i,column],3), x, y,
                          gp=gpar(col=ifelse(round(charac_rna[i,column],3)>thd_rna[column],
                                             'white','black')))
              },)))
pdf(paste0(folder,"/rna.pdf"))
draw(ht_rna[[1]]+ht_rna[[2]]+ht_rna[[3]]+ht_rna[[4]]+ht_rna[[5]],
     column_title = "RNA data",
     column_title_gp = gpar(fontsize = 16))
dev.off()

## ----
## Plot DNAm heatmap
## ----
col_fun_met = sapply(seq(ncol(charac_met)-2),function(x)
  colorRamp2(seq(from=min(charac_met[,x]),
                 to=max(charac_met[,x]),
                 length.out=9), RColorBrewer::brewer.pal(name="Purples",n=9)))
names(col_fun_met) = colnames(charac_met)[seq(ncol(charac_met)-2)]

names_ht_met = c("nFeatures","nCellTypes","Mean R2","Phenotypic volume","Kurtosis","Technology")
thd_met = c(6e5,8,.9,-35)
ht_met = list(Heatmap(charac_met %>% select(nFeatures),
                      col = col_fun_met$nFeatures,
                      name = "nFeatures",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_names_side = "left",
                      row_split = charac_met$data_type,
                      row_title = c("in silico"," in vitro","in vivo"),
                      row_title_gp = gpar(fontface='EUC', fontsize=10, fontfamily='HersheySerif',fill=c("firebrick","navy","green4"),col=c("black","grey90","black")),
                      row_order = c(order(rownames(charac_met)[charac_met$data_type=="silico"]),
                                    order(rownames(charac_met)[charac_met$data_type=="vitro"])+sum(charac_met$data_type=="silico"),
                                    order(rownames(charac_met)[charac_met$data_type=="vivo"])+sum(charac_met$data_type%in%c("silico","vitro"))),
                      column_names_side = 'bottom',
                      border = T,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(round(charac_met[i,j],3), x, y,
                                  gp=gpar(col=ifelse(round(charac_met[i,j],3)>thd_met[1],
                                                     'white','black')))
                      },))
ht_met = c(ht_met,lapply(seq(2,ncol(charac_met)-2), function(column)
  Heatmap(charac_met %>% select(colnames(charac_met)[column]),
          col = col_fun_met[[column]],
          name = names_ht_met[column],
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_order = c(order(rownames(charac_met)[charac_met$data_type=="silico"]),
                        order(rownames(charac_met)[charac_met$data_type=="vitro"])+sum(charac_met$data_type=="silico"),
                        order(rownames(charac_met)[charac_met$data_type=="vivo"])+sum(charac_met$data_type%in%c("silico","vitro"))),
          row_names_gp = gpar(fontsize=0),
          column_names_side = 'bottom',
          row_split = charac_met$data_type,
          border = T,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(charac_met[i,column],3), x, y,
                      gp=gpar(col=ifelse(round(charac_met[i,column],3)>thd_met[column],
                                         'white','black')))
          },)))
ht_met = c(ht_met,
           Heatmap(charac_met %>% select(Technology),
                   col = c("#DBDBDB","#6DC0E0","#92673C"),
                   name = "Technology",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_split = charac_met$data_type,
                   row_order = c(order(rownames(charac_met)[charac_met$data_type=="silico"]),
                                 order(rownames(charac_met)[charac_met$data_type=="vitro"])+sum(charac_met$data_type=="silico"),
                                 order(rownames(charac_met)[charac_met$data_type=="vivo"])+sum(charac_met$data_type%in%c("silico","vitro"))),
                   row_names_gp = gpar(fontsize=0),
                   column_names_side = 'bottom',
                   border = T,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(charac_met[i,7], x, y)
                   },))
pdf(paste0(folder,"/dnam.pdf"))
draw(ht_met[[1]]+ht_met[[2]]+ht_met[[3]]+ht_met[[4]]+ht_met[[6]],
     column_title = "DNAm data",
     column_title_gp = gpar(fontsize = 16))
dev.off()
