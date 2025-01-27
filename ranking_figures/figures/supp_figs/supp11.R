## ----
## Set parameters, put your own
## ----
ref_path = "../../../data/references/"
dataset_rna = c('BrCL1','BrCL2','BlCL','PaCL1','PaCL2','BrMIX','PaMIX','BlREAL1','SkREAL')
dataset_dnam = c('BrCL1','LuCL','PaCL1','PaCL2','BlMIX','PaMIX','BlREAL2','BlREAL3')
folder = "supp11"

## ----
## Load libraries
## ----
library(dplyr)
library(circlize)
library(ComplexHeatmap)

## ----
## Load RNA reference profiles
## ----
# silico
names_silico = grep("CL", dataset_rna, value = T)
data_silico = sapply(names_silico, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_silico = lapply(data_silico, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("rna",x,value=T))}
})

# vitro
names_vitro = grep("MIX", dataset_rna, value = T)
data_vitro = sapply(names_vitro, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_vitro = lapply(data_vitro, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("rna",x,value=T))}
})

# vivo
names_vivo = grep("REAL", dataset_rna, value = T)
data_vivo = sapply(names_vivo, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_vivo = lapply(data_vivo, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("rna",x,value=T))}
})  

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

rownames(charac_rna) = c(names_silico,names_vitro,names_vivo)
saveRDS(charac_rna,paste0(folder,"/charac_rna.rds"))

## ----
## Load DNAm reference profiles
## ----
# silico
names_silico = grep("CL", dataset_dnam, value = T)
data_silico = sapply(names_silico, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_silico = lapply(data_silico, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("dnam",x,value=T))}
})

# vitro
names_vitro = grep("MIX", dataset_dnam, value = T)
data_vitro = sapply(names_vitro, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_vitro = lapply(data_vitro, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("dnam",x,value=T))}
})

# vivo
names_vivo = grep("REAL", dataset_dnam, value = T)
data_vivo = sapply(names_vivo, function(x)
  paste0(ref_path,list.files(ref_path,pattern=x)))
data_vivo = lapply(data_vivo, function(x) {
  if (length(x)==1) {readRDS(x)}
  else {readRDS(grep("dnam",x,value=T))}
})  

## ----
## Compute DNAm charac
## ----
data_all_dnam = c(data_silico,data_vitro,data_vivo)
charac_dnam = data.frame("nFeatures" = sapply(data_all_dnam, nrow),
                        "nCellTypes" = sapply(data_all_dnam, ncol),
                        "Mean R2" = pbapply::pbsapply(data_all_dnam, function(x) {
                          mean(sapply(seq(ncol(x)-1), function(y)
                            mean(sapply(seq(y+1,ncol(x)), function(z) cor(x[,y],x[,z])))))}),
                        "Phenotypic volume" = sapply(data_all_dnam, function(x) {
                          eign = eigen(cov(t(x[TOAST::findRefinx(as.matrix(x), nmarker=ncol(x)-1),])), symmetric = T, only.values = T)$values
                          log(prod(eign[eign!=0])/(ncol(x)-1))
                        }),
                        "Kurtosis" = pbapply::pbsapply(data_all_dnam, function(x) mean(apply(x,2,moments::kurtosis))),
                        "data_type" = c(rep("silico",length(data_silico)),
                                        rep("vitro",length(data_vitro)),
                                        rep("vivo",length(data_vivo))),
                        "Technology"=c('EPIC','800k','450k','800k',
                                       '450k','800k',
                                       '450k','450k'))
rownames(charac_dnam) = c(names_silico,names_vitro,names_vivo)
saveRDS(charac_dnam,paste0(folder,"/charac_dnam.rds"))

## ----
## Plot RNA heatmap (supp figure 11 panel B)
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
## Plot DNAm heatmap (supp figure 11 panel A)
## ----
col_fun_dnam = sapply(seq(ncol(charac_dnam)-2),function(x)
  colorRamp2(seq(from=min(charac_dnam[,x]),
                 to=max(charac_dnam[,x]),
                 length.out=9), RColorBrewer::brewer.pal(name="Purples",n=9)))
names(col_fun_dnam) = colnames(charac_dnam)[seq(ncol(charac_dnam)-2)]

names_ht_dnam = c("nFeatures","nCellTypes","Mean R2","Phenotypic volume","Kurtosis","Technology")
thd_dnam = c(6e5,8,.9,-35)
ht_dnam = list(Heatmap(charac_dnam %>% select(nFeatures),
                      col = col_fun_dnam$nFeatures,
                      name = "nFeatures",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_names_side = "left",
                      row_split = charac_dnam$data_type,
                      row_title = c("in silico"," in vitro","in vivo"),
                      row_title_gp = gpar(fontface='EUC', fontsize=10, fontfamily='HersheySerif',fill=c("firebrick","navy","green4"),col=c("black","grey90","black")),
                      row_order = c(order(rownames(charac_dnam)[charac_dnam$data_type=="silico"]),
                                    order(rownames(charac_dnam)[charac_dnam$data_type=="vitro"])+sum(charac_dnam$data_type=="silico"),
                                    order(rownames(charac_dnam)[charac_dnam$data_type=="vivo"])+sum(charac_dnam$data_type%in%c("silico","vitro"))),
                      column_names_side = 'bottom',
                      border = T,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(round(charac_dnam[i,j],3), x, y,
                                  gp=gpar(col=ifelse(round(charac_dnam[i,j],3)>thd_dnam[1],
                                                     'white','black')))
                      },))
ht_dnam = c(ht_dnam,lapply(seq(2,ncol(charac_dnam)-2), function(column)
  Heatmap(charac_dnam %>% select(colnames(charac_dnam)[column]),
          col = col_fun_dnam[[column]],
          name = names_ht_dnam[column],
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_order = c(order(rownames(charac_dnam)[charac_dnam$data_type=="silico"]),
                        order(rownames(charac_dnam)[charac_dnam$data_type=="vitro"])+sum(charac_dnam$data_type=="silico"),
                        order(rownames(charac_dnam)[charac_dnam$data_type=="vivo"])+sum(charac_dnam$data_type%in%c("silico","vitro"))),
          row_names_gp = gpar(fontsize=0),
          column_names_side = 'bottom',
          row_split = charac_dnam$data_type,
          border = T,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(charac_dnam[i,column],3), x, y,
                      gp=gpar(col=ifelse(round(charac_dnam[i,column],3)>thd_dnam[column],
                                         'white','black')))
          },)))
ht_dnam = c(ht_dnam,
           Heatmap(charac_dnam %>% select(Technology),
                   col = c("#DBDBDB","#6DC0E0","#92673C"),
                   name = "Technology",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_split = charac_dnam$data_type,
                   row_order = c(order(rownames(charac_dnam)[charac_dnam$data_type=="silico"]),
                                 order(rownames(charac_dnam)[charac_dnam$data_type=="vitro"])+sum(charac_dnam$data_type=="silico"),
                                 order(rownames(charac_dnam)[charac_dnam$data_type=="vivo"])+sum(charac_dnam$data_type%in%c("silico","vitro"))),
                   row_names_gp = gpar(fontsize=0),
                   column_names_side = 'bottom',
                   border = T,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(charac_dnam[i,7], x, y)
                   },))
pdf(paste0(folder,"/dnam.pdf"))
draw(ht_dnam[[1]]+ht_dnam[[2]]+ht_dnam[[3]]+ht_dnam[[4]]+ht_dnam[[6]],
     column_title = "DNAm data",
     column_title_gp = gpar(fontsize = 16))
dev.off()
