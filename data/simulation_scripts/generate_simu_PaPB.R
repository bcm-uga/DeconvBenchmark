set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs & functions
## ----
source("generic_functions_sim.R")
generate_simu_pb <- function(data_sc, n_samples, alpha, varCrit=10) {
  celltypes = sort(unique(data_sclabel_full))
  num.cell.types <- length(celltypes)
  proportions <- generate_proportions((n_samples+30), celltypes, alpha, varCrit) #prop x sample (+ extra more samples in case the rounded prop selects no patients)
  
  data = matrix(nrow=nrow(data_sc),ncol=n_samples)
  sample_idx = 1
  sample_keep = c()
  for (i in seq(n_samples)) {
    prop = round(100*proportions[,sample_idx])
    # restrict to patients that have the correct number of cells
    # as advised in https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03292-w#Sec11
    patient_subset = unique(data_sc$orig.ident)[sapply(unique(data_sc$orig.ident), function(x)
      all(table(data_sc$orig.ident,data_sc$label_full)[x,]>prop))]
    while (length(patient_subset)==0) {
      print("prop not fitting with cell numbers in any patients")
      sample_idx = sample_idx+1
      prop = round(100*proportions[,sample_idx])
      patient_subset = unique(data_sc$orig.ident)[sapply(unique(data_sc$orig.ident), function(x)
        all(table(data_sc$orig.ident,data_sc$label_full)[x,]>prop))]
    }
    sample_keep = c(sample_keep,sample_idx)
    sample_idx = sample_idx+1
    # pick as many cells from each type from a single patient as needed and sum
    patient = sample(patient_subset, 1)
    sc_subset = subset(data_sc, subset = orig.ident %in% patient)
    sc_subset$cell = paste0("cell",seq(ncol(sc_subset)))
    cells <- unlist(mapply(function(prop,type) sample(sc_subset$cell[sc_subset$label_full==type], prop, replace = TRUE),
           prop[prop>0], names(prop)[prop>0], SIMPLIFY = F))
    mat = as.matrix(sc_subset@assays$RNA@counts[,sc_subset$cell %in% cells])
    colnames(mat) = unique(sc_subset$cell[sc_subset$cell %in% cells])
    data[,i] = rowSums(mat[,cells])
  }
  proportion_rounded = round(100*proportions[celltypes,sample_keep])/100
  colnames(data) = colnames(proportion_rounded)
  
  return(list(D_rna_sim=data, Amat=proportion_rounded))
}

## ----
## Load data using your own path and cell types
## ----
print("-> Loading data...")
sc_data = readRDS("~/projects/datashare/pancreas_deconv_data/results/20210624/02.1_output/peng_pool_labelled.rds")
cells_keep = c("B cells","T cells","Endothelial","Fibroblasts","Macrophages","Cancer basal","Cancer classical")
sc_data$label_full[sc_data$label_full == "Endothelial cells"] = "Endothelial"
sc_data$label_full[sc_data$label_full == "Fibro/Stellate"] = "Fibroblasts"
sc_data$label_full[sc_data$label_full == "Mono.Macro"] = "Macrophages"
sc_data$label_full[sc_data$label_full == "Basal"] = "Cancer basal"
sc_data$label_full[sc_data$label_full == "Classic"] = "Cancer classical"
sc_data <- subset(sc_data, subset = label_full %in% cells_keep)

T_all <- readRDS("../references/PaCL2.rds")
# Merge T cells and remove neutrophils
T_rna = data.frame(T_all[,c("B cells", "Cancer basal", "Cancer classical", "Endothelial", "Fibroblasts", "Macrophages")])
T_tcell = data.frame("T cells" = rowMeans(T_all[,c("CD4 T cells","CD8 T cells")]))
T_rna = cbind(T_rna,T_tcell)
T_rna <- T_rna[!apply(T_rna, 1, function(x) all(x==0)),]

# Common genes between bulk and sc
common_genes = intersect(rownames(sc_data),rownames(T_rna))
sc_data = sc_data[common_genes,]
T_rna = T_rna[common_genes,]

saveRDS(T_rna, file="../references/PaPB.rds")

## ----
## Fix parameters for the simulation
## ----
n_rep = 10
n_samples = 120
varCrit = 3
alpha = c(0.01, 0.07, 0.15, 0.46, 0.01, 0.30)
p = .1

## ----
## Generate simulations stored in a folder called 'simulations'
## ----
print(paste0("-> Generate simu..."))
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  print(paste0("Simu n°",i))
  data_simu_tot <- generate_simu_pb(sc_data, n_samples=n_samples, alpha=alpha, varCrit=params$varCrit)
  Amat <- data_simu_tot$Amat
  # rescale to same library size as Tref
  scale_fact = round(mean(colSums(T_rna)) / mean(colSums(data_simu_tot$D_rna_sim)))
  data_simu_tot$D_rna_sim = data_simu_tot$*D_rna_sim * scale_fact
  saveRDS(list(D_rna_sim = data_simu_tot$D_rna_sim,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_PaPB_",sim_txt,i,".rds"))
}
