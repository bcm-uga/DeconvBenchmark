set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs
## ----
lot="Cobos"
source(paste0("../../src/0simu_dirichlet_functions_",lot,".R"))

## ----
## load data
## ----
print("-> Loading data...")
T_rna <- readRDS("../../../datashare/Cobos23/T_raw.rds")
T_rna <- T_rna[,sort(colnames(T_rna))]
row0 = apply(T_rna, 1, function(x) all(x==0))
T_rna <- T_rna[!row0,]
saveRDS(T_rna, file=paste0("simulations/rna/",today,"_",lot,"_T_rna_ref.rds"))
# K cell line
# Mesenchymal stem cells
# T cells
# K cell line
# K cell line
# monocytes

## ----
## generate simulated set
## ----
n_samples=120
params = data.frame(varCrit=10,sd_rna=1)
print(paste0("-> Generate training set..."))
n_rep=10
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  if (!file.exists(paste0("simulations/rna/",today,"_",lot,"_",sim_txt,i,".rds"))) {
    print(paste0("Simu n°",i))
    data_simu_clean0 <- generate_simu_set(T_rna, n_samples=n_samples, varCrit=params$varCrit)
    A <- data_simu_clean0[[2]]
    data_simu_clean <- data_simu_clean0[[1]]
    data_simu <- generate_simu_noise(data_simu_clean,
                                   p=0.1,sd_rna=params$sd_rna) # data raw
    saveRDS(list(D_rna_sim=data_simu,
                 A_ref=A), file=paste0("simulations/rna/",today,"_",lot,"_",sim_txt,i,".rds"))
  }
}

