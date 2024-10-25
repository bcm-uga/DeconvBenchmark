
set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs
## ----
lot="He"
source(paste0("../../src/0simu_dirichlet_functions_",lot,".R"))

## ----
## load data
## ----
print("-> Loading data...")
T_met <- readRDS("../../../datashare/He22/T_met_raw.rds")
row0 = apply(T_met, 1, function(x) all(x==0))
T_met <- T_met[!row0,]
saveRDS(T_met, file=paste0("simulations/met/",today,"_",lot,"_T_met_ref.rds"))
# K cell line
# Normal cells
# Imm cells

## ----
## generate simulated set
## ----
n_samples=120
params = data.frame(varCrit=10,sd_met=3)
print(paste0("-> Generate training set..."))
n_rep=10
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  if (!file.exists(paste0("simulations/rna/",today,"_",lot,"_",sim_txt,i,".rds"))) {
    print(paste0("Simu n°",i))
    data_simu_clean0 <- generate_simu_set(T_met, n_samples=n_samples, varCrit=params$varCrit)
    A <- data_simu_clean0[[2]]
    data_simu_clean <- data_simu_clean0[[1]]
    data_simu <- generate_simu_noise(data_simu_clean,
                                   p=0.1,sd_met=params$sd_met) # data raw
    saveRDS(list(D_met_sim=data_simu,
                 A_ref=A), file=paste0("simulations/met/",today,"_",lot,"_",sim_txt,i,".rds"))
  }
}
